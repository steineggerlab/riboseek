#include "CommandDeclarations.h"
#include "Debug.h"
#include "MMseqsMPI.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "LinearFoldBridge.h"
#ifdef HAVE_INFERNAL_BRIDGE
#include "infernal/InfernalBridge.h"
#endif

#include <array>
#include <algorithm>
#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Native generatecm command:
// maps DB -> DB covariance-model generation to MMseqs profile construction.
int generatecm(int argc, const char **argv, const Command &command) {
    return result2profile(argc, argv, command);
}

namespace {

struct AlnSeq {
    std::string id;
    std::string aln;
};

static inline std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])) != 0) {
        ++b;
    }
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) != 0) {
        --e;
    }
    return s.substr(b, e - b);
}

static inline char normalizeBase(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (c == 'T') {
        return 'U';
    }
    return c;
}

static inline bool isGap(char c) {
    return c == '-' || c == '.';
}

static inline int baseToIdx(char c) {
    c = normalizeBase(c);
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'U') return 3;
    return -1;
}

static bool parseFastaAlignment(const std::string &path, std::vector<AlnSeq> &seqs) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        return false;
    }
    std::string line;
    AlnSeq cur;
    bool sawHeader = false;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            sawHeader = true;
            if (!cur.id.empty()) {
                seqs.push_back(cur);
                cur = AlnSeq();
            }
            cur.id = trim(line.substr(1));
            continue;
        }
        if (!sawHeader) {
            return false;
        }
        std::string s = trim(line);
        for (size_t i = 0; i < s.size(); ++i) {
            const char c = s[i];
            if (std::isspace(static_cast<unsigned char>(c)) != 0) {
                continue;
            }
            cur.aln.push_back(c);
        }
    }
    if (!cur.id.empty()) {
        seqs.push_back(cur);
    }
    return !seqs.empty();
}

static bool parseStockholmAlignment(const std::string &path, std::vector<AlnSeq> &seqs, std::string &ssCons, std::string &rfCons) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        return false;
    }
    std::string first;
    if (!std::getline(in, first)) {
        return false;
    }
    first = trim(first);
    if (first.find("# STOCKHOLM") != 0) {
        return false;
    }

    std::map<std::string, std::string> acc;
    std::vector<std::string> order;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line == "//") {
            break;
        }
        if (line[0] == '#') {
            if (line.find("#=GC") == 0) {
                std::stringstream ss(line);
                std::string t0, t1, t2;
                ss >> t0 >> t1 >> t2;
                if (t1 == "SS_cons" && !t2.empty()) {
                    ssCons += t2;
                } else if (t1 == "RF" && !t2.empty()) {
                    rfCons += t2;
                }
            }
            continue;
        }
        std::stringstream ss(line);
        std::string name, chunk;
        ss >> name >> chunk;
        if (name.empty() || chunk.empty()) {
            continue;
        }
        if (acc.find(name) == acc.end()) {
            order.push_back(name);
        }
        acc[name] += chunk;
    }

    for (size_t i = 0; i < order.size(); ++i) {
        AlnSeq s;
        s.id = order[i];
        s.aln = acc[order[i]];
        seqs.push_back(s);
    }
    return !seqs.empty();
}

static void ensureValidAlignment(const std::vector<AlnSeq> &seqs) {
    if (seqs.empty()) {
        Debug(Debug::ERROR) << "cmbuild: alignment has no sequences\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t L = seqs[0].aln.size();
    if (L == 0) {
        Debug(Debug::ERROR) << "cmbuild: alignment has zero length\n";
        EXIT(EXIT_FAILURE);
    }
    for (size_t i = 1; i < seqs.size(); ++i) {
        if (seqs[i].aln.size() != L) {
            Debug(Debug::ERROR) << "cmbuild: alignment sequences must have equal length\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

static std::vector<int> parseSsPairs(const std::string &ss) {
    std::vector<int> pair(ss.size(), -1);
    std::vector<int> st;
    st.reserve(ss.size());
    for (size_t i = 0; i < ss.size(); ++i) {
        const char c = ss[i];
        if (c == '(' || c == '<') {
            st.push_back(static_cast<int>(i));
        } else if (c == ')' || c == '>') {
            if (st.empty()) {
                Debug(Debug::WARNING) << "cmbuild: unbalanced SS_cons, ignoring structure constraints\n";
                return std::vector<int>();
            }
            const int j = st.back();
            st.pop_back();
            pair[static_cast<size_t>(j)] = static_cast<int>(i);
            pair[i] = j;
        } else if (c == '.' || c == ':' || c == '_' || c == '-' || c == ',') {
            continue;
        } else {
            Debug(Debug::WARNING) << "cmbuild: unsupported SS_cons symbol '" << c
                                  << "', treating as unpaired\n";
        }
    }
    if (!st.empty()) {
        Debug(Debug::WARNING) << "cmbuild: unbalanced SS_cons, ignoring structure constraints\n";
        return std::vector<int>();
    }
    return pair;
}

static std::string buildConsensusForColumns(const std::vector<AlnSeq> &seqs,
                                            const std::vector<int> &keepCols) {
    std::string cons;
    cons.reserve(keepCols.size());
    for (size_t k = 0; k < keepCols.size(); ++k) {
        int cnt[4] = {0, 0, 0, 0};
        for (size_t s = 0; s < seqs.size(); ++s) {
            const int bi = baseToIdx(seqs[s].aln[static_cast<size_t>(keepCols[k])]);
            if (bi >= 0) {
                cnt[bi] += 1;
            }
        }
        int bestIdx = 0;
        int bestCnt = cnt[0];
        for (int b = 1; b < 4; ++b) {
            if (cnt[b] > bestCnt) {
                bestCnt = cnt[b];
                bestIdx = b;
            }
        }
        static const char bases[4] = {'A', 'C', 'G', 'U'};
        cons.push_back(bases[bestIdx]);
    }
    return cons;
}

static bool predictRawPairsLinearFold(const std::vector<AlnSeq> &seqs,
                                      const std::vector<int> &keepRaw,
                                      size_t Lraw,
                                      std::vector<int> &rawPairs,
                                      std::string &ssConsOut) {
    rawPairs.assign(Lraw, -1);
    ssConsOut.assign(Lraw, '.');
    if (keepRaw.empty()) {
        return false;
    }
    const std::string consensus = buildConsensusForColumns(seqs, keepRaw);
    if (consensus.size() < 2) {
        return false;
    }
    std::string pred;
    double score = 0.0;
    if (!linearFoldPredictDotBracket(consensus, pred, &score)) {
        return false;
    }
    if (pred.size() != keepRaw.size()) {
        return false;
    }
    for (size_t i = 0; i < keepRaw.size(); ++i) {
        ssConsOut[static_cast<size_t>(keepRaw[i])] = pred[i];
    }
    const std::vector<int> kp = parseSsPairs(pred);
    if (kp.empty()) {
        return false;
    }
    for (size_t i = 0; i < kp.size(); ++i) {
        const int j = kp[i];
        if (j < 0) {
            continue;
        }
        const int rawI = keepRaw[i];
        const int rawJ = keepRaw[static_cast<size_t>(j)];
        rawPairs[static_cast<size_t>(rawI)] = rawJ;
    }
    return true;
}

static std::string buildStockholmText(const std::string &id,
                                      const std::vector<AlnSeq> &seqs,
                                      const std::string &ssCons) {
    std::ostringstream out;
    out << "# STOCKHOLM 1.0\n";
    out << "#=GF ID " << (id.empty() ? "mmseqs_model" : id) << "\n";
    for (size_t i = 0; i < seqs.size(); ++i) {
        out << seqs[i].id << " " << seqs[i].aln << "\n";
    }
    if (!ssCons.empty()) {
        out << "#=GC SS_cons " << ssCons << "\n";
    }
    out << "//\n";
    return out.str();
}

static std::string ntName(const std::string &prefix, int l, int r) {
    std::stringstream ss;
    ss << prefix << "_" << l << "_" << r;
    return ss.str();
}

struct MidEstimates {
    std::vector<double> matchOcc; // 1-based
    std::vector<double> delOcc;   // 1-based
    std::vector<double> tMM;      // 1-based edge p->p+1
    std::vector<double> tMD;      // 1-based edge p->p+1
    std::vector<double> tDM;      // 1-based edge p->p+1
    std::vector<double> tDD;      // 1-based edge p->p+1
    std::map< std::pair<int, int>, double > pairMatchOcc;
    std::map< std::pair<int, int>, double > pairDelOcc;
};

static void emitUnaryRules(std::ostream &out,
                           const std::string &lhs,
                           const double cnt[4],
                           const double bg[4],
                           double pc) {
    const double unaryWeight = 1.35;
    double s = 0.0;
    for (int b = 0; b < 4; ++b) {
        s += cnt[b] + pc;
    }
    static const char bases[4] = {'A', 'C', 'G', 'U'};
    for (int b = 0; b < 4; ++b) {
        const double p = (cnt[b] + pc) / s;
        const double q = std::max(bg[b], 1e-12);
        out << "UNARY " << lhs << " " << bases[b] << " " << (unaryWeight * std::log2(p / q)) << "\n";
    }
}

static void emitPairRules(std::ostream &out,
                          const std::string &lhs,
                          const std::string &innerOrDash,
                          const double cnt[16],
                          const double bg[4],
                          double pc,
                          double transLogp) {
    const double pairWeight = 2.30;
    double s = 0.0;
    for (int k = 0; k < 16; ++k) {
        s += cnt[k] + pc;
    }
    static const char bases[4] = {'A', 'C', 'G', 'U'};
    for (int l = 0; l < 4; ++l) {
        for (int r = 0; r < 4; ++r) {
            const int idx = l * 4 + r;
            const double p = (cnt[idx] + pc) / s;
            const double q = std::max(bg[l] * bg[r], 1e-12);
            out << "PAIR " << lhs << " " << innerOrDash << " " << bases[l] << " " << bases[r] << " " << (pairWeight * std::log2(p / q) + transLogp) << "\n";
        }
    }
}

static void buildInterval(std::ostream &out,
                          int l,
                          int r,
                          const std::vector< std::array<double, 4> > &singleCnt,
                          const std::map< std::pair<int, int>, std::array<double, 16> > &pairCnt,
                          const std::vector<int> &pairMap,
                          const MidEstimates &mid,
                          const double bg[4],
                          double nEffSeq,
                          const std::string &prefix,
                          double pc,
                          bool allowPairSideIndel,
                          std::vector< std::vector<unsigned char> > &seen) {
    if (l > r || seen[static_cast<size_t>(l)][static_cast<size_t>(r)] != 0) {
        return;
    }
    seen[static_cast<size_t>(l)][static_cast<size_t>(r)] = 1;

    const std::string lhs = ntName(prefix, l, r);
    if (l == r) {
        emitUnaryRules(out, lhs, singleCnt[static_cast<size_t>(l)].data(), bg, pc);
        return;
    }

    const int p = pairMap.empty() ? -1 : pairMap[static_cast<size_t>(l)];
    if (p >= l && p <= r) {
        if (p == r) {
            std::map< std::pair<int, int>, std::array<double, 16> >::const_iterator it = pairCnt.find(std::make_pair(l, r));
            std::array<double, 16> empty = {};
            const std::array<double, 16> &cnt = (it == pairCnt.end()) ? empty : it->second;
            const std::map< std::pair<int, int>, double >::const_iterator mit = mid.pairMatchOcc.find(std::make_pair(l, r));
            const double pairOcc = (mit != mid.pairMatchOcc.end()) ? mit->second : 0.5;
            const double pairTrans = std::log2(std::max(pairOcc, 1e-6));
            if (l + 1 <= r - 1) {
                const std::string inner = ntName(prefix, l + 1, r - 1);
                const std::string leftNt = ntName(prefix, l, l);
                const std::string rightNt = ntName(prefix, r, r);
                buildInterval(out, l + 1, r - 1, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
                if (allowPairSideIndel) {
                    buildInterval(out, l, l, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
                    buildInterval(out, r, r, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
                }
                emitPairRules(out, lhs, inner, cnt.data(), bg, pc, pairTrans);
                // Allow deleting an entire paired consensus column pair with a strong penalty.
                // This approximates Infernal-like delete flexibility while keeping structure preference.
                const std::map< std::pair<int, int>, double >::const_iterator dit = mid.pairDelOcc.find(std::make_pair(l, r));
                const double pDelPair = (dit != mid.pairDelOcc.end()) ? dit->second : 0.01;
                out << "EPS " << lhs << " " << inner << " " << std::log2(std::max(pDelPair, 1e-6)) << "\n";
                if (allowPairSideIndel) {
                    // One-sided paired-column indel flexibility (strongly penalized).
                    out << "BINARY " << lhs << " " << leftNt << " " << inner << " -5.5\n";
                    out << "BINARY " << lhs << " " << inner << " " << rightNt << " -5.5\n";
                }
            } else {
                emitPairRules(out, lhs, "-", cnt.data(), bg, pc, pairTrans);
                if (allowPairSideIndel) {
                    // Adjacent pair fallback: allow one side to remain (strongly penalized).
                    const std::string leftNt = ntName(prefix, l, l);
                    const std::string rightNt = ntName(prefix, r, r);
                    buildInterval(out, l, l, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
                    buildInterval(out, r, r, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
                    out << "EPS " << lhs << " " << leftNt << " -6.0\n";
                    out << "EPS " << lhs << " " << rightNt << " -6.0\n";
                }
            }
            return;
        }
        const std::string left = ntName(prefix, l, p);
        const std::string right = ntName(prefix, p + 1, r);
        buildInterval(out, l, p, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
        buildInterval(out, p + 1, r, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
        out << "BINARY " << lhs << " " << left << " " << right << " -0.3\n";
        return;
    }

    const std::string left = ntName(prefix, l, l);
    const std::string right = ntName(prefix, l + 1, r);
    buildInterval(out, l, l, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
    buildInterval(out, l + 1, r, singleCnt, pairCnt, pairMap, mid, bg, nEffSeq, prefix, pc, allowPairSideIndel, seen);
    const double pMatch = (l < static_cast<int>(mid.matchOcc.size())) ? mid.matchOcc[static_cast<size_t>(l)] : 0.5;
    const double pDel = (l < static_cast<int>(mid.delOcc.size())) ? mid.delOcc[static_cast<size_t>(l)] : 0.1;
    double pKeep = pMatch;
    double pSkip = pDel;
    if (l < static_cast<int>(mid.tMM.size())) {
        pKeep *= std::max(mid.tMM[static_cast<size_t>(l)], 1e-6);
    }
    if (l < static_cast<int>(mid.tDM.size())) {
        pSkip *= std::max(mid.tDM[static_cast<size_t>(l)] + mid.tDD[static_cast<size_t>(l)], 1e-6);
    }
    out << "BINARY " << lhs << " " << left << " " << right << " " << std::log2(std::max(pKeep, 1e-8)) << "\n";
    out << "EPS " << lhs << " " << right << " " << std::log2(std::max(pSkip, 1e-8)) << "\n";
}

static double percentile(std::vector<double> v, double q) {
    if (v.empty()) {
        return 0.0;
    }
    if (q <= 0.0) {
        return *std::min_element(v.begin(), v.end());
    }
    if (q >= 1.0) {
        return *std::max_element(v.begin(), v.end());
    }
    const size_t k = static_cast<size_t>(q * static_cast<double>(v.size() - 1));
    std::nth_element(v.begin(), v.begin() + static_cast<ptrdiff_t>(k), v.end());
    return v[k];
}

} // namespace

// Native cmbuild command:
// Build a native CM grammar from aligned RNA sequences (Stockholm or aligned FASTA).
int cmbuild(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (MMseqsMPI::isMaster() == false) {
        return EXIT_SUCCESS;
    }

    std::vector<AlnSeq> seqs;
    std::string ssCons;
    std::string rfCons;
    if (!parseStockholmAlignment(par.db1, seqs, ssCons, rfCons)) {
        seqs.clear();
        ssCons.clear();
        rfCons.clear();
        if (!parseFastaAlignment(par.db1, seqs)) {
            Debug(Debug::ERROR) << "cmbuild: failed reading input alignment as Stockholm or FASTA: " << par.db1 << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    ensureValidAlignment(seqs);

    const size_t Lraw = seqs[0].aln.size();
    std::vector<double> seqWeight(seqs.size(), 0.0);
    {
        // Henikoff position-based sequence weights on raw alignment columns.
        for (size_t c = 0; c < Lraw; ++c) {
            int counts[4] = {0, 0, 0, 0};
            int distinct = 0;
            for (size_t s = 0; s < seqs.size(); ++s) {
                const int bi = baseToIdx(seqs[s].aln[c]);
                if (bi >= 0) {
                    if (counts[bi] == 0) {
                        distinct += 1;
                    }
                    counts[bi] += 1;
                }
            }
            if (distinct == 0) {
                continue;
            }
            for (size_t s = 0; s < seqs.size(); ++s) {
                const int bi = baseToIdx(seqs[s].aln[c]);
                if (bi >= 0 && counts[bi] > 0) {
                    seqWeight[s] += 1.0 / (static_cast<double>(distinct) * static_cast<double>(counts[bi]));
                }
            }
        }
        double sumW = 0.0;
        for (size_t s = 0; s < seqWeight.size(); ++s) {
            sumW += seqWeight[s];
        }
        if (sumW <= 0.0) {
            for (size_t s = 0; s < seqWeight.size(); ++s) {
                seqWeight[s] = 1.0;
            }
        } else {
            const double scale = static_cast<double>(seqs.size()) / sumW;
            for (size_t s = 0; s < seqWeight.size(); ++s) {
                seqWeight[s] *= scale;
            }
        }
    }
    double nEffSeq = 0.0;
    for (size_t s = 0; s < seqWeight.size(); ++s) {
        nEffSeq += seqWeight[s];
    }
    std::vector< std::vector<int8_t> > alnCode(seqs.size(), std::vector<int8_t>(Lraw, static_cast<int8_t>(-1)));
    for (size_t si = 0; si < seqs.size(); ++si) {
        for (size_t c = 0; c < Lraw; ++c) {
            alnCode[si][c] = static_cast<int8_t>(baseToIdx(seqs[si].aln[c]));
        }
    }
    // Infernal-like fragment inference (--fragthresh default 0.5):
    // terminal gaps of fragment sequences are treated as missing data.
    std::vector<unsigned char> isFrag(seqs.size(), 0);
    std::vector<int> fragFirst(seqs.size(), -1);
    std::vector<int> fragLast(seqs.size(), -1);
    for (size_t si = 0; si < seqs.size(); ++si) {
        int first = -1, last = -1;
        for (size_t c = 0; c < Lraw; ++c) {
            if (!isGap(seqs[si].aln[c])) {
                if (first < 0) first = static_cast<int>(c);
                last = static_cast<int>(c);
            }
        }
        fragFirst[si] = first;
        fragLast[si] = last;
        if (first >= 0 && last >= first) {
            const int span = last - first + 1;
            if (static_cast<double>(span) <= 0.5 * static_cast<double>(Lraw)) {
                isFrag[si] = 1;
            }
        }
    }
    const auto isMissingTerminal = [&](size_t si, int rawPos) -> bool {
        if (isFrag[si] == 0) return false;
        const int f = fragFirst[si];
        const int l = fragLast[si];
        if (f < 0 || l < 0) return true;
        return (rawPos < f || rawPos > l);
    };
    std::vector<int> keepRaw;
    keepRaw.reserve(Lraw);
    if (!rfCons.empty() && rfCons.size() == Lraw) {
        for (size_t c = 0; c < Lraw; ++c) {
            const char rc = rfCons[c];
            if (rc != '.' && rc != '-' && rc != '_' && rc != '~') {
                keepRaw.push_back(static_cast<int>(c));
            }
        }
    } else {
        // Infernal-like column selection: terminal gaps of inferred fragments
        // are missing data and do not count against nongap fraction.
        for (size_t c = 0; c < Lraw; ++c) {
            size_t valid = 0;
            size_t nongap = 0;
            for (size_t s = 0; s < seqs.size(); ++s) {
                if (isMissingTerminal(s, static_cast<int>(c))) {
                    continue;
                }
                ++valid;
                if (!isGap(seqs[s].aln[c])) {
                    ++nongap;
                }
            }
            if (valid == 0) {
                continue;
            }
            const size_t minNongap = (valid + 1) / 2;
            if (nongap >= minNongap) {
                keepRaw.push_back(static_cast<int>(c));
            }
        }
    }
    if (keepRaw.empty()) {
        Debug(Debug::ERROR) << "cmbuild: all alignment columns are gaps\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<int> rawPairs;
    bool haveStructure = false;
    if (!ssCons.empty()) {
        if (ssCons.size() == Lraw) {
            rawPairs = parseSsPairs(ssCons);
            haveStructure = !rawPairs.empty();
            if (!haveStructure) {
                Debug(Debug::WARNING) << "cmbuild: SS_cons parse failed, trying LinearFold consensus prediction\n";
            }
        } else {
            Debug(Debug::WARNING) << "cmbuild: SS_cons length mismatch, trying LinearFold consensus prediction\n";
        }
    }
    if (!haveStructure) {
        std::string predSsCons;
        if (predictRawPairsLinearFold(seqs, keepRaw, Lraw, rawPairs, predSsCons)) {
            ssCons = predSsCons;
            haveStructure = true;
            Debug(Debug::INFO) << "cmbuild: predicted SS_cons using native LinearFold on consensus sequence\n";
        } else {
            rawPairs.clear();
            Debug(Debug::WARNING) << "cmbuild: no usable structure constraints (SS_cons missing/invalid and LinearFold prediction failed)\n";
        }
    }
#ifdef HAVE_INFERNAL_BRIDGE
    if (InfernalBridge::isConfigured()) {
        std::string cmText;
        std::string infernalErr;
        const std::string stoText = buildStockholmText("mmseqs_model", seqs, ssCons);
        if (InfernalBridge::buildCmFromStockholmText(stoText, cmText, infernalErr)) {
            DBWriter writer(par.db2.c_str(),
                            par.db2Index.c_str(),
                            1,
                            par.compressed,
                            Parameters::DBTYPE_GENERIC_DB);
            writer.open();
            writer.writeData(cmText.c_str(), cmText.size(), 0, 0);
            writer.close();
            Debug(Debug::INFO) << "cmbuild: wrote Infernal-backed CM through CMake bridge to " << par.db2 << "\n";
            return EXIT_SUCCESS;
        }
        Debug(Debug::WARNING) << "cmbuild: Infernal bridge failed (" << infernalErr
                              << "), falling back to native MMseqs builder\n";
    }
#endif
    struct BuiltModel {
        int N;
        std::vector<int> keep;
        std::vector<int> pairMap; // 1-based
        std::vector< std::array<double, 4> > singleCnt; // 1-based
        std::map< std::pair<int, int>, std::array<double, 16> > pairCnt;
    };
    const auto buildModel = [&](const std::vector<int> &keepCols) -> BuiltModel {
        BuiltModel bm;
        bm.keep = keepCols;
        bm.N = static_cast<int>(keepCols.size());
        bm.pairMap.assign(static_cast<size_t>(bm.N) + 1, -1);
        std::vector<int> rawToNew(Lraw, -1);
        for (size_t i = 0; i < keepCols.size(); ++i) {
            rawToNew[static_cast<size_t>(keepCols[i])] = static_cast<int>(i) + 1;
        }
        if (!rawPairs.empty()) {
            for (size_t oi = 0; oi < Lraw; ++oi) {
                const int oj = rawPairs[oi];
                if (oj < 0) {
                    continue;
                }
                if (oi < static_cast<size_t>(oj)) {
                    const int ni = rawToNew[oi];
                    const int nj = rawToNew[static_cast<size_t>(oj)];
                    if (ni > 0 && nj > 0) {
                        bm.pairMap[static_cast<size_t>(ni)] = nj;
                        bm.pairMap[static_cast<size_t>(nj)] = ni;
                    }
                }
            }
        }
        bm.singleCnt.assign(static_cast<size_t>(bm.N) + 1, std::array<double, 4>());
        for (size_t si = 0; si < seqs.size(); ++si) {
            const double w = seqWeight[si];
            for (int ni = 1; ni <= bm.N; ++ni) {
                const int raw = bm.keep[static_cast<size_t>(ni - 1)];
                const int bi = alnCode[si][static_cast<size_t>(raw)];
                if (bi >= 0) {
                    bm.singleCnt[static_cast<size_t>(ni)][static_cast<size_t>(bi)] += w;
                }
            }
            for (int ni = 1; ni <= bm.N; ++ni) {
                const int nj = bm.pairMap[static_cast<size_t>(ni)];
                if (nj <= ni) {
                    continue;
                }
                const int rawI = bm.keep[static_cast<size_t>(ni - 1)];
                const int rawJ = bm.keep[static_cast<size_t>(nj - 1)];
                const int li = alnCode[si][static_cast<size_t>(rawI)];
                const int ri = alnCode[si][static_cast<size_t>(rawJ)];
                if (li >= 0 && ri >= 0) {
                    bm.pairCnt[std::make_pair(ni, nj)][static_cast<size_t>(li * 4 + ri)] += w;
                }
            }
        }
        return bm;
    };

    BuiltModel mainModel = buildModel(keepRaw);
    std::vector<int> rawToMain(Lraw, -1);
    for (size_t i = 0; i < keepRaw.size(); ++i) {
        rawToMain[static_cast<size_t>(keepRaw[i])] = static_cast<int>(i) + 1;
    }

    const auto estimateMid = [&](const BuiltModel &bm) -> MidEstimates {
        MidEstimates mid;
        mid.matchOcc.assign(static_cast<size_t>(bm.N) + 1, 0.5);
        mid.delOcc.assign(static_cast<size_t>(bm.N) + 1, 0.1);
        mid.tMM.assign(static_cast<size_t>(bm.N) + 1, 0.7);
        mid.tMD.assign(static_cast<size_t>(bm.N) + 1, 0.3);
        mid.tDM.assign(static_cast<size_t>(bm.N) + 1, 0.7);
        mid.tDD.assign(static_cast<size_t>(bm.N) + 1, 0.3);

        for (int p = 1; p <= bm.N; ++p) {
            const std::array<double, 4> &cnt = bm.singleCnt[static_cast<size_t>(p)];
            const double obs = cnt[0] + cnt[1] + cnt[2] + cnt[3];
            const double occ = std::max(0.01, std::min(0.99, (obs + 1.0) / (nEffSeq + 2.0)));
            const double del = std::max(0.01, std::min(0.95, 1.0 - occ));
            mid.matchOcc[static_cast<size_t>(p)] = occ;
            mid.delOcc[static_cast<size_t>(p)] = del;
        }

        for (int p = 1; p < bm.N; ++p) {
            double cMM = 1.0, cMD = 1.0, cDM = 1.0, cDD = 1.0;
            const int rawP = bm.keep[static_cast<size_t>(p - 1)];
            const int rawQ = bm.keep[static_cast<size_t>(p)];
            for (size_t si = 0; si < seqs.size(); ++si) {
                const double w = seqWeight[si];
                if (isMissingTerminal(si, rawP) || isMissingTerminal(si, rawQ)) {
                    continue;
                }
                const bool curM = (alnCode[si][static_cast<size_t>(rawP)] >= 0);
                const bool nxtM = (alnCode[si][static_cast<size_t>(rawQ)] >= 0);
                if (curM && nxtM) cMM += w;
                else if (curM && !nxtM) cMD += w;
                else if (!curM && nxtM) cDM += w;
                else cDD += w;
            }
            const double rowM = cMM + cMD;
            const double rowD = cDM + cDD;
            mid.tMM[static_cast<size_t>(p)] = std::max(0.01, std::min(0.99, cMM / std::max(rowM, 1e-9)));
            mid.tMD[static_cast<size_t>(p)] = std::max(0.01, std::min(0.99, cMD / std::max(rowM, 1e-9)));
            mid.tDM[static_cast<size_t>(p)] = std::max(0.01, std::min(0.99, cDM / std::max(rowD, 1e-9)));
            mid.tDD[static_cast<size_t>(p)] = std::max(0.01, std::min(0.99, cDD / std::max(rowD, 1e-9)));
        }

        for (std::map< std::pair<int, int>, std::array<double, 16> >::const_iterator it = bm.pairCnt.begin();
             it != bm.pairCnt.end(); ++it) {
            double obs = 0.0;
            for (int k = 0; k < 16; ++k) {
                obs += it->second[static_cast<size_t>(k)];
            }
            const double occ = std::max(0.01, std::min(0.99, (obs + 1.0) / (nEffSeq + 2.0)));
            const double del = std::max(0.01, std::min(0.95, 1.0 - occ));
            mid.pairMatchOcc[it->first] = occ;
            mid.pairDelOcc[it->first] = del;
        }
        return mid;
    };

    MidEstimates midMain = estimateMid(mainModel);

    // Build additional one-column-deletion variants for top gappy unpaired columns.
    // This increases indel flexibility and improves score calibration on shifted homologs.
    std::vector<BuiltModel> delModels;
    std::vector<std::string> delPrefixes;
    std::vector<int> insertionSitesMain; // 1-based consensus positions in main model
    if (!rawPairs.empty() && keepRaw.size() > 1) {
        struct DropCand {
            int raw;
            double gapFrac;
        };
        std::vector<DropCand> cands;
        for (size_t k = 0; k < keepRaw.size(); ++k) {
            const int raw = keepRaw[k];
            if (rawPairs[static_cast<size_t>(raw)] >= 0) {
                continue;
            }
            size_t gaps = 0;
            for (size_t s = 0; s < seqs.size(); ++s) {
                if (isGap(seqs[s].aln[static_cast<size_t>(raw)])) {
                    ++gaps;
                }
            }
            const double gapFrac = static_cast<double>(gaps) / static_cast<double>(seqs.size());
            if (gapFrac >= 0.10) {
                DropCand dc;
                dc.raw = raw;
                dc.gapFrac = gapFrac;
                cands.push_back(dc);
            }
        }
        std::sort(cands.begin(), cands.end(), [](const DropCand &a, const DropCand &b) {
            if (a.gapFrac != b.gapFrac) {
                return a.gapFrac > b.gapFrac;
            }
            return a.raw < b.raw;
        });
        const size_t maxInsertionSites = std::min<size_t>(12, cands.size());
        for (size_t si = 0; si < maxInsertionSites; ++si) {
            const int pos = rawToMain[static_cast<size_t>(cands[si].raw)];
            if (pos > 0 && pos < mainModel.N) {
                insertionSitesMain.push_back(pos);
            }
        }
        std::sort(insertionSitesMain.begin(), insertionSitesMain.end());
        insertionSitesMain.erase(std::unique(insertionSitesMain.begin(), insertionSitesMain.end()), insertionSitesMain.end());

        const size_t maxDelVariants = std::min<size_t>(3, cands.size());
        for (size_t vi = 0; vi < maxDelVariants; ++vi) {
            const int dropRaw = cands[vi].raw;
            std::vector<int> keepAlt;
            keepAlt.reserve(keepRaw.size() - 1);
            for (size_t i = 0; i < keepRaw.size(); ++i) {
                if (keepRaw[i] != dropRaw) {
                    keepAlt.push_back(keepRaw[i]);
                }
            }
            if (keepAlt.empty()) {
                continue;
            }
            BuiltModel dm = buildModel(keepAlt);
            delModels.push_back(dm);
            char pref[16];
            std::snprintf(pref, sizeof(pref), "D%zu", vi + 1);
            delPrefixes.push_back(std::string(pref));
            Debug(Debug::INFO) << "cmbuild: added deletion variant " << delPrefixes.back()
                               << " removing unpaired column raw=" << (dropRaw + 1)
                               << " gapFraction=" << cands[vi].gapFrac
                               << " altModelLen=" << dm.N << "\n";
        }
    }

    std::ostringstream out;

    const double pseudo = 0.5;
    double bgCnt[4] = {1.0, 1.0, 1.0, 1.0};
    for (size_t si = 0; si < seqs.size(); ++si) {
        const double w = seqWeight[si];
        for (size_t k = 0; k < keepRaw.size(); ++k) {
            const int bi = alnCode[si][static_cast<size_t>(keepRaw[k])];
            if (bi >= 0) {
                bgCnt[bi] += w;
            }
        }
    }
    const double bgSum = bgCnt[0] + bgCnt[1] + bgCnt[2] + bgCnt[3];
    double bg[4] = {
        bgCnt[0] / bgSum,
        bgCnt[1] / bgSum,
        bgCnt[2] / bgSum,
        bgCnt[3] / bgSum
    };
    const auto scoreKeptSequence = [&](const std::vector<int> &obs) -> double {
        double s = 0.0;
        for (int p = 1; p <= mainModel.N; ++p) {
            const int q = mainModel.pairMap[static_cast<size_t>(p)];
            if (q > p) {
                const int li = obs[static_cast<size_t>(p - 1)];
                const int ri = obs[static_cast<size_t>(q - 1)];
                if (li >= 0 && ri >= 0) {
                    std::map< std::pair<int, int>, std::array<double, 16> >::const_iterator it = mainModel.pairCnt.find(std::make_pair(p, q));
                    std::array<double, 16> empty = {};
                    const std::array<double, 16> &cnt = (it == mainModel.pairCnt.end()) ? empty : it->second;
                    double den = 0.0;
                    for (int k = 0; k < 16; ++k) {
                        den += cnt[static_cast<size_t>(k)] + pseudo;
                    }
                    const double pij = (cnt[static_cast<size_t>(li * 4 + ri)] + pseudo) / den;
                    const double qij = std::max(bg[li] * bg[ri], 1e-12);
                    s += std::log2(pij / qij);
                } else {
                    s += -2.0;
                }
            } else if (q < p) {
                continue;
            } else {
                const int bi = obs[static_cast<size_t>(p - 1)];
                if (bi >= 0) {
                    const std::array<double, 4> &cnt = mainModel.singleCnt[static_cast<size_t>(p)];
                    double den = 0.0;
                    for (int b = 0; b < 4; ++b) {
                        den += cnt[static_cast<size_t>(b)] + pseudo;
                    }
                    const double pi = (cnt[static_cast<size_t>(bi)] + pseudo) / den;
                    const double qi = std::max(bg[bi], 1e-12);
                    s += std::log2(pi / qi);
                } else {
                    s += -1.0;
                }
            }
        }
        return s;
    };

    std::vector<double> trainScores;
    trainScores.reserve(seqs.size());
    std::vector<int> obs(static_cast<size_t>(mainModel.N), -1);
    for (size_t si = 0; si < seqs.size(); ++si) {
        for (int p = 1; p <= mainModel.N; ++p) {
            const int raw = mainModel.keep[static_cast<size_t>(p - 1)];
            obs[static_cast<size_t>(p - 1)] = alnCode[si][static_cast<size_t>(raw)];
        }
        trainScores.push_back(scoreKeptSequence(obs));
    }
    std::vector<double> posScores;
    std::vector<double> negScores;
    posScores.reserve(trainScores.size());
    negScores.reserve(trainScores.size());
    for (size_t i = 0; i < trainScores.size(); ++i) {
        if (trainScores[i] >= 0.0) {
            posScores.push_back(trainScores[i]);
        } else {
            negScores.push_back(trainScores[i]);
        }
    }
    if (negScores.empty() && !trainScores.empty()) {
        uint32_t rng = 0x12345678u;
        const int decoys = std::min(static_cast<int>(trainScores.size()), 512);
        negScores.reserve(static_cast<size_t>(decoys));
        for (int d = 0; d < decoys; ++d) {
            for (int p = 1; p <= mainModel.N; ++p) {
                rng = 1664525u * rng + 1013904223u;
                const uint32_t r = rng >> 16;
                double c = (static_cast<double>(r % 100000u) + 0.5) / 100000.0;
                int bi = 0;
                if (c < bg[0]) {
                    bi = 0;
                } else if (c < bg[0] + bg[1]) {
                    bi = 1;
                } else if (c < bg[0] + bg[1] + bg[2]) {
                    bi = 2;
                } else {
                    bi = 3;
                }
                obs[static_cast<size_t>(p - 1)] = bi;
            }
            negScores.push_back(scoreKeptSequence(obs));
        }
        Debug(Debug::INFO) << "cmbuild: synthesized " << negScores.size()
                           << " decoy negatives for CAL_NEG fit\n";
    }
    // Emit scores in bit units, then apply a conservative global slope so
    // native scores are on an Infernal-like bit scale (zero intercept).
    // ---- Emit Infernal CM format (MATL-only topology) ----
    // This routes through runInfernalExactScan in cmsearch, giving consistent
    // scoring with infernal-built --noss CMs.
    const int N = mainModel.N;
    // State layout:
    //   ROOT:   S(0), IL(1), IR(2)
    //   MATL_p: ML(3+3*(p-1)), D(4+3*(p-1)), IL(5+3*(p-1))  for p=1..N
    //   END:    E(3+3*N)
    const int totalStates = 3 * N + 4;
    // Infernal's configured W is typically closer to ~2*CLEN for structured
    // models; too-small W truncates alignments and depresses scores.
    const int W = std::max(N + N / 2 + 10, 2 * N + 13);
    const int sEnd = 3 + 3 * N;

    // For paired positions, marginalize pairCnt to get single-base emit counts.
    std::vector< std::array<double, 4> > emitCnt(static_cast<size_t>(N) + 1);
    for (int p = 1; p <= N; ++p) {
        emitCnt[static_cast<size_t>(p)] = mainModel.singleCnt[static_cast<size_t>(p)];
        const int q = mainModel.pairMap[static_cast<size_t>(p)];
        if (q > p) {
            // Left side of pair: marginalize over right base
            std::map< std::pair<int,int>, std::array<double,16> >::const_iterator it =
                mainModel.pairCnt.find(std::make_pair(p, q));
            if (it != mainModel.pairCnt.end()) {
                for (int b = 0; b < 4; ++b) {
                    double m = 0.0;
                    for (int r = 0; r < 4; ++r) m += it->second[static_cast<size_t>(b * 4 + r)];
                    emitCnt[static_cast<size_t>(p)][static_cast<size_t>(b)] = m;
                }
            }
        } else if (q >= 0 && q < p) {
            // Right side of pair: marginalize over left base
            std::map< std::pair<int,int>, std::array<double,16> >::const_iterator it =
                mainModel.pairCnt.find(std::make_pair(q, p));
            if (it != mainModel.pairCnt.end()) {
                for (int b = 0; b < 4; ++b) {
                    double m = 0.0;
                    for (int l = 0; l < 4; ++l) m += it->second[static_cast<size_t>(l * 4 + b)];
                    emitCnt[static_cast<size_t>(p)][static_cast<size_t>(b)] = m;
                }
            }
        }
    }
    const auto posteriorMeanMixDirichlet4 = [](const std::array<double, 4> &cnt) -> std::array<double, 4> {
        static const double q[10] = {
            0.081706, 0.104534, 0.048944, 0.064111, 0.085266,
            0.045348, 0.100949, 0.108835, 0.234814, 0.125493
        };
        static const double alpha[10][4] = {
            {0.963855, 3.273863, 0.444739, 1.958731},
            {0.589011, 0.648423, 0.360672, 5.771004},
            {2.609834, 0.127100, 1.180559, 0.134264},
            {1.259286, 0.659029, 4.874613, 0.882126},
            {4.664219, 0.628128, 0.448894, 0.661556},
            {0.250974, 9.700414, 0.206184, 0.338607},
            {0.178455, 0.049385, 7.914643, 0.100802},
            {23.818220, 0.064454, 0.119891, 0.101866},
            {2.980233, 1.817786, 1.818483, 3.042635},
            {0.024428, 0.064315, 0.008054, 0.107062}
        };
        double logw[10];
        double maxLogw = -1e300;
        double csum = 0.0;
        for (int i = 0; i < 4; ++i) {
            csum += cnt[static_cast<size_t>(i)];
        }
        for (int m = 0; m < 10; ++m) {
            double asum = 0.0;
            for (int i = 0; i < 4; ++i) {
                asum += alpha[m][i];
            }
            double lw = std::log(std::max(q[m], 1e-300)) + std::lgamma(asum) - std::lgamma(asum + csum);
            for (int i = 0; i < 4; ++i) {
                lw += std::lgamma(alpha[m][i] + cnt[static_cast<size_t>(i)]) - std::lgamma(alpha[m][i]);
            }
            logw[m] = lw;
            if (lw > maxLogw) {
                maxLogw = lw;
            }
        }
        double wsum = 0.0;
        double wnorm[10];
        for (int m = 0; m < 10; ++m) {
            wnorm[m] = std::exp(logw[m] - maxLogw);
            wsum += wnorm[m];
        }
        wsum = std::max(wsum, 1e-300);
        std::array<double, 4> p = {};
        for (int m = 0; m < 10; ++m) {
            const double wm = wnorm[m] / wsum;
            double asum = 0.0;
            for (int i = 0; i < 4; ++i) {
                asum += alpha[m][i];
            }
            const double den = std::max(csum + asum, 1e-12);
            for (int i = 0; i < 4; ++i) {
                p[static_cast<size_t>(i)] += wm * ((cnt[static_cast<size_t>(i)] + alpha[m][i]) / den);
            }
        }
        return p;
    };
    const auto posteriorMeanMixDirichlet16 = [](const std::array<double, 16> &cnt) -> std::array<double, 16> {
        static const double q[10] = {
            0.016584, 0.000948, 0.185395, 0.082929, 0.039651,
            0.141227, 0.132571, 0.249417, 0.140727, 0.010551
        };
        static const double alpha[10][16] = {
            {0.142252, 0.180113, 0.153776, 0.222524, 0.539721, 0.170380, 0.230123, 0.190004, 0.583682, 0.222992, 0.134277, 0.187596, 12172.002267, 0.546735, 14.841962, 17.271555},
            {2.547410, 14.293143, 0.015263, 7.761130, 0.029915, 3.493007, 14.049507, 1.480341, 3.754643, 7.140983, 4.733217, 44.190624, 2.758417, 1.687945, 2.421882, 2.724272},
            {0.054512, 0.067070, 0.054506, 1.210822, 0.119647, 0.030366, 3.188992, 0.076098, 0.060153, 1.426299, 0.042134, 0.362385, 2.308941, 0.048026, 0.695604, 0.146166},
            {0.481661, 0.414811, 0.419836, 3.024237, 0.421853, 0.232594, 3.637964, 0.328914, 0.400575, 2.647559, 0.269173, 1.022533, 3.376215, 0.380735, 1.397263, 0.695235},
            {0.145102, 0.122876, 3.107999, 0.099093, 10.564395, 4.450523, 15500.159054, 0.049506, 0.032312, 0.368096, 4.375012, 0.111267, 0.904329, 0.074819, 0.376995, 0.093335},
            {0.016163, 0.040913, 0.014116, 0.527169, 0.003200, 0.001437, 0.074551, 0.013706, 0.019149, 0.413953, 0.012037, 0.268400, 0.078554, 0.005712, 0.020960, 0.037299},
            {0.004230, 0.045568, 0.000699, 0.258190, 0.001391, 0.020574, 0.073793, 0.001111, 0.017598, 7.014687, 0.015465, 0.189706, 0.057044, 0.014820, 0.013400, 0.001169},
            {0.008027, 0.006602, 0.012884, 0.089652, 0.040423, 0.011659, 0.789524, 0.021433, 0.011990, 0.091424, 0.011034, 0.019953, 0.389278, 0.009006, 0.198688, 0.027512},
            {0.068663, 0.176455, 0.077881, 2.165192, 0.035566, 0.051544, 1.087382, 0.048265, 0.057469, 5.631915, 0.048459, 0.906756, 0.904423, 0.086167, 0.270333, 0.159528},
            {0.478576, 0.402540, 18.466281, 16947.982248, 0.389092, 0.386664, 0.619656, 20.908826, 0.375696, 4.605442, 0.396373, 13.623423, 0.513956, 0.363145, 0.606193, 18.301915}
        };
        double logw[10];
        double maxLogw = -1e300;
        double csum = 0.0;
        for (int i = 0; i < 16; ++i) csum += cnt[static_cast<size_t>(i)];
        for (int m = 0; m < 10; ++m) {
            double asum = 0.0;
            for (int i = 0; i < 16; ++i) asum += alpha[m][i];
            double lw = std::log(std::max(q[m], 1e-300)) + std::lgamma(asum) - std::lgamma(asum + csum);
            for (int i = 0; i < 16; ++i) {
                lw += std::lgamma(alpha[m][i] + cnt[static_cast<size_t>(i)]) - std::lgamma(alpha[m][i]);
            }
            logw[m] = lw;
            if (lw > maxLogw) maxLogw = lw;
        }
        double wsum = 0.0;
        double wnorm[10];
        for (int m = 0; m < 10; ++m) {
            wnorm[m] = std::exp(logw[m] - maxLogw);
            wsum += wnorm[m];
        }
        wsum = std::max(wsum, 1e-300);
        std::array<double, 16> p = {};
        for (int m = 0; m < 10; ++m) {
            const double wm = wnorm[m] / wsum;
            double asum = 0.0;
            for (int i = 0; i < 16; ++i) asum += alpha[m][i];
            const double den = std::max(csum + asum, 1e-12);
            for (int i = 0; i < 16; ++i) {
                p[static_cast<size_t>(i)] += wm * ((cnt[static_cast<size_t>(i)] + alpha[m][i]) / den);
            }
        }
        return p;
    };
    const int nbps = [&]() {
        int c = 0;
        for (int p = 1; p <= N; ++p) if (mainModel.pairMap[static_cast<size_t>(p)] > p) c++;
        return c;
    }();
    const auto setTargetRelEnt = [&](int clen, int numBp) -> double {
        const double esigma = 45.0;
        const double reTarget = (numBp > 0) ? 0.59 : 0.38; // DEFAULT_ETARGET / DEFAULT_ETARGET_HMMFILTER
        const double et = (esigma - std::log2(2.0 / (static_cast<double>(clen) * static_cast<double>(clen + 1)))) / static_cast<double>(clen);
        return std::max(et, reTarget);
    };
    const auto meanMatchRelEntAtScale = [&](double scale) -> double {
        double re = 0.0;
        for (int p = 1; p <= N; ++p) {
            std::array<double, 4> sc = {};
            for (int b = 0; b < 4; ++b) sc[static_cast<size_t>(b)] = scale * emitCnt[static_cast<size_t>(p)][static_cast<size_t>(b)];
            const std::array<double, 4> post = posteriorMeanMixDirichlet4(sc);
            double r = 0.0;
            for (int b = 0; b < 4; ++b) {
                const double pb = std::max(post[static_cast<size_t>(b)], 1e-12);
                const double qb = std::max(bg[static_cast<size_t>(b)], 1e-12);
                r += pb * std::log2(pb / qb);
            }
            re += r;
        }
        return (N > 0) ? (re / static_cast<double>(N)) : 0.0;
    };
    // Infernal --eent style effective sequence number: rescale weighted counts so
    // posterior mean match relative entropy reaches a length-dependent target.
    const double neffMin = 0.1;
    const double neffMax = std::max(static_cast<double>(seqs.size()), neffMin);
    const double targetRe = setTargetRelEnt(N, nbps);
    double lo = neffMin / static_cast<double>(seqs.size());
    double hi = neffMax / static_cast<double>(seqs.size());
    auto fx = [&](double scale) { return meanMatchRelEntAtScale(scale) - targetRe; };
    double fLo = fx(lo);
    double fHi = fx(hi);
    double effScale = hi;
    if (fLo >= 0.0) {
        effScale = lo;
    } else if (fHi <= 0.0) {
        effScale = hi;
    } else {
        for (int it = 0; it < 40; ++it) {
            const double mid = 0.5 * (lo + hi);
            const double fMid = fx(mid);
            if (fMid > 0.0) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        effScale = 0.5 * (lo + hi);
    }
    // Infernal default prior: MATL singlet emission prior == insert emission prior.
    std::vector< std::array<double, 4> > emitProb(static_cast<size_t>(N) + 1);
    for (int p = 1; p <= N; ++p) {
        std::array<double, 4> sc = {};
        for (int b = 0; b < 4; ++b) sc[static_cast<size_t>(b)] = effScale * emitCnt[static_cast<size_t>(p)][static_cast<size_t>(b)];
        emitProb[static_cast<size_t>(p)] = posteriorMeanMixDirichlet4(sc);
    }
    std::vector<double> pairMiHalfBonus(static_cast<size_t>(N) + 1, 0.0);
    if (nEffSeq >= 8.0) {
        for (std::map< std::pair<int, int>, std::array<double, 16> >::const_iterator it = mainModel.pairCnt.begin();
             it != mainModel.pairCnt.end(); ++it) {
            const int p = it->first.first;
            const int q = it->first.second;
            if (p < 1 || p > N || q < 1 || q > N) continue;
            std::array<double, 16> sc = {};
            for (int k = 0; k < 16; ++k) sc[static_cast<size_t>(k)] = effScale * it->second[static_cast<size_t>(k)];
            const std::array<double, 16> post = posteriorMeanMixDirichlet16(sc);
            double pi[4] = {0.0, 0.0, 0.0, 0.0};
            double pj[4] = {0.0, 0.0, 0.0, 0.0};
            for (int a = 0; a < 4; ++a) {
                for (int b = 0; b < 4; ++b) {
                    const double pij = post[static_cast<size_t>(a * 4 + b)];
                    pi[a] += pij;
                    pj[b] += pij;
                }
            }
            double mi = 0.0;
            for (int a = 0; a < 4; ++a) {
                for (int b = 0; b < 4; ++b) {
                    const double pij = std::max(post[static_cast<size_t>(a * 4 + b)], 1e-12);
                    const double qij = std::max(pi[a] * pj[b], 1e-12);
                    mi += pij * std::log2(pij / qij);
                }
            }
            const double bonus = 2.0 * mi;
            pairMiHalfBonus[static_cast<size_t>(p)] = std::max(pairMiHalfBonus[static_cast<size_t>(p)], bonus);
            pairMiHalfBonus[static_cast<size_t>(q)] = std::max(pairMiHalfBonus[static_cast<size_t>(q)], bonus);
        }
    }

    // Compute log2-odds emit score for base b at position p using Infernal-like
    // mixture-Dirichlet posterior mean parameters.
    const auto emitScore = [&](int p, int b) -> double {
        const double prob = std::max(emitProb[static_cast<size_t>(p)][static_cast<size_t>(b)], 1e-12);
        const double bgp = std::max(bg[static_cast<size_t>(b)], 1e-12);
        return std::log2(prob / bgp) + pairMiHalfBonus[static_cast<size_t>(p)];
    };

    // Format a log2 score; "*" for -inf.
    const auto fsc = [](double v) -> std::string {
        if (!std::isfinite(v) || v < -999.0) { return "*"; }
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.3f", v);
        return std::string(buf);
    };

    // Infernal default prior alphas for MATL-only transitions (next node type MATL/END).
    const std::array<double, 4> ALPHA_ROOT_S_MATL  = {0.028518011579, 0.024705844026, 1.464047470747, 0.074164509948};
    const std::array<double, 4> ALPHA_ROOT_IL_MATL = {0.250101882938, 0.155728904821, 0.370945030932, 0.027811408475};
    const std::array<double, 3> ALPHA_ROOT_IR_MATL = {0.601223387577, 0.939499051719, 0.092516097691};
    const std::array<double, 3> ALPHA_MATL_ML_MATL = {0.015185708311, 1.809432933023, 0.038601480352};
    const std::array<double, 3> ALPHA_MATL_D_MATL  = {0.005679808868, 0.127365862719, 0.277086556814};
    const std::array<double, 3> ALPHA_MATL_IL_MATL = {0.601223387577, 0.939499051719, 0.092516097691};
    const std::array<double, 2> ALPHA_MATL_ML_END  = {0.009635966745, 1.220143960207};
    const std::array<double, 2> ALPHA_MATL_D_END   = {0.019509171372, 6.781321301695};
    const std::array<double, 2> ALPHA_MATL_IL_END  = {0.264643213319, 0.671462565227};

    // State line helper (10 fixed fields + transitions + emissions).
    // Format: TYPE  v  p1  p2  cfirst  cnum  x  dmin  dmax  y
    const auto stateLine = [&](const char *type, int v, int p1, int p2,
                                int cfirst, int cnum, int dmin, int dmax) -> std::string {
        char buf[128];
        std::snprintf(buf, sizeof(buf),
                      "%6s %5d %5d %d %5d %5d %5d %5d %5d %5d",
                      type, v, p1, p2, cfirst, cnum, 0, dmin, dmax, W + 1);
        return std::string(buf);
    };

    // Write Infernal CM header
    out << "INFERNAL1/a [MMseqs2-cmbuild]\n";
    out << "NAME     mmseqs_model\n";
    out << "STATES   " << totalStates << "\n";
    out << "CLEN     " << N << "\n";
    out << "W        " << W << "\n";
    out << "ALPH     RNA\n";
    out << "NULL     0.000  0.000  0.000  0.000\n";
    out << "CM\n";

    // ROOT node
    out << "                                             [ ROOT    0 ]      -      -\n";
    const auto logPost4 = [](const std::array<double, 4> &cnt, const std::array<double, 4> &alpha) -> std::array<double, 4> {
        std::array<double, 4> lp = {};
        double den = 0.0;
        for (int i = 0; i < 4; ++i) den += cnt[static_cast<size_t>(i)] + alpha[static_cast<size_t>(i)];
        den = std::max(den, 1e-12);
        for (int i = 0; i < 4; ++i) {
            const double p = (cnt[static_cast<size_t>(i)] + alpha[static_cast<size_t>(i)]) / den;
            lp[static_cast<size_t>(i)] = std::log2(std::max(p, 1e-12));
        }
        return lp;
    };
    const auto logPost3 = [](const std::array<double, 3> &cnt, const std::array<double, 3> &alpha) -> std::array<double, 3> {
        std::array<double, 3> lp = {};
        double den = 0.0;
        for (int i = 0; i < 3; ++i) den += cnt[static_cast<size_t>(i)] + alpha[static_cast<size_t>(i)];
        den = std::max(den, 1e-12);
        for (int i = 0; i < 3; ++i) {
            const double p = (cnt[static_cast<size_t>(i)] + alpha[static_cast<size_t>(i)]) / den;
            lp[static_cast<size_t>(i)] = std::log2(std::max(p, 1e-12));
        }
        return lp;
    };
    const auto logPost2 = [](const std::array<double, 2> &cnt, const std::array<double, 2> &alpha) -> std::array<double, 2> {
        std::array<double, 2> lp = {};
        double den = 0.0;
        for (int i = 0; i < 2; ++i) den += cnt[static_cast<size_t>(i)] + alpha[static_cast<size_t>(i)];
        den = std::max(den, 1e-12);
        for (int i = 0; i < 2; ++i) {
            const double p = (cnt[static_cast<size_t>(i)] + alpha[static_cast<size_t>(i)]) / den;
            lp[static_cast<size_t>(i)] = std::log2(std::max(p, 1e-12));
        }
        return lp;
    };
    std::array<double, 4> rootSCnt = {0.0, 0.0, 0.0, 0.0};   // S: IL, IR, ML, D
    std::array<double, 4> rootILCnt = {0.0, 0.0, 0.0, 0.0};  // IL: IL, IR, ML, D
    std::array<double, 3> rootIRCnt = {0.0, 0.0, 0.0};       // IR: IR, ML, D
    if (N > 0) {
        const int raw0 = mainModel.keep[0];
        for (size_t si = 0; si < seqs.size(); ++si) {
            const double w = effScale * seqWeight[si];
            if (isMissingTerminal(si, raw0)) continue;
            const bool m0 = (alnCode[si][static_cast<size_t>(raw0)] >= 0);
            if (m0) rootSCnt[2] += w;
            else    rootSCnt[3] += w;
        }
    }
    const std::array<double, 4> rootSLog = logPost4(rootSCnt, ALPHA_ROOT_S_MATL);
    const std::array<double, 4> rootILLog = logPost4(rootILCnt, ALPHA_ROOT_IL_MATL);
    const std::array<double, 3> rootIRLog = logPost3(rootIRCnt, ALPHA_ROOT_IR_MATL);
    // S(0): children IL(1), IR(2), ML(3), D(4)
    out << stateLine("S", 0, -1, 0, 1, 4, 1, W)
        << "  " << fsc(rootSLog[0]) << "  " << fsc(rootSLog[1])
        << "  " << fsc(rootSLog[2]) << "  " << fsc(rootSLog[3]) << "\n";
    // IL(1): self-loop + IR + ML + D; emit uniform
    out << stateLine("IL", 1, 1, 2, 1, 4, N, W)
        << "  " << fsc(rootILLog[0]) << "  " << fsc(rootILLog[1])
        << "  " << fsc(rootILLog[2]) << "  " << fsc(rootILLog[3])
        << "  0.000  0.000  0.000  0.000\n";
    // IR(2): self-loop + ML + D; emit uniform
    out << stateLine("IR", 2, 2, 3, 2, 3, N, W)
        << "  " << fsc(rootIRLog[0]) << "  " << fsc(rootIRLog[1]) << "  " << fsc(rootIRLog[2])
        << "  0.000  0.000  0.000  0.000\n";

    // MATL nodes p=1..N
    for (int p = 1; p <= N; ++p) {
        const int sML = 3 + 3 * (p - 1);
        const int sD  = sML + 1;
        const int sIL = sML + 2;
        const bool isLast = (p == N);

        char emitBuf[80];
        std::snprintf(emitBuf, sizeof(emitBuf), "  %.3f  %.3f  %.3f  %.3f",
                      emitScore(p, 0), emitScore(p, 1),
                      emitScore(p, 2), emitScore(p, 3));

        out << "                                             [ MATL   " << p << " ]      " << p << "      -\n";

        if (!isLast) {
            const int nIL = sIL;       // self; but next MATL's IL is sML + 2 same as sIL when iterating
            std::array<double, 3> cML = {0.0, 0.0, 0.0}; // IL, ML, D
            std::array<double, 3> cD  = {0.0, 0.0, 0.0}; // IL, ML, D
            std::array<double, 3> cIL = {0.0, 0.0, 0.0}; // IL, ML, D
            const int rawP = mainModel.keep[static_cast<size_t>(p - 1)];
            const int rawQ = mainModel.keep[static_cast<size_t>(p)];
            for (size_t si = 0; si < seqs.size(); ++si) {
                const double w = effScale * seqWeight[si];
                if (isMissingTerminal(si, rawP) || isMissingTerminal(si, rawQ)) continue;
                const bool curM = (alnCode[si][static_cast<size_t>(rawP)] >= 0);
                const bool nxtM = (alnCode[si][static_cast<size_t>(rawQ)] >= 0);
                if (curM && nxtM) cML[1] += w;
                else if (curM && !nxtM) cML[2] += w;
                if (!curM && nxtM) cD[1] += w;
                else if (!curM && !nxtM) cD[2] += w;
            }
            const std::array<double, 3> lML = logPost3(cML, ALPHA_MATL_ML_MATL);
            const std::array<double, 3> lD  = logPost3(cD, ALPHA_MATL_D_MATL);
            const std::array<double, 3> lIL = logPost3(cIL, ALPHA_MATL_IL_MATL);
            // ML
            out << stateLine("ML", sML, sML - 1, 3, nIL, 3, 1, W)
                << "  " << fsc(lML[0]) << "  " << fsc(lML[1]) << "  " << fsc(lML[2])
                << emitBuf << "\n";
            // D
            out << stateLine("D", sD, sML - 1, 3, nIL, 3, 0, W)
                << "  " << fsc(lD[0]) << "  " << fsc(lD[1]) << "  " << fsc(lD[2]) << "\n";
            // IL: self-loop
            out << stateLine("IL", sIL, sIL, 3, nIL, 3, 1, W)
                << "  " << fsc(lIL[0]) << "  " << fsc(lIL[1]) << "  " << fsc(lIL[2])
                << "  0.000  0.000  0.000  0.000\n";
        } else {
            // Last position: transitions to IL_last and E
            std::array<double, 2> cML = {0.0, 0.0}; // IL, E
            std::array<double, 2> cD  = {0.0, 0.0}; // IL, E
            std::array<double, 2> cIL = {0.0, 0.0}; // IL, E
            if (N > 0) {
                const int rawL = mainModel.keep[static_cast<size_t>(N - 1)];
                for (size_t si = 0; si < seqs.size(); ++si) {
                    const double w = effScale * seqWeight[si];
                    if (isMissingTerminal(si, rawL)) continue;
                    const bool m = (alnCode[si][static_cast<size_t>(rawL)] >= 0);
                    if (m) cML[1] += w;
                    else   cD[1] += w;
                }
            }
            const std::array<double, 2> lML = logPost2(cML, ALPHA_MATL_ML_END);
            const std::array<double, 2> lD  = logPost2(cD, ALPHA_MATL_D_END);
            const std::array<double, 2> lIL = logPost2(cIL, ALPHA_MATL_IL_END);
            // ML: to IL(self) then E
            out << stateLine("ML", sML, sML - 1, 3, sIL, 2, 1, 1)
                << "  " << fsc(lML[0]) << "  " << fsc(lML[1]) << emitBuf << "\n";
            // D: to IL(self) then E
            out << stateLine("D", sD, sML - 1, 3, sIL, 2, 0, 0)
                << "  " << fsc(lD[0]) << "  " << fsc(lD[1]) << "\n";
            // IL: self-loop then E
            out << stateLine("IL", sIL, sIL, 3, sIL, 2, 1, W)
                << "  " << fsc(lIL[0]) << "  " << fsc(lIL[1])
                << "  0.000  0.000  0.000  0.000\n";
        }
    }

    // END node
    out << "                                             [ END   " << (N + 1) << " ]      -      -\n";
    out << stateLine("E", sEnd, sEnd - 1, 3, -1, 0, 0, 0) << "\n";
    out << "//\n";

    const std::string cmText = out.str();
    DBWriter writer(par.db2.c_str(),
                    par.db2Index.c_str(),
                    1,
                    par.compressed,
                    Parameters::DBTYPE_GENERIC_DB);
    writer.open();
    writer.writeData(cmText.c_str(), cmText.size(), 0, 0);
    writer.close();
    Debug(Debug::INFO) << "cmbuild: wrote Infernal CM (MATL topology, CLEN=" << N
                       << " STATES=" << totalStates << ") to " << par.db2 << "\n";
    return EXIT_SUCCESS;
}
