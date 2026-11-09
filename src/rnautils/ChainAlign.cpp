// chainalign: blastn-style ungapped seed-and-chain alignment of prefilter pairs (prototype for the segmentation search).
//
// For every (query, target) pair of a prefilter db:
//   1. word hits: every exact match of W consecutive letters between query and target (query words hashed once),
//   2. ungapped X-drop extension of each hit along its diagonal (hits inside an earlier extension on the same diagonal
//      are skipped) -> HSPs,
//   3. chaining: the best colinear chain of HSPs, by dynamic programming over HSPs sorted by query start. A link j -> i
//      trims i where it overlaps j, then bridges the space between them with an ungapped stretch on the shorter side
//      plus one gap of the diagonal shift (cost gap open + gap extend x length; the gap at either end of the bridge,
//      whichever scores higher), so the chain score approximates a gapped alignment score,
//   4. optionally (--chain-extend X > 0) gapped X-drop extension outward from both chain ends (affine-gap dynamic
//      programming row by row, cells more than X below the best score pruned, as BLAST's gapped extension), so the
//      alignment runs on through stretches with no word hit, where the chain alone stops,
//   5. the chain is written as one alignment (the bridges' gaps in its backtrace) when it passes -e / -c / --cov-mode,
//      E-value from the chain score with the same gapped statistics as align.
// Sequences are taken as stored (dinucleotide-encoded, amino acid typed dbs). Raw scores in matrix units (bit factor 2).
// Options: --chain-word word length, --chain-xdrop ungapped X-drop, --chain-min-hsp minimum HSP score,
// --chain-bridge longest bridge in letters, --chain-max-hsp HSPs kept per pair (best first), --chain-extend end
// extension X-drop, --chain-extend-eval only chains this significant are extended (most chains are chance ones).
// Gaps cost open + extend x (length - 1) in the extension, as in align.
#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "EvalueComputation.h"
#include "Alignment.h"
#include "Matcher.h"
#include "QueryMatcher.h"
#include "StripedSmithWaterman.h"
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {
struct Hsp {
    int qs, ts, len, score;
};

// Gapped X-drop extension of A[0..) vs B[0..) from the cell before both (the chain end, already scored): affine gaps
// (Gotoh), H not clamped at 0, cells below best - X pruned, each row spanning only its live columns.  Returns the best
// score gained and, in ops, the alignment from the origin to the best cell ('M', 'I' = A letter, 'D' = B letter);
// aLen / bLen: letters of A and B it covers.
struct XdropExtender {
    enum : int { NEG = -1000000000 };
    std::vector<int> H0, F0, H1, F1;
    std::vector<uint8_t> tb;               // per cell: bits 0-1 H source (0 diag, 1 E, 2 F), bit 2 E extended, bit 3 F
    std::vector<size_t> rowOff;
    std::vector<int> rowLo, rowHi;

    int extend(const unsigned char *A, const short *biasA, int la, const unsigned char *B, int lb, short **M, int open,
               int ext, int X, std::string &ops, int &aLen, int &bLen) {
        ops.clear();
        aLen = bLen = 0;
        if (la == 0 || lb == 0) {
            return 0;
        }
        tb.clear();
        rowOff.clear();
        rowLo.clear();
        rowHi.clear();
        // row 0: horizontal gap from the origin
        int best = 0, bi = 0, bj = 0;
        H0.assign(1, 0);
        F0.assign(1, NEG);
        int lo = 0, hi = 0;
        rowOff.push_back(0);
        tb.push_back(0);
        {
            int e = NEG;
            for (int j = 1; j <= lb; j++) {
                const int hl = H0[j - 1];
                const int eo = hl == NEG ? NEG : hl - open, ee = e == NEG ? NEG : e - ext;
                e = std::max(eo, ee);
                if (e < best - X) {
                    break;
                }
                H0.push_back(e);
                F0.push_back(NEG);
                tb.push_back(static_cast<uint8_t>(1 | (ee >= eo && ee != NEG ? 4 : 0)));
                hi = j;
            }
        }
        rowLo.push_back(lo);
        rowHi.push_back(hi);
        for (int i = 1; i <= la; i++) {
            // previous row: columns [lo, hi] in H0 / F0 (index j - lo)
            const int plo = lo, phi = hi;
            int nlo = -1, nhi = -1;
            H1.clear();
            F1.clear();
            const size_t off = tb.size();
            int e = NEG;
            int j = plo;
            for (;; j++) {
                if (j > lb) {
                    break;
                }
                const bool inPrev = j <= phi;
                const int hUp = inPrev ? H0[j - plo] : NEG, fUp = inPrev ? F0[j - plo] : NEG;
                const int hDiag = (j - 1 >= plo && j - 1 <= phi) ? H0[j - 1 - plo] : NEG;
                const int fo = hUp == NEG ? NEG : hUp - open, fe = fUp == NEG ? NEG : fUp - ext;
                const int f = std::max(fo, fe);
                const int hLeft = H1.empty() ? NEG : H1.back();
                const int eo = hLeft == NEG ? NEG : hLeft - open, ee = e == NEG ? NEG : e - ext;
                e = std::max(eo, ee);
                const int d = (hDiag == NEG || j == 0) ? NEG : hDiag + M[A[i - 1]][B[j - 1]] + biasA[i - 1];
                int h = d;
                uint8_t bits = 0;
                if (e > h) {
                    h = e;
                    bits = 1;
                }
                if (f > h) {
                    h = f;
                    bits = 2;
                }
                if (ee >= eo && ee != NEG) {
                    bits |= 4;
                }
                if (fe >= fo && fe != NEG) {
                    bits |= 8;
                }
                if (h < best - X || h == NEG) {
                    h = NEG;
                }
                if (j > phi + 1 && h == NEG) {
                    break;                 // past the previous row only the horizontal gap can live
                }
                H1.push_back(h);
                F1.push_back(h == NEG ? NEG : f);
                tb.push_back(bits);
                if (h != NEG) {
                    if (nlo < 0) {
                        nlo = j;
                    }
                    nhi = j;
                    if (h > best) {
                        best = h;
                        bi = i;
                        bj = j;
                    }
                }
            }
            if (nlo < 0) {
                break;
            }
            // keep columns [nlo, nhi] of this row
            const int cur0 = plo;
            H0.assign(H1.begin() + (nlo - cur0), H1.begin() + (nhi - cur0) + 1);
            F0.assign(F1.begin() + (nlo - cur0), F1.begin() + (nhi - cur0) + 1);
            rowOff.push_back(off);
            rowLo.push_back(cur0);             // tb row i starts at column cur0
            rowHi.push_back(cur0 + static_cast<int>(H1.size()) - 1);
            lo = nlo;
            hi = nhi;
            if (tb.size() > (64u << 20)) {
                break;
            }
        }
        // traceback from the best cell
        int i = bi, j = bj, state = 0;
        while (i > 0 || j > 0) {
            const uint8_t bits = tb[rowOff[i] + (j - rowLo[i])];
            if (state == 0) {
                const int src = bits & 3;
                if (src == 0) {
                    ops.push_back('M');
                    i--;
                    j--;
                    continue;
                }
                state = src;
            }
            if (state == 1) {
                ops.push_back('D');
                state = (bits & 4) ? 1 : 0;
                j--;
            } else {
                ops.push_back('I');
                state = (bits & 8) ? 2 : 0;
                i--;
            }
        }
        std::reverse(ops.begin(), ops.end());
        aLen = bi;
        bLen = bj;
        return best;
    }
};
}

int chainalign(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> tdbr(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    tdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> prefdbr(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                   DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    prefdbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    if (Parameters::isEqualDbtype(qdbr.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) == false ||
        Parameters::isEqualDbtype(tdbr.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) == false) {
        Debug(Debug::ERROR) << "chainalign needs dinucleotide-encoded (amino acid typed) query and target databases\n";
        return EXIT_FAILURE;
    }

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    short **const M = subMat.subMatrix;
    const int gapOpen = par.gapOpen.values.aminoacid();
    const int gapExtend = par.gapExtend.values.aminoacid();
    EvalueComputation evaluer(tdbr.getAminoAcidDBSize(), &subMat, gapOpen, gapExtend);
    const int W = par.chainWord;
    const int XU = par.chainXdrop;
    const int minHsp = par.chainMinHsp;
    const int maxBridge = par.chainBridge;
    const size_t maxHsp = static_cast<size_t>(par.chainMaxHsp);
    const int XG = par.chainExtend;
    const double XE = par.chainExtendEval;
    Debug(Debug::INFO) << "word " << W << ", ungapped X-drop " << XU << ", min HSP " << minHsp << ", bridge "
                       << maxBridge << ", HSPs per pair " << maxHsp << ", end extension X-drop " << XG << " (chains with E <= " << XE << ")\n";
    const uint64_t wordMask = (W == 12) ? ~0ULL : ((1ULL << (5 * W)) - 1);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    size_t nPairs = 0, nNoHsp = 0, nPassed = 0, nHsp = 0, nChained = 0, nExtLetters = 0;
    Debug::Progress progress(prefdbr.getSize());
#pragma omp parallel reduction(+: nPairs, nNoHsp, nPassed, nHsp, nChained, nExtLetters)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<unsigned char> q, t;
        std::vector<uint64_t> qCode;          // word code ending at each query position
        std::vector<int> head, next;          // query word hash: bucket -> last position, position -> previous
        std::vector<int> diagEnd;             // target position up to which a diagonal has been extended
        std::vector<int> touched;
        std::vector<Hsp> hsps;
        std::vector<int> best, prev, trimOf, bridgeQ, gapFirst;
        char buffer[1024 + 32768 * 4];
        std::string out, backtrace, btF, btB;
        XdropExtender xd;
        std::vector<unsigned char> ra, rb;
        std::vector<short> qb, rqb;             // query composition-bias correction per position (0 when off)
        std::vector<float> fbias;

        // score of the ungapped stretch q[qs..qs+n) vs t[ts..ts+n)
        auto stretch = [&](int qs, int ts, int n) {
            int s = 0;
            for (int k = 0; k < n; k++) {
                s += M[q[qs + k]][t[ts + k]] + qb[qs + k];
            }
            return s;
        };

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < prefdbr.getSize(); i++) {
            progress.updateProgress();
            const unsigned int qKey = prefdbr.getDbKey(i);
            const size_t qId = qdbr.getId(qKey);
            if (qId == UINT_MAX) {
                Debug(Debug::ERROR) << "Query " << qKey << " is not in " << par.db1 << "\n";
                EXIT(EXIT_FAILURE);
            }
            const char *qSeq = qdbr.getData(qId, thread_idx);
            const int qLen = static_cast<int>(qdbr.getSeqLen(qId));
            q.resize(qLen);
            for (int k = 0; k < qLen; k++) {
                q[k] = subMat.aa2num[static_cast<int>(qSeq[k])];
            }
            // the query's local composition-bias correction, rounded as align does (added to every score of a query
            // letter: word extensions, bridges, end extension), so scores and E-values follow align's
            qb.assign(qLen, 0);
            if (par.compBiasCorrection) {
                fbias.resize(qLen);
                SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, q.data(), qLen, fbias.data(),
                                                              par.compBiasCorrectionScale);
                for (int k = 0; k < qLen; k++) {
                    qb[k] = static_cast<short>(fbias[k] < 0.0f ? fbias[k] - 0.5f : fbias[k] + 0.5f);
                }
            }
            // hash the query's words
            size_t hsize = 16;
            while (hsize < 2 * static_cast<size_t>(qLen)) {
                hsize <<= 1;
            }
            head.assign(hsize, -1);
            next.assign(qLen, -1);
            qCode.assign(qLen, 0);
            uint64_t code = 0;
            for (int k = 0; k < qLen; k++) {
                code = ((code << 5) | q[k]) & wordMask;
                qCode[k] = code;
                if (k >= W - 1) {
                    const size_t b = (code * 0x9E3779B97F4A7C15ULL) >> 32 & (hsize - 1);
                    next[k] = head[b];
                    head[b] = k;
                }
            }
            out.clear();
            char *data = prefdbr.getData(i, thread_idx);
            while (*data != '\0') {
                hit_t h = QueryMatcher::parsePrefilterHit(data);
                data = Util::skipLine(data);
                const size_t tId = tdbr.getId(h.seqId);
                if (tId == UINT_MAX) {
                    continue;
                }
                nPairs++;
                const char *tSeq = tdbr.getData(tId, thread_idx);
                const int tLen = static_cast<int>(tdbr.getSeqLen(tId));
                if (Util::canBeCovered(par.covThr, par.covMode, static_cast<float>(qLen), static_cast<float>(tLen)) == false) {
                    continue;
                }
                t.resize(tLen);
                for (int k = 0; k < tLen; k++) {
                    t[k] = subMat.aa2num[static_cast<int>(tSeq[k])];
                }
                // 1-2. word hits -> ungapped X-drop extensions (diagonal index qp - tp + tLen)
                if (diagEnd.size() < static_cast<size_t>(qLen + tLen + 1)) {
                    diagEnd.resize(qLen + tLen + 1, -1);
                }
                hsps.clear();
                code = 0;
                for (int tp = 0; tp < tLen; tp++) {
                    code = ((code << 5) | t[tp]) & wordMask;
                    if (tp < W - 1) {
                        continue;
                    }
                    const size_t b = (code * 0x9E3779B97F4A7C15ULL) >> 32 & (hsize - 1);
                    for (int qp = head[b]; qp >= 0; qp = next[qp]) {
                        if (qCode[qp] != code) {
                            continue;
                        }
                        const int d = qp - tp + tLen;
                        if (diagEnd[d] > tp) {
                            continue;
                        }
                        if (diagEnd[d] < 0) {
                            touched.push_back(d);
                        }
                        // word occupies q[qp-W+1..qp], t[tp-W+1..tp]; extend right then left
                        int score = stretch(qp - W + 1, tp - W + 1, W);
                        int s = score, bestS = score, r = 0, bestR = 0;
                        while (qp + r + 1 < qLen && tp + r + 1 < tLen) {
                            r++;
                            s += M[q[qp + r]][t[tp + r]] + qb[qp + r];
                            if (s > bestS) {
                                bestS = s;
                                bestR = r;
                            } else if (bestS - s > XU) {
                                break;
                            }
                        }
                        s = bestS;
                        int l = 0, bestL = 0;
                        const int q0 = qp - W + 1, t0 = tp - W + 1;
                        while (q0 - l - 1 >= 0 && t0 - l - 1 >= 0) {
                            l++;
                            s += M[q[q0 - l]][t[t0 - l]] + qb[q0 - l];
                            if (s > bestS) {
                                bestS = s;
                                bestL = l;
                            } else if (bestS - s > XU) {
                                break;
                            }
                        }
                        diagEnd[d] = tp + bestR + 1;
                        if (bestS >= minHsp) {
                            hsps.push_back(Hsp{q0 - bestL, t0 - bestL, W + bestL + bestR, bestS});
                        }
                    }
                }
                for (int d : touched) {
                    diagEnd[d] = -1;
                }
                touched.clear();
                if (hsps.empty()) {
                    nNoHsp++;
                    continue;
                }
                nHsp += hsps.size();
                if (hsps.size() > maxHsp) {
                    std::partial_sort(hsps.begin(), hsps.begin() + maxHsp, hsps.end(),
                                      [](const Hsp &a, const Hsp &b) { return a.score > b.score; });
                    hsps.resize(maxHsp);
                }
                std::sort(hsps.begin(), hsps.end(), [](const Hsp &a, const Hsp &b) {
                    return a.qs < b.qs || (a.qs == b.qs && a.ts < b.ts);
                });
                // 3. chaining
                const int n = static_cast<int>(hsps.size());
                best.assign(n, 0);
                prev.assign(n, -1);
                trimOf.assign(n, 0);
                bridgeQ.assign(n, 0);
                gapFirst.assign(n, 0);
                int top = 0;
                for (int a = 0; a < n; a++) {
                    const Hsp &hi = hsps[a];
                    best[a] = hi.score;
                    for (int j = 0; j < a; j++) {
                        const Hsp &hj = hsps[j];
                        const int qe = hj.qs + hj.len, te = hj.ts + hj.len;   // exclusive ends of j
                        if (hi.qs + hi.len <= qe || hi.ts + hi.len <= te || hi.ts <= hj.ts) {
                            continue;
                        }
                        const int ov = std::max(0, std::max(qe - hi.qs, te - hi.ts));
                        if (ov >= hi.len) {
                            continue;
                        }
                        const int qs = hi.qs + ov, ts = hi.ts + ov;
                        const int gq = qs - qe, gt = ts - te;
                        if (gq > maxBridge || gt > maxBridge) {
                            continue;
                        }
                        const int iScore = ov ? hi.score - stretch(hi.qs, hi.ts, ov) : hi.score;
                        const int m = std::min(gq, gt), g = std::abs(gq - gt);
                        const int gapCost = g ? gapOpen + gapExtend * g : 0;
                        // bridge: m aligned letters, then the gap (gap first: after j; else: before i)
                        const int sAfter = stretch(qe, te, m);                  // M right after j, gap before i
                        const int sBefore = stretch(qs - m, ts - m, m);         // gap right after j, M before i
                        const int bridge = std::max(sAfter, sBefore) - gapCost;
                        const int cand = best[j] + bridge + iScore;
                        if (cand > best[a]) {
                            best[a] = cand;
                            prev[a] = j;
                            trimOf[a] = ov;
                            bridgeQ[a] = m;
                            gapFirst[a] = sBefore > sAfter;
                        }
                    }
                    if (best[a] > best[top]) {
                        top = a;
                    }
                }
                // 4. the chain as one alignment
                std::vector<int> chain;
                for (int a = top; a >= 0; a = prev[a]) {
                    chain.push_back(a);
                }
                std::reverse(chain.begin(), chain.end());
                nChained += chain.size();
                backtrace.clear();
                int ident = 0;
                auto addM = [&](int qs, int ts, int len) {
                    for (int k = 0; k < len; k++) {
                        ident += (q[qs + k] == t[ts + k]);
                    }
                    backtrace.append(len, 'M');
                };
                const Hsp &h0 = hsps[chain[0]];
                int qe = h0.qs + h0.len, te = h0.ts + h0.len;
                addM(h0.qs, h0.ts, h0.len);
                for (size_t c = 1; c < chain.size(); c++) {
                    const int a = chain[c];
                    const Hsp &hi = hsps[a];
                    const int qs = hi.qs + trimOf[a], ts = hi.ts + trimOf[a], len = hi.len - trimOf[a];
                    const int gq = qs - qe, gt = ts - te, m = bridgeQ[a];
                    const char gapCh = gq > gt ? 'I' : 'D';
                    const int g = std::abs(gq - gt);
                    if (gapFirst[a]) {
                        backtrace.append(g, gapCh);
                        addM(qs - m, ts - m, m);
                    } else {
                        addM(qe, te, m);
                        backtrace.append(g, gapCh);
                    }
                    addM(qs, ts, len);
                    qe = qs + len;
                    te = ts + len;
                }
                int rawScore = best[top];
                int qs0 = h0.qs, ts0 = h0.ts;
                if (XG > 0 && evaluer.computeEvalue(rawScore, qLen) <= XE) {
                    // 4. gapped X-drop extension outward from the chain's last and first columns
                    int la, lb;
                    const int gF = xd.extend(q.data() + qe, qb.data() + qe, qLen - qe, t.data() + te, tLen - te, M, gapOpen, gapExtend,
                                             XG, btF, la, lb);
                    if (gF > 0) {
                        int qi = qe, ti = te;
                        for (char op : btF) {
                            if (op == 'M') {
                                ident += (q[qi] == t[ti]);
                                qi++;
                                ti++;
                            } else if (op == 'I') {
                                qi++;
                            } else {
                                ti++;
                            }
                        }
                        rawScore += gF;
                        backtrace.append(btF);
                        nExtLetters += btF.size();
                        qe += la;
                        te += lb;
                    }
                    ra.assign(q.rend() - qs0, q.rend());
                    rqb.assign(qb.rend() - qs0, qb.rend());
                    rb.assign(t.rend() - ts0, t.rend());
                    const int gB = xd.extend(ra.data(), rqb.data(), qs0, rb.data(), ts0, M, gapOpen, gapExtend, XG, btB, la, lb);
                    if (gB > 0) {
                        // btB runs outward (right to left in the sequences): reverse it into sequence order
                        std::reverse(btB.begin(), btB.end());
                        int qi = qs0 - la, ti = ts0 - lb;
                        for (char op : btB) {
                            if (op == 'M') {
                                ident += (q[qi] == t[ti]);
                                qi++;
                                ti++;
                            } else if (op == 'I') {
                                qi++;
                            } else {
                                ti++;
                            }
                        }
                        rawScore += gB;
                        backtrace.insert(0, btB);
                        nExtLetters += btB.size();
                        qs0 -= la;
                        ts0 -= lb;
                    }
                }
                const double evalue = evaluer.computeEvalue(rawScore, qLen);
                const int bits = static_cast<int>(evaluer.computeBitScore(rawScore) + 0.5);
                const unsigned int alnLen = static_cast<unsigned int>(backtrace.size());
                const float qCov = SmithWaterman::computeCov(qs0, qe - 1, qLen);
                const float tCov = SmithWaterman::computeCov(ts0, te - 1, tLen);
                Matcher::result_t res(h.seqId, bits, qCov, tCov, static_cast<float>(ident) / alnLen, evalue, alnLen,
                                      qs0, qe - 1, qLen, ts0, te - 1, tLen, backtrace);
                const bool isIdentity = (qKey == h.seqId);
                if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr,
                                             par.covMode, par.covThr) == false) {
                    continue;
                }
                nPassed++;
                size_t len = Matcher::resultToBuffer(buffer, res, par.addBacktrace, true);
                out.append(buffer, len);
            }
            writer.writeData(out.c_str(), out.size(), qKey, thread_idx);
        }
    }
    writer.close();
    Debug(Debug::INFO) << nPairs << " pairs, " << nNoHsp << " without an HSP, " << nHsp << " HSPs, " << nChained
                       << " HSPs in the best chains, " << nExtLetters << " columns added by end extension, " << nPassed
                       << " alignments written\n";
    prefdbr.close();
    tdbr.close();
    qdbr.close();
    return EXIT_SUCCESS;
}
