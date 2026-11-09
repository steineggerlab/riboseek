// envalign: localized alignment of long query units.
//
// Seeds only localize, Smith-Waterman aligns. The query units are cut into non-overlapping seeding windows
// beforehand (splitsequence --sequence-overlap 0, whose headers carry the source key and offset of each window), the
// windows are prefiltered against whole target units, and this module
//   1. maps the window hits to query-unit diagonals (D = diagonal + window offset),
//   2. chains the hits of each (query unit, target) pair on a consistent diagonal (|dD| <= --chain-band, query gap
//      <= --chain-gap) into envelopes; envelopes that overlap on both sequences once padded become one,
//   3. runs the same Smith-Waterman as align inside each envelope plus --env-margin on every side, on sub-ranges of
//      the sequences in memory, and doubles the margin while the alignment reaches a clipped envelope edge (at most
//      --env-max-iter times),
//   4. writes the alignments in unit coordinates (alignment db keyed by the query unit, as align); -c / --cov-mode
//      apply to unit lengths, before (pairs that cannot be covered are skipped) and after alignment.
// A target of at most --local-min residues is always taken whole. E-values use the whole target database size, as in
// align, so they are comparable with the alignments of the short queries that go through prefilter + align directly.
#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Orf.h"
#include "QueryMatcher.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "Util.h"

#include <algorithm>
#include <climits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

struct Hit {
    int qa, qb;      // query window [qa, qb) in unit coordinates
    int diag;        // query-unit diagonal: qpos - tpos
};

struct Chain {
    int qa, qb, lastDiag, ta, tb;
};

struct Envelope {
    unsigned int tKey;
    int cq0, cq1, ct0, ct1;   // chain core
};

// chain the hits of one (query, target) pair, sorted by window start
std::vector<Chain> chainHits(std::vector<Hit> &hits, int band, int gap) {
    std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
        return std::tie(a.qa, a.qb, a.diag) < std::tie(b.qa, b.qb, b.diag);
    });
    std::vector<Chain> chains;
    for (size_t i = 0; i < hits.size(); i++) {
        const Hit &h = hits[i];
        bool joined = false;
        for (size_t c = 0; c < chains.size(); c++) {
            Chain &x = chains[c];
            if (std::abs(h.diag - x.lastDiag) <= band && h.qa - x.qb <= gap) {
                x.qb = std::max(x.qb, h.qb);
                x.lastDiag = h.diag;
                x.ta = std::min(x.ta, h.qa - h.diag);
                x.tb = std::max(x.tb, h.qb - h.diag);
                joined = true;
                break;
            }
        }
        if (joined == false) {
            Chain x = {h.qa, h.qb, h.diag, h.qa - h.diag, h.qb - h.diag};
            chains.push_back(x);
        }
    }
    return chains;
}

void envelopeOf(const Envelope &e, int margin, int qLen, int tLen, int localMin,
                int &eq0, int &eq1, int &et0, int &et1) {
    eq0 = std::max(0, e.cq0 - margin);
    eq1 = std::min(qLen, e.cq1 + margin);
    if (tLen <= localMin) {
        et0 = 0;
        et1 = tLen;
    } else {
        et0 = std::max(0, e.ct0 - margin);
        et1 = std::min(tLen, e.ct1 + margin);
    }
}

}

int envalign(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> tdbr(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    tdbr.open(DBReader<unsigned int>::NOSORT);
    std::string winHeader = par.db3 + "_h";
    std::string winHeaderIndex = par.db3 + "_h.index";
    DBReader<unsigned int> whdr(winHeader.c_str(), winHeaderIndex.c_str(), par.threads,
                                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    whdr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> prefdbr(par.db4.c_str(), par.db4Index.c_str(), par.threads,
                                   DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    prefdbr.open(DBReader<unsigned int>::NOSORT);

    const int querySeqType = qdbr.getDbtype();
    const int targetSeqType = tdbr.getDbtype();
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS) == false ||
        Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_AMINO_ACIDS) == false) {
        Debug(Debug::ERROR) << "envalign needs dinucleotide-encoded (amino acid typed) query and target databases\n";
        return EXIT_FAILURE;
    }

    // window -> (query unit key, offset)
    struct Window { size_t prefId; int from, len; };
    std::unordered_map<unsigned int, std::vector<Window> > windowsOf;   // unit key -> its windows
    for (size_t i = 0; i < prefdbr.getSize(); i++) {
        unsigned int wKey = prefdbr.getDbKey(i);
        size_t hId = whdr.getId(wKey);
        if (hId == UINT_MAX) {
            Debug(Debug::ERROR) << "Window " << wKey << " has no header in " << winHeader << "\n";
            return EXIT_FAILURE;
        }
        Orf::SequenceLocation loc = Orf::parseOrfHeader(whdr.getData(hId, 0));
        if (loc.id == UINT_MAX) {
            Debug(Debug::ERROR) << "Window " << wKey << " header carries no source location (splitsequence output?)\n";
            return EXIT_FAILURE;
        }
        Window w = {i, static_cast<int>(loc.from), static_cast<int>(loc.to - loc.from + 1)};
        windowsOf[loc.id].push_back(w);
    }
    std::vector<unsigned int> units;
    units.reserve(windowsOf.size());
    for (auto &kv : windowsOf) {
        units.push_back(kv.first);
    }
    std::sort(units.begin(), units.end());

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    const int gapOpen = par.gapOpen.values.aminoacid();
    const int gapExtend = par.gapExtend.values.aminoacid();
    EvalueComputation evaluer(tdbr.getAminoAcidDBSize(), &subMat, gapOpen, gapExtend);
    const unsigned int swMode = Matcher::SCORE_COV_SEQID;
    const size_t maxLen = std::max(qdbr.getMaxSeqLen(), tdbr.getMaxSeqLen()) + 1;

    DBWriter writer(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    size_t nEnvelopes = 0, nAligned = 0, nWidened = 0, nBroken = 0, nCapped = 0;
    Debug::Progress progress(units.size());
#pragma omp parallel reduction(+: nEnvelopes, nAligned, nWidened, nBroken, nCapped)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence qSeq(maxLen, querySeqType, &subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(maxLen, targetSeqType, &subMat, 0, false, par.compBiasCorrection);
        Matcher matcher(querySeqType, maxLen, &subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale,
                        gapOpen, gapExtend, par.correlationScoreWeight, par.zdrop);
        char buffer[1024 + 32768 * 4];
        std::string out;
        std::vector<Matcher::result_t> results;

#pragma omp for schedule(dynamic, 1)
        for (size_t u = 0; u < units.size(); u++) {
            progress.updateProgress();
            const unsigned int qKey = units[u];
            const size_t qId = qdbr.getId(qKey);
            if (qId == UINT_MAX) {
                Debug(Debug::ERROR) << "Query unit " << qKey << " is not in " << par.db1 << "\n";
                EXIT(EXIT_FAILURE);
            }
            const char *qData = qdbr.getData(qId, thread_idx);
            const int qLen = static_cast<int>(qdbr.getSeqLen(qId));

            // 1. window hits -> unit diagonals, grouped by target
            std::map<unsigned int, std::vector<Hit> > byTarget;
            const std::vector<Window> &wins = windowsOf.at(qKey);
            for (size_t w = 0; w < wins.size(); w++) {
                const int off = wins[w].from;
                char *data = prefdbr.getData(wins[w].prefId, thread_idx);
                while (*data != '\0') {
                    hit_t h = QueryMatcher::parsePrefilterHit(data);
                    Hit x;
                    x.qa = off;
                    x.qb = std::min(qLen, off + wins[w].len);
                    // the prefilter keeps diagonals as 16 bit (qpos - tpos of the window); a real diagonal lies in
                    // (-tLen, window length), so values at or above the window length wrapped around
                    int d = static_cast<short>(h.diagonal);
                    if (d >= wins[w].len) {
                        d -= 65536;
                    }
                    x.diag = d + off;
                    byTarget[h.seqId].push_back(x);
                    data = Util::skipLine(data);
                }
            }

            // 2. chains -> envelopes (merged when they overlap on both sequences once padded)
            std::vector<Envelope> envs;
            for (auto &kv : byTarget) {
                const size_t tId = tdbr.getId(kv.first);
                if (tId == UINT_MAX) {
                    continue;
                }
                const int tLen = static_cast<int>(tdbr.getSeqLen(tId));
                // -c / --cov-mode on unit lengths, as align: skip pairs that cannot reach the coverage
                if (Util::canBeCovered(par.covThr, par.covMode, static_cast<float>(qLen), static_cast<float>(tLen)) == false) {
                    continue;
                }
                std::vector<Chain> chains = chainHits(kv.second, par.chainBand, par.chainGap);
                std::vector<Envelope> pairEnvs;
                for (size_t c = 0; c < chains.size(); c++) {
                    Envelope e = {kv.first, chains[c].qa, chains[c].qb,
                                  std::max(0, chains[c].ta), std::min(tLen, chains[c].tb)};
                    if (e.ct1 <= e.ct0) {
                        continue;
                    }
                    bool merged = false;
                    for (size_t p = 0; p < pairEnvs.size(); p++) {
                        Envelope &m = pairEnvs[p];
                        if (std::min(m.cq1, e.cq1) + 2 * par.envMargin > std::max(m.cq0, e.cq0) &&
                            std::min(m.ct1, e.ct1) + 2 * par.envMargin > std::max(m.ct0, e.ct0)) {
                            m.cq0 = std::min(m.cq0, e.cq0); m.cq1 = std::max(m.cq1, e.cq1);
                            m.ct0 = std::min(m.ct0, e.ct0); m.ct1 = std::max(m.ct1, e.ct1);
                            merged = true;
                            break;
                        }
                    }
                    if (merged == false) {
                        pairEnvs.push_back(e);
                    }
                }
                envs.insert(envs.end(), pairEnvs.begin(), pairEnvs.end());
            }
            nEnvelopes += envs.size();

            // 3. Smith-Waterman inside each envelope, margin doubled while the alignment touches a clipped edge
            results.clear();
            std::set<std::tuple<unsigned int, int, int, int, int> > seen;
            for (size_t i = 0; i < envs.size(); i++) {
                const Envelope &e = envs[i];
                const size_t tId = tdbr.getId(e.tKey);
                const char *tData = tdbr.getData(tId, thread_idx);
                const int tLen = static_cast<int>(tdbr.getSeqLen(tId));
                int margin = par.envMargin;
                int eq0 = 0, eq1 = 0, et0 = 0, et1 = 0;
                envelopeOf(e, margin, qLen, tLen, par.localMin, eq0, eq1, et0, et1);
                bool have = false;
                Matcher::result_t best;
                for (int it = 0; it <= par.envMaxIter; it++) {
                    qSeq.mapSequence(qId, qKey, qData + eq0, eq1 - eq0);
                    matcher.initQuery(&qSeq);
                    tSeq.mapSequence(tId, e.tKey, tData + et0, et1 - et0);
                    Matcher::result_t res = matcher.getSWResult(&tSeq, INT_MAX, false, par.covMode, 0.0f, par.evalThr,
                                                                swMode, par.seqIdMode, false);
                    if (res.eval > par.evalThr) {
                        break;              // rejected before the start is computed (start stays -1)
                    }
                    if (res.qStartPos < 0 || res.dbStartPos < 0) {
                        nBroken++;          // passing score but no start: 16-bit score saturation
                        break;
                    }
                    if (res.seqId < par.seqIdThr - 1e-6 || res.alnLength < static_cast<unsigned int>(par.alnLenThr)) {
                        break;
                    }
                    best = res;
                    have = true;
                    const bool edge = (res.qStartPos == 0 && eq0 > 0) || (res.qEndPos == eq1 - eq0 - 1 && eq1 < qLen) ||
                                      (res.dbStartPos == 0 && et0 > 0) || (res.dbEndPos == et1 - et0 - 1 && et1 < tLen);
                    if (edge == false) {
                        break;
                    }
                    int n0, n1, n2, n3;
                    envelopeOf(e, margin * 2, qLen, tLen, par.localMin, n0, n1, n2, n3);
                    if (n0 == eq0 && n1 == eq1 && n2 == et0 && n3 == et1) {
                        break;
                    }
                    if (it == par.envMaxIter) {
                        nCapped++;
                        break;
                    }
                    margin *= 2;
                    eq0 = n0; eq1 = n1; et0 = n2; et1 = n3;
                    nWidened++;
                }
                if (have == false) {
                    continue;
                }
                // 4. unit coordinates
                best.qStartPos += eq0; best.qEndPos += eq0; best.qLen = qLen;
                best.dbStartPos += et0; best.dbEndPos += et0; best.dbLen = tLen;
                best.qcov = SmithWaterman::computeCov(best.qStartPos, best.qEndPos, qLen);
                best.dbcov = SmithWaterman::computeCov(best.dbStartPos, best.dbEndPos, tLen);
                if (Util::hasCoverage(par.covThr, par.covMode, best.qcov, best.dbcov) == false) {
                    continue;
                }
                std::tuple<unsigned int, int, int, int, int> key(best.dbKey, best.qStartPos, best.qEndPos,
                                                                 best.dbStartPos, best.dbEndPos);
                if (seen.insert(key).second == false) {
                    continue;       // two envelopes widened onto the same alignment
                }
                results.push_back(best);
                nAligned++;
            }
            std::sort(results.begin(), results.end(), Matcher::compareHits);
            out.clear();
            for (size_t i = 0; i < results.size(); i++) {
                size_t len = Matcher::resultToBuffer(buffer, results[i], par.addBacktrace);
                out.append(buffer, len);
            }
            writer.writeData(out.c_str(), out.length(), qKey, thread_idx);
        }
    }
    writer.close();
    Debug(Debug::INFO) << units.size() << " query units, " << nEnvelopes << " envelopes, " << nAligned
                       << " alignments, " << nWidened << " widenings, " << nCapped << " still at an edge after "
                       << par.envMaxIter << " widenings, " << nBroken << " saturated (start -1) skipped\n";
    prefdbr.close();
    whdr.close();
    tdbr.close();
    qdbr.close();
    return EXIT_SUCCESS;
}
