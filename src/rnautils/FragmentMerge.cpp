#include "FragmentMerge.h"
#include "RnaRealign.h"

#include "Matcher.h"
#include "DBReader.h"

#include <algorithm>
#include <string>
#include <vector>

// Splitting a long target with `splitsequence` makes every chunk its own target, so an
// alignment that straddles a chunk boundary is reported as two (or more) colinear
// fragments. Once offsetalignment has mapped the coordinates back onto the original
// sequence those fragments become adjacent, and the whole alignment can be recovered by
// re-running Smith-Waterman once over a bounded window of the ORIGINAL target -- a window
// that spans the former boundary, so nothing is cut any more.
//
// Groups that are already a single, well-covered alignment are passed through untouched,
// so the cost is one windowed SW per genuinely fragmented locus.

namespace {

// Two fragments belong to the same alignment when they are adjacent (or slightly
// overlapping) in BOTH query and target coordinates. The tolerance absorbs the few bases
// each chunk's own local alignment trims at its edge.
// TODO: these two parameters need tuning
const int JOIN_TOLERANCE = 50;
// A lone hit covering less than this fraction of the query may have been cut by a boundary.
const float FULL_COVERAGE = 0.95f;

void mergeFragments(std::vector<Matcher::result_t> &results, unsigned int queryKey,
                    DBReader<unsigned int> *qSource, DBReader<unsigned int> *tSource,
                    const unsigned char *splitEntries, size_t splitEntriesSize,
                    int threadIdx) {
    if (qSource == NULL || tSource == NULL || results.size() == 0 || splitEntries == NULL) {
        return;
    }
    // Only entries that splitsequence actually divided can carry a chunk boundary.
    const auto entryWasSplit = [&](unsigned int dbKey) {
        return (size_t)dbKey < splitEntriesSize && splitEntries[dbKey] != 0;
    };

    std::stable_sort(results.begin(), results.end(),
                     [](const Matcher::result_t &a, const Matcher::result_t &b) {
        if (a.dbKey != b.dbKey) return a.dbKey < b.dbKey;
        const bool ar = a.qStartPos > a.qEndPos, br = b.qStartPos > b.qEndPos;
        if (ar != br) return ar < br;
        return std::min(a.dbStartPos, a.dbEndPos) < std::min(b.dbStartPos, b.dbEndPos);
    });

    std::vector<Matcher::result_t> out;
    out.reserve(results.size());

    size_t i = 0;
    while (i < results.size()) {
        const bool reverse = results[i].qStartPos > results[i].qEndPos;
        int tLo = std::min(results[i].dbStartPos, results[i].dbEndPos);
        int tHi = std::max(results[i].dbStartPos, results[i].dbEndPos);
        int qLo = std::min(results[i].qStartPos, results[i].qEndPos);
        int qHi = std::max(results[i].qStartPos, results[i].qEndPos);
        size_t j = i + 1;
        while (j < results.size()
               && results[j].dbKey == results[i].dbKey
               && (results[j].qStartPos > results[j].qEndPos) == reverse
               && std::min(results[j].dbStartPos, results[j].dbEndPos) <= tHi + JOIN_TOLERANCE
               && std::min(results[j].qStartPos, results[j].qEndPos) <= qHi + JOIN_TOLERANCE) {
            tLo = std::min(tLo, std::min(results[j].dbStartPos, results[j].dbEndPos));
            tHi = std::max(tHi, std::max(results[j].dbStartPos, results[j].dbEndPos));
            qLo = std::min(qLo, std::min(results[j].qStartPos, results[j].qEndPos));
            qHi = std::max(qHi, std::max(results[j].qStartPos, results[j].qEndPos));
            j++;
        }

        const unsigned int qLen = results[i].qLen;
        const float qcov = qLen ? static_cast<float>(qHi - qLo + 1) / static_cast<float>(qLen) : 1.0f;
        const bool looksFragmented = entryWasSplit(results[i].dbKey)
                                     && (((j - i) > 1) || (qcov < FULL_COVERAGE));

        bool replaced = false;
        if (looksFragmented) {
            const size_t tId = tSource->getId(results[i].dbKey);
            const size_t qId = qSource->getId(queryKey);
            if (tId != UINT_MAX && qId != UINT_MAX) {
                const char *tData = tSource->getData(tId, threadIdx);
                const int tLen = static_cast<int>(tSource->getSeqLen(tId));
                const char *qData = qSource->getData(qId, threadIdx);
                const unsigned int qDataLen = qSource->getSeqLen(qId);
                // window the ORIGINAL target so the realignment can cross the old boundary
                const int pad = static_cast<int>(qDataLen);
                const int wStart = std::max(0, tLo - pad);
                const int wEnd = std::min(tLen, tHi + pad + 1);
                const int wLen = wEnd - wStart;
                if (wLen > 0 && (unsigned int)wLen <= rnaRealignMaxWindow()) {
                    RnaRealignHit hit;
                    std::string backtrace;
                    if (rnaRealignWindow(qData, qDataLen, reverse, tData + wStart, wLen,
                                         queryKey, results[i].dbKey,
                                         tSource->getAminoAcidDBSize(), threadIdx,
                                         hit, backtrace)) {
                        // The realignment now scores the query as a profile, exactly as
                        // the align module does, so its score and E-value are on the same
                        // scale as every other hit and can be adopted directly.
                        Matcher::result_t m = results[i];   // keep orf position fields
                        m.score = hit.score;
                        m.eval = hit.eval;
                        m.qcov = hit.qcov;
                        m.dbcov = hit.dbcov;
                        m.seqId = hit.seqId;
                        m.alnLength = hit.alnLength;
                        m.qStartPos = hit.qStartPos;
                        m.qEndPos = hit.qEndPos;
                        m.dbStartPos = hit.dbStartPos + wStart;
                        m.dbEndPos = hit.dbEndPos + wStart;
                        m.qLen = qDataLen;
                        m.dbLen = tLen;
                        // offsetalignment writes with compress=false, so the backtrace it
                        // hands on must already be in run-length form ("1542M"); convertalis
                        // parses those counts and would otherwise compute alnLen 0.
                        m.backtrace = Matcher::compressAlignment(backtrace);
                        out.push_back(m);
                        replaced = true;
                    }
                }
            }
        }
        if (replaced == false) {
            for (size_t k = i; k < j; k++) {
                out.push_back(results[k]);
            }
        }
        i = j;
    }
    // Grouping needed position order; the rest of the pipeline expects hits ordered by
    // significance. Restore it here rather than relying on offsetalignment's mergeQuery
    // path to re-sort, which does not run when --merge-query 0.
    std::stable_sort(out.begin(), out.end(), Matcher::compareHits);
    results.swap(out);
}

}  // namespace

void registerFragmentMerger() {
    Matcher::registerFragmentMerger(mergeFragments);
}
