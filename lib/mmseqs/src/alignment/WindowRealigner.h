#ifndef WINDOW_REALIGNER_H
#define WINDOW_REALIGNER_H

#include "Matcher.h"
#include "Sequence.h"
#include "DBReader.h"
#include "Parameters.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cstdlib>

// Realigns the hits of result2msa and result2profile whose record has no backtrace.
// A record with positions (ALN_RES_WITHOUT_BT_COL_CNT columns) is realigned inside the target
// window those positions span, on the strand they give: a reverse hit (qStartPos > qEndPos) is
// aligned against the reverse complement of the window. The result is in the record convention,
// with target positions relative to the window: a reverse hit has qStartPos > qEndPos and
// dbStartPos < dbEndPos, and its backtrace runs from qStartPos down (computeMSA flips it back).
// The aligner always runs as Smith-Waterman: a nucleotide query type would select the banded
// nucleotide aligner, which cannot score dinucleotide letters.
// Built on first use and grown to the longest query or window seen, so modules pay for it only
// when a result has records without backtraces.
class WindowRealigner {
public:
    WindowRealigner(int queryDbtype, BaseMatrix *subMat, EvalueComputation *evaluer,
                    bool compBiasCorrection, float compBiasCorrectionScale, int gapOpen, int gapExtend, int zdrop)
        : subMat(subMat), evaluer(evaluer), compBiasCorrection(compBiasCorrection),
          compBiasCorrectionScale(compBiasCorrectionScale), gapOpen(gapOpen), gapExtend(gapExtend), zdrop(zdrop),
          alignType(Parameters::isEqualDbtype(queryDbtype, Parameters::DBTYPE_NUCLEOTIDES)
                    ? DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS,
                                                                DBReader<unsigned int>::getExtendedDbtype(queryDbtype))
                    : queryDbtype),
          capacity(0), matcher(NULL), window(NULL), query(NULL), active(NULL), queryReady(false) {}

    ~WindowRealigner() {
        delete matcher;
        delete window;
        delete query;
    }

    // The query of the hits that follow: center is the module's mapped query, data its raw entry.
    // A query of another type is mapped here, right away: data may be a per-thread buffer that
    // the next read from the same database overwrites.
    void newQuery(Sequence *center, size_t id, unsigned int key, const char *data, unsigned int len) {
        if (center->getSequenceType() == alignType) {
            active = center;
        } else {
            if (query == NULL) {
                query = new Sequence(std::max(len, 2u), alignType, subMat, 0, false, compBiasCorrection);
            }
            query->mapSequence(id, key, data, len);
            active = query;
        }
        queryReady = false;
    }

    // win holds the window's letters as the module mapped them for the hit's strand
    // (mapSequenceReverse for a reverse hit); tLen is the whole target's length.
    Matcher::result_t align(const unsigned char *win, unsigned int winLen, bool reverse,
                            unsigned int dbKey, unsigned int tLen) {
        ensureCapacity(std::max(winLen, static_cast<unsigned int>(active->L)));
        if (queryReady == false) {
            matcher->initQuery(active);
            queryReady = true;
        }
#ifdef RIBOSEEK
        if (reverse) {
            const unsigned char *revcomp = subMat->num2revcompnum;
            for (unsigned int j = 0; j < winLen; j++) {
                window->numSequence[j] = revcomp[win[winLen - 1 - j]];
            }
        } else
#endif
        {
            reverse = false;
            std::copy(win, win + winLen, window->numSequence);
        }
        window->L = static_cast<int>(winLen);
        Matcher::result_t r = matcher->getSWResult(window, INT_MAX, false, 0, 0.0, FLT_MAX,
                                                   Matcher::SCORE_COV_SEQID, 0, false);
        if (r.dbStartPos >= 0 && r.dbEndPos >= 0) {
            if (reverse) {
                // reverse-complement position j is window position winLen - 1 - j
                const int dbStart = r.dbStartPos;
                r.dbStartPos = static_cast<int>(winLen) - 1 - r.dbEndPos;
                r.dbEndPos = static_cast<int>(winLen) - 1 - dbStart;
                std::swap(r.qStartPos, r.qEndPos);
                std::reverse(r.backtrace.begin(), r.backtrace.end());
            }
            r.dbcov = static_cast<float>(std::abs(r.dbEndPos - r.dbStartPos) + 1) / static_cast<float>(std::max(tLen, 1u));
        }
        r.dbKey = dbKey;
        r.dbLen = tLen;
        return r;
    }

private:
    void ensureCapacity(unsigned int len) {
        if (matcher != NULL && len <= capacity) {
            return;
        }
        capacity = std::max(len, 2 * capacity);
        delete matcher;
        delete window;
        matcher = new Matcher(alignType, capacity, subMat, evaluer, compBiasCorrection, compBiasCorrectionScale,
                              gapOpen, gapExtend, 0.0, zdrop);
        window = new Sequence(capacity, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, false);
        queryReady = false;
    }

    BaseMatrix *subMat;
    EvalueComputation *evaluer;
    bool compBiasCorrection;
    float compBiasCorrectionScale;
    int gapOpen;
    int gapExtend;
    int zdrop;
    const int alignType;

    unsigned int capacity;
    Matcher *matcher;
    Sequence *window;
    Sequence *query;    // the query mapped as alignType, when the module's query has another type
    Sequence *active;   // what the aligner is initialised with: the module's query or query
    bool queryReady;
};

#endif
