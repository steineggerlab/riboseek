#include "RnaRealign.h"

#include "LocalParameters.h"
#include "DBReader.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "Sequence.h"
#include "RnaMatcher.h"

#include <algorithm>
#include <mutex>
#include <vector>

namespace {

const unsigned int MAX_WINDOW = 64 * 1024;

struct Context {
    SubstitutionMatrix *subMat;
    EvalueComputation *evaluer;
    int gapOpen;
    int gapExtend;
    int querySeqType;
    int targetSeqType;
    bool compBiasCorrection;
    float compBiasCorrectionScale;
    int covMode;
    float covThr;
    double evalThr;
    int seqIdMode;
};

struct ThreadLocal {
    Sequence *qSeq;
    Sequence *wSeq;
    RnaMatcher *matcher;
};

Context *ctx = NULL;
std::mutex ctxMutex;
bool ctxReady = false;
std::vector<ThreadLocal> threadLocals;

// Built on first use: the hook is installed at start-up, long before parameters are
// parsed, and the E-value statistics need the real target database size.
void ensureContext(size_t targetTotalResidues) {
    if (ctxReady) {
        return;
    }
    std::lock_guard<std::mutex> lock(ctxMutex);
    if (ctxReady) {
        return;
    }
    LocalParameters &par = LocalParameters::getLocalInstance();
    ctx = new Context();
    ctx->subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, par.scoreBias);
    ctx->gapOpen = par.gapOpen.values.aminoacid();
    ctx->gapExtend = par.gapExtend.values.aminoacid();
    ctx->evaluer = new EvalueComputation(targetTotalResidues ? targetTotalResidues : 1,
                                         ctx->subMat, ctx->gapOpen, ctx->gapExtend, true, false);
    // Match the align module: the query is scored as a PROFILE. splitstrand types the
    // query DB as HMM_PROFILE | SRC_SEQUENCE while leaving raw characters on disk, and
    // Sequence::mapSequence then builds the profile on the fly via
    // dinucMapProfile -> dinucBuildProfileFromSequence. Declaring the same type here
    // reproduces that, so ssw_init takes the profile path with the identical scale
    // instead of falling back to the raw substitution matrix.
    const uint16_t dinuc = LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
    ctx->querySeqType = DBReader<unsigned int>::setExtendedDbtype(
        Parameters::DBTYPE_HMM_PROFILE, dinuc | Parameters::DBTYPE_EXTENDED_SRC_SEQUENCE);
    ctx->targetSeqType = DBReader<unsigned int>::setExtendedDbtype(
        Parameters::DBTYPE_AMINO_ACIDS, dinuc);
    ctx->compBiasCorrection = par.compBiasCorrection;
    ctx->compBiasCorrectionScale = par.compBiasCorrectionScale;
    ctx->covMode = par.covMode;
    ctx->covThr = par.covThr;
    ctx->evalThr = par.evalThr;
    ctx->seqIdMode = par.seqIdMode;
    size_t nThreads = par.threads > 0 ? (size_t)par.threads : 1;
    threadLocals.assign(nThreads, ThreadLocal());
    for (size_t i = 0; i < nThreads; i++) {
        threadLocals[i].qSeq = NULL;
        threadLocals[i].wSeq = NULL;
        threadLocals[i].matcher = NULL;
    }
    ctxReady = true;
}

char complementBase(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': case 'U': case 'u': return 'A';
        default: return 'N';
    }
}

}  // namespace

unsigned int rnaRealignMaxWindow() {
    return MAX_WINDOW;
}

bool rnaRealignWindow(const char *querySeq, unsigned int queryLen, bool reverse,
                      const char *targetWindow, unsigned int windowLen,
                      unsigned int queryKey, unsigned int targetKey,
                      size_t targetTotalResidues, int threadIdx,
                      RnaRealignHit &out, std::string &backtrace) {
    if (querySeq == NULL || targetWindow == NULL || queryLen == 0 || windowLen == 0
        || windowLen > MAX_WINDOW || queryLen > MAX_WINDOW) {
        return false;
    }
    ensureContext(targetTotalResidues);
    if (ctx == NULL || (size_t)threadIdx >= threadLocals.size()) {
        return false;
    }
    ThreadLocal &tl = threadLocals[threadIdx];
    if (tl.matcher == NULL) {
        tl.qSeq = new Sequence(MAX_WINDOW, ctx->querySeqType, ctx->subMat, 0, false, ctx->compBiasCorrection);
        tl.wSeq = new Sequence(MAX_WINDOW, ctx->targetSeqType, ctx->subMat, 0, false, ctx->compBiasCorrection);
        tl.matcher = new RnaMatcher(ctx->querySeqType, MAX_WINDOW, ctx->subMat, ctx->evaluer,
                                    ctx->compBiasCorrection, ctx->compBiasCorrectionScale,
                                    ctx->gapOpen, ctx->gapExtend, 0.0f, 0);
    }

    std::string qRev;
    const char *qUse = querySeq;
    if (reverse) {
        qRev.assign(queryLen, 'N');
        for (unsigned int i = 0; i < queryLen; i++) {
            qRev[i] = complementBase(querySeq[queryLen - 1 - i]);
        }
        qUse = qRev.c_str();
    }

    // dbKey 0 (even) on purpose: dinucBuildProfileFromSequence reverses the profile for
    // odd keys, and the reverse strand is already handled by revcomp'ing qUse above.
    tl.qSeq->mapSequence(0, 0, qUse, queryLen);
    tl.wSeq->mapSequence(0, targetKey, targetWindow, windowLen);
    tl.matcher->initQuery(tl.qSeq);
    RnaMatcher::result_t r = tl.matcher->getSWResult(tl.wSeq, 0, false, ctx->covMode, ctx->covThr,
                                                     ctx->evalThr, RnaMatcher::SCORE_COV_SEQID,
                                                     ctx->seqIdMode, false, false, false);
    if (r.dbStartPos < 0 || r.qStartPos < 0) {
        return false;
    }
    out.dbStartPos = r.dbStartPos;
    out.dbEndPos = r.dbEndPos;
    if (reverse) {
        // hand back forward-query coordinates (descending, as the reverse-strand convention)
        out.qStartPos = static_cast<int>(queryLen) - 1 - r.qStartPos;
        out.qEndPos = static_cast<int>(queryLen) - 1 - r.qEndPos;
    } else {
        out.qStartPos = r.qStartPos;
        out.qEndPos = r.qEndPos;
    }
    out.alnLength = r.alnLength;
    out.score = r.score;
    out.qcov = r.qcov;
    out.dbcov = r.dbcov;
    out.seqId = r.seqId;
    out.eval = r.eval;
    backtrace = r.backtrace;
    return true;
}
