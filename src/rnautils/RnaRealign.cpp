#include "RnaRealign.h"

#include "LocalParameters.h"
#include "DBReader.h"
#include "SubstitutionMatrix.h"
#include "RnaEvalueComputation.h"
#include "Sequence.h"
#include "RnaMatcher.h"
#include "TargetGcStat.h"
#include "EvalueModel.h"
#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <cstdlib>
#include <mutex>
#include <vector>

namespace {

const unsigned int MAX_WINDOW = 64 * 1024;

struct Context {
    SubstitutionMatrix *subMat;
    RnaEvalueComputation *evaluer;
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
    bool profileQuery;   // the realigned query is a genuine profile (iterations after the first)
    bool strandDb;       // queries come from the iteration's strand-split DB (exact rnaalign query)
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
// parsed, and the E-value model needs the real target database.
void ensureContext(DBReader<unsigned int> *targetDb, int queryDbType, int threadIdx) {
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
    const size_t targetTotalResidues = targetDb->getAminoAcidDBSize();
    ctx->evaluer = new RnaEvalueComputation(targetTotalResidues ? targetTotalResidues : 1,
                                         ctx->subMat, ctx->gapOpen, ctx->gapExtend, true, false);
    // The query is the one rnaalign aligned in this iteration: a genuine profile (strand-split
    // result2profile output, base type HMM_PROFILE) or the sequence. It gets the same per-query
    // model and the same kind of target GC histogram as rnaalign -- of the ORIGINAL target, which
    // is what the realigned hit is reported against (2026-10-20).
    const uint16_t dinuc = LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
    const uint16_t strandSplit = LocalParameters::DBTYPE_EXTENDED_STRAND_SPLIT;
    ctx->profileQuery = queryDbType >= 0 && Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE)
                        && (DBReader<unsigned int>::getExtendedDbtype(queryDbType) & Parameters::DBTYPE_EXTENDED_SRC_SEQUENCE) == 0;
    ctx->strandDb = queryDbType >= 0;
    ctx->evaluer->setProfileQuery(ctx->profileQuery);
    bool precomputed = false;
    const std::string targetName = targetDb->getDataFileName();
    TargetGcStat st = TargetGcStat::forTarget(targetName, targetDb, (double)targetTotalResidues, precomputed, threadIdx);
    if (st.empty()) {
        Debug(Debug::ERROR) << "Cannot compute E-values for re-aligned fragments: the target database "
                            << targetName << " has no A/C/G/U residues to build its GC statistic from.\n";
        EXIT(EXIT_FAILURE);
    }
    ctx->evaluer->setTargetGcHistogram(st.deriveHistogram(4, EvalueModel::R0_TARGET_LEN, (double)targetTotalResidues));
    Debug(Debug::INFO) << "E-value for re-aligned fragments: per-query " << (ctx->profileQuery ? "profile" : "sequence")
                       << " model (query " << (ctx->strandDb ? "as aligned in this iteration" : "sequence, legacy path")
                       << "), target GC histogram "
                       << (precomputed ? "from " : "sampled from ") << st.nSeq << " sequences of " << targetName
                       << (precomputed ? " (precomputed)" : "") << "\n";
    // Match the align module's query type exactly (the dbtype the workflow gives rnaalign's query):
    //  - first iteration, a sequence: HMM_PROFILE | strand-split | dinucleotide | SRC_SEQUENCE. Raw
    //    characters on disk; Sequence::mapSequence builds the profile on the fly
    //    (dinucMapProfile -> dinucBuildProfileFromSequence), reversed for an odd (reverse) key;
    //  - later iterations, a profile: HMM_PROFILE | strand-split; the odd key holds splitstrand's
    //    reverse-complemented profile, and ssw_init turns its features forward-equivalent.
    // The legacy path (no strand DB) keeps its old type: sequence, reverse built from the text.
    if (ctx->profileQuery) {
        ctx->querySeqType = DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_HMM_PROFILE, strandSplit);
    } else if (ctx->strandDb) {
        ctx->querySeqType = DBReader<unsigned int>::setExtendedDbtype(
            Parameters::DBTYPE_HMM_PROFILE, strandSplit | dinuc | Parameters::DBTYPE_EXTENDED_SRC_SEQUENCE);
    } else {
        ctx->querySeqType = DBReader<unsigned int>::setExtendedDbtype(
            Parameters::DBTYPE_HMM_PROFILE, dinuc | Parameters::DBTYPE_EXTENDED_SRC_SEQUENCE);
    }
    // A padded (GPU) target already holds dinucleotide letters: read them as they are, the way
    // rnaalign does, instead of encoding them a second time
    if (DBReader<unsigned int>::getExtendedDbtype(targetDb->getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU) {
        ctx->targetSeqType = DBReader<unsigned int>::encodedGpuDbtype(targetDb->getDbtype());
    } else {
        ctx->targetSeqType = DBReader<unsigned int>::setExtendedDbtype(
            Parameters::DBTYPE_AMINO_ACIDS, dinuc);
    }
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

bool rnaRealignWindow(const char *queryData, unsigned int queryLen, int queryDbType,
                      unsigned int queryEntryKey, bool reverse,
                      const char *targetWindow, unsigned int windowLen,
                      unsigned int queryKey, unsigned int targetKey,
                      DBReader<unsigned int> *targetDb, int threadIdx,
                      RnaRealignHit &out, std::string &backtrace) {
    if (queryData == NULL || targetWindow == NULL || queryLen == 0 || windowLen == 0
        || windowLen > MAX_WINDOW || queryLen > MAX_WINDOW) {
        return false;
    }
    ensureContext(targetDb, queryDbType, threadIdx);
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

    if (ctx->strandDb) {
        // the strand copy rnaalign aligned, under its own key (parity = strand)
        tl.qSeq->mapSequence(0, queryEntryKey, queryData, queryLen);
    } else {
        std::string qRev;
        const char *qUse = queryData;
        if (reverse) {
            qRev.assign(queryLen, 'N');
            for (unsigned int i = 0; i < queryLen; i++) {
                qRev[i] = complementBase(queryData[queryLen - 1 - i]);
            }
            qUse = qRev.c_str();
        }
        // dbKey 0 (even) on purpose: dinucBuildProfileFromSequence reverses the profile for
        // odd keys, and the reverse strand is already handled by revcomp'ing qUse above.
        tl.qSeq->mapSequence(0, 0, qUse, queryLen);
    }
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
