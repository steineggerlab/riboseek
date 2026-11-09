// xdropalign: blastn-style seed-and-extend alignment of prefilter pairs (prototype for the segmentation search).
//
// For every (query, target) pair of a prefilter db (kmermatcher pairs: target key, score, k-mer diagonal):
//   1. ungapped extension along the k-mer diagonal (BlockAligner::ungappedAlign, the best local segment of the
//      diagonal, as linclust's align2clust),
//   2. a seed inside that segment: the three identical consecutive letters closest to its middle,
//   3. banded gapped X-drop extension from the seed in both directions (BlockAligner::bandedalign, block aligner),
//      which stops where the score drops X below its best, e.g. at a long insertion,
//   4. the alignment (with backtrace) is written when it passes -e / -c / --cov-mode / --min-seq-id, as align.
// The cost of a pair is about its alignment length x the band, not query x target length. Sequences are taken as
// stored (dinucleotide-encoded, amino acid typed dbs). X-drop: RIBOSEEK_XDROP (raw score), default as align2clust
// (32 x gap extend + gap open). RIBOSEEK_XD_TRIGGER: gapped extension only when the ungapped segment's E-value is at
// most this (blastn's trigger; default: always). RIBOSEEK_XD_SAMPLE=n: write 'ungapped-E gapped-E' of the pairs of
// every n-th query to stderr.
#include "LocalParameters.h"
#include "BlockAligner.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "EvalueComputation.h"
#include "Alignment.h"
#include "Matcher.h"
#include "QueryMatcher.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <cstdlib>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

int xdropalign(int argc, const char **argv, const Command &command) {
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
        Debug(Debug::ERROR) << "xdropalign needs dinucleotide-encoded (amino acid typed) query and target databases\n";
        return EXIT_FAILURE;
    }

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(subMat);
    const int gapOpen = par.gapOpen.values.aminoacid();
    const int gapExtend = par.gapExtend.values.aminoacid();
    EvalueComputation evaluer(tdbr.getAminoAcidDBSize(), &subMat, gapOpen, gapExtend);
    const char *xe = getenv("RIBOSEEK_XDROP");
    const int xDrop = xe ? atoi(xe) : 32 * gapExtend + gapOpen;
    const char *te = getenv("RIBOSEEK_XD_TRIGGER");
    const double trigger = te ? atof(te) : -1.0;
    const char *se = getenv("RIBOSEEK_XD_SAMPLE");
    const size_t sample = se ? strtoul(se, NULL, 10) : 0;
    Debug(Debug::INFO) << "X-drop " << xDrop << ", trigger " << trigger << "\n";
    const size_t maxLen = std::max(qdbr.getMaxSeqLen(), tdbr.getMaxSeqLen()) + 1;

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    size_t nPairs = 0, nNoSeed = 0, nPassed = 0, nTrig = 0;
    Debug::Progress progress(prefdbr.getSize());
#pragma omp parallel reduction(+: nPairs, nNoSeed, nPassed, nTrig)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence query(maxLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection);
        Sequence target(maxLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection);
        BlockAligner aligner(Parameters::DBTYPE_AMINO_ACIDS, maxLen, &subMat, &fastMatrix, &evaluer,
                             par.compBiasCorrection, par.compBiasCorrectionScale, -gapOpen, -gapExtend);
        char buffer[1024 + 32768 * 4];
        std::string out, backtrace;

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < prefdbr.getSize(); i++) {
            progress.updateProgress();
            const unsigned int qKey = prefdbr.getDbKey(i);
            const size_t qId = qdbr.getId(qKey);
            if (qId == UINT_MAX) {
                Debug(Debug::ERROR) << "Query " << qKey << " is not in " << par.db1 << "\n";
                EXIT(EXIT_FAILURE);
            }
            char *qSeq = qdbr.getData(qId, thread_idx);
            const int qLen = static_cast<int>(qdbr.getSeqLen(qId));
            query.mapSequence(qId, qKey, qSeq, qLen);
            aligner.initQuery(&query);
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
                char *tSeq = tdbr.getData(tId, thread_idx);
                const int tLen = static_cast<int>(tdbr.getSeqLen(tId));
                if (Util::canBeCovered(par.covThr, par.covMode, static_cast<float>(qLen), static_cast<float>(tLen)) == false) {
                    continue;
                }
                target.mapSequence(tId, h.seqId, tSeq, tLen);
                // 1. ungapped extension on the k-mer diagonal
                BlockAligner::UngappedAln_res u = aligner.ungappedAlign(&target, h.diagonal);
                if (u.qStart < 0 || u.tStart < 0 || u.alnLen < 3) {
                    nNoSeed++;
                    continue;
                }
                if (trigger >= 0 && u.eval > trigger) {
                    nTrig++;
                    continue;
                }
                // 2. seed: three identical consecutive letters, closest to the middle of the ungapped segment
                const int mid = u.alnLen / 2;
                int seed = -1;
                for (int off = 0; off <= u.alnLen && seed < 0; off++) {
                    for (int s = -1; s <= 1 && seed < 0; s += 2) {
                        const int b = mid + s * off;
                        if (b < 0 || b + 2 >= u.alnLen) {
                            continue;
                        }
                        const int qp = u.qStart + b, tp = u.tStart + b;
                        if (qSeq[qp] == tSeq[tp] && qSeq[qp + 1] == tSeq[tp + 1] && qSeq[qp + 2] == tSeq[tp + 2]) {
                            seed = b + 1;
                        }
                        if (off == 0) {
                            break;
                        }
                    }
                }
                if (seed < 0) {
                    nNoSeed++;
                    continue;
                }
                // 3. banded gapped X-drop extension in both directions
                backtrace.clear();
                s_align a = aligner.bandedalign(&target, u.qStart + seed, u.tStart + seed, backtrace, xDrop, 0.0f, 0);
                if (sample && qKey % sample == 0) {
                    fprintf(stderr, "XDS\t%g\t%g\n", u.eval, (a.evalue < 0 || backtrace.empty()) ? 1e300 : a.evalue);
                }
                if (a.evalue < 0 || backtrace.empty()) {
                    continue;
                }
                const unsigned int alnLen = static_cast<unsigned int>(backtrace.size());
                const float seqId = alnLen ? static_cast<float>(a.identicalAACnt) / alnLen : 0.0f;
                Matcher::result_t res(h.seqId, a.score1, a.qCov, a.tCov, seqId, a.evalue, alnLen,
                                      a.qStartPos1, a.qEndPos1, qLen, a.dbStartPos1, a.dbEndPos1, tLen, backtrace);
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
    Debug(Debug::INFO) << nPairs << " pairs, " << nNoSeed << " without a seed, " << nTrig << " below the trigger, "
                       << nPassed << " alignments written\n";
    delete[] fastMatrix.matrix;
    delete[] fastMatrix.matrixData;
    prefdbr.close();
    tdbr.close();
    qdbr.close();
    return EXIT_SUCCESS;
}
