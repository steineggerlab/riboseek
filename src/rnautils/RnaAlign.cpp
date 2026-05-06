#include "LocalParameters.h"
#include "LocalCommandDeclarations.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "QueryMatcher.h"
#include "FileUtil.h"
#include "Sequence.h"
#include "Util.h"
#include "FastSort.h"
#include "RnaMatcher.h"
#include "DinucleotideMapping.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>

static bool checkCriteria(RnaMatcher::result_t &res, bool isIdentity, double evalThr, double seqIdThr,
                          int alnLenThr, int covMode, float covThr) {
    const bool evalOk = (res.eval <= evalThr);
    const bool seqIdOK = (res.seqId >= seqIdThr);
    const bool covOK = Util::hasCoverage(covThr, covMode, res.qcov, res.dbcov);
    const bool alnLenOK = Util::hasAlignmentLength(alnLenThr, res.alnLength);
    return isIdentity || (evalOk && seqIdOK && covOK && alnLenOK);
}


int rnaalign(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE,
        "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id",
        NULL, 0);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    // RNA always searches both strands
    const bool bothStrands = true;
    const size_t targetDbSize = par.dbSize;
    float covThr = par.covThr;
    const int covMode = par.covMode;
    const int seqIdMode = par.seqIdMode;
    const double evalThr = par.evalThr;
    const double seqIdThr = par.seqIdThr;
    const int alnLenThr = par.alnLenThr;
    const bool includeIdentity = par.includeIdentity;
    bool addBacktrace = par.addBacktrace;
    const float scoreBias = par.scoreBias;
    const unsigned int threads = static_cast<unsigned int>(par.threads);
    const unsigned int compressed = par.compressed;
    const size_t maxSeqLen = par.maxSeqLen;
    const bool compBiasCorrection = par.compBiasCorrection;
    const float compBiasCorrectionScale = par.compBiasCorrectionScale;
    const unsigned int maxAccept = static_cast<unsigned int>(par.maxAccept);
    const unsigned int maxReject = static_cast<unsigned int>(par.maxRejected);
    const float correlationScoreWeight = par.correlationScoreWeight;
    const int altAlignment = par.altAlignment;  // # of alternative alignments per target

    // Realignment parameters
    const bool realign = par.realign;
    const float realignScoreBias = par.realignScoreBias;
    const int realignMaxSeqs = par.realignMaxSeqs;
    float realignCov = 0.0f;
    unsigned int realignSwMode = RnaMatcher::SCORE_COV;
    if (realign) {
        realignCov = covThr;
        covThr = 0.0f;
        if (addBacktrace == false) {
            addBacktrace = true;
        }
    }

    unsigned int alignmentMode = par.alignmentMode;
    if (alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED) {
        Debug(Debug::ERROR) << "Use rescorediagonal for ungapped alignment mode.\n";
        return EXIT_FAILURE;
    }
    if (addBacktrace) {
        alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }

    // Compute realignSwMode from the full alignment mode (after addBacktrace upgrade)
    // This ensures realignment computes backtrace when needed
    if (realign) {
        unsigned int realignAlnMode = std::max(alignmentMode, (unsigned int)Parameters::ALIGNMENT_MODE_SCORE_COV);
        switch (realignAlnMode) {
            case Parameters::ALIGNMENT_MODE_SCORE_COV:
                realignSwMode = RnaMatcher::SCORE_COV;
                break;
            case Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID:
                realignSwMode = RnaMatcher::SCORE_COV_SEQID;
                break;
            default:
                realignSwMode = RnaMatcher::SCORE_COV;
                break;
        }
    }

    unsigned int swMode;
    switch (alignmentMode) {
        case Parameters::ALIGNMENT_MODE_SCORE_COV:
            swMode = RnaMatcher::SCORE_COV;
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID:
            swMode = RnaMatcher::SCORE_COV_SEQID;
            break;
        case Parameters::ALIGNMENT_MODE_FAST_AUTO:
            if (covThr > 0.0 && seqIdThr == 0.0) {
                swMode = RnaMatcher::SCORE_COV;
            } else if (covThr > 0.0 && seqIdThr > 0.0) {
                swMode = RnaMatcher::SCORE_COV_SEQID;
            } else {
                swMode = RnaMatcher::SCORE_ONLY;
            }
            break;
        default:
            swMode = RnaMatcher::SCORE_ONLY;
            break;
    }

    // When realign is enabled, initial alignment uses SCORE_ONLY (fast pass),
    // then realignment checks coverage with biased matrix
    if (realign) {
        swMode = RnaMatcher::SCORE_ONLY;
    }

    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    IndexReader *tDbrIdx = new IndexReader(par.db2, par.threads,
                              extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                              (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    DBReader<unsigned int> *tdbr = tDbrIdx->sequenceReader;
    int targetSeqType = tdbr->getDbtype();

    bool sameQTDB = (par.db2.compare(par.db1) == 0);
    IndexReader *qDbrIdx = tDbrIdx;
    DBReader<unsigned int> *qdbr = tdbr;
    int querySeqType = targetSeqType;
    if (!sameQTDB) {
        qDbrIdx = new IndexReader(par.db1, par.threads,
                                  extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                                  (touch) ? IndexReader::PRELOAD_INDEX : 0);
        qdbr = qDbrIdx->sequenceReader;
        querySeqType = qdbr->getDbtype();
    }

    // RNA always uses dinucleotide pair encoding.
    // Force DINUCLEOTIDE flag on both query and target so dinucMapSequence activates.
    // Keep the base types as-is (AMINO_ACIDS for sequences, HMM_PROFILE for profiles).
    {
        uint16_t tExt = DBReader<unsigned int>::getExtendedDbtype(targetSeqType);
        tExt |= LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
        targetSeqType = DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS, tExt);
    }
    {
        uint16_t qExt = DBReader<unsigned int>::getExtendedDbtype(querySeqType);
        qExt |= LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
        if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE)) {
            querySeqType = DBReader<unsigned int>::setExtendedDbtype(querySeqType, qExt);
        } else {
            // Sequence input: keep as AMINO_ACIDS + DINUCLEOTIDE.
            // ssw_init will use the substitution matrix directly (non-profile path).
            querySeqType = DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS, qExt);
        }
    }

    DBReader<unsigned int> prefdbr(par.db3.c_str(), par.db3Index.c_str(), threads,
                                    DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    prefdbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    bool reversePrefilterResult = Parameters::isEqualDbtype(prefdbr.getDbtype(), Parameters::DBTYPE_PREFILTER_REV_RES);

    // Use dinucleotide substitution matrix for RNA
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, scoreBias);
    int gapOpen = par.gapOpen.values.aminoacid();
    int gapExtend = par.gapExtend.values.aminoacid();

    // Create biased substitution matrix for realignment
    BaseMatrix *realign_m = NULL;
    if (realign && realignScoreBias != 0.0f) {
        realign_m = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, scoreBias + realignScoreBias);
    }

    size_t effectiveDbSize = targetDbSize == 0 ? tdbr->getAminoAcidDBSize() : targetDbSize;
    EvalueComputation evaluer(effectiveDbSize, &subMat, gapOpen, gapExtend,
                              bothStrands, /*useRnaCorrection=*/false);

    Debug(Debug::INFO) << "Query database size: " << qdbr->getSize() << " type: " << Parameters::getDbTypeName(querySeqType) << "\n";
    Debug(Debug::INFO) << "Target database size: " << tdbr->getSize() << " type: " << Parameters::getDbTypeName(targetSeqType) << "\n";

    int dbtype = Parameters::DBTYPE_ALIGNMENT_RES;
    dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype, DBReader<unsigned int>::getExtendedDbtype(prefdbr.getDbtype()));
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), threads, compressed, dbtype);
    dbw.open();

    const size_t dbSize = prefdbr.getSize();
    if (dbSize == 0) {
        dbw.close();
        prefdbr.close();
        if (!sameQTDB) delete qDbrIdx;
        delete tDbrIdx;
        return EXIT_SUCCESS;
    }

    size_t totalMemory = Util::getTotalSystemMemory();
    size_t flushSize = 1000000;
    if (totalMemory > prefdbr.getTotalDataSize()) {
        flushSize = dbSize;
    }
    size_t iterations = static_cast<size_t>(ceil(static_cast<double>(dbSize) / static_cast<double>(flushSize)));

    size_t alignmentsNum = 0;
    size_t totalPassedNum = 0;

    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(dbSize - (i * flushSize), flushSize);
        Debug::Progress progress(bucketSize);

#pragma omp parallel num_threads(threads)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            std::string alnResultsOutString;
            alnResultsOutString.reserve(1024 * 1024);
            char buffer[1024 + 32768 * 4];

            Sequence qSeq(maxSeqLen, querySeqType, &subMat, 0, false, compBiasCorrection);
            Sequence dbSeq(maxSeqLen, targetSeqType, &subMat, 0, false, compBiasCorrection);
            const size_t maxMatcherSeqLen = std::max(tdbr->getMaxSeqLen(), qdbr->getMaxSeqLen());

            std::vector<RnaMatcher::result_t> swResults;
            swResults.reserve(300);
            RnaMatcher matcher(querySeqType, maxMatcherSeqLen, &subMat, &evaluer,
                               compBiasCorrection, compBiasCorrectionScale,
                               gapOpen, gapExtend, correlationScoreWeight, 0);

            std::vector<RnaMatcher::result_t> swRealignResults;
            RnaMatcher *realigner = NULL;
            if (realign) {
                swRealignResults.reserve(300);
                realigner = &matcher;
                if (realign_m != NULL) {
                    realigner = new RnaMatcher(querySeqType, maxMatcherSeqLen, realign_m, &evaluer,
                                               compBiasCorrection, compBiasCorrectionScale,
                                               gapOpen, gapExtend, 0.0, 0);
                }
            }

            const char *words[10];

#pragma omp for schedule(dynamic, 5) reduction(+: alignmentsNum, totalPassedNum)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();

                char *data, *origData;
                data = origData = prefdbr.getData(id, thread_idx);
                unsigned int queryDbKey = prefdbr.getDbKey(id);
                size_t origQueryLen = 0;
                // Derive reverse from dbKey: even = forward, odd = reverse
                bool reverse = (queryDbKey % 2 == 1);

                if (*data != '\0') {
                    size_t qId = qdbr->getId(queryDbKey);
                    char *querySeqData = qdbr->getData(qId, thread_idx);
                    if (querySeqData == NULL) {
                        Debug(Debug::ERROR) << "Query sequence " << queryDbKey
                                            << " is required in the prefiltering, but is not contained in the query sequence database.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    size_t queryLen = qdbr->getSeqLen(qId);
                    origQueryLen = queryLen;

                    qSeq.mapSequence(qId, queryDbKey, querySeqData, queryLen);
                    matcher.initQuery(&qSeq);
                }

                size_t passedNum = 0;
                unsigned int rejected = 0;
                while (*data != '\0' && passedNum < maxAccept && rejected < maxReject) {
                    Util::parseKey(data, buffer);
                    const unsigned int dbKey = (unsigned int)strtoul(buffer, NULL, 10);
                    size_t elements = Util::getWordsOfLine(data, words, 10);

                    short diagonal = 0;
                    bool isReverse = false;
                    if (elements == 3) {
                        hit_t hit = QueryMatcher::parsePrefilterHit(data);
                        isReverse = reversePrefilterResult && (hit.prefScore < 0);
                        diagonal = static_cast<short>(hit.diagonal);
                    }
                    data = Util::skipLine(data);

                    size_t dbId = tdbr->getId(dbKey);
                    char *dbSeqData = tdbr->getData(dbId, thread_idx);
                    if (dbSeqData == NULL) {
                        Debug(Debug::ERROR) << "Sequence " << dbKey
                                            << " is required in the prefiltering, but is not contained in the target sequence database!\n";
                        EXIT(EXIT_FAILURE);
                    }
                    dbSeq.mapSequence(dbId, dbKey, dbSeqData, tdbr->getSeqLen(dbId));
                    if (reverse) {
                        dinucEncodeReverse(&dbSeq);
                    }

                    if (Util::canBeCovered(covThr, covMode, static_cast<float>(origQueryLen), static_cast<float>(dbSeq.L)) == false) {
                        rejected++;
                        continue;
                    }

                    const bool isIdentity = (queryDbKey == dbKey && (includeIdentity || sameQTDB));

                    RnaMatcher::result_t res = matcher.getSWResult(&dbSeq, static_cast<int>(diagonal), isReverse,
                                                                    covMode, covThr, evalThr, swMode, seqIdMode,
                                                                    isIdentity, false, reverse);
                    alignmentsNum++;

                    if (isIdentity) {
                        res.qcov = 1.0f;
                        res.dbcov = 1.0f;
                        res.seqId = 1.0f;
                    }

                    // Skip results where the aligner returned early without computing start positions
                    // When realign is active, initial pass uses SCORE_ONLY (no start pos),
                    // realignment will recompute full alignment boundaries
                    if (addBacktrace && !realign && res.qStartPos == -1) {
                        rejected++;
                        continue;
                    }

                    if (checkCriteria(res, isIdentity, evalThr, seqIdThr, alnLenThr, covMode, covThr)) {
                        swResults.emplace_back(res);
                        passedNum++;
                        totalPassedNum++;
                        rejected = 0;

                        // Alternative alignments — mask matched TARGET region with X
                        // and re-run the matcher up to altAlignment times to find
                        // additional non-overlapping hits in the same target.
                        // Mirrors lib/mmseqs/src/alignment/Alignment.cpp:590-599.
                        if (altAlignment > 0 && !isIdentity
                            && res.dbStartPos >= 0 && res.dbEndPos > res.dbStartPos) {
                            const unsigned char xIndex = subMat.aa2num[static_cast<int>('X')];
                            // Mask the just-matched region of dbSeq.numSequence
                            const int dbLen = dbSeq.L;
                            int p0 = std::max(0, res.dbStartPos);
                            int p1 = std::min(dbLen, res.dbEndPos);
                            for (int p = p0; p < p1; ++p) {
                                dbSeq.numSequence[p] = xIndex;
                            }
                            for (int altIdx = 0; altIdx < altAlignment; ++altIdx) {
                                RnaMatcher::result_t altRes = matcher.getSWResult(
                                    &dbSeq, INT_MAX, isReverse, covMode, covThr,
                                    evalThr, swMode, seqIdMode, isIdentity, false, reverse);
                                if (!checkCriteria(altRes, isIdentity, evalThr,
                                                   seqIdThr, alnLenThr, covMode, covThr)) {
                                    break;
                                }
                                if (altRes.dbStartPos < 0 || altRes.dbEndPos <= altRes.dbStartPos) {
                                    break;
                                }
                                swResults.emplace_back(altRes);
                                passedNum++;
                                totalPassedNum++;
                                int ap0 = std::max(0, altRes.dbStartPos);
                                int ap1 = std::min(dbLen, altRes.dbEndPos);
                                for (int p = ap0; p < ap1; ++p) {
                                    dbSeq.numSequence[p] = xIndex;
                                }
                            }
                        }
                    } else {
                        rejected++;
                    }
                }

                if (swResults.size() > 1) {
                    SORT_SERIAL(swResults.begin(), swResults.end(), RnaMatcher::compareHits);
                }

                std::vector<RnaMatcher::result_t> *returnRes = &swResults;
                if (realign && *origData != '\0') {
                    // For profile queries, remap profile_for_alignment using consensus + biased matrix
                    if (Parameters::isEqualDbtype(qSeq.getSequenceType(), Parameters::DBTYPE_HMM_PROFILE) && realign_m != NULL) {
                        int8_t *profileForAlignment = const_cast<int8_t*>(qSeq.getAlignmentProfile());
                        for (int pos = 0; pos < qSeq.L; pos++) {
                            unsigned char queryLetter = qSeq.numSequence[pos];
                            for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                                if (reverse) {
                                    profileForAlignment[aa_num * qSeq.L + pos] = realign_m->subMatrix[queryLetter][realign_m->num2revcompnum[aa_num]];
                                } else {
                                    profileForAlignment[aa_num * qSeq.L + pos] = realign_m->subMatrix[queryLetter][aa_num];
                                }
                            }
                        }
                        // Zero out the X row
                        if (realign_m->alphabetSize - Sequence::PROFILE_AA_SIZE != 0) {
                            memset(&profileForAlignment[(realign_m->alphabetSize - 1) * qSeq.L], 0, qSeq.L);
                        }
                    }
                    realigner->initQuery(&qSeq);
                    int realignAccepted = 0;
                    for (size_t result = 0; result < swResults.size() && realignAccepted < realignMaxSeqs; result++) {
                        size_t dbId = tdbr->getId(swResults[result].dbKey);
                        char *dbSeqData = tdbr->getData(dbId, thread_idx);
                        if (dbSeqData == NULL) {
                            Debug(Debug::ERROR) << "Sequence " << swResults[result].dbKey
                                                << " is required in the prefiltering, but is not contained in the target sequence database!\n";
                            EXIT(EXIT_FAILURE);
                        }
                        dbSeq.mapSequence(dbId, swResults[result].dbKey, dbSeqData, tdbr->getSeqLen(dbId));
                        if (reverse) {
                            dinucEncodeReverse(&dbSeq);
                        }

                        const bool isIdentity = (queryDbKey == swResults[result].dbKey && (includeIdentity || sameQTDB));
                        RnaMatcher::result_t res = realigner->getSWResult(&dbSeq, INT_MAX, false, covMode, realignCov, FLT_MAX,
                                                                           realignSwMode, seqIdMode, isIdentity, false, reverse);

                        const bool covOK = Util::hasCoverage(realignCov, covMode, res.qcov, res.dbcov);
                        if (covOK || isIdentity) {
                            res.score = swResults[result].score;
                            res.eval = swResults[result].eval;
                            swRealignResults.emplace_back(res);
                            realignAccepted++;
                        }
                    }

                    if (swRealignResults.size() > 1) {
                        SORT_SERIAL(swRealignResults.begin(), swRealignResults.end(), RnaMatcher::compareHits);
                    }

                    returnRes = &swRealignResults;
                }

                for (size_t result = 0; result < returnRes->size(); result++) {
                    size_t len = RnaMatcher::resultToBuffer(buffer, (*returnRes)[result], addBacktrace);
                    alnResultsOutString.append(buffer, len);
                }

                dbw.writeData(alnResultsOutString.c_str(), alnResultsOutString.length(), queryDbKey, thread_idx);
                alnResultsOutString.clear();
                swResults.clear();
                swRealignResults.clear();
            }

            if (realigner != NULL && realigner != &matcher) {
                delete realigner;
            }

            if (i != (iterations - 1)) {
#pragma omp barrier
                if (thread_idx == 0) {
                    prefdbr.remapData();
                }
#pragma omp barrier
            }
        }
    }

    if (realign_m != NULL) {
        delete realign_m;
    }

    dbw.close();
    prefdbr.close();
    if (!sameQTDB) {
        delete qDbrIdx;
    }
    delete tDbrIdx;

    Debug(Debug::INFO) << alignmentsNum << " alignments calculated\n";
    Debug(Debug::INFO) << totalPassedNum << " sequence pairs passed the thresholds";
    if (alignmentsNum > 0) {
        Debug(Debug::INFO) << " (" << ((float)totalPassedNum / (float)alignmentsNum) << " of overall calculated)";
    }
    Debug(Debug::INFO) << "\n";
    if (dbSize > 0) {
        size_t hits = totalPassedNum / dbSize;
        size_t hits_rest = totalPassedNum % dbSize;
        float hits_f = ((float)hits) + ((float)hits_rest) / (float)dbSize;
        Debug(Debug::INFO) << hits_f << " hits per query sequence\n";
    }

    return EXIT_SUCCESS;
}
