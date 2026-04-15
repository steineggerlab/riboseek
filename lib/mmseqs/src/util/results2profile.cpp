#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "tantan.h"
#include "IndexReader.h"
#include "Masker.h"

#ifdef OPENMP
#include <omp.h>
#endif

int results2profile(int argc, const char **argv, const Command &command, bool returnAlnRes) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    // default for result2profile to filter MSA
    par.filterMsa = 1;
    if (returnAlnRes) {
        par.PARAM_INCLUDE_IDENTITY.description = "keep the query (representative) sequence";
        par.PARAM_INCLUDE_IDENTITY.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_MAX_SEQ_ID.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_QID.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_QSC.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_COV.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_NDIFF.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.parseParameters(argc, argv, command, false, 0, 0);
    par.evalProfile = (par.evalThr < par.evalProfile || returnAlnRes) ? par.evalThr : par.evalProfile;
    par.printParameters(command.cmd, argc, argv, *command.params);

    std::vector<std::string> qid_str_vec = Util::split(par.qid, ",");
    std::vector<int> qid_vec;
    for (size_t qid_idx = 0; qid_idx < qid_str_vec.size(); qid_idx++) {
        float qid_float = strtod(qid_str_vec[qid_idx].c_str(), NULL);
        qid_vec.push_back(static_cast<int>(qid_float*100));
    }
    std::sort(qid_vec.begin(), qid_vec.end());

    DBReader<unsigned int> *tDbr = NULL;
    IndexReader *tDbrIdx = NULL;
    
    // vector of targetdbs and resultdbs
    std::vector<std::string> targetDbPaths, resultDbPaths;
    // Make sure that the length of par.filenames is even
    if (par.filenames.size() % 2 != 0) {
        Debug(Debug::ERROR) << "Internal error: DB paths not provided appropriately, please check the input again\n";
        return EXIT_FAILURE;
    }
    size_t i = 1;
    for (; i < par.filenames.size() / 2; i++) {
        targetDbPaths.push_back(par.filenames[i]);
    }
    for (; i < par.filenames.size() - 1; i++) {
        resultDbPaths.push_back(par.filenames[i]);
    }

    std::vector<DBReader<unsigned int>*> resultDbrs;
    DBReader<unsigned int> *resultDbr = NULL;
    std::vector<std::vector<unsigned int>> resultDbIds;
    resultDbIds.reserve(resultDbPaths.size());

    for (size_t i = 0; i < resultDbPaths.size(); i++) {
        resultDbr = new DBReader<unsigned int>(resultDbPaths[i].c_str(), (resultDbPaths[i] + ".index").c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        resultDbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
        resultDbrs.push_back(resultDbr);
    }
    size_t dbFrom = 0;
    size_t dbSize = 0;
#ifdef HAVE_MPI
    resultDbrs[0].decomposeDomainByAminoAcid(MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom + dbSize << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.filenames[par.filenames.size() - 1], par.filenames[par.filenames.size() - 1] + ".index", MMseqsMPI::rank);
#else
    dbSize = resultDbrs[0]->getSize(); // getSize should be same for all result DBs
    std::pair<std::string, std::string> tmpOutput = std::make_pair(par.filenames[par.filenames.size() - 1], par.filenames[par.filenames.size() - 1] + ".index");
#endif

    // reserve dbSize in each of the resultDbIds
    for (size_t i = 0; i < resultDbrs.size(); i++) {
        resultDbIds.push_back(std::vector<unsigned int>());
        resultDbIds[i].reserve(dbSize);
        for (size_t j = 0; j < dbSize; j++) {
            resultDbIds[i].push_back(UINT_MAX);
        }
    }

    // save the resultDb ids to the queryKeys
    for (size_t i = 0; i < resultDbrs.size(); i++) {
        DBReader<unsigned int> *currentDbr = resultDbrs[i];
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            unsigned int queryKey = currentDbr->getDbKey(id);
            resultDbIds[i][queryKey] = id;
        }
    }

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, dbSize), (size_t)1);
#endif

    // vector of DBReader
    std::vector<DBReader<unsigned int>*> tDbrs;
    std::vector<IndexReader*> tDbrIdxs;
    bool templateDBIsIndex = false;
    bool needSrcIndex = false;
    int targetSeqType = -1;
    int targetDbtype = FileUtil::parseDbType(par.db2.c_str());
    if (Parameters::isEqualDbtype(targetDbtype, Parameters::DBTYPE_INDEX_DB)) {
        uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.filenames[par.filenames.size() - 2].c_str()));
        needSrcIndex = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
        bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
        tDbrIdx = new IndexReader(par.db2, par.threads,
                                  needSrcIndex ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                                  (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        tDbr = tDbrIdx->sequenceReader;
        templateDBIsIndex = true;
        targetSeqType = tDbr->getDbtype();
    }

    if (templateDBIsIndex == false) {
        // Iterate through targetDbPaths
        for (size_t i = 0; i < targetDbPaths.size(); i++) {
            // tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            tDbr = new DBReader<unsigned int>(targetDbPaths[i].c_str(), (targetDbPaths[i] + ".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            tDbr->open(DBReader<unsigned int>::NOSORT);
            tDbrs.push_back(tDbr);
        }
        targetSeqType = tDbr->getDbtype();
    }

    DBReader<unsigned int> *qDbr = NULL;
    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;
    if (!sameDatabase) {
        qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        qDbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            qDbr->readMmapedDataInMemory();
        }
    } else {
        qDbr = tDbr;
    }
    unsigned int tDbrMax = 0;
    for (size_t i = 0; i < tDbrs.size(); i++) {
        if (tDbrs[i]->getMaxSeqLen() > tDbrMax) {
            tDbrMax = tDbrs[i]->getMaxSeqLen();
        }
    }
    const unsigned int maxSequenceLength = std::max(tDbrMax, qDbr->getMaxSeqLen());

    // qDbr->readMmapedDataInMemory();
    // make sure to touch target after query, so if there is not enough memory for the query, at least the targets
    // might have had enough space left to be residung in the page cache
    if (sameDatabase == false && templateDBIsIndex == false && par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        // tDbr->readMmapedDataInMemory();
        for (size_t i = 0; i < tDbrs.size(); i++) {
            tDbrs[i]->readMmapedDataInMemory();
        }
    }

    int type = Parameters::DBTYPE_HMM_PROFILE;
    const int writePlain = par.profileOutputMode == 1;
    if (par.profileOutputMode == 1) {
        type = Parameters::DBTYPE_OMIT_FILE;
        par.compressed = false;
    } else if (returnAlnRes) {
        type = Parameters::DBTYPE_ALIGNMENT_RES;
        if (needSrcIndex) {
            type = DBReader<unsigned int>::setExtendedDbtype(type, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
        }
    } else if (par.pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC) {
        type = DBReader<unsigned int>::setExtendedDbtype(type, Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS);
    }
    DBWriter resultWriter(tmpOutput.first.c_str(), tmpOutput.second.c_str(), localThreads, par.compressed, type);
    resultWriter.open();

    // + 1 for query
    size_t maxSetSize = 0;
    for (size_t i = 0; i < tDbrs.size(); i++) {
        maxSetSize += resultDbrs[i]->maxCount('\n') + 1;
    }
    // size_t maxSetSize = resultReader.maxCount('\n') + 1;

    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);
    size_t autoDbSize = 0;
    for (size_t i = 0; i < tDbrs.size(); i++) {
        autoDbSize += tDbrs[i]->getAminoAcidDBSize();
    }
    size_t effectiveDbSize = par.dbSize == 0 ? autoDbSize : par.dbSize;
    if (par.strand == 2) {
        effectiveDbSize *= 2;
    }
    EvalueComputation evalueComputation(effectiveDbSize, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());

    if (qDbr->getDbtype() == -1 || targetSeqType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(qDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database\n";
        return EXIT_FAILURE;
    }

    Debug(Debug::INFO) << "Query database size: " << qDbr->getSize() << " type: " << qDbr->getDbTypeName() << "\n";
    size_t dbSizeSum = 0;
    for (size_t i = 0; i < tDbrs.size(); i++) {
        Debug(Debug::INFO) << "Target database " << i << " size: " << tDbrs[i]->getSize() << " type: " << tDbrs[i]->getDbTypeName() << "\n";
        dbSizeSum += tDbrs[i]->getSize();
    }
    Debug(Debug::INFO) << "Total target database size: " << dbSizeSum << " type: " << Parameters::getDbTypeName(targetSeqType) << "\n";

    const bool isFiltering = par.filterMsa != 0 || returnAlnRes;
    Debug::Progress progress(dbSize - dbFrom);
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        Matcher matcher(qDbr->getDbtype(), maxSequenceLength, &subMat, &evalueComputation,
                        par.compBiasCorrection, par.compBiasCorrectionScale,
                        par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 0.0, par.zdrop);
        Masker masker(subMat);
        MultipleAlignment aligner(maxSequenceLength, &subMat);
        PSSMCalculator calculator(
            &subMat, maxSequenceLength, maxSetSize, par.pcmode, par.pca, par.pcb
#ifdef GAP_POS_SCORING
            , par.gapOpen.values.aminoacid()
            , par.gapPseudoCount
#endif
        );
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
        Sequence centerSequence(maxSequenceLength, qDbr->getDbtype(), &subMat, 0, false, par.compBiasCorrection);
        Sequence edgeSequence(maxSequenceLength, targetSeqType, &subMat, 0, false, false);

        char dbKey[255];
        const char *entry[255];
        char buffer[1024 + 32768*4];
        float * pNullBuffer = new float[maxSequenceLength + 1];

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::vector<std::vector<unsigned char>> seqSet;
        seqSet.reserve(300);

        DBReader<unsigned int> *localResultReader = NULL;
        DBReader<unsigned int> *localTDbr = NULL;

        std::string result;
        result.reserve((maxSequenceLength + 1) * Sequence::PROFILE_READIN_SIZE);

#pragma omp for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();

            unsigned int queryKey = resultDbrs[0]->getDbKey(id);
            size_t queryId = qDbr->getId(queryKey);
            if (queryId == UINT_MAX) {
                Debug(Debug::WARNING) << "Invalid query sequence " << queryKey << "\n";
                continue;
            }
            centerSequence.mapSequence(queryId, queryKey, qDbr->getData(queryId, thread_idx), qDbr->getSeqLen(queryId));

            // Iterate through the target and result DBs
            for (size_t i = 0; i < tDbrs.size(); i++) {
                localTDbr = tDbrs[i];
                localResultReader = resultDbrs[i];
                bool isQueryInit = false;
                char *data = localResultReader->getData(resultDbIds[i][queryKey], thread_idx);
                while (*data != '\0') {
                    Util::parseKey(data, dbKey);
                    const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                    // in the same database case, we have the query repeated
                    if (key == queryKey && sameDatabase == true) {
                        if(returnAlnRes && par.includeIdentity){
                            Matcher::result_t res = Matcher::parseAlignmentRecord(data);
                            size_t len = Matcher::resultToBuffer(buffer, res, true);
                            result.append(buffer, len);
                        }

                        data = Util::skipLine(data);
                        continue;
                    }

                    const size_t columns = Util::getWordsOfLine(data, entry, 255);
                    float evalue = 0.0;
                    if (returnAlnRes == false && columns >= 4) {
                        evalue = strtod(entry[3], NULL);
                    }

                    if (returnAlnRes == true || evalue < par.evalProfile) {
                        const size_t edgeId = localTDbr->getId(key);
                        if (edgeId == UINT_MAX) {
                            Debug(Debug::ERROR) << "Sequence " << key << " does not exist in target sequence database\n";
                            EXIT(EXIT_FAILURE);
                        }
                        
                        if (columns > Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                            alnResults.emplace_back(Matcher::parseAlignmentRecord(data));
                        } else {
                            // Recompute if not all the backtraces are present
                            if (isQueryInit == false) {
                                matcher.initQuery(&centerSequence);
                                isQueryInit = true;
                            }
                            alnResults.emplace_back(matcher.getSWResult(&edgeSequence, INT_MAX, false, 0, 0.0, FLT_MAX, Matcher::SCORE_COV_SEQID, 0, false));
                        }
                        // Check if it is on the reverse strand or not
                        bool reverse = false;
                        if (alnResults.back().qStartPos > alnResults.back().qEndPos) {
                            reverse = true;
                        }
                        edgeSequence.mapSequence(edgeId, key, localTDbr->getData(edgeId, thread_idx), localTDbr->getSeqLen(edgeId), false, NULL, reverse, false, 1.0f, localTDbr->isPadded());
                        seqSet.emplace_back(std::vector<unsigned char>(edgeSequence.numSequence, edgeSequence.numSequence + edgeSequence.L));
                    }
                    data = Util::skipLine(data);
                }
            }

            // Recompute if not all the backtraces are present
            MultipleAlignment::MSAResult res = aligner.computeMSA(&centerSequence, seqSet, alnResults, true);

            // do not count query
            size_t filteredSetSize = (isFiltering == true)  ?
                                     filter.filter(res, alnResults, (int)(par.covMSAThr * 100), qid_vec, par.qsc, (int)(par.filterMaxSeqId * 100), par.Ndiff, par.filterMinEnable)
                                     :
                                     res.setSize;
             //MultipleAlignment::print(res, &subMat);

            if (returnAlnRes) {
                for (size_t i = 0; i < (filteredSetSize - 1); ++i) {
                    size_t len = Matcher::resultToBuffer(buffer, alnResults[i], true);
                    result.append(buffer, len);
                }
            } else {
                for (size_t pos = 0; pos < res.centerLength; pos++) {
                    if (res.msaSequence[0][pos] == MultipleAlignment::GAP) {
                        Debug(Debug::ERROR) << "Error in computePSSMFromMSA. First sequence of MSA is not allowed to contain gaps.\n";
                        EXIT(EXIT_FAILURE);
                    }
                }

                PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(filteredSetSize, res.centerLength,
                                                                                (const char **) res.msaSequence,
#ifdef GAP_POS_SCORING
                                                                                alnResults,
#endif
                                                                                par.wg, 0.0);
                if (writePlain) {
                    result.clear();
                    result.append("Query profile of sequence ");
                    result.append(SSTR(queryKey));
                    result.push_back('\n');
                    calculator.profileToString(result, res.centerLength);
                } else {                                                                
                    if (par.compBiasCorrection == true){
                        SubstitutionMatrix::calcGlobalAaBiasCorrection(&subMat, pssmRes.pssm, pNullBuffer,
                                                                    Sequence::PROFILE_AA_SIZE,
                                                                    res.centerLength);
                    }
                    if (par.maskProfile == true) {
                        masker.maskPssm(centerSequence, par.maskProb, pssmRes);
                    }
                    pssmRes.toBuffer(centerSequence, subMat, result);
                } 
            }
            resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx, writePlain == false);
            result.clear();
            alnResults.clear();

            MultipleAlignment::deleteMSA(&res);
            seqSet.clear();
        }
        delete[] pNullBuffer;
    }
    resultWriter.close(returnAlnRes == false || writePlain == true);
    if (writePlain) {
        FileUtil::remove(par.db4Index.c_str());
    }
    // resultReader.close();
    for (size_t i = 0; i < resultDbrs.size(); i++) {
        resultDbrs[i]->close();
        delete resultDbrs[i];
    }

    if (!sameDatabase) {
        qDbr->close();
        delete qDbr;
    }
    if (tDbrIdx == NULL) {
        for (size_t i = 0; i < tDbrs.size(); i++) {
            tDbrs[i]->close();
            delete tDbrs[i];
        }
        // tDbr->close();
        // delete tDbr;
    } else {
        delete tDbrIdx;
    }

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // master reduces results
    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int procs = 0; procs < MMseqsMPI::numProc; procs++) {
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(par.filenames[par.filenames.size() - 1], par.filenames[par.filenames.size() - 1] + ".index", procs);
            splitFiles.push_back(std::make_pair(tmpFile.first, tmpFile.second));

        }
        DBWriter::mergeResults(par.filenames[par.filenames.size() - 1], par.filenames[par.filenames.size() - 1] + ".index", splitFiles);
    }
#endif

    if (MMseqsMPI::isMaster() && returnAlnRes == false) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.filenames[par.filenames.size() - 1], DBFiles::SEQUENCE_ANCILLARY);
    }

    return EXIT_SUCCESS;
}

int results2profile(int argc, const char **argv, const Command &command) {
    return results2profile(argc, argv, command, false);
}

int filterresults(int argc, const char **argv, const Command &command) {
    return results2profile(argc, argv, command, true);
}

