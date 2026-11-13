//
// Created by Martin Steinegger on 17.09.18.
//

#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
#include "NucleotideMatrix.h"
#include "FastSort.h"
#include "SubstitutionMatrixProfileStates.h"
#include "IndexReader.h"
#include "QueryMatcherTaxonomyHook.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <chrono>
#include <thread>

#ifdef OPENMP
#include <omp.h>
#endif
// #define HAVE_CUDA 1
#ifdef HAVE_CUDA
#include "GpuUtil.h"
#include "Alignment.h"
#include <signal.h>
#include <algorithm>
#include <cuda_runtime.h>
#endif

#ifdef HAVE_CUDA

volatile sig_atomic_t keepRunningClient = 1;
void intHandlerClient(int) {
    keepRunningClient = 0;
}

void runFilterOnGpu(Parameters & par, BaseMatrix * subMat,
                    DBReader<unsigned int> * qdbr, DBReader<unsigned int> * tdbr,
                    bool sameDB, DBWriter & resultWriter, EvalueComputation * evaluer,
                    QueryMatcherTaxonomyHook *taxonomyHook){
    Debug::Progress progress(qdbr->getSize());
    const int querySeqType = qdbr->getDbtype();
    Sequence qSeq(par.maxSeqLen, querySeqType, subMat, 0, false, par.compBiasCorrection);

    std::vector<Marv::Result> results;
    results.reserve(par.maxResListLen);
    std::vector<hit_t> shortResults;
    std::vector<Matcher::result_t> resultsAln;

    size_t profileBufferLength = par.maxSeqLen;
    int8_t* profile = NULL;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) == false) {
        profile = (int8_t*)malloc(subMat->alphabetSize * profileBufferLength * sizeof(int8_t));
    }

    std::string resultBuffer;
    resultBuffer.reserve(262144);
    char buffer[1024+32768];

    size_t compBufferSize = (par.maxSeqLen + 1) * sizeof(float);
    float *compositionBias = NULL;
    if (par.compBiasCorrection == true) {
        compositionBias = (float*)malloc(compBufferSize);
        memset(compositionBias, 0, compBufferSize);
    }

    std::string hash = "";
    if (par.gpuServer != 0) {
        hash = GPUSharedMemory::getShmHash(par.db2, par.gpuServerNamespace);
        std::string path = "/dev/shm/" + hash;
        // Debug(Debug::WARNING) << path << "\n";
        int waitTimeout = par.gpuServerWaitTimeout;
        std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
        bool statusPrinted = false;
        while (true) {
            size_t shmSize = FileUtil::getFileSize(path);
            // server is ready once the shm file exists and is not 0 byte large
            if (shmSize != (size_t)-1 && shmSize > 0) {
                break;
            }

            if (waitTimeout == 0) {
                Debug(Debug::ERROR) 
                    << "gpuserver for database " << par.db2 << " not found.\n"
                    << "Please start gpuserver with the same CUDA_VISIBLE_DEVICES\n";
                EXIT(EXIT_FAILURE);
            }

            if (waitTimeout > 0) {
                if (statusPrinted == false) {
                    Debug(Debug::WARNING) << "Waiting for `gpuserver`\n";
                    statusPrinted = true;
                } else {
                    Debug(Debug::WARNING) << ".";
                }
                std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();
                if (elapsed >= waitTimeout) {
                    Debug(Debug::ERROR)
                        << "\ngpuserver for database " << par.db2 << " not found after " << elapsed <<  "seconds.\n"
                        << "Please start gpuserver with the same CUDA_VISIBLE_DEVICES\n";
                    EXIT(EXIT_FAILURE);
                }
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
        if (waitTimeout > 0 && statusPrinted) {
            Debug(Debug::INFO) << "\n";
        }
    }

    size_t* offsetData = NULL;
    int32_t* lengthData = NULL;
    std::vector<size_t> offsets;
    std::vector<int32_t> lengths;
    GPUSharedMemory* layout = NULL;
    if (hash.empty()) {
        offsets.reserve(tdbr->getSize() + 1);
        lengths.reserve(tdbr->getSize());
        for (size_t id = 0; id < tdbr->getSize(); id++) {
            offsets.emplace_back(tdbr->getIndex()[id].offset);
            lengths.emplace_back(tdbr->getIndex()[id].length - 2);
        }
        offsets.emplace_back(offsets.back() + lengths.back());
        offsetData = offsets.data();
        lengthData = lengths.data();
    } else {
        layout = GPUSharedMemory::openSharedMemory(hash);
    }

    const bool serverMode = par.gpuServer;
    Marv* marv = NULL;

    // Database chunking variables for large databases
    std::vector<std::pair<size_t, size_t>> dbChunks; // (startIdx, endIdx) for each chunk
    std::vector<void*> dbHandles;
    size_t numChunks = 1;

    if (serverMode == 0) {
       if (offsetData == NULL || lengthData == NULL) {
           Debug(Debug::ERROR) << "Invalid GPU database\n";
           EXIT(EXIT_FAILURE);
       }
        int32_t maxTargetLength = lengths.back();
        Marv::AlignmentType type = (par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED) ?
                Marv::AlignmentType::GAPLESS_SMITH_WATERMAN : Marv::AlignmentType::GAPLESS;
        // marv = new Marv(tdbr->getSize(), subMat->alphabetSize, maxTargetLength,
        //                 par.maxResListLen, type);
        // void* h = marv->loadDb(
        //     tdbr->getDataForFile(0), offsetData, lengthData, tdbr->getDataSizeForFile(0)
        // );
        // marv->setDb(h);
        // Get total available GPU memory across all devices
        std::vector<int> deviceIds = Marv::getDeviceIds();
        size_t totalGpuMem = 0;
        size_t totalFreeMem = 0;
        for (int deviceId : deviceIds) {
            cudaSetDevice(deviceId);
            size_t freeMem = 0, totalMem = 0;
            cudaMemGetInfo(&freeMem, &totalMem);
            totalGpuMem += totalMem;
            totalFreeMem += freeMem;
        }
        cudaSetDevice(deviceIds[0]); // Reset to first device
        
        // Estimate database size (sequence data + overhead for GPU structures)
        size_t dbDataSize = tdbr->getDataSizeForFile(0);
        // GPU overhead: ~1x for database data
        size_t estimatedGpuUsage = dbDataSize;
        
        // Use 80% of available free memory as usable, leave headroom for kernels and temp storage
        size_t usableGpuMem = (size_t)(totalFreeMem * 0.8);
        
        if (estimatedGpuUsage > usableGpuMem) {
            // Database is too large, need to split into chunks
            numChunks = (estimatedGpuUsage + usableGpuMem - 1) / usableGpuMem;
            size_t seqsPerChunk = (tdbr->getSize() + numChunks - 1) / numChunks;
            
            Debug(Debug::INFO) << "Database size (" << (dbDataSize / (1024*1024)) << " MB) exceeds available GPU memory ("
                               << (usableGpuMem / (1024*1024)) << " MB). Splitting into " << numChunks << " chunks.\n";
            
            for (size_t chunk = 0; chunk < numChunks; chunk++) {
                size_t startIdx = 0;
                size_t domainSize = 0;
                tdbr->decomposeDomainByAminoAcid(chunk, numChunks, &startIdx, &domainSize);
                if (startIdx < tdbr->getSize()) {
                    dbChunks.push_back({startIdx, startIdx + domainSize});
                }
            }
            numChunks = dbChunks.size();
        } else {
            // Database fits in GPU memory, single chunk
            dbChunks.push_back({0, tdbr->getSize()});
        }
        
        // Load database
        if (numChunks == 1) {
            // Single chunk - create Marv and load entire database
            marv = new Marv(tdbr->getSize(), subMat->alphabetSize, maxTargetLength,
                            par.maxResListLen, type);
            void* h = marv->loadDb(
                tdbr->getDataForFile(0), offsetData, lengthData, tdbr->getDataSizeForFile(0)
            );
            marv->setDb(h);
            dbHandles.push_back(h);
        } else {
            // Multiple chunks - Marv instances will be created on-demand during the query loop
            Debug(Debug::INFO) << "Using chunked database loading with " << numChunks << " chunks\n";
        }
    } else if (layout == NULL) {
       Debug(Debug::ERROR) << "No GPU server shared memory connection\n";
       EXIT(EXIT_FAILURE);
    } else {
        struct sigaction act;
        // Set up the handler for SIGINT and SIGTERM
        memset(&act, 0, sizeof(act));
        act.sa_handler = intHandlerClient;
        sigaction(SIGINT, &act, NULL);
        sigaction(SIGTERM, &act, NULL);
        dbChunks.push_back({0, tdbr->getSize()}); // Server mode uses single "chunk"
    }

    // For multi-chunk mode: store accumulated results per query across all chunks
    // Map from query id to vector of results from all chunks
    std::vector<std::vector<Marv::Result>> allQueryResults;
    if (numChunks > 1) {
        allQueryResults.resize(qdbr->getSize());
        for (size_t i = 0; i < qdbr->getSize(); i++) {
            allQueryResults[i].reserve(numChunks * par.maxResListLen);
        }
    }

    // Outer loop: iterate through database chunks
    for (size_t chunkIdx = 0; chunkIdx < numChunks; chunkIdx++) {
        if (!keepRunningClient) {
            break;
        }

        Marv* currentMarv = marv;
        size_t chunkStart = dbChunks[chunkIdx].first;
        size_t chunkEnd = dbChunks[chunkIdx].second;
        size_t chunkSize = chunkEnd - chunkStart;
        
        // For multi-chunk mode, create and load chunk-specific Marv instance
        std::vector<size_t> chunkOffsets;
        std::vector<int32_t> chunkLengths;
        if (numChunks > 1) {
            Debug(Debug::INFO) << "Processing chunk " << (chunkIdx + 1) << "/" << numChunks 
                               << " (sequences " << chunkStart << "-" << chunkEnd << ")\n";
            
            // Prepare chunk-specific offsets and lengths
            chunkOffsets.resize(chunkSize + 1);
            chunkLengths.resize(chunkSize);
            size_t baseOffset = offsets[chunkStart];
            for (size_t i = 0; i < chunkSize; i++) {
                chunkOffsets[i] = offsets[chunkStart + i] - baseOffset;
                chunkLengths[i] = lengths[chunkStart + i];
            }
            chunkOffsets[chunkSize] = offsets[chunkEnd] - baseOffset;
            
            // Calculate chunk data size and find max length for this chunk
            size_t chunkDataSize = chunkOffsets[chunkSize];
            int32_t chunkMaxLength = *std::max_element(chunkLengths.begin(), chunkLengths.end());
            
            // Create Marv instance for this chunk
            Marv::AlignmentType type = (par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED) ?
                    Marv::AlignmentType::GAPLESS_SMITH_WATERMAN : Marv::AlignmentType::GAPLESS;
            currentMarv = new Marv(chunkSize, subMat->alphabetSize, chunkMaxLength,
                                   par.maxResListLen, type);
            
            // Load chunk database
            void* chunkHandle = currentMarv->loadDb(
                tdbr->getDataForFile(0) + baseOffset,
                chunkOffsets.data(),
                chunkLengths.data(),
                chunkDataSize
            );
            currentMarv->setDb(chunkHandle);
        }

        // Inner loop: iterate through all queries for this chunk
        for (size_t id = 0; id < qdbr->getSize(); id++) {
            if (!keepRunningClient) {
                break;
            }
            size_t queryKey = qdbr->getDbKey(id);
            unsigned int querySeqLen = qdbr->getSeqLen(id);
            char *querySeqData = qdbr->getData(id, 0);
            qSeq.mapSequence(id, queryKey, querySeqData, querySeqLen);
            if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE)) {
                profile = qSeq.profile_for_alignment;
            } else {
                if ((size_t)qSeq.L >= profileBufferLength) {
                    profileBufferLength = (size_t)qSeq.L * 1.5;
                    profile = (int8_t*)realloc(profile, subMat->alphabetSize * profileBufferLength * sizeof(int8_t));
                }
                if (compositionBias != NULL) {
                    if ((size_t)qSeq.L >= compBufferSize) {
                        compBufferSize = (size_t)qSeq.L * 1.5 * sizeof(float);
                        compositionBias = (float*)realloc(compositionBias, compBufferSize);
                        // memset(compositionBias, 0, compBufferSize);
                    }
                    SubstitutionMatrix::calcLocalAaBiasCorrection(subMat, qSeq.numSequence, qSeq.L, compositionBias, par.compBiasCorrectionScale);
                }
                for (size_t j = 0; j < (size_t)subMat->alphabetSize; ++j) {
                    for (size_t i = 0; i < (size_t)qSeq.L; ++i) {
                        short bias = 0;
                        if (compositionBias != NULL) {
                            bias = static_cast<short>((compositionBias[i] < 0.0) ? (compositionBias[i] - 0.5) : (compositionBias[i] + 0.5));
                        }
                        profile[j * qSeq.L  + i] = subMat->subMatrix[j][qSeq.numSequence[i]] + bias;
                    }
                }
            }
            Marv::Stats stats;
            stats.results = 0;
            if (serverMode == 0) {
                stats = currentMarv->scan(reinterpret_cast<const char *>(qSeq.numSequence), qSeq.L, profile, results.data());
                if (numChunks > 1) {
                    for (size_t i = 0; i < stats.results; i++) {
                        Marv::Result res = results[i];
                        res.id = res.id + chunkStart; // Convert local ID to global ID
                        allQueryResults[id].push_back(res);
                    }
                }
            } else {
                bool claimed = false;
                while (!claimed) {
                    if (layout->serverExit.load(std::memory_order_acquire) == true) {
                        // server has shut down
                        Debug(Debug::ERROR) << "GPU server has unexpectedly shut down\n";
                        EXIT(EXIT_FAILURE);
                    }
                    if (keepRunningClient == false) {
                        EXIT(EXIT_FAILURE);
                    }

                    int expected = GPUSharedMemory::IDLE;
                    int desired = GPUSharedMemory::RESERVED;
                    if (layout->state.compare_exchange_strong(expected, desired, std::memory_order_acq_rel)) {
                        // Debug(Debug::ERROR) << "switch to reserved\n";
                        claimed = true;
                        memcpy(layout->getQueryPtr(), qSeq.numSequence, qSeq.L);
                        memcpy(layout->getProfilePtr(), profile, subMat->alphabetSize * qSeq.L);
                        layout->queryLen = qSeq.L;
                        std::atomic_thread_fence(std::memory_order_release);
                        // Debug(Debug::ERROR) << "switch to ready\n";
                        layout->state.store(GPUSharedMemory::READY, std::memory_order_release);

                        while (true) {
                            if (layout->serverExit.load(std::memory_order_acquire) == true) {
                                Debug(Debug::ERROR) << "GPU server has unexpectedly shut down\n";
                                EXIT(EXIT_FAILURE);
                            }

                            if (layout->state.load(std::memory_order_acquire) == GPUSharedMemory::DONE) {
                                break;
                            } else {
                                std::this_thread::yield();
                            }
                        }

                        std::atomic_thread_fence(std::memory_order_acquire);
                        memcpy(results.data(), layout->getResultsPtr(), layout->resultLen * sizeof(Marv::Result));
                        stats.results = layout->resultLen;
                        // Debug(Debug::ERROR) << "switch to idle\n";
                        layout->state.store(GPUSharedMemory::IDLE, std::memory_order_release);
                        if (keepRunningClient == false) {
                            EXIT(EXIT_FAILURE);
                        }
                    } else {
                        std::this_thread::yield();
                    }
                }
            }

            // For single chunk or server mode: process results immediately
            if (numChunks == 1 || serverMode != 0) {
                if (keepRunningClient == false) {
                    EXIT(EXIT_FAILURE);
                }

                for(size_t i = 0; i < stats.results; i++){
                    unsigned int targetKey = tdbr->getDbKey(results[i].id);
                    int score = results[i].score;
                    if(taxonomyHook != NULL){
                        TaxID currTax = taxonomyHook->taxonomyMapping->lookup(targetKey);
                        if (taxonomyHook->expression[0]->isAncestor(currTax) == false) {
                            continue;
                        }
                    }
                    bool hasDiagScore = (score > par.minDiagScoreThr);
                    const bool isIdentity = (queryKey == targetKey && (par.includeIdentity || sameDB))? true : false;
                    
                    if (isIdentity || hasDiagScore) {
                        if(par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED){
                            Matcher::result_t res;
                            res.dbKey = targetKey;
                            res.eval = evaluer->computeEvalue(score, qSeq.L);
                            res.dbEndPos = results[i].dbEndPos;
                            res.dbLen = tdbr->getSeqLen(results[i].id);
                            res.qEndPos =  results[i].qEndPos;
                            res.qLen = qSeq.L;
                            unsigned int qAlnLen = std::max(static_cast<unsigned int>(res.qEndPos), static_cast<unsigned int>(1));
                            unsigned int dbAlnLen = std::max(static_cast<unsigned int>(res.dbEndPos), static_cast<unsigned int>(1));
                            res.seqId = Matcher::estimateSeqIdByScorePerCol(score, qAlnLen, dbAlnLen);
                            res.qcov = SmithWaterman::computeCov(0, res.qEndPos, res.qLen );
                            res.dbcov = SmithWaterman::computeCov(0, res.dbEndPos, res.dbLen );
                            res.score = evaluer->computeBitScore(score);
                            if(Alignment::checkCriteria(res, isIdentity, par.evalThr,  par.seqIdThr,  par.alnLenThr,  par.covMode,  par.covThr)){
                                resultsAln.emplace_back(res);
                            }
                        } else {
                            hit_t hit;
                            hit.seqId = targetKey;
                            hit.prefScore = score;
                            hit.diagonal = 0;
                            shortResults.emplace_back(hit);
                        }
                    }
                }
                if(par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED) {
                    SORT_PARALLEL(resultsAln.begin(), resultsAln.end(), Matcher::compareHits);
                    size_t maxSeqs = std::min(par.maxResListLen, resultsAln.size());
                    for (size_t i = 0; i < maxSeqs; ++i) {
                        size_t len = Matcher::resultToBuffer(buffer, resultsAln[i], false);
                        resultBuffer.append(buffer, len);
                    }
                }else{
                    SORT_PARALLEL(shortResults.begin(), shortResults.end(), hit_t::compareHitsByScoreAndId);
                    size_t maxSeqs = std::min(par.maxResListLen, shortResults.size());
                    for (size_t i = 0; i < maxSeqs; ++i) {
                        size_t len = QueryMatcher::prefilterHitToBuffer(buffer, shortResults[i]);
                        resultBuffer.append(buffer, len);
                    }


                }

                resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, 0);
                resultBuffer.clear();
                shortResults.clear();
                resultsAln.clear();
                progress.updateProgress();
            }
        } // end query loop

        // Cleanup chunk Marv instance (only for multi-chunk mode)
        if (numChunks > 1 && currentMarv != NULL) {
            delete currentMarv;
        }
    } // end chunk loop
    
    // For multi-chunk mode: merge results from all chunks and write output
    if (numChunks > 1 && serverMode == 0) {
        Debug(Debug::INFO) << "Merging results from all chunks...\n";
        for (size_t id = 0; id < qdbr->getSize(); id++) {
            if (!keepRunningClient) {
                break;
            }
            size_t queryKey = qdbr->getDbKey(id);
            unsigned int querySeqLen = qdbr->getSeqLen(id);
            
            // Sort accumulated results by score (descending)
            std::vector<Marv::Result>& queryResults = allQueryResults[id];
            std::sort(queryResults.begin(), queryResults.end(), 
                [](const Marv::Result& a, const Marv::Result& b) {
                    return a.score > b.score;
                });
            
            // Keep top N results
            size_t numResults = std::min(queryResults.size(), (size_t)par.maxResListLen);
            
            for(size_t i = 0; i < numResults; i++){
                unsigned int targetKey = tdbr->getDbKey(queryResults[i].id);
                int score = queryResults[i].score;
                if(taxonomyHook != NULL){
                    TaxID currTax = taxonomyHook->taxonomyMapping->lookup(targetKey);
                    if (taxonomyHook->expression[0]->isAncestor(currTax) == false) {
                        continue;
                    }
                }
                bool hasDiagScore = (score > par.minDiagScoreThr);
                const bool isIdentity = (queryKey == targetKey && (par.includeIdentity || sameDB))? true : false;
                
                if (isIdentity || hasDiagScore) {
                    if(par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED){
                        Matcher::result_t res;
                        res.dbKey = targetKey;
                        res.eval = evaluer->computeEvalue(score, querySeqLen);
                        res.dbEndPos = queryResults[i].dbEndPos;
                        res.dbLen = tdbr->getSeqLen(queryResults[i].id);
                        res.qEndPos = queryResults[i].qEndPos;
                        res.qLen = querySeqLen;
                        unsigned int qAlnLen = std::max(static_cast<unsigned int>(res.qEndPos), static_cast<unsigned int>(1));
                        unsigned int dbAlnLen = std::max(static_cast<unsigned int>(res.dbEndPos), static_cast<unsigned int>(1));
                        res.seqId = Matcher::estimateSeqIdByScorePerCol(score, qAlnLen, dbAlnLen);
                        res.qcov = SmithWaterman::computeCov(0, res.qEndPos, res.qLen );
                        res.dbcov = SmithWaterman::computeCov(0, res.dbEndPos, res.dbLen );
                        res.score = evaluer->computeBitScore(score);
                        if(Alignment::checkCriteria(res, isIdentity, par.evalThr,  par.seqIdThr,  par.alnLenThr,  par.covMode,  par.covThr)){
                            resultsAln.emplace_back(res);
                        }
                    } else {
                        hit_t hit;
                        hit.seqId = targetKey;
                        hit.prefScore = score;
                        hit.diagonal = 0;
                        shortResults.emplace_back(hit);
                    }
                }
            }
            if(par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED) {
                SORT_PARALLEL(resultsAln.begin(), resultsAln.end(), Matcher::compareHits);
                size_t maxSeqs = std::min(par.maxResListLen, resultsAln.size());
                for (size_t i = 0; i < maxSeqs; ++i) {
                    size_t len = Matcher::resultToBuffer(buffer, resultsAln[i], false);
                    resultBuffer.append(buffer, len);
                }
            }else{
                SORT_PARALLEL(shortResults.begin(), shortResults.end(), hit_t::compareHitsByScoreAndId);
                size_t maxSeqs = std::min(par.maxResListLen, shortResults.size());
                for (size_t i = 0; i < maxSeqs; ++i) {
                    size_t len = QueryMatcher::prefilterHitToBuffer(buffer, shortResults[i]);
                    resultBuffer.append(buffer, len);
                }
            }

            resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, 0);
            resultBuffer.clear();
            shortResults.clear();
            resultsAln.clear();
            progress.updateProgress();
        }
    }

    if (marv != NULL) {
        delete marv;
    } else {
        GPUSharedMemory::unmap(layout);
    }

    if (compositionBias != NULL) {
        free(compositionBias);
    }
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) == false) {
        free(profile);
    }
}
#endif

void runFilterOnCpu(Parameters & par, BaseMatrix * subMat, int8_t * tinySubMat,
                    DBReader<unsigned int> * qdbr, DBReader<unsigned int> * tdbr,
                    SequenceLookup * sequenceLookup, bool sameDB, DBWriter & resultWriter, EvalueComputation * evaluer,
                    QueryMatcherTaxonomyHook *taxonomyHook, int alignmentMode){
    std::vector<hit_t> shortResults;
    shortResults.reserve(tdbr->getSize()/2);
    Debug::Progress progress(qdbr->getSize());
    const int targetSeqType = tdbr->getDbtype();
    const int querySeqType = qdbr->getDbtype();
#ifdef OPENMP
    omp_set_nested(1);
#endif

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        char buffer[1024+32768];
        std::vector<hit_t> threadShortResults;
        Sequence qSeq(par.maxSeqLen, querySeqType, subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, targetSeqType, subMat, 0, false, par.compBiasCorrection);
        SmithWaterman aligner(par.maxSeqLen, subMat->alphabetSize,
                              par.compBiasCorrection, par.compBiasCorrectionScale, NULL);

        std::string resultBuffer;
        resultBuffer.reserve(262144);
        for (size_t id = 0; id < qdbr->getSize(); id++) {
            char *querySeqData = qdbr->getData(id, thread_idx);
            size_t queryKey = qdbr->getDbKey(id);
            unsigned int querySeqLen = qdbr->getSeqLen(id);

            qSeq.mapSequence(id, queryKey, querySeqData, querySeqLen);
//            qSeq.printProfileStatePSSM();
            if(Parameters::isEqualDbtype(qSeq.getSeqType(), Parameters::DBTYPE_HMM_PROFILE) ){
                aligner.ssw_init(&qSeq, qSeq.getAlignmentProfile(), subMat);
            }else{
                aligner.ssw_init(&qSeq, tinySubMat, subMat);
            }
#pragma omp for schedule(static) nowait
            for (size_t tId = 0; tId < tdbr->getSize(); tId++) {
                unsigned int targetKey = tdbr->getDbKey(tId);
                if(taxonomyHook != NULL){
                    TaxID currTax = taxonomyHook->taxonomyMapping->lookup(targetKey);
                    if (taxonomyHook->expression[thread_idx]->isAncestor(currTax) == false) {
                        continue;
                    }
                }

                const bool isIdentity = (queryKey == targetKey && (par.includeIdentity || sameDB))? true : false;
                if(sequenceLookup == NULL){
                    char * targetSeq = tdbr->getData(tId, thread_idx);
                    unsigned int targetSeqLen = tdbr->getSeqLen(tId);
                    tSeq.mapSequence(tId, targetKey, targetSeq, targetSeqLen);
                    // mask numSequence
                    unsigned char xChar = subMat->aa2num[static_cast<int>('X')];
                    for (int i = 0; i < tSeq.L; i++) {
                        tSeq.numSequence[i] = ((targetSeq[i] >= 32 && targetSeq[i] <= 52) || targetSeq[i] >= 97)  ? xChar : tSeq.numSequence[i];
                    }
                }else{
                    tSeq.mapSequence(tId, targetKey, sequenceLookup->getSequence(tId));
                }
                float queryLength = qSeq.L;
                float targetLength = tSeq.L;
                if(Util::canBeCovered(par.covThr, par.covMode, queryLength, targetLength)==false){
                    continue;
                }

                bool hasEvalue = true;
                int score;
                if (alignmentMode == 0) {
                    score = aligner.ungapped_alignment(tSeq.numSequence, tSeq.L);
                } else {
                    std::string backtrace;
                    s_align res;
                    if (isIdentity) {
                        res = aligner.scoreIdentical(
                                tSeq.numSequence, tSeq.L, evaluer, Matcher::SCORE_ONLY, backtrace
                        );
                    } else {
                        res = aligner.ssw_align(
                                tSeq.numSequence,
                                tSeq.L,
                                backtrace,
                                par.gapOpen.values.aminoacid(),
                                par.gapExtend.values.aminoacid(),
                                Matcher::SCORE_ONLY,
                                par.evalThr,
                                evaluer,
                                par.covMode,
                                par.covThr,
                                par.correlationScoreWeight,
                                qSeq.L / 2
                        );
                    }
                    score = res.score1;
                    // check if evalThr != inf
                    double evalue = 0.0;
                    if (par.evalThr < std::numeric_limits<double>::max()) {
                        evalue = evaluer->computeEvalue(score, qSeq.L);
                    }
                    hasEvalue = (evalue <= par.evalThr);
                }
                bool hasDiagScore = (score > par.minDiagScoreThr);
                if (isIdentity || (hasDiagScore && hasEvalue)) {
                    hit_t hit;
                    hit.seqId = targetKey;
                    hit.prefScore = score;
                    hit.diagonal = 0;
                    threadShortResults.emplace_back(hit);
                }
            }
#pragma omp critical
            {
                shortResults.insert(shortResults.end(), threadShortResults.begin(), threadShortResults.end());
                threadShortResults.clear();
            }
#pragma omp barrier
#pragma omp master
            {
                SORT_PARALLEL(shortResults.begin(), shortResults.end(), hit_t::compareHitsByScoreAndId);
                size_t maxSeqs = std::min(par.maxResListLen, shortResults.size());
                for (size_t i = 0; i < maxSeqs; ++i) {
                    size_t len = QueryMatcher::prefilterHitToBuffer(buffer, shortResults[i]);
                    resultBuffer.append(buffer, len);
                }

                resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, 0);
                resultBuffer.clear();
                shortResults.clear();
                progress.updateProgress();
            }
#pragma omp barrier
        }
    }
}

int prefilterInternal(int argc, const char **argv, const Command &command, int mode) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    int outputDbtype = (par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED)
                      ? Parameters::DBTYPE_ALIGNMENT_RES : Parameters::DBTYPE_PREFILTER_RES;
    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, outputDbtype);
    resultWriter.open();
    bool sameDB = (par.db2.compare(par.db1) == 0);
    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader tDbrIdx(par.db2, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0 );
    IndexReader * qDbrIdx = NULL;
    DBReader<unsigned int> * qdbr = NULL;
    DBReader<unsigned int> * tdbr = tDbrIdx.sequenceReader;

    if (par.gpu == true) {
        const bool isGpuDb = DBReader<unsigned int>::getExtendedDbtype(tdbr->getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU;
        if (isGpuDb == false) {
            Debug(Debug::ERROR) << "Database " << FileUtil::baseName(par.db2) << " is not a valid GPU database\n" 
                                << "Please call: makepaddedseqdb " << FileUtil::baseName(par.db2) << " " << FileUtil::baseName(par.db2) << "_pad\n";
            EXIT(EXIT_FAILURE);
        }
    }

    const int targetSeqType = tdbr->getDbtype();
    int querySeqType;
    if (sameDB == true) {
        qDbrIdx = &tDbrIdx;
        qdbr = tdbr;
        querySeqType = targetSeqType;
    } else {
        // open the sequence, prefiltering and output databases
        qDbrIdx = new IndexReader(par.db1, par.threads, IndexReader::SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0);
        qdbr = qDbrIdx->sequenceReader;
        querySeqType = qdbr->getDbtype();
    }

    SequenceLookup * sequenceLookup = NULL;
    if(Parameters::isEqualDbtype(tDbrIdx.getDbtype(), Parameters::DBTYPE_INDEX_DB)){
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tDbrIdx.index);
        if(data.splits == 1){
            sequenceLookup = PrefilteringIndexReader::getSequenceLookup(0, tDbrIdx.index, par.preloadMode);
        }
    }
    BaseMatrix *subMat;
    EvalueComputation * evaluer;
    int8_t * tinySubMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
        evaluer = new EvalueComputation(tdbr->getAminoAcidDBSize(), subMat, par.gapOpen.values.nucleotide(), par.gapExtend.values.nucleotide());
        tinySubMat = new int8_t[subMat->alphabetSize*subMat->alphabetSize];
        for (int i = 0; i < subMat->alphabetSize; i++) {
            for (int j = 0; j < subMat->alphabetSize; j++) {
                tinySubMat[i*subMat->alphabetSize + j] = subMat->subMatrix[i][j];
            }
        }
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        evaluer = new EvalueComputation(tdbr->getAminoAcidDBSize(), subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
        tinySubMat = new int8_t[subMat->alphabetSize*subMat->alphabetSize];
        for (int i = 0; i < subMat->alphabetSize; i++) {
            for (int j = 0; j < subMat->alphabetSize; j++) {
                tinySubMat[i*subMat->alphabetSize + j] = subMat->subMatrix[i][j];
            }
        }
    }


    QueryMatcherTaxonomyHook * taxonomyHook = NULL;
    if(par.PARAM_TAXON_LIST.wasSet){
        taxonomyHook = new QueryMatcherTaxonomyHook(par.db2, tdbr, par.taxonList, par.threads);
    }
    if(par.gpu){
#ifdef HAVE_CUDA
        runFilterOnGpu(par, subMat, qdbr, tdbr, sameDB,
                       resultWriter, evaluer, taxonomyHook);
#else
        Debug(Debug::ERROR) << "MMseqs2 was compiled without CUDA support\n";
        EXIT(EXIT_FAILURE);
#endif
    }else{
        runFilterOnCpu(par, subMat, tinySubMat, qdbr, tdbr, sequenceLookup, sameDB,
                   resultWriter, evaluer, taxonomyHook,  mode);
    }

    resultWriter.close();

    if(taxonomyHook != NULL){
        delete taxonomyHook;
    }

    if (sequenceLookup != NULL) {
        delete sequenceLookup;
    }

    if(sameDB == false){
        delete qDbrIdx;
    }

    delete [] tinySubMat;
    delete subMat;
    delete evaluer;

    return 0;
}

int ungappedprefilter(int argc, const char **argv, const Command &command) {
    return prefilterInternal(argc, argv, command, 0);
}

int gappedprefilter(int argc, const char **argv, const Command &command) {
    return prefilterInternal(argc, argv, command, 1);
}
