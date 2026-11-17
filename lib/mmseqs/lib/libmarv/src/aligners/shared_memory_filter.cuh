#ifndef LIBMARV_SHARED_MEMORY_FILTER_CUH
#define LIBMARV_SHARED_MEMORY_FILTER_CUH

#include "../cuda_hip_rename.h"
#include "../util.cuh"
#include "../cuda_errorcheck.cuh"
#include "../types.hpp"

#include <utility>
#include <vector>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN

    // Shared memory of the PSSM a kernel holds, as computed at its launch. The tuned configs assume
    // 21-row protein PSSMs; with more rows (MARV_PSSM_ROWS) some no longer fit on the device.
    template<class PssmScoreType, class SmemScoreType, int groupsize, int numRegs>
    constexpr int gaplessSharedPssmBytes(){
        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
        return sizeof(SharedPSSM_singletile<SmemScoreType, MARV_PSSM_ROWS, numColumnsPSSM>);
    }

    template<class PssmScoreType, int groupsize, int numRegs>
    constexpr int smithWatermanSharedPssmBytes(){
        return sizeof(SharedPSSM_singletile<PssmScoreType, MARV_PSSM_ROWS, groupsize * numRegs>);
    }

    // Drops the aligners whose PSSM exceeds the device's opt-in shared memory limit (they would fail to
    // launch), keeping the performance hints, which are indexed like the aligners, in step.
    template<class AlignerVector>
    void dropAlignersExceedingSharedMemory(AlignerVector& aligners, std::vector<int>& hints){
        int deviceId = 0;
        CUDACHECK(cudaGetDevice(&deviceId));
        int maxSharedMemory = 0;
        CUDACHECK(cudaDeviceGetAttribute(&maxSharedMemory, cudaDevAttrMaxSharedMemoryPerBlockOptin, deviceId));
        const bool withHints = hints.size() == aligners.size();
        size_t kept = 0;
        for(size_t i = 0; i < aligners.size(); i++){
            if(aligners[i]->getSharedMemoryBytes() > maxSharedMemory){
                continue;
            }
            if(withHints){
                hints[kept] = hints[i];
            }
            aligners[kept++] = std::move(aligners[i]);
        }
        aligners.resize(kept);
        if(withHints){
            hints.resize(kept);
        }
    }

LIBMARV_NAMESPACE_WITH_NESTING_END

#endif
