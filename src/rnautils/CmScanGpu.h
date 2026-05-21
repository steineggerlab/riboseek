#ifndef RIBOSEEK_CMSCAN_GPU_H
#define RIBOSEEK_CMSCAN_GPU_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <array>

struct CmFastDeckGpuState {
    int type = 0;
    int dmin = 0;
    int dmax = -1;
    int bLeft = -1;
    int bRight = -1;
    int dConsume = 0;
    int niShift = 0;
    int emitSize = 0;
    int emitOffset = -1;
    int trCount = 0;
    size_t trOff = 0;
    uint8_t consumeMask = 0;
    float endSc = 0.0f;
    float null2Agg[4] = {0.0f, 0.0f, 0.0f, 0.0f};
};

struct CmFastDeckGpuBatchJob {
    int N = 0;
    int minSpan = 1;
    int maxSpan = 0;
    int forcedI = -1;
    int forcedD = -1;
};

struct CmFastDeckGpuBatchResult {
    int found = 0;
    int traceOverflow = 0;
    int bestI = 1;
    int bestD = 0;
    int minUsed = 0;
    int maxUsed = 0;
    int traceLen = 0;
    int bestMode = 'J';
    double bestSc = 0.0;
    double bestNull3Corr = 0.0;
    int obsCount[4] = {0, 0, 0, 0};
    double modelAggRaw[4] = {0.0, 0.0, 0.0, 0.0};
};

#ifdef HAVE_CUDA
static constexpr size_t CM_FASTDECK_GPU_TBCELL_BYTES = sizeof(int) * 3;

bool cmFastDeckGpuGetMemoryInfo(size_t *freeBytesOut, size_t *totalBytesOut);

bool runInfernalExactScanFastDeckGpu(
    int N,
    int M,
    int iStride,
    size_t stateStride,
    size_t cells,
    const int8_t *seqCode,
    const std::vector<CmFastDeckGpuState> &states,
    const std::vector<int> &activeStates,
    const std::vector<size_t> &stateBase,
    const std::vector<int> &bSplitBegByVD,
    const std::vector<int> &bSplitEndByVD,
    const std::vector<uint16_t> &trDst,
    const std::vector<float> &trSc,
    const std::vector<float> &emitData,
    bool hasLocalCfg,
    float elSelf,
    float *hostVitOut,
    std::string *errorOut);

bool runInfernalExactScanFastDeckGpuBatch(
    int M,
    int rootState,
    int maxN,
    int batchSize,
    int iStride,
    size_t stateStride,
    size_t cellsPerJob,
    const std::vector<CmFastDeckGpuBatchJob> &jobs,
    const std::vector<int8_t> &packedSeqCode,
    const std::vector<int> &packedPrefA,
    const std::vector<int> &packedPrefC,
    const std::vector<int> &packedPrefG,
    const std::vector<int> &packedPrefU,
    const std::vector<double> &log2Int,
    const std::vector<CmFastDeckGpuState> &states,
    const std::vector<int> &activeStates,
    const std::vector<size_t> &stateBase,
    const std::vector<int> &bSplitBegByVD,
    const std::vector<int> &bSplitEndByVD,
    const std::vector<uint16_t> &trDst,
    const std::vector<float> &trSc,
    const std::vector<float> &emitData,
    bool truncModesEnabled,
    bool null3Enabled,
    double log2Omega3,
    const std::array<double, 4> &log2Null,
    int traceCapacity,
    std::vector<CmFastDeckGpuBatchResult> *hostResultsOut,
    std::vector<uint16_t> *hostTraceStatesOut,
    std::string *errorOut);
#else
static constexpr size_t CM_FASTDECK_GPU_TBCELL_BYTES = sizeof(int) * 3;

inline bool cmFastDeckGpuGetMemoryInfo(size_t *, size_t *) {
    return false;
}

inline bool runInfernalExactScanFastDeckGpu(
    int,
    int,
    int,
    size_t,
    size_t,
    const int8_t *,
    const std::vector<CmFastDeckGpuState> &,
    const std::vector<int> &,
    const std::vector<size_t> &,
    const std::vector<int> &,
    const std::vector<int> &,
    const std::vector<uint16_t> &,
    const std::vector<float> &,
    const std::vector<float> &,
    bool,
    float,
    float *,
    std::string *) {
    return false;
}

inline bool runInfernalExactScanFastDeckGpuBatch(
    int,
    int,
    int,
    int,
    int,
    size_t,
    size_t,
    const std::vector<CmFastDeckGpuBatchJob> &,
    const std::vector<int8_t> &,
    const std::vector<int> &,
    const std::vector<int> &,
    const std::vector<int> &,
    const std::vector<int> &,
    const std::vector<double> &,
    const std::vector<CmFastDeckGpuState> &,
    const std::vector<int> &,
    const std::vector<size_t> &,
    const std::vector<int> &,
    const std::vector<int> &,
    const std::vector<uint16_t> &,
    const std::vector<float> &,
    const std::vector<float> &,
    bool,
    bool,
    double,
    const std::array<double, 4> &,
    int,
    std::vector<CmFastDeckGpuBatchResult> *,
    std::vector<uint16_t> *,
    std::string *) {
    return false;
}
#endif

#endif
