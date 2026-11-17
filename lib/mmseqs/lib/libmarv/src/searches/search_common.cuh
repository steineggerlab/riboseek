#ifndef LIBMARV_SEARCH_COMMON_CUH
#define LIBMARV_SEARCH_COMMON_CUH

#include "../cuda_hip_rename.h"

#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
#endif

#include <memory>
#include <vector>

#include "../aligners/score_only_alignment_interfaces.cuh"
#include "../aligners/score_with_endpos_alignment_interfaces.cuh"

#include "../config.hpp"
#include "../cuda_errorcheck.cuh"

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

struct BenchmarkStats{
    double seconds{};
    double gcups{};
};

struct SearchResult{
    std::vector<int> scores{};
    std::vector<ReferenceIdT> referenceIds{};
    std::vector<AlignmentEndPosition> endPositions{};
    std::unique_ptr<BenchmarkStats> stats{};
};

struct IsOverflowScore{
    int threshold;

    __host__ __device__
    IsOverflowScore(int t) : threshold(t){}

    __host__ __device__
    bool operator()(int score){
        return score >= threshold;
    }
};

struct IsIndexToOverflowScore{
    int threshold;
    const int* scores;

    __host__ __device__
    IsIndexToOverflowScore(int t, const int* s) : threshold(t), scores(s){}

    __host__ __device__
    bool operator()(int index){
        return scores[index] >= threshold;
    }
};


struct OverflowSettings{
    enum class AlignerType{
        ScoreOnly,
        ScoreEndpos
    };

    bool checkForOverflow{};
    AlignerType alignerType{};
    cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*> overflowAlignersScoreOnly{};
    cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> overflowAlignersScoreWithEndpos{};


    static OverflowSettings withoutOverflowCorrection(){
        OverflowSettings result;
        result.checkForOverflow = false;
        return result;
    }

    //re-compute all alignments with a score >= threshold
    static OverflowSettings withOverflowCorrection(cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*> aligners){
        OverflowSettings result;
        result.checkForOverflow = true;
        result.overflowAlignersScoreOnly = aligners;
        result.alignerType = AlignerType::ScoreOnly;
        return result;
    }

    //re-compute all alignments with a score >= threshold
    static OverflowSettings withOverflowCorrection(cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> aligners){
        OverflowSettings result;
        result.checkForOverflow = true;
        result.overflowAlignersScoreWithEndpos = aligners;
        result.alignerType = AlignerType::ScoreEndpos;
        return result;
    }
};

struct MaskingOptions{
    MaskingOptions() = default;
    MaskingOptions(int threshold) : MaskingOptions(threshold, 20){}
    MaskingOptions(int threshold, char letter) : maskingLetter(letter), maskingThreshold(threshold){}

    char maskingLetter = 20;
    int maskingThreshold = 6; //mask runs of length >= maskingThreshold
};




    // Gathers one GPU's top results into the master GPU's buffer. Copies between two GPUs have been
    // seen to deliver wrong data silently on some multi-GPU systems even with peer access reported as
    // available, so results of another device go through host memory (riboseek). The caller has made
    // srcDevice current; it is current again on return.
    inline void copyResultsToMasterDevice(void* masterDst, const void* src, size_t bytes,
                                          int srcDevice, int masterDevice,
                                          cudaStream_t srcStream, cudaStream_t masterStream){
        if(srcDevice == masterDevice){
            CUDACHECK(cudaMemcpyAsync(masterDst, src, bytes, cudaMemcpyDeviceToDevice, srcStream));
            return;
        }
        std::vector<char> staging(bytes);
        CUDACHECK(cudaStreamSynchronize(srcStream));
        CUDACHECK(cudaMemcpy(staging.data(), src, bytes, cudaMemcpyDeviceToHost));
        CUDACHECK(cudaSetDevice(masterDevice));
        CUDACHECK(cudaStreamSynchronize(masterStream));
        CUDACHECK(cudaMemcpy(masterDst, staging.data(), bytes, cudaMemcpyHostToDevice));
        CUDACHECK(cudaSetDevice(srcDevice));
    }

LIBMARV_NAMESPACE_END

#endif