// Batched pairwise Smith-Waterman scores for sparse (query, target) pairs (see pairsw.h).
//
// One cooperative group of G threads aligns one pair; groups pull pairs from a shared counter.  Thread k of a group
// holds R consecutive query columns (W = G * R per tile) and walks the target rows with a lag of k rows (wavefront):
// at step s it computes row s - k, taking H / E of the column left of its strip from thread k - 1 (shuffle) and
// keeping H / F of its own columns in registers.  Queries longer than W run tile by tile; the last thread of a tile
// stores H / E of the tile's last column per row in global memory for the next tile (MULTI kernels only).
// Recurrence (StripedSmithWaterman): E[c] = max(E[c-1] - ge, H[c-1] - go), F[j] = max(F[j-1] - ge, H[j-1] - go),
// H = max(0, H[j-1][c-1] + s, E, F); int32 throughout.
#include "pairsw.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include <cooperative_groups.h>
#include <cuda_runtime.h>

namespace cg = cooperative_groups;

#ifndef PAIRSW_DPX
#define PAIRSW_DPX 1
#endif

#define PSW_CHECK(x)                                                                                              \
    do {                                                                                                          \
        cudaError_t err_ = (x);                                                                                   \
        if (err_ != cudaSuccess) {                                                                                \
            fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(err_), __FILE__, __LINE__);            \
            exit(EXIT_FAILURE);                                                                                   \
        }                                                                                                         \
    } while (0)

namespace {

constexpr int NEG = -(1 << 28);
constexpr int PADL = 31;          // query padding letter: its column scores PADS
constexpr int PADS = -100000;
constexpr int BLOCK = 256;
constexpr int R = 16;             // query columns per thread
constexpr int PF = 8;             // target rows prefetched ahead
constexpr int NCLASS = 8;         // kernel classes: group size G = 4, 8, 16, 32 (tile width 16 G), single tile / tiles
constexpr int TBUCKETS = 4096;    // target length buckets (16 nt) for ordering within a class
constexpr size_t TEMP_BYTES = 256ull << 20;

template <int G, bool MULTI>
__global__ void __launch_bounds__(BLOCK) pairswKernel(const uint8_t *__restrict__ letters, const int8_t *__restrict__ bias,
                                                      const PairSW::Pair *__restrict__ pairs, const uint32_t *__restrict__ order,
                                                      uint32_t n, uint32_t *counter, int32_t *__restrict__ scores,
                                                      int2 *temp, size_t tempStride, size_t nGroups,
                                                      const int8_t *__restrict__ subMat,
                                                      const uint8_t *__restrict__ lmapG, int alph, int go, int ge) {
    __shared__ int smat[32 * 32];  // [target letter][query letter]
    __shared__ uint8_t lmap[256];  // raw byte -> letter
    for (int i = threadIdx.x; i < 256; i += blockDim.x) lmap[i] = lmapG[i];
    for (int i = threadIdx.x; i < 32 * 32; i += blockDim.x) {
        const int t = i / 32, q = i % 32;
        smat[i] = (q == PADL) ? PADS : (t < alph && q < alph) ? int(subMat[t * alph + q]) : PADS;
    }
    __syncthreads();
    auto tile = cg::tiled_partition<G>(cg::this_thread_block());
    const int lane = tile.thread_rank();
    const size_t gid = (size_t(blockIdx.x) * blockDim.x + threadIdx.x) / G;
    // the grid is rounded up to whole blocks: groups past nGroups have no tile buffer (temp holds nGroups strides)
    if (gid >= nGroups) return;
    int2 *tmp = MULTI ? temp + gid * tempStride : nullptr;
    constexpr int W = G * R;
    while (true) {
        uint32_t k = 0;
        if (lane == 0) k = atomicAdd(counter, 1u);
        k = tile.shfl(k, 0);
        if (k >= n) break;
        const uint32_t p = order[k];
        const PairSW::Pair P = pairs[p];
        const uint8_t *q = letters + P.qOff;
        const int8_t *qb = bias + P.qOff;
        const uint8_t *t = letters + P.tOff;
        const int qL = int(P.qLen), tL = int(P.tLen);
        const int nTiles = MULTI ? (qL + W - 1) / W : 1;
        int best = 0;
        for (int tt = 0; tt < nTiles; tt++) {
            int ql[R], b[R], H[R], F[R];
#pragma unroll
            for (int r = 0; r < R; r++) {
                const int c = tt * W + lane * R + r;
                ql[r] = (c < qL) ? int(lmap[q[c]]) : PADL;
                b[r] = (c < qL) ? int(qb[c]) : 0;
                H[r] = 0;
                F[r] = NEG;
            }
            int outH = 0, outE = NEG, lastIn = 0;
            // prefetch: the next PF target letters (and, for lane 0 of later tiles, boundary values) in a shift
            // register, so the global loads of a row are issued PF steps before the row is computed
            int tq[PF];
            int2 bq[PF];
#pragma unroll
            for (int i = 0; i < PF; i++) {
                tq[i] = (i < tL) ? int(t[i]) : 0;
                bq[i] = (MULTI && tt > 0 && lane == 0 && i < tL) ? tmp[i] : make_int2(0, NEG);
            }
            tile.sync();
            for (int s = 0; s < tL + G - 1; s++) {
                int inH = tile.shfl_up(outH, 1), inE = tile.shfl_up(outE, 1);
                const int j = s - lane;
                if (lane == 0) { inH = bq[0].x; inE = bq[0].y; }
                if (j >= 0 && j < tL) {
                    const int *srow = smat + int(lmap[tq[0]]) * 32;
                    // shift the prefetch registers and issue the load of row j + PF
#pragma unroll
                    for (int i = 0; i + 1 < PF; i++) { tq[i] = tq[i + 1]; bq[i] = bq[i + 1]; }
                    const int jn = j + PF;
                    tq[PF - 1] = (jn < tL) ? int(t[jn]) : 0;
                    bq[PF - 1] = (MULTI && tt > 0 && lane == 0 && jn < tL) ? tmp[jn] : make_int2(0, NEG);
                    int hd = (j == 0) ? 0 : lastIn, hl = inH, e = inE;
#pragma unroll
                    for (int r = 0; r < R; r++) {
#if PAIRSW_DPX
                        // DPX: fused add-max / three-way max with ReLU (hardware on sm_90+)
                        const int f = __viaddmax_s32(F[r], -ge, H[r] - go);
                        e = __viaddmax_s32(e, -ge, hl - go);
                        const int h = __vimax3_s32_relu(hd + srow[ql[r]] + b[r], e, f);
#else
                        const int f = max(F[r] - ge, H[r] - go);
                        e = max(e - ge, hl - go);
                        const int h = max(max(0, hd + srow[ql[r]] + b[r]), max(e, f));
#endif
                        hd = H[r];
                        H[r] = h;
                        F[r] = f;
                        hl = h;
                        best = max(best, h);
                    }
                    outH = hl;
                    outE = e;
                    lastIn = inH;
                    if (MULTI && lane == G - 1 && tt + 1 < nTiles) tmp[j] = make_int2(outH, outE);
                }
            }
            tile.sync();
        }
        // reduce the best score over the group
#pragma unroll
        for (int o = G / 2; o > 0; o /= 2) best = max(best, tile.shfl_down(best, o));
        if (lane == 0) scores[p] = best;
    }
}

// kernel class of a pair: the group size whose tiling pads the query least (cost = padded columns x (rows + lag),
// x 1.1 for tiling); class = 2 x log2(G / 4) + (tiles > 1).  One column past a tile width used to double the work
// (257 nt in a 512-wide tile: 142 instead of 538 GCUPS).
int classOf(uint32_t qLen, uint32_t tLen) {
    static const double pen = getenv("PAIRSW_TILE_PEN") ? atof(getenv("PAIRSW_TILE_PEN")) : 0.1;
    int bestC = 0;
    double bestCost = 0;
    for (int gi = 3; gi >= 0; gi--) {  // ties: the larger group (fewer tiles)
        const uint32_t W = 16u * (4u << gi);
        const uint32_t tiles = (qLen + W - 1) / W;
        const double cost = double(tiles) * W * (double(tLen) + (4 << gi) - 1) * (1.0 + pen * (tiles - 1));
        if (gi == 3 || cost < bestCost) { bestCost = cost; bestC = 2 * gi + (tiles > 1); }
    }
    return bestC;
}

template <int G, bool MULTI>
void launch(const PairSW &, cudaStream_t st, int numSMs, size_t groupsCap, const uint8_t *dL, const int8_t *dB,
            const PairSW::Pair *dP, const uint32_t *dO, uint32_t n, uint32_t *dC, int32_t *dS, int2 *dT,
            size_t stride, const int8_t *dM, const uint8_t *dMap, int alph, int go, int ge) {
    if (n == 0) return;
    int perSM = 0;
    PSW_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&perSM, pairswKernel<G, MULTI>, BLOCK, 0));
    size_t groups = std::min<size_t>(size_t(std::max(perSM, 1)) * numSMs * (BLOCK / G), n);
    groups = std::min(groups, groupsCap);
    const unsigned blocks = unsigned((groups * G + BLOCK - 1) / BLOCK);
    PSW_CHECK(cudaMemsetAsync(dC, 0, sizeof(uint32_t), st));
    pairswKernel<G, MULTI><<<blocks, BLOCK, 0, st>>>(dL, dB, dP, dO, n, dC, dS, dT, stride, groups, dM, dMap, alph, go, ge);
    PSW_CHECK(cudaGetLastError());
}

}  // namespace

int PairSW::deviceCount() {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess) return 0;
    return n;
}

size_t PairSW::deviceBytes(size_t maxLetters, size_t maxPairs) {
    return maxLetters * 2 + maxPairs * (sizeof(Pair) + 4 + 4) + TEMP_BYTES + 32 * 32 + 64;
}

PairSW::PairSW(const int8_t *subMat, int alphabetSize, int gapOpen, int gapExtend, size_t maxLetters, size_t maxPairs,
               const uint8_t *letterMap)
    : alphabetSize(alphabetSize), gapOpen(gapOpen), gapExtend(gapExtend), maxLetters(maxLetters), maxPairs(maxPairs),
      tempBytes(TEMP_BYTES) {
    if (alphabetSize > PADL) { fprintf(stderr, "PairSW: alphabet larger than 31 letters\n"); exit(EXIT_FAILURE); }
    cudaStream_t st;
    PSW_CHECK(cudaStreamCreate(&st));
    stream = st;
    PSW_CHECK(cudaMalloc(&dLetters, maxLetters));
    PSW_CHECK(cudaMalloc(&dBias, maxLetters));
    PSW_CHECK(cudaMalloc(&dPairs, maxPairs * sizeof(Pair)));
    PSW_CHECK(cudaMalloc(&dScores, maxPairs * sizeof(int32_t)));
    PSW_CHECK(cudaMalloc(&dOrder, maxPairs * sizeof(uint32_t)));
    PSW_CHECK(cudaMalloc(&dCounter, 64));
    PSW_CHECK(cudaMalloc(&dTemp, tempBytes));
    void *m;
    PSW_CHECK(cudaMalloc(&m, size_t(alphabetSize) * alphabetSize));
    PSW_CHECK(cudaMemcpy(m, subMat, size_t(alphabetSize) * alphabetSize, cudaMemcpyHostToDevice));
    dMat = m;
    uint8_t ident[256];
    for (int i = 0; i < 256; i++) ident[i] = uint8_t(i < 32 ? i : PADL);
    if (letterMap)
        for (int i = 0; i < 256; i++) ident[i] = letterMap[i] < alphabetSize ? letterMap[i] : uint8_t(PADL);
    PSW_CHECK(cudaMalloc(&dMap, 256));
    PSW_CHECK(cudaMemcpy(dMap, ident, 256, cudaMemcpyHostToDevice));
}

PairSW::~PairSW() {
    cudaFree(dLetters); cudaFree(dBias); cudaFree(dPairs); cudaFree(dScores); cudaFree(dOrder); cudaFree(dCounter);
    cudaFree(dTemp); cudaFree(dMat); cudaFree(dMap);
    cudaStreamDestroy(static_cast<cudaStream_t>(stream));
}

void PairSW::score(const uint8_t *letters, const int8_t *bias, size_t nLetters, const Pair *pairs, size_t nPairs,
                   int32_t *scores) {
    if (nLetters > maxLetters || nPairs > maxPairs) { fprintf(stderr, "PairSW: batch larger than allocated\n"); exit(EXIT_FAILURE); }
    if (nPairs == 0) return;
    cudaStream_t st = static_cast<cudaStream_t>(stream);
    // order: by query length class, then target length descending (similar pairs share warps; long ones first)
    std::vector<uint32_t> cnt(size_t(NCLASS) * TBUCKETS + 1, 0), order(nPairs);
    auto key = [&](const Pair &p) {
        const uint32_t tb = std::min<uint32_t>(p.tLen >> 4, TBUCKETS - 1);
        return size_t(classOf(p.qLen, p.tLen)) * TBUCKETS + (TBUCKETS - 1 - tb);
    };
    for (size_t i = 0; i < nPairs; i++) cnt[key(pairs[i]) + 1]++;
    for (size_t i = 1; i < cnt.size(); i++) cnt[i] += cnt[i - 1];
    size_t start[NCLASS + 1];
    for (int c = 0; c <= NCLASS; c++) start[c] = cnt[size_t(c) * TBUCKETS];
    uint32_t maxT[NCLASS] = {0};
    for (size_t i = 0; i < nPairs; i++) {
        order[cnt[key(pairs[i])]++] = uint32_t(i);
        const int c = classOf(pairs[i].qLen, pairs[i].tLen);
        maxT[c] = std::max(maxT[c], pairs[i].tLen);
    }
    PSW_CHECK(cudaMemcpyAsync(dLetters, letters, nLetters, cudaMemcpyHostToDevice, st));
    PSW_CHECK(cudaMemcpyAsync(dBias, bias, nLetters, cudaMemcpyHostToDevice, st));
    PSW_CHECK(cudaMemcpyAsync(dPairs, pairs, nPairs * sizeof(Pair), cudaMemcpyHostToDevice, st));
    PSW_CHECK(cudaMemcpyAsync(dOrder, order.data(), nPairs * sizeof(uint32_t), cudaMemcpyHostToDevice, st));
    int dev = 0, numSMs = 0;
    PSW_CHECK(cudaGetDevice(&dev));
    PSW_CHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, dev));
    auto L = static_cast<const uint8_t *>(dLetters);
    auto B = static_cast<const int8_t *>(dBias);
    auto P = static_cast<const Pair *>(dPairs);
    auto O = static_cast<const uint32_t *>(dOrder);
    auto C = static_cast<uint32_t *>(dCounter);
    auto S = static_cast<int32_t *>(dScores);
    auto T = static_cast<int2 *>(dTemp);
    auto M = static_cast<const int8_t *>(dMat);
    auto MP = static_cast<const uint8_t *>(dMap);
    const size_t NOCAP = ~size_t(0);
    static const bool prof = getenv("PAIRSW_PROFILE") != nullptr;
    cudaEvent_t ev[NCLASS + 2];
    if (prof) for (auto &e : ev) PSW_CHECK(cudaEventCreate(&e));
    if (prof) PSW_CHECK(cudaEventRecord(ev[0], st));
    auto n = [&](int c) { return uint32_t(start[c + 1] - start[c]); };
    // classes run one after another on the stream (the counter and the tile buffer are shared)
    auto cap = [&](int c) { return std::max<size_t>(1, tempBytes / (size_t(std::max<uint32_t>(maxT[c], 1)) * sizeof(int2))); };
    auto stride = [&](int c) { return size_t(std::max<uint32_t>(maxT[c], 1)); };
#define PSW_LAUNCH(c, G, MULTI)                                                                                        \
    launch<G, MULTI>(*this, st, numSMs, MULTI ? cap(c) : NOCAP, L, B, P, O + start[c], n(c), C, S, T,                 \
                     MULTI ? stride(c) : 0, M, MP, alphabetSize, gapOpen, gapExtend);                                  \
    if (prof) PSW_CHECK(cudaEventRecord(ev[c + 1], st));
    PSW_LAUNCH(0, 4, false) PSW_LAUNCH(1, 4, true) PSW_LAUNCH(2, 8, false) PSW_LAUNCH(3, 8, true)
    PSW_LAUNCH(4, 16, false) PSW_LAUNCH(5, 16, true) PSW_LAUNCH(6, 32, false) PSW_LAUNCH(7, 32, true)
#undef PSW_LAUNCH
    PSW_CHECK(cudaMemcpyAsync(scores, dScores, nPairs * sizeof(int32_t), cudaMemcpyDeviceToHost, st));
    PSW_CHECK(cudaStreamSynchronize(st));
    if (prof) {
        double cl[NCLASS] = {0};
        for (size_t i = 0; i < nPairs; i++) cl[classOf(pairs[i].qLen, pairs[i].tLen)] += double(pairs[i].qLen) * pairs[i].tLen;
        fprintf(stderr, "PairSW batch: %zu pairs;", nPairs);
        for (int c = 0; c < NCLASS; c++) {
            float ms = 0;
            PSW_CHECK(cudaEventElapsedTime(&ms, ev[c], ev[c + 1]));
            fprintf(stderr, " class %d: %u pairs %.0f ms %.0f GCUPS;", c, n(c), ms, ms > 0 ? cl[c] / ms / 1e6 : 0.0);
        }
        fprintf(stderr, "\n");
        for (auto &e : ev) cudaEventDestroy(e);
    }
}
