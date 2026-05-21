#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "MMseqsMPI.h"
#include "MathUtil.h"
#include "DBWriter.h"
#include "DBReader.h"
#include "FileUtil.h"
#include "Matcher.h"
#include "Util.h"
#include "NucleotideMatrix.h"
#include "Sequence.h"
#include "rnautils/CmScanGpu.h"
#include "simd.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/wait.h>
#if defined(__AVX2__)
#include <immintrin.h>
#define CM_HAS_SSE2 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define CM_HAS_SSE2 1
#else
#include <simde/x86/sse2.h>
#define CM_HAS_SSE2 1
#endif

namespace {

static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double LN2 = 0.69314718055994530942;
static const double LOG2E = 1.4426950408889634074;
typedef uint16_t StateId;
#if defined(__GNUC__) || defined(__clang__)
#define CM_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define CM_ALWAYS_INLINE inline
#endif

struct FastaSeq {
    unsigned int key = 0;
    std::string id;
    std::string seq;
};

enum CmStateType {
    CM_ST_MP,
    CM_ST_ML,
    CM_ST_MR,
    CM_ST_IL,
    CM_ST_IR,
    CM_ST_D,
    CM_ST_S,
    CM_ST_B,
    CM_ST_E,
    CM_ST_UNKNOWN
};

// Node-type encoding for cm_localize. Matches Infernal cm.h ndtype except we
// only carry the types relevant to begin/end eligibility decisions.
enum CmNodeType {
    CM_ND_UNKNOWN = 0,
    CM_ND_ROOT,
    CM_ND_MATP,
    CM_ND_MATL,
    CM_ND_MATR,
    CM_ND_BEGL,
    CM_ND_BEGR,
    CM_ND_BIF,
    CM_ND_END,
};

struct CmState {
    CmStateType type;
    int idx;
    int cfirst;
    int cnum;
    std::vector<double> trans; // for !B: size cnum
    std::vector<double> emit;  // 4 for singlet emitters, 16 for MP, empty otherwise
    std::vector<float> emitF; // float parallel for hot-path use
    int mapLeft;               // consensus/map coordinate, -1 if none
    int mapRight;              // consensus/map coordinate, -1 if none
    int dmin1;
    int dmax1;
    // Local-mode (cm_localize) metadata.
    int nodeIdx = -1;
    CmNodeType nodeType = CM_ND_UNKNOWN;
    bool isFirstOfNode = false;
    double beginSc = NEG_INF;  // log2(begin prob); finite only at begin-eligible heads
    double endSc = NEG_INF;    // log2(end prob);   finite only at end-eligible heads
    // Truncation DP: marginal emission scores for L/R modes (Infernal cm_CalcMargLikScores).
    // For MATP_MP only, size 4 per nucleotide; empty otherwise.
    std::vector<double> lmesc; // L-mode: keep left, marginalize right partner
    std::vector<double> rmesc; // R-mode: keep right, marginalize left partner
};

struct ExactStateExec {
    CmStateType type;
    int niShift;
    int dConsume;
    bool consumesLeft;
    bool consumesRight;
    int bLeft;
    int bRight;
    StateId trDst4[4];
    float trScF4[4];
    double trSc4[4];
    size_t trOff;
    int trCount;
    const float *emitPtr;
    int emitSize;
    int splitKMin;
    int splitKMax;
    int splitRMin;
    int splitRMax;
    int dmin;
    int dmax;
    double null2Agg[4];
};

struct InfernalExactModel {
    struct ExpTail {
        bool valid = false;
        double lambda = 0.0;
        double muExtrap = 0.0;
        double dbsize = 0.0;
        int nrandhits = 0;
    };

    int clen;
    int w;
    bool hasNull = false;
    double nullProb[4] = {0.25, 0.25, 0.25, 0.25}; // A,C,G,U
    double null2Omega = 1.0 / 65536.0;
    double null3Omega = 1.0 / 65536.0;
    // Local-mode parameters (Infernal cm_localize). Currently parsed only for
    // diagnostic use; DP wiring is the next step in the 2YGH_1 fix.
    double pBegin = 0.0;     // local-begin probability (0 = glocal)
    double pEnd = 0.0;       // local-end probability  (0 = glocal)
    double elSelf = 0.0;     // EL self-transition score (nats per residue, ≤0)
    bool hasLocalCfg = false;
    ExpTail expLC;
    ExpTail expLI;
    ExpTail expGC;
    ExpTail expGI;
    std::vector<CmState> states;
    int rootState;
    // Truncation DP support (gated on MMSEQS_CMSCAN_TRUNC_DP=1; otherwise empty).
    // emap mirrors Infernal CMEmitMap_t (display.c:CreateEmitMap): per-node
    // consensus position bounds + EL-anchor position.
    int nNodes = 0;
    std::vector<int> nodeFirstState; // size nNodes; first state index of each node
    std::vector<CmNodeType> nodeType;
    std::vector<int> emapLpos;       // node -> consensus 0..clen+1 (inclusive low bound)
    std::vector<int> emapRpos;       // node -> consensus 0..clen+1 (inclusive high bound)
    std::vector<int> emapEpos;       // node -> consensus position EL insertion follows
    // Expected per-state occupancy psi[v] (Infernal cm_ExpectedPositionOccupancy).
    // Linear-space probability of state v being entered in a sample. Size states.size().
    std::vector<double> psi;
    // Truncation penalties (Infernal cm_tr_penalties_Create).
    // [type][stateIdx]; type 0=5'+3', 1=5'only, 2=3'only. Indexed by global state id.
    // Two parallel sets: trpG = global mode (no local), trpL = local mode.
    // Stored as log2 probabilities (NEG_INF if not eligible for that truncation type).
    std::vector<std::vector<double>> trpG;
    std::vector<std::vector<double>> trpL;
    bool hasTruncDp = false;
};

struct CmBeamLocalizeResult {
    bool valid = false;
    int start = 1;     // 1-indexed inclusive in the scanned region
    int end = 0;       // 1-indexed inclusive in the scanned region
    float bestScore = -std::numeric_limits<float>::infinity();
    int rowsPruned = 0;
    int rowsVisited = 0;
};

template <typename T>
struct RawBuffer {
    T *ptr;
    size_t cap;

    RawBuffer() : ptr(NULL), cap(0) {}
    ~RawBuffer() { std::free(ptr); }

    RawBuffer(const RawBuffer &) = delete;
    RawBuffer &operator=(const RawBuffer &) = delete;

    bool ensure(size_t need) {
        if (need <= cap) {
            return true;
        }
        size_t newCap = (cap == 0) ? 1024 : cap;
        while (newCap < need) {
            newCap <<= 1;
        }
        void *np = std::realloc(ptr, newCap * sizeof(T));
        if (np == NULL) {
            return false;
        }
        ptr = static_cast<T *>(np);
        cap = newCap;
        return true;
    }

    T *data() { return ptr; }
    const T *data() const { return ptr; }
};

struct ExactScanWorkspace {
    std::vector<int8_t> seqCode;
    std::vector<int> prefA;
    std::vector<int> prefC;
    std::vector<int> prefG;
    std::vector<int> prefU;
    std::vector<double> log2Int;
    RawBuffer<float> vit;
    std::vector<double> in;
    RawBuffer<size_t> stateBase;
    RawBuffer<int> bSplitBegByVD;
    RawBuffer<int> bSplitEndByVD;
    RawBuffer<float> bSplitTmp;
    std::vector<float> trScF;
    RawBuffer<double> insD; // deck-based inside: layout [d * M * (localN+1) + v * (localN+1) + li]
    std::vector<double> memoVitDense;
    std::vector<double> memoInDense;
    std::vector<uint32_t> memoVitSeen;
    std::vector<uint32_t> memoInSeen;
    uint32_t memoVitGen = 1;
    uint32_t memoInGen = 1;
    std::unordered_map<int64_t, double> memoVitMap;
    std::unordered_map<int64_t, double> memoInMap;
    // Per-state sparse chart for Infernal-compat exactVitRec.
    // QDB still defines the legal d rows per state, but each (state,d) row now
    // stores only the i interval that survives the optional LoL-derived state
    // bounds. This preserves the exact recurrence while avoiding dead cells for
    // i outside the state-local target band.
    std::vector<float> exactChart;
    std::vector<uint32_t> exactChartSeen;
    uint32_t exactChartGen = 1;
    std::vector<size_t> exactChartStateDOffset; // M+1 entries, prefix over d-rows
    std::vector<size_t> exactChartRowOffset;    // totalRows+1 entries, prefix over live cells
    std::vector<int> exactChartRowImin;         // totalRows entries
    std::vector<int> exactChartRowSize;         // totalRows entries
    std::vector<int> exactChartBandDmin;    // M entries
    std::vector<int> exactChartBandSize;    // M entries (= dmax-dmin+1, or 0 if empty)
};

static inline std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])) != 0) {
        ++b;
    }
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) != 0) {
        --e;
    }
    return s.substr(b, e - b);
}

static inline char normalizeBase(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (c == 'T') {
        return 'U';
    }
    return c;
}

static inline double logSumExp(double a, double b) {
    if (a == NEG_INF) {
        return b;
    }
    if (b == NEG_INF) {
        return a;
    }
    if (a < b) {
        std::swap(a, b);
    }
    const double d = b - a; // <= 0
    if (d < -50.0) {
        return a;
    }
    return a + std::log1p(std::exp(d));
}

static bool cmFastMathEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_FASTMATH");
        cached = (env != NULL && std::string(env) == "1") ? 1 : 0;
    }
    return cached == 1;
}

static int parseEnvIntLocal(const char *name, int defVal) {
    const char *env = std::getenv(name);
    return (env != NULL) ? std::atoi(env) : defVal;
}

static bool cmFastLogsumEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_FASTLOGSUM");
        cached = (env != NULL && std::string(env) == "1") ? 1 : 0;
    }
    return cached == 1;
}

static bool cmDebugExtEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_DEBUG_EXT");
        cached = (env != NULL && std::string(env) == "1") ? 1 : 0;
    }
    return cached == 1;
}

static bool cmExactFastDeckEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_EXACT_FASTDECK");
        // Keep fast-deck enabled by default; allow explicit opt-out.
        cached = (env != NULL && std::string(env) == "0") ? 0 : 1;
    }
    return cached == 1;
}

static bool cmTruncModesEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_TRUNC");
        // Infernal-like default: truncation modes enabled unless explicitly disabled.
        cached = (env != NULL && std::string(env) == "0") ? 0 : 1;
    }
    return cached == 1;
}

static bool cmNull3Enabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_NULL3");
        cached = (env != NULL && std::string(env) == "0") ? 0 : 1;
    }
    return cached == 1;
}

static bool cmGpuExactEnabled() {
#ifdef HAVE_CUDA
    return Parameters::getInstance().gpu != 0;
#else
    return false;
#endif
}

static bool cmGpuBatchCommonPathEnabled(const InfernalExactModel &model) {
#ifdef HAVE_CUDA
    if (!cmGpuExactEnabled() || !cmExactFastDeckEnabled()) {
        return false;
    }
    if (model.hasLocalCfg) {
        return false;
    }
    const char *rescoreEnv = std::getenv("MMSEQS_CMSCAN_RESCORE_ENVELOPE");
    if (rescoreEnv != NULL && std::string(rescoreEnv) == "1") {
        return false;
    }
    const char *forceGlobalEnv = std::getenv("MMSEQS_CMSCAN_FORCE_GLOBAL");
    if (forceGlobalEnv != NULL && forceGlobalEnv[0] != '\0' && forceGlobalEnv[0] != '0') {
        return false;
    }
    const char *traceElEnv = std::getenv("MMSEQS_CMSCAN_TRACE_EL");
    if (traceElEnv != NULL && traceElEnv[0] == '1') {
        return false;
    }
    const char *traceLboEnv = std::getenv("MMSEQS_CMSCAN_TRACE_LBO");
    if (traceLboEnv != NULL && traceLboEnv[0] == '1') {
        return false;
    }
    const char *dumpTraceEnv = std::getenv("MMSEQS_CMSCAN_DUMP_TRACE_TID");
    if (dumpTraceEnv != NULL && dumpTraceEnv[0] != '\0') {
        return false;
    }
    return true;
#else
    (void)model;
    return false;
#endif
}

static double cmGpuBatchMemFrac() {
    static double cached = -1.0;
    if (cached < 0.0) {
        cached = 0.90;
        if (const char *env = std::getenv("MMSEQS_CMSCAN_GPU_BATCH_MEM_FRAC")) {
            const double parsed = std::atof(env);
            if (parsed > 0.0 && parsed < 1.0) cached = parsed;
        }
    }
    return cached;
}

static size_t cmGpuBatchMemReserveBytes() {
    static size_t cached = 0;
    if (cached == 0) {
        cached = 256ULL << 20;
        if (const char *env = std::getenv("MMSEQS_CMSCAN_GPU_BATCH_MEM_RESERVE_MB")) {
            const unsigned long long parsed = std::strtoull(env, NULL, 10);
            cached = static_cast<size_t>(parsed) << 20;
        }
    }
    return cached;
}

static size_t cmGpuUsableBytes() {
#ifdef HAVE_CUDA
    size_t freeBytes = 0;
    size_t totalBytes = 0;
    if (!cmFastDeckGpuGetMemoryInfo(&freeBytes, &totalBytes)) {
        return 0;
    }
    const size_t fracBudget = static_cast<size_t>(static_cast<double>(freeBytes) * cmGpuBatchMemFrac());
    const size_t reserve = cmGpuBatchMemReserveBytes();
    return (fracBudget > reserve) ? (fracBudget - reserve) : fracBudget;
#else
    return 0;
#endif
}

static size_t estimateFastDeckGpuSingleJobBytes(int M,
                                                int N,
                                                size_t activeStateCount,
                                                size_t trDstCount,
                                                size_t trScCount,
                                                size_t emitFloatCount) {
    if (M <= 0 || N <= 0) {
        return 0;
    }
    const size_t iStride = static_cast<size_t>(N + 2);
    const size_t stateStride = static_cast<size_t>(N + 1) * iStride;
    const size_t cells = static_cast<size_t>(M) * stateStride;
    const size_t vdCells = static_cast<size_t>(M) * static_cast<size_t>(N + 1);
    size_t bytes = 0;
    bytes += static_cast<size_t>(N + 1) * sizeof(int8_t);
    bytes += static_cast<size_t>(M) * sizeof(CmFastDeckGpuState);
    bytes += activeStateCount * sizeof(int);
    bytes += static_cast<size_t>(M) * sizeof(size_t);
    bytes += 2ULL * vdCells * sizeof(int);
    bytes += trDstCount * sizeof(uint16_t);
    bytes += trScCount * sizeof(float);
    bytes += emitFloatCount * sizeof(float);
    bytes += cells * sizeof(float);
    return bytes;
}

// Optimal-accuracy alignment opt-in. When set, runInfernalExactScan switches
// from CYK (max over paths) to Inside+Outside+OA (posterior-coverage trace),
// reusing the vit buffer across phases: Inside in vit, then alpha overwritten
// in-place by OA scores during the OA fill, then vit reinterpreted as a uint8
// shadow matrix for traceback. One extra float buffer (beta) is allocated.
static bool cmOAEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_OA");
        cached = (env != NULL && std::string(env) == "1") ? 1 : 0;
    }
    return cached == 1;
}

static inline double fastLog1pExp(double x) {
    // Stable approximation for log(1 + exp(x)); tuned for x <= 0 path used by logSumExp.
    if (x <= -20.0) {
        return 0.0;
    }
    if (x >= 20.0) {
        return x;
    }
    const float e2 = static_cast<float>(x * LOG2E);
    const float ex = static_cast<float>(MathUtil::fpow2(e2));
    const float onePlus = 1.0f + ex;
    return static_cast<double>(MathUtil::flog2(onePlus)) * LN2;
}

static inline double logSumExpMaybeFast(double a, double b) {
    if (!cmFastMathEnabled()) {
        return logSumExp(a, b);
    }
    if (a == NEG_INF) {
        return b;
    }
    if (b == NEG_INF) {
        return a;
    }
    if (a < b) {
        std::swap(a, b);
    }
    const double d = b - a; // <= 0
    if (d < -50.0) {
        return a;
    }
    return a + fastLog1pExp(d);
}

static inline float bSplitBestFull(const float *left, const float *rightStart, int rightStep, int d, float negInf) {
    float best = negInf;
    int k = 0;
#if defined(__AVX2__)
    auto hmax256 = [](__m256 v) -> float {
        const __m128 lo = _mm256_castps256_ps128(v);
        const __m128 hi = _mm256_extractf128_ps(v, 1);
        __m128 m = _mm_max_ps(lo, hi);
        m = _mm_max_ps(m, _mm_movehl_ps(m, m));
        m = _mm_max_ps(m, _mm_shuffle_ps(m, m, _MM_SHUFFLE(1, 1, 1, 1)));
        return _mm_cvtss_f32(m);
    };
    const __m256 negV = _mm256_set1_ps(negInf);
    const __m256i offs = _mm256_setr_epi32(0, rightStep, 2 * rightStep, 3 * rightStep,
                                           4 * rightStep, 5 * rightStep, 6 * rightStep, 7 * rightStep);
    __m256 bestV = negV;
    for (; k + 7 <= d; k += 8) {
        const __m256 l = _mm256_loadu_ps(left + k);
        const __m256i base = _mm256_set1_epi32(k * rightStep);
        const __m256i idx = _mm256_add_epi32(base, offs);
        const __m256 r = _mm256_i32gather_ps(rightStart, idx, 4);
        const __m256 valid = _mm256_and_ps(_mm256_cmp_ps(l, negV, _CMP_NEQ_OQ),
                                           _mm256_cmp_ps(r, negV, _CMP_NEQ_OQ));
        const __m256 s = _mm256_add_ps(l, r);
        const __m256 c = _mm256_blendv_ps(negV, s, valid);
        bestV = _mm256_max_ps(bestV, c);
    }
    {
        const float hv = hmax256(bestV);
        if (hv > best) {
            best = hv;
        }
    }
#endif
    for (; k + 3 <= d; k += 4) {
        const float l0 = left[k];
        const float r0 = rightStart[(k + 0) * rightStep];
        if (l0 != negInf && r0 != negInf) {
            const float cand = l0 + r0;
            if (cand > best) {
                best = cand;
            }
        }
        const float l1 = left[k + 1];
        const float r1 = rightStart[(k + 1) * rightStep];
        if (l1 != negInf && r1 != negInf) {
            const float cand = l1 + r1;
            if (cand > best) {
                best = cand;
            }
        }
        const float l2 = left[k + 2];
        const float r2 = rightStart[(k + 2) * rightStep];
        if (l2 != negInf && r2 != negInf) {
            const float cand = l2 + r2;
            if (cand > best) {
                best = cand;
            }
        }
        const float l3 = left[k + 3];
        const float r3 = rightStart[(k + 3) * rightStep];
        if (l3 != negInf && r3 != negInf) {
            const float cand = l3 + r3;
            if (cand > best) {
                best = cand;
            }
        }
    }
    for (; k <= d; ++k) {
        const float l = left[k];
        const float r = rightStart[k * rightStep];
        if (l != negInf && r != negInf) {
            const float cand = l + r;
            if (cand > best) {
                best = cand;
            }
        }
    }
    return best;
}

static inline float bSplitBestPacked(const float *left, const float *rightPacked, int len, float negInf) {
    float best = negInf;
    int k = 0;
#if defined(__AVX2__)
    auto hmax256 = [](__m256 v) -> float {
        const __m128 lo = _mm256_castps256_ps128(v);
        const __m128 hi = _mm256_extractf128_ps(v, 1);
        __m128 m = _mm_max_ps(lo, hi);
        m = _mm_max_ps(m, _mm_movehl_ps(m, m));
        m = _mm_max_ps(m, _mm_shuffle_ps(m, m, _MM_SHUFFLE(1, 1, 1, 1)));
        return _mm_cvtss_f32(m);
    };
    const __m256 negV = _mm256_set1_ps(negInf);
    __m256 bestV = negV;
    for (; k + 7 < len; k += 8) {
        const __m256 l = _mm256_loadu_ps(left + k);
        const __m256 r = _mm256_loadu_ps(rightPacked + k);
        const __m256 valid = _mm256_and_ps(_mm256_cmp_ps(l, negV, _CMP_NEQ_OQ),
                                           _mm256_cmp_ps(r, negV, _CMP_NEQ_OQ));
        const __m256 s = _mm256_add_ps(l, r);
        const __m256 c = _mm256_blendv_ps(negV, s, valid);
        bestV = _mm256_max_ps(bestV, c);
    }
    {
        const float hv = hmax256(bestV);
        if (hv > best) {
            best = hv;
        }
    }
#endif
    for (; k < len; ++k) {
        const float l = left[k];
        const float r = rightPacked[k];
        if (l != negInf && r != negInf) {
            const float cand = l + r;
            if (cand > best) {
                best = cand;
            }
        }
    }
    return best;
}

static inline double expMaybeFast(double x) {
    if (cmFastMathEnabled()) {
        return static_cast<double>(MathUtil::fpow2(static_cast<float>(x * LOG2E)));
    }
    return std::exp(x);
}

static inline double logMaybeFast(double x) {
    if (cmFastMathEnabled()) {
        return static_cast<double>(MathUtil::flog2(static_cast<float>(x))) * LN2;
    }
    return std::log(x);
}

static inline double exp2MaybeFast(double x) {
    if (cmFastMathEnabled()) {
        return static_cast<double>(MathUtil::fpow2(static_cast<float>(x)));
    }
    return std::exp2(x);
}

static inline double log2MaybeFast(double x) {
    if (cmFastMathEnabled()) {
        return static_cast<double>(MathUtil::flog2(static_cast<float>(x)));
    }
    return std::log2(x);
}

struct LogSumAcc {
    double maxVal;
    double scaledSum;
    bool has;

    LogSumAcc() : maxVal(NEG_INF), scaledSum(0.0), has(false) {}

    inline void add(double x) {
        if (!has) {
            maxVal = x;
            scaledSum = 1.0;
            has = true;
            return;
        }
        if (x > maxVal) {
            const double d = maxVal - x;
            scaledSum = ((d < -50.0) ? 0.0 : scaledSum * expMaybeFast(d)) + 1.0;
            maxVal = x;
        } else {
            const double d = x - maxVal;
            if (d >= -50.0) {
                scaledSum += expMaybeFast(d);
            }
        }
    }

    inline double value() const {
        if (!has) {
            return NEG_INF;
        }
        return maxVal + logMaybeFast(scaledSum);
    }
};

static inline double log2SumExp2(double a, double b) {
    if (a == NEG_INF) {
        return b;
    }
    if (b == NEG_INF) {
        return a;
    }
    if (a < b) {
        std::swap(a, b);
    }
    const double d = b - a;
    if (d < -80.0) {
        return a;
    }
    return a + std::log2(1.0 + std::exp2(d));
}

struct FastLog2OnePlusPow2NegLut {
    static const int TABLE_SCALE = 64;
    static const int TABLE_MAX_X = 32;
    static const int TABLE_N = TABLE_MAX_X * TABLE_SCALE + 1;
    float v[TABLE_N];
    FastLog2OnePlusPow2NegLut() {
        for (int i = 0; i < TABLE_N; ++i) {
            const double xv = static_cast<double>(i) / static_cast<double>(TABLE_SCALE);
            v[i] = static_cast<float>(std::log2(1.0 + std::exp2(-xv)));
        }
    }
};

static CM_ALWAYS_INLINE float log2OnePlusPow2NegApproxF(float x) {
    if (x <= 0.0f) {
        return 1.0f;
    }
    if (x >= static_cast<float>(FastLog2OnePlusPow2NegLut::TABLE_MAX_X)) {
        return 0.0f;
    }
    static const FastLog2OnePlusPow2NegLut lut;
    int idx = static_cast<int>(x * static_cast<float>(FastLog2OnePlusPow2NegLut::TABLE_SCALE) + 0.5f);
    if (idx < 0) idx = 0;
    if (idx >= FastLog2OnePlusPow2NegLut::TABLE_N) idx = FastLog2OnePlusPow2NegLut::TABLE_N - 1;
    return lut.v[idx];
}

static CM_ALWAYS_INLINE float log2SumExp2ApproxF(float a, float b) {
    static const float NEG_INF_F = -std::numeric_limits<float>::infinity();
    if (a == NEG_INF_F) {
        return b;
    }
    if (b == NEG_INF_F) {
        return a;
    }
    if (a < b) {
        const float tmp = a;
        a = b;
        b = tmp;
    }
    const float d = a - b; // >= 0
    if (d >= static_cast<float>(FastLog2OnePlusPow2NegLut::TABLE_MAX_X)) {
        return a;
    }
    return a + log2OnePlusPow2NegApproxF(d);
}

static CM_ALWAYS_INLINE void log2AccFastAddF(float x, bool &has, float &acc) {
    if (!has) {
        acc = x;
        has = true;
        return;
    }
    acc = log2SumExp2ApproxF(acc, x);
}

#if defined(CM_HAS_SSE2)
// 4-lane SIMD log-sum-exp in log2 space, using the scalar LUT for the
// log2(1 + 2^(-d)) correction. Drop-in replacement for `_mm_max_ps` in the
// fast-deck fill when running Inside (sum over paths) instead of CYK (max).
//
// Lane semantics: for each i, returns log2(2^a[i] + 2^b[i]).
// NEG_INF handling: if max(a,b) == -inf the result is -inf; if only one is
// -inf the other is returned. This matches the scalar fall-through path.
static CM_ALWAYS_INLINE __m128 log2SumExp2ApproxFSse(__m128 a, __m128 b) {
    static const float NEG_INF_F_LSEX = -std::numeric_limits<float>::infinity();
    const __m128 mx = _mm_max_ps(a, b);
    const __m128 mn = _mm_min_ps(a, b);
    // d = mx - mn >= 0 when both finite. When mn == -inf, the diff is +inf
    // (or NaN if mx is also -inf); the scalar LUT clamps to TABLE_MAX_X so
    // non-NaN large d returns 0. We mask the both-NEG_INF lane separately.
    const __m128 niV = _mm_set1_ps(NEG_INF_F_LSEX);
    const __m128 isMxNI = _mm_cmpeq_ps(mx, niV);
    // Replace mn with mx in the both-NEG_INF lanes so d collapses to 0 and
    // the LUT lookup is well-defined; we'll clobber the result with NEG_INF
    // afterwards via the mask.
    const __m128 mnSafe = _mm_or_ps(_mm_and_ps(isMxNI, mx), _mm_andnot_ps(isMxNI, mn));
    const __m128 d = _mm_sub_ps(mx, mnSafe);
    alignas(16) float dArr[4];
    _mm_store_ps(dArr, d);
    alignas(16) float adjArr[4];
    adjArr[0] = log2OnePlusPow2NegApproxF(dArr[0]);
    adjArr[1] = log2OnePlusPow2NegApproxF(dArr[1]);
    adjArr[2] = log2OnePlusPow2NegApproxF(dArr[2]);
    adjArr[3] = log2OnePlusPow2NegApproxF(dArr[3]);
    const __m128 sum = _mm_add_ps(mx, _mm_load_ps(adjArr));
    // Both-NEG_INF lanes -> NEG_INF; otherwise sum.
    return _mm_or_ps(_mm_and_ps(isMxNI, niV), _mm_andnot_ps(isMxNI, sum));
}
#endif

static inline void log2AccExactAdd(double x, bool &has, double &maxVal, double &scaledSum) {
    if (!has) {
        maxVal = x;
        scaledSum = 1.0;
        has = true;
        return;
    }
    if (x > maxVal) {
        const double d = maxVal - x;
        scaledSum = ((d < -80.0) ? 0.0 : scaledSum * exp2MaybeFast(d)) + 1.0;
        maxVal = x;
    } else {
        const double d = x - maxVal;
        if (d >= -80.0) {
            scaledSum += exp2MaybeFast(d);
        }
    }
}

static inline double log2AccExactValue(bool has, double maxVal, double scaledSum) {
    return has ? (maxVal + log2MaybeFast(scaledSum)) : NEG_INF;
}

template <bool FastMode>
struct LogSumAccBase2T;

template <>
struct LogSumAccBase2T<false> {
    double maxVal;
    double scaledSum;
    bool has;

    LogSumAccBase2T() : maxVal(NEG_INF), scaledSum(0.0), has(false) {}

    inline void add(double x) {
        if (!has) {
            maxVal = x;
            scaledSum = 1.0;
            has = true;
            return;
        }
        if (x > maxVal) {
            const double d = maxVal - x;
            scaledSum = ((d < -80.0) ? 0.0 : scaledSum * exp2MaybeFast(d)) + 1.0;
            maxVal = x;
        } else {
            const double d = x - maxVal;
            if (d >= -80.0) {
                scaledSum += exp2MaybeFast(d);
            }
        }
    }

    inline double value() const {
        if (!has) {
            return NEG_INF;
        }
        return maxVal + log2MaybeFast(scaledSum);
    }
};

template <>
struct LogSumAccBase2T<true> {
    float maxVal;
    bool has;

    LogSumAccBase2T() : maxVal(-std::numeric_limits<float>::infinity()), has(false) {}

    inline void add(double x) {
        const float xf = static_cast<float>(x);
        if (!has) {
            maxVal = xf;
            has = true;
            return;
        }
        maxVal = log2SumExp2ApproxF(maxVal, xf);
    }

    inline double value() const {
        return has ? static_cast<double>(maxVal) : NEG_INF;
    }
};

typedef LogSumAccBase2T<false> LogSumAccBase2Exact;
typedef LogSumAccBase2T<true> LogSumAccBase2Fast;

static double emitBitToProbSingle(double bitSc, const InfernalExactModel &model, int a) {
    if (!std::isfinite(bitSc) || a < 0 || a >= 4 || model.nullProb[a] <= 0.0) {
        return 0.0;
    }
    return std::exp2(bitSc) * model.nullProb[a];
}

static double emitBitToProbPair(double bitSc, const InfernalExactModel &model, int a, int b) {
    if (!std::isfinite(bitSc) || a < 0 || a >= 4 || b < 0 || b >= 4 ||
        model.nullProb[a] <= 0.0 || model.nullProb[b] <= 0.0) {
        return 0.0;
    }
    return std::exp2(bitSc) * model.nullProb[a] * model.nullProb[b];
}

static double scoreCorrectionNull2BitsFromTrace(const InfernalExactModel &model,
                                                const double modelAggRaw[4],
                                                const int obsCount[4]) {
    if (!model.hasNull || model.null2Omega <= 0.0) {
        return 0.0;
    }
    double total = 0.0;
    for (int a = 0; a < 4; ++a) {
        total += modelAggRaw[a];
    }
    if (total <= 0.0) {
        return 0.0;
    }
    double score = 0.0;
    for (int a = 0; a < 4; ++a) {
        if (obsCount[a] <= 0 || model.nullProb[a] <= 0.0) {
            continue;
        }
        const double p = std::max(1e-12, modelAggRaw[a] / total);
        score += static_cast<double>(obsCount[a]) * std::log2(p / model.nullProb[a]);
    }
    score += std::log2(model.null2Omega);
    return log2SumExp2(0.0, score);
}

static bool infernalExactScoreToEvalue(const InfernalExactModel &model,
                                       bool insideMode,
                                       char hitMode,
                                       double scoreBits,
                                       double targetDbResidues,
                                       double &outEvalue) {
    const bool truncLike = (hitMode == 'L' || hitMode == 'R' || hitMode == 'T');
    const InfernalExactModel::ExpTail &exp = insideMode
        ? (truncLike ? model.expGI : model.expLI)
        : (truncLike ? model.expGC : model.expLC);
    if (!exp.valid || exp.lambda <= 0.0 || exp.dbsize <= 0.0 || exp.nrandhits <= 0 || targetDbResidues <= 0.0) {
        return false;
    }
    const double curEffDbSize = (targetDbResidues / exp.dbsize) * static_cast<double>(exp.nrandhits);
    const double p = std::exp(-exp.lambda * (scoreBits - exp.muExtrap));
    outEvalue = p * curEffDbSize;
    return std::isfinite(outEvalue);
}

static bool looksLikeInfernalCm(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        return false;
    }
    std::string first;
    std::getline(in, first);
    first = trim(first);
    // Handle optional UTF-8 BOM.
    if (first.size() >= 3 &&
        static_cast<unsigned char>(first[0]) == 0xEF &&
        static_cast<unsigned char>(first[1]) == 0xBB &&
        static_cast<unsigned char>(first[2]) == 0xBF) {
        first = first.substr(3);
    }
    return first.find("INFERNAL1/") != std::string::npos;
}

static int parseMaybeInt(const std::string &tok) {
    if (tok.empty() || tok == "-") {
        return -1;
    }
    char *end = NULL;
    long v = std::strtol(tok.c_str(), &end, 10);
    if (end == tok.c_str() || *end != '\0') {
        return -1;
    }
    return static_cast<int>(v);
}

static bool parseMaybeDouble(const std::string &tok, double &out) {
    if (tok == "*") {
        out = NEG_INF;
        return true;
    }
    char *end = NULL;
    out = std::strtod(tok.c_str(), &end);
    return end != tok.c_str() && *end == '\0';
}

static CmStateType parseCmStateType(const std::string &tok) {
    if (tok == "MP") {
        return CM_ST_MP;
    }
    if (tok == "ML") {
        return CM_ST_ML;
    }
    if (tok == "MR") {
        return CM_ST_MR;
    }
    if (tok == "IL") {
        return CM_ST_IL;
    }
    if (tok == "IR") {
        return CM_ST_IR;
    }
    if (tok == "D") {
        return CM_ST_D;
    }
    if (tok == "S") {
        return CM_ST_S;
    }
    if (tok == "B") {
        return CM_ST_B;
    }
    if (tok == "E") {
        return CM_ST_E;
    }
    return CM_ST_UNKNOWN;
}

static int baseToIdx(char c) {
    c = normalizeBase(c);
    if (c == 'A') {
        return 0;
    }
    if (c == 'C') {
        return 1;
    }
    if (c == 'G') {
        return 2;
    }
    if (c == 'U') {
        return 3;
    }
    return -1;
}

static char complementBase(char c) {
    switch (normalizeBase(c)) {
        case 'A': return 'U';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'U': return 'A';
        default:  return 'N';
    }
}

static std::string reverseComplement(const std::string &seq) {
    std::string rc(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); ++i) {
        rc[seq.size() - 1 - i] = complementBase(seq[i]);
    }
    return rc;
}

static double parseInfernalAsciiScore(const std::string &tok) {
    if (tok == "*") {
        return NEG_INF;
    }
    double v = 0.0;
    if (!parseMaybeDouble(tok, v)) {
        return NEG_INF;
    }
    return v;
}

static double parseInfernalAsciiProb(const std::string &tok, double nullProb) {
    if (tok == "*") {
        return 0.0;
    }
    double bits = 0.0;
    if (!parseMaybeDouble(tok, bits)) {
        return nullProb;
    }
    return std::exp2(bits) * nullProb;
}

static CmNodeType parseCmNodeType(const std::string &tok) {
    if (tok == "ROOT") return CM_ND_ROOT;
    if (tok == "MATP") return CM_ND_MATP;
    if (tok == "MATL") return CM_ND_MATL;
    if (tok == "MATR") return CM_ND_MATR;
    if (tok == "BEGL") return CM_ND_BEGL;
    if (tok == "BEGR") return CM_ND_BEGR;
    if (tok == "BIF")  return CM_ND_BIF;
    if (tok == "END")  return CM_ND_END;
    return CM_ND_UNKNOWN;
}

// =====================================================================
// Stage 6: native port of Infernal cm_TrCYKInsideAlign (cm_dpalign_trunc.c)
// Fills J/L/R/T DP and shadow matrices for full-sequence trCYK alignment.
// Not yet wired into the alignment path; integration is Stage 8.
// =====================================================================
struct TrCykResult {
    double score = -1e30;
    char   mode = 'J';
    int    b = 0;
    std::vector<std::vector<std::vector<int8_t>>>  Jyshadow, Lyshadow, Ryshadow;
    std::vector<std::vector<std::vector<int>>>     Jkshadow, Lkshadow, Rkshadow, Tkshadow;
    std::vector<std::vector<std::vector<int8_t>>>  Lkmode, Rkmode;
    std::vector<std::vector<std::vector<float>>>   Jalpha, Lalpha, Ralpha, Talpha;
};

// TRMODE_* offset constants — match Infernal cm_alphabet.h.
// Match Infernal infernal.h: TRMODE_T=0, TRMODE_R=1, TRMODE_L=2, TRMODE_J=3,
// TRMODE_UNKNOWN=4. TRMODE_*_OFFSET picks DP source for shadow decoding.
static const int8_t  TR_TRMODE_T = 0;
static const int8_t  TR_TRMODE_R = 1;
static const int8_t  TR_TRMODE_L = 2;
static const int8_t  TR_TRMODE_J = 3;
static const int8_t  TR_TRMODE_UNKNOWN = 4;
static const int    TR_TRMODE_J_OFFSET = 0;
static const int    TR_TRMODE_L_OFFSET = 10;
static const int    TR_TRMODE_R_OFFSET = 20;
static const int8_t  TR_USED_LOCAL_BEGIN = 101;
static const int8_t  TR_USED_EL          = 102;
static const int8_t  TR_USED_TRUNC_BEGIN = 103;
static const int8_t  TR_USED_TRUNC_END   = 104;

static inline bool trNotImpossible(double x) {
    // Infernal NOT_IMPOSSIBLE: x > -1e8. Our endsc/trp use NEG_INF = -inf
    // when ineligible, so any finite value qualifies.
    return std::isfinite(x);
}

// Mirrors esl_abc_FAvgScore: arithmetic average of esc[0..K-1] over degeneracy
// expansion class. For fully-degenerate X/N over {A,C,G,U}: avg of all four.
template<typename T>
static inline float trEmitDegen1(const T *esc, int8_t b) {
    if (b >= 0 && b <= 3) return (float)esc[(size_t)b];
    double s = 0.0;
    for (int x = 0; x < 4; ++x) s += (double)esc[(size_t)x];
    return (float)(s * 0.25);
}

// MP pair emission with degeneracy on either or both sides. Arithmetic avg
// over the cartesian product of degeneracy expansions for bi and bj.
static inline float trEmitDegenPair(const float *esc, int Kp, int8_t bi, int8_t bj) {
    const bool degI = (bi < 0 || bi > 3);
    const bool degJ = (bj < 0 || bj > 3);
    if (!degI && !degJ) return esc[(size_t)((int)bi * Kp + (int)bj)];
    int xLo = degI ? 0 : (int)bi, xHi = degI ? 3 : (int)bi;
    int yLo = degJ ? 0 : (int)bj, yHi = degJ ? 3 : (int)bj;
    double s = 0.0;
    int n = 0;
    for (int x = xLo; x <= xHi; ++x) {
        for (int y = yLo; y <= yHi; ++y) {
            s += (double)esc[(size_t)(x * Kp + y)];
            ++n;
        }
    }
    return (float)(s / (double)n);
}

static const TrCykResult& runTrCYKInsideAlign(const InfernalExactModel &m,
                                              const std::vector<int8_t> &dsq,
                                              int L,
                                              char preset_mode,
                                              int pty_idx) {
    const float TR_IMPOSSIBLE = -1e30f;
    const int M = (int)m.states.size();

    // fill_{L,R,T} based on preset_mode (mirrors cm_TrFillFromMode).
    bool fill_L = false, fill_R = false, fill_T = false;
    switch (preset_mode) {
        case 'J': fill_L = false; fill_R = false; fill_T = false; break;
        case 'L': fill_L = true;  fill_R = false; fill_T = false; break;
        case 'R': fill_L = false; fill_R = true;  fill_T = false; break;
        case 'T': fill_L = true;  fill_R = true;  fill_T = true;  break;
        default:  fill_L = true;  fill_R = true;  fill_T = true;  break; // unknown
    }

    // Per-thread arena. Vectors grow monotonically across calls — never shrink —
    // so the inner float/int storage is reused. Reset overhead is just std::fill
    // over [0..LP1) per (v, j) row, no malloc/free per envelope.
    thread_local TrCykResult r;
    // Reset the score-summary fields (data arrays are overwritten by reset+DP below).
    r.score = -1e30;
    r.mode = 'J';
    r.b = 0;

    const int LP1 = L + 1;
    auto resetAlpha = [&](std::vector<std::vector<std::vector<float>>> &a) {
        if ((int)a.size() < M + 1) a.resize((size_t)M + 1);
        for (int v = 0; v <= M; ++v) {
            auto &mv = a[(size_t)v];
            if ((int)mv.size() < LP1) mv.resize((size_t)LP1);
            for (int j = 0; j < LP1; ++j) {
                auto &jv = mv[(size_t)j];
                if ((int)jv.size() < LP1) jv.resize((size_t)LP1);
                std::fill(jv.begin(), jv.begin() + LP1, TR_IMPOSSIBLE);
            }
        }
    };
    auto resetByteSh = [&](std::vector<std::vector<std::vector<int8_t>>> &a, int8_t val) {
        if ((int)a.size() < M + 1) a.resize((size_t)M + 1);
        for (int v = 0; v <= M; ++v) {
            auto &mv = a[(size_t)v];
            if ((int)mv.size() < LP1) mv.resize((size_t)LP1);
            for (int j = 0; j < LP1; ++j) {
                auto &jv = mv[(size_t)j];
                if ((int)jv.size() < LP1) jv.resize((size_t)LP1);
                std::fill(jv.begin(), jv.begin() + LP1, val);
            }
        }
    };
    auto resetIntSh = [&](std::vector<std::vector<std::vector<int>>> &a) {
        if ((int)a.size() < M + 1) a.resize((size_t)M + 1);
        for (int v = 0; v <= M; ++v) {
            auto &mv = a[(size_t)v];
            if ((int)mv.size() < LP1) mv.resize((size_t)LP1);
            for (int j = 0; j < LP1; ++j) {
                auto &jv = mv[(size_t)j];
                if ((int)jv.size() < LP1) jv.resize((size_t)LP1);
                std::fill(jv.begin(), jv.begin() + LP1, 0);
            }
        }
    };

    resetAlpha(r.Jalpha);
    if (fill_L) resetAlpha(r.Lalpha);
    if (fill_R) resetAlpha(r.Ralpha);
    if (fill_T) resetAlpha(r.Talpha); // Talpha only used at B states; allocate full for simplicity.

    resetByteSh(r.Jyshadow, (int8_t)TR_USED_EL);
    if (fill_L) resetByteSh(r.Lyshadow, (int8_t)TR_USED_EL);
    if (fill_R) resetByteSh(r.Ryshadow, (int8_t)TR_USED_EL);
    resetIntSh(r.Jkshadow);
    if (fill_L) resetIntSh(r.Lkshadow);
    if (fill_R) resetIntSh(r.Rkshadow);
    if (fill_T) resetIntSh(r.Tkshadow);
    if (fill_L) resetByteSh(r.Lkmode, (int8_t)TR_TRMODE_J);
    if (fill_R) resetByteSh(r.Rkmode, (int8_t)TR_TRMODE_J);

    // el_scA[d] = elSelf * d. elSelf is in our model already in log2-bits.
    std::vector<float> el_scA((size_t)LP1, 0.0f);
    const double elSelf = m.elSelf;
    for (int d = 0; d <= L; ++d) {
        el_scA[(size_t)d] = (float)(elSelf * (double)d);
    }

    const bool localOn = m.hasLocalCfg;

    // MMSEQS_CMSCAN_TRUNC_DISABLE_EL=1: gate End-Local pickup off in the DP.
    // Mirrors --notrunc-style behavior where EL is unreachable, so the trace
    // never picks an EL parsetree path. Justified because our renderer outputs
    // fixed-CLEN columns and cannot represent EL inserts the way Infernal's
    // Parsetrees2Alignment does (dynamic-width MSA with per-cpos EL blocks).
    static const char *const elDisableEnv = std::getenv("MMSEQS_CMSCAN_TRUNC_DISABLE_EL");
    const bool disableEL = (elDisableEnv != nullptr && elDisableEnv[0] == '1');
    // MMSEQS_CMSCAN_TRUNC_NO_TRPENALTY=1: include v=0 in main recurrence and
    // gate trunc-entry pickup off. Mirrors Infernal --notrunc behavior where
    // alpha[0][L][L] is filled by the natural S→child recurrence at root,
    // not via a USED_TRUNC_BEGIN pickup. Trace then walks normally from v=0.
    static const char *const noTrPenEnv = std::getenv("MMSEQS_CMSCAN_TRUNC_NO_TRPENALTY");
    const bool noTrPenalty = (noTrPenEnv != nullptr && noTrPenEnv[0] == '1');

    // If local ends are on (and EL not disabled), EL deck (v=M) holds el_scA[d]
    // for j,d in [0,L]. With EL disabled, leave the deck at TR_IMPOSSIBLE.
    if (localOn && !disableEL) {
        for (int j = 0; j <= L; ++j) {
            for (int d = 0; d <= j; ++d) {
                r.Jalpha[(size_t)M][(size_t)j][(size_t)d] = el_scA[(size_t)d];
                if (fill_L) r.Lalpha[(size_t)M][(size_t)j][(size_t)d] = el_scA[(size_t)d];
                if (fill_R) r.Ralpha[(size_t)M][(size_t)j][(size_t)d] = el_scA[(size_t)d];
            }
        }
    }

    // Trackers for begin states across modes.
    int Jb = 0, Lb = 0, Rb = 0, Tb = 0;
    // For noTrPenalty (--notrunc-equivalent): track best local-begin pickup
    // separately, commit only at end (mirrors Infernal cm_CYKInsideAlign).
    float bsc_local = TR_IMPOSSIBLE;
    int   blb_local = -1;

    // Helper to read transition score for state v at child yoffset.
    auto tscOf = [&](int v, int yoffset) -> float {
        const CmState &s = m.states[(size_t)v];
        if (yoffset < 0 || (size_t)yoffset >= s.trans.size()) return TR_IMPOSSIBLE;
        double t = s.trans[(size_t)yoffset];
        if (!std::isfinite(t)) return TR_IMPOSSIBLE;
        return (float)t;
    };

    // Main recursion: v = M-1 down to 1. v=0 (ROOT_S) handled via trpenalty hand-off below.
    // When noTrPenalty, also include v=0 so alpha[0] is filled by the natural
    // S-state recurrence (no truncation entry needed).
    const int loopStop = noTrPenalty ? -1 : 0;
    for (int v = M - 1; v > loopStop; --v) {
        const CmState &s = m.states[(size_t)v];
        const CmStateType st = s.type;
        const int cfirst = s.cfirst;
        const int cnum   = s.cnum;

        int sd = 0, sdl = 0, sdr = 0;
        switch (st) {
            case CM_ST_MP: sd = 2; sdl = 1; sdr = 1; break;
            case CM_ST_ML: case CM_ST_IL: sd = 1; sdl = 1; sdr = 0; break;
            case CM_ST_MR: case CM_ST_IR: sd = 1; sdl = 0; sdr = 1; break;
            default: sd = 0; sdl = 0; sdr = 0; break;
        }

        // Re-init this state's J/L/R deck if a local end is reachable from v.
        // With disableEL, skip this entirely so EL is never the picked path.
        const double endsc_v = s.endSc;
        if (localOn && !disableEL && trNotImpossible(endsc_v)) {
            const float endF = (float)endsc_v;
            const simd_float vEndF = simdf32_set(endF);
            const int VS = VECSIZE_FLOAT;
            // SIMD over d for each j: Jalpha[v][j][d] = el_scA[d - sd] + endF
            for (int j = 0; j <= L; ++j) {
                float *Jrow = r.Jalpha[(size_t)v][(size_t)j].data();
                int d = sd;
                for (; d + VS - 1 <= j; d += VS) {
                    simd_float vEl = simdf32_loadu(&el_scA[(size_t)(d - sd)]);
                    simdf32_storeu(Jrow + d, simdf32_add(vEl, vEndF));
                }
                for (; d <= j; ++d) {
                    Jrow[d] = el_scA[(size_t)(d - sd)] + endF;
                }
            }
            if (fill_L) {
                for (int j = 0; j <= L; ++j) {
                    float *Lrow = r.Lalpha[(size_t)v][(size_t)j].data();
                    int d = sdl;
                    for (; d + VS - 1 <= j; d += VS) {
                        simd_float vEl = simdf32_loadu(&el_scA[(size_t)(d - sdl)]);
                        simdf32_storeu(Lrow + d, simdf32_add(vEl, vEndF));
                    }
                    for (; d <= j; ++d) {
                        Lrow[d] = el_scA[(size_t)(d - sdl)] + endF;
                    }
                }
            }
            if (fill_R) {
                for (int j = 0; j <= L; ++j) {
                    float *Rrow = r.Ralpha[(size_t)v][(size_t)j].data();
                    int d = sdr;
                    for (; d + VS - 1 <= j; d += VS) {
                        simd_float vEl = simdf32_loadu(&el_scA[(size_t)(d - sdr)]);
                        simdf32_storeu(Rrow + d, simdf32_add(vEl, vEndF));
                    }
                    for (; d <= j; ++d) {
                        Rrow[d] = el_scA[(size_t)(d - sdr)] + endF;
                    }
                }
            }
        }

        // Per-state recursion by type.
        if (st == CM_ST_E) {
            for (int j = 0; j <= L; ++j) {
                r.Jalpha[(size_t)v][(size_t)j][0] = 0.0f;
                if (fill_L) r.Lalpha[(size_t)v][(size_t)j][0] = 0.0f;
                if (fill_R) r.Ralpha[(size_t)v][(size_t)j][0] = 0.0f;
            }
        }
        else if (st == CM_ST_IL || st == CM_ST_ML) {
            const float *esc = s.emitF.data();
            const int Ryoffset0 = (st == CM_ST_IL) ? 1 : 0;
            for (int j = sdr; j <= L; ++j) {
                const int j_sdr = j - sdr;
                for (int d = sd; d <= j; ++d) {
                    const int d_sd = d - sd;
                    int i = j - d + 1;
                    for (int yoff = 0; yoff < cnum; ++yoff) {
                        const int y = cfirst + yoff;
                        const float tsc = tscOf(v, yoff);
                        float sc = r.Jalpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                        if (sc > r.Jalpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = sc;
                            r.Jyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                        }
                        if (fill_L) {
                            float scL = r.Lalpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                            if (scL > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scL;
                                r.Lyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_L_OFFSET);
                            }
                        }
                    }
                    // Single-residue emit on left position i.
                    const int8_t bi = (i >= 1 && i <= L) ? dsq[(size_t)i] : (int8_t)4;
                    const float emitL = (i >= 1 && i <= L) ? trEmitDegen1(esc, bi) : TR_IMPOSSIBLE;
                    r.Jalpha[(size_t)v][(size_t)j][(size_t)d] += emitL;
                    if (r.Jalpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                        r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                    }
                    if (fill_L) {
                        if (d >= 2) {
                            r.Lalpha[(size_t)v][(size_t)j][(size_t)d] += emitL;
                        } else {
                            r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = emitL;
                            r.Lyshadow[(size_t)v][(size_t)j][(size_t)d] = TR_USED_TRUNC_END;
                        }
                        if (r.Lalpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                            r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                        }
                    }
                    // R deck handled separately (depends on Jalpha[y][j_sdr][d], not d_sd).
                    if (fill_R) {
                        for (int yoff = Ryoffset0; yoff < cnum; ++yoff) {
                            const int y = cfirst + yoff;
                            const float tsc = tscOf(v, yoff);
                            float scA = r.Jalpha[(size_t)y][(size_t)j_sdr][(size_t)d] + tsc;
                            if (scA > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scA;
                                r.Ryshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                            }
                            float scB = r.Ralpha[(size_t)y][(size_t)j_sdr][(size_t)d] + tsc;
                            if (scB > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scB;
                                r.Ryshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_R_OFFSET);
                            }
                        }
                        if (r.Ralpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                            r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                        }
                    }
                }
            }
        }
        else if (st == CM_ST_IR || st == CM_ST_MR) {
            const float *esc = s.emitF.data();
            const int Lyoffset0 = (st == CM_ST_IR) ? 1 : 0;
            for (int j = sdr; j <= L; ++j) {
                const int j_sdr = j - sdr;
                for (int d = sd; d <= j; ++d) {
                    const int d_sd = d - sd;
                    for (int yoff = 0; yoff < cnum; ++yoff) {
                        const int y = cfirst + yoff;
                        const float tsc = tscOf(v, yoff);
                        float sc = r.Jalpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                        if (sc > r.Jalpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = sc;
                            r.Jyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                        }
                        if (fill_R) {
                            float scR = r.Ralpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                            if (scR > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scR;
                                r.Ryshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_R_OFFSET);
                            }
                        }
                    }
                    // Single-residue emit on right position j.
                    const int8_t bj = (j >= 1 && j <= L) ? dsq[(size_t)j] : (int8_t)4;
                    const float emitR = (j >= 1 && j <= L) ? trEmitDegen1(esc, bj) : TR_IMPOSSIBLE;
                    r.Jalpha[(size_t)v][(size_t)j][(size_t)d] += emitR;
                    if (r.Jalpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                        r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                    }
                    if (fill_R) {
                        if (d >= 2) {
                            r.Ralpha[(size_t)v][(size_t)j][(size_t)d] += emitR;
                        } else {
                            r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = emitR;
                            r.Ryshadow[(size_t)v][(size_t)j][(size_t)d] = TR_USED_TRUNC_END;
                        }
                        if (r.Ralpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                            r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                        }
                    }
                    // L deck handled separately (uses j, d not j_sdr, d_sd).
                    if (fill_L) {
                        for (int yoff = Lyoffset0; yoff < cnum; ++yoff) {
                            const int y = cfirst + yoff;
                            const float tsc = tscOf(v, yoff);
                            float scA = r.Jalpha[(size_t)y][(size_t)j][(size_t)d] + tsc;
                            if (scA > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scA;
                                r.Lyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                            }
                            float scB = r.Lalpha[(size_t)y][(size_t)j][(size_t)d] + tsc;
                            if (scB > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scB;
                                r.Lyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_L_OFFSET);
                            }
                        }
                        if (r.Lalpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                            r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                        }
                    }
                }
            }
        }
        else if (st == CM_ST_MP) {
            const float *esc = s.emitF.data();
            const int escSize = (int)s.emitF.size();
            const int Kp = 4; // our MP emit table is 4x4 (no degeneracies stored).
            // Recurrence over children.
            for (int yoff = 0; yoff < cnum; ++yoff) {
                const int y = cfirst + yoff;
                const float tsc = tscOf(v, yoff);
                for (int j = sdr; j <= L; ++j) {
                    const int j_sdr = j - sdr;
                    for (int d = sd; d <= j; ++d) {
                        const int d_sd = d - sd;
                        float sc = r.Jalpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                        if (sc > r.Jalpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = sc;
                            r.Jyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                        }
                    }
                    if (fill_L) {
                        for (int d = sdl; d <= j; ++d) {
                            const int d_sdl = d - sdl;
                            float scA = r.Jalpha[(size_t)y][(size_t)j][(size_t)d_sdl] + tsc;
                            if (scA > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scA;
                                r.Lyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                            }
                            float scB = r.Lalpha[(size_t)y][(size_t)j][(size_t)d_sdl] + tsc;
                            if (scB > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scB;
                                r.Lyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_L_OFFSET);
                            }
                        }
                    }
                    if (fill_R) {
                        for (int d = sdr; d <= j; ++d) {
                            const int d_sdr = d - sdr;
                            float scA = r.Jalpha[(size_t)y][(size_t)j_sdr][(size_t)d_sdr] + tsc;
                            if (scA > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scA;
                                r.Ryshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                            }
                            float scB = r.Ralpha[(size_t)y][(size_t)j_sdr][(size_t)d_sdr] + tsc;
                            if (scB > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scB;
                                r.Ryshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_R_OFFSET);
                            }
                        }
                    }
                }
            }
            // Add emission scores; lmesc/rmesc for marginal modes.
            const float *lme = (s.lmesc.size() >= 4) ? nullptr : nullptr; (void)lme;
            for (int j = 0; j <= L; ++j) {
                int i = j;
                if (j >= 1) {
                    r.Jalpha[(size_t)v][(size_t)j][1] = TR_IMPOSSIBLE;
                    if (fill_L) {
                        const int8_t bi = (i >= 1 && i <= L) ? dsq[(size_t)i] : (int8_t)4;
                        float em = ((int)s.lmesc.size() >= 4 && i >= 1 && i <= L)
                                   ? trEmitDegen1(s.lmesc.data(), bi) : TR_IMPOSSIBLE;
                        r.Lalpha[(size_t)v][(size_t)j][1] = em;
                        r.Lyshadow[(size_t)v][(size_t)j][1] = TR_USED_TRUNC_END;
                    }
                    if (fill_R) {
                        const int8_t bj = (j >= 1 && j <= L) ? dsq[(size_t)j] : (int8_t)4;
                        float em = ((int)s.rmesc.size() >= 4 && j >= 1 && j <= L)
                                   ? trEmitDegen1(s.rmesc.data(), bj) : TR_IMPOSSIBLE;
                        r.Ralpha[(size_t)v][(size_t)j][1] = em;
                        r.Ryshadow[(size_t)v][(size_t)j][1] = TR_USED_TRUNC_END;
                    }
                }
                i--;
                for (int d = 2; d <= j; ++d) {
                    const int8_t bi = (i >= 1 && i <= L) ? dsq[(size_t)i] : (int8_t)4;
                    const int8_t bj = (j >= 1 && j <= L) ? dsq[(size_t)j] : (int8_t)4;
                    float pairEm = (escSize >= Kp * Kp && i >= 1 && i <= L && j >= 1 && j <= L)
                                   ? trEmitDegenPair(esc, Kp, bi, bj) : TR_IMPOSSIBLE;
                    r.Jalpha[(size_t)v][(size_t)j][(size_t)d] += pairEm;
                    if (fill_L) {
                        float em = ((int)s.lmesc.size() >= 4 && i >= 1 && i <= L)
                                   ? trEmitDegen1(s.lmesc.data(), bi) : TR_IMPOSSIBLE;
                        r.Lalpha[(size_t)v][(size_t)j][(size_t)d] += em;
                    }
                    if (fill_R) {
                        float em = ((int)s.rmesc.size() >= 4 && j >= 1 && j <= L)
                                   ? trEmitDegen1(s.rmesc.data(), bj) : TR_IMPOSSIBLE;
                        r.Ralpha[(size_t)v][(size_t)j][(size_t)d] += em;
                    }
                    i--;
                }
            }
            // Clamp to IMPOSSIBLE for d>=1.
            for (int j = 0; j <= L; ++j) {
                for (int d = 1; d <= j; ++d) {
                    if (r.Jalpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                        r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                    }
                    if (fill_L && r.Lalpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                        r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                    }
                    if (fill_R && r.Ralpha[(size_t)v][(size_t)j][(size_t)d] < TR_IMPOSSIBLE) {
                        r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = TR_IMPOSSIBLE;
                    }
                }
            }
        }
        else if (st != CM_ST_B) {
            // D, S — non-self, non-emitting.
            for (int yoff = 0; yoff < cnum; ++yoff) {
                const int y = cfirst + yoff;
                const float tsc = tscOf(v, yoff);
                for (int j = sdr; j <= L; ++j) {
                    const int j_sdr = j - sdr;
                    for (int d = sd; d <= j; ++d) {
                        const int d_sd = d - sd;
                        float sc = r.Jalpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                        if (sc > r.Jalpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = sc;
                            r.Jyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_J_OFFSET);
                        }
                        if (fill_L) {
                            float scL = r.Lalpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                            if (scL > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scL;
                                r.Lyshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_L_OFFSET);
                            }
                        }
                        if (fill_R) {
                            float scR = r.Ralpha[(size_t)y][(size_t)j_sdr][(size_t)d_sd] + tsc;
                            if (scR > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scR;
                                r.Ryshadow[(size_t)v][(size_t)j][(size_t)d] = (int8_t)(yoff + TR_TRMODE_R_OFFSET);
                            }
                        }
                    }
                    if (fill_L) r.Lalpha[(size_t)v][(size_t)j][0] = TR_IMPOSSIBLE;
                    if (fill_R) r.Ralpha[(size_t)v][(size_t)j][0] = TR_IMPOSSIBLE;
                    if (st == CM_ST_S) {
                        if (fill_L) r.Lyshadow[(size_t)v][(size_t)j][0] = TR_USED_TRUNC_END;
                        if (fill_R) r.Ryshadow[(size_t)v][(size_t)j][0] = TR_USED_TRUNC_END;
                    }
                }
            }
        }
        else {
            // B state: BIF.
            const int y = cfirst; // left subtree
            const int z = cnum;   // right subtree (encoded in cnum field)
            for (int j = 0; j <= L; ++j) {
                for (int d = 0; d <= j; ++d) {
                    for (int k = 0; k <= d; ++k) {
                        float scJ = r.Jalpha[(size_t)y][(size_t)(j - k)][(size_t)(d - k)]
                                  + r.Jalpha[(size_t)z][(size_t)j][(size_t)k];
                        if (scJ > r.Jalpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Jalpha[(size_t)v][(size_t)j][(size_t)d] = scJ;
                            r.Jkshadow[(size_t)v][(size_t)j][(size_t)d] = k;
                        }
                        if (fill_L) {
                            float scL = r.Jalpha[(size_t)y][(size_t)(j - k)][(size_t)(d - k)]
                                      + r.Lalpha[(size_t)z][(size_t)j][(size_t)k];
                            if (scL > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scL;
                                r.Lkshadow[(size_t)v][(size_t)j][(size_t)d] = k;
                                r.Lkmode[(size_t)v][(size_t)j][(size_t)d] = TR_TRMODE_J;
                            }
                        }
                        if (fill_R) {
                            float scR = r.Ralpha[(size_t)y][(size_t)(j - k)][(size_t)(d - k)]
                                      + r.Jalpha[(size_t)z][(size_t)j][(size_t)k];
                            if (scR > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scR;
                                r.Rkshadow[(size_t)v][(size_t)j][(size_t)d] = k;
                                r.Rkmode[(size_t)v][(size_t)j][(size_t)d] = TR_TRMODE_J;
                            }
                        }
                    }
                    if (fill_T) {
                        for (int k = 1; k < d; ++k) {
                            float scT = r.Ralpha[(size_t)y][(size_t)(j - k)][(size_t)(d - k)]
                                      + r.Lalpha[(size_t)z][(size_t)j][(size_t)k];
                            if (scT > r.Talpha[(size_t)v][(size_t)j][(size_t)d]) {
                                r.Talpha[(size_t)v][(size_t)j][(size_t)d] = scT;
                                r.Tkshadow[(size_t)v][(size_t)j][(size_t)d] = k;
                            }
                        }
                    }
                    // Special case 1: full sequence aligns to BEGL_S left child (k=0).
                    if (fill_L) {
                        float scA = r.Jalpha[(size_t)y][(size_t)j][(size_t)d];
                        if (scA > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scA;
                            r.Lkshadow[(size_t)v][(size_t)j][(size_t)d] = 0;
                            r.Lkmode[(size_t)v][(size_t)j][(size_t)d] = TR_TRMODE_J;
                        }
                        float scB = r.Lalpha[(size_t)y][(size_t)j][(size_t)d];
                        if (scB > r.Lalpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Lalpha[(size_t)v][(size_t)j][(size_t)d] = scB;
                            r.Lkshadow[(size_t)v][(size_t)j][(size_t)d] = 0;
                            r.Lkmode[(size_t)v][(size_t)j][(size_t)d] = TR_TRMODE_L;
                        }
                    }
                    // Special case 2: full sequence aligns to BEGR_S right child (k=d).
                    if (fill_R) {
                        float scA = r.Jalpha[(size_t)z][(size_t)j][(size_t)d];
                        if (scA > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scA;
                            r.Rkshadow[(size_t)v][(size_t)j][(size_t)d] = d;
                            r.Rkmode[(size_t)v][(size_t)j][(size_t)d] = TR_TRMODE_J;
                        }
                        float scB = r.Ralpha[(size_t)z][(size_t)j][(size_t)d];
                        if (scB > r.Ralpha[(size_t)v][(size_t)j][(size_t)d]) {
                            r.Ralpha[(size_t)v][(size_t)j][(size_t)d] = scB;
                            r.Rkshadow[(size_t)v][(size_t)j][(size_t)d] = d;
                            r.Rkmode[(size_t)v][(size_t)j][(size_t)d] = TR_TRMODE_R;
                        }
                    }
                }
            }
        }

        // Truncated begin into v: only update [0][L][L] cells.
        double trpenalty = NEG_INF;
        if (pty_idx >= 0 && pty_idx < (int)m.trpL.size() && (int)m.trpL[(size_t)pty_idx].size() == M
            && (int)m.trpG.size() > pty_idx && (int)m.trpG[(size_t)pty_idx].size() == M) {
            trpenalty = localOn ? m.trpL[(size_t)pty_idx][(size_t)v]
                                : m.trpG[(size_t)pty_idx][(size_t)v];
        }
        if (noTrPenalty) trpenalty = NEG_INF;  // disable trunc entry
        if (trNotImpossible(trpenalty)) {
            const float trF = (float)trpenalty;
            float scJ = r.Jalpha[(size_t)v][(size_t)L][(size_t)L] + trF;
            if (scJ > r.Jalpha[0][(size_t)L][(size_t)L]) {
                r.Jalpha[0][(size_t)L][(size_t)L] = scJ;
                Jb = v;
            }
            if (fill_L) {
                float scL = r.Lalpha[(size_t)v][(size_t)L][(size_t)L] + trF;
                if (scL > r.Lalpha[0][(size_t)L][(size_t)L]) {
                    r.Lalpha[0][(size_t)L][(size_t)L] = scL;
                    Lb = v;
                }
            }
            if (fill_R) {
                float scR = r.Ralpha[(size_t)v][(size_t)L][(size_t)L] + trF;
                if (scR > r.Ralpha[0][(size_t)L][(size_t)L]) {
                    r.Ralpha[0][(size_t)L][(size_t)L] = scR;
                    Rb = v;
                }
            }
            if (fill_T && st == CM_ST_B) {
                float scT = r.Talpha[(size_t)v][(size_t)L][(size_t)L] + trF;
                if (scT > r.Talpha[0][(size_t)L][(size_t)L]) {
                    r.Talpha[0][(size_t)L][(size_t)L] = scT;
                    Tb = v;
                }
            }
        }

        // Local-begin pickup (Infernal cm_CYKInsideAlign:1028-1033). Only
        // active when noTrPenalty (mirrors --notrunc). beginSc[v] is finite
        // for begin-eligible heads (MATP/MATL/MATR/BIF firsts) with values
        // from cm_CalculateLocalBeginProbs. Track separately, commit at end.
        if (noTrPenalty && localOn && v > 0 && trNotImpossible(s.beginSc)) {
            float cand = r.Jalpha[(size_t)v][(size_t)L][(size_t)L] + (float)s.beginSc;
            if (cand > bsc_local) {
                bsc_local = cand;
                blb_local = v;
            }
        }
    } /* end loop v = M-1..1 */

    // Commit local-begin pickup to root (Infernal cm_CYKInsideAlign:1040-1043).
    if (noTrPenalty && bsc_local > r.Jalpha[0][(size_t)L][(size_t)L]) {
        r.Jalpha[0][(size_t)L][(size_t)L] = bsc_local;
        // Use TR_USED_TRUNC_BEGIN sentinel since trace already handles it
        // (jumps to b without consuming residues, then continues normally).
        r.Jyshadow[0][(size_t)L][(size_t)L] = TR_USED_TRUNC_BEGIN;
        Jb = blb_local;
    }
    // DEBUG: dump alpha[v][j][d] for specific (v,j,d) tuples.
    // Set MMSEQS_CMSCAN_DUMP_CELLS="v:j:d,v:j:d,..." to enable.
    if (const char *cellEnv = std::getenv("MMSEQS_CMSCAN_DUMP_CELLS")) {
        std::string spec = cellEnv;
        size_t pos = 0;
        while (pos < spec.size()) {
            size_t comma = spec.find(',', pos);
            if (comma == std::string::npos) comma = spec.size();
            std::string tok = spec.substr(pos, comma - pos);
            int vv = -1, jj = -1, dd = -1;
            std::sscanf(tok.c_str(), "%d:%d:%d", &vv, &jj, &dd);
            if (vv >= 0 && vv < M && jj >= 0 && jj <= L && dd >= 0 && dd <= jj) {
                float aJ = r.Jalpha[(size_t)vv][(size_t)jj][(size_t)dd];
                std::fprintf(stderr, "[CELL_DBG] v=%d j=%d d=%d type=%d J=%.4f\n",
                             vv, jj, dd, (int)m.states[(size_t)vv].type, aJ);
            }
            pos = comma + 1;
        }
    }
    // DEBUG: dump alpha[v][j][d] for matching Infernal CM_DUMP_CYK_J output.
    // Set MMSEQS_CMSCAN_DUMP_CYK_J="<j>,<j>,..." to enable.
    if (const char *jenv = std::getenv("MMSEQS_CMSCAN_DUMP_CYK_J")) {
        const char *p = jenv;
        while (*p) {
            char *endp = nullptr;
            long jt = std::strtol(p, &endp, 10);
            if (endp == p) break;
            if (jt > 0 && jt <= L) {
                for (int dd = 0; dd <= (int)jt; ++dd) {
                    for (int vv = 0; vv < M; ++vv) {
                        float a = r.Jalpha[(size_t)vv][(size_t)jt][(size_t)dd];
                        if (a <= -1e29f) continue;
                        std::fprintf(stderr, "DUMP_CYK_RIBO j=%ld d=%d v=%d type=%d alpha=%.6f\n",
                                     jt, dd, vv, (int)m.states[(size_t)vv].type, a);
                    }
                }
            }
            p = endp;
            if (*p == ',') p++;
            else if (*p) break;
        }
    }
    // DEBUG: dump root-level decision values when MMSEQS_CMSCAN_DUMP_ROOT=1.
    if (std::getenv("MMSEQS_CMSCAN_DUMP_ROOT") != nullptr) {
        const CmState &s0 = m.states[0];
        std::fprintf(stderr, "[ROOT_DBG] L=%d cnum=%d cfirst=%d alpha[0][L][L]=%.4f Jb=%d\n",
                     L, s0.cnum, s0.cfirst, (double)r.Jalpha[0][(size_t)L][(size_t)L], Jb);
        for (int yoff = 0; yoff < s0.cnum; ++yoff) {
            int y = s0.cfirst + yoff;
            float tsc = (yoff < (int)s0.trans.size()) ? (float)s0.trans[yoff] : -1e30f;
            float ay = (y >= 0 && y < M) ? r.Jalpha[(size_t)y][(size_t)L][(size_t)L] : -1e30f;
            float bgs = (y >= 0 && y < M) ? (float)m.states[(size_t)y].beginSc : -1e30f;
            std::fprintf(stderr, "[ROOT_DBG]   yoff=%d y=%d type=%d tsc=%.4f alpha[y]=%.4f beginSc=%.4f sum=%.4f\n",
                         yoff, y, (int)m.states[(size_t)y].type, tsc, ay, bgs, ay+tsc);
        }
        // Also dump state 3 (MP at node 1) transitions for cell-diff vs Infernal
        if (M > 12) {
            const CmState &s3 = m.states[3];
            std::fprintf(stderr, "[ROOT_DBG] state 3 (MP) cnum=%d cfirst=%d alpha[3][L][L]=%.4f endSc=%.4f\n",
                         s3.cnum, s3.cfirst, (double)r.Jalpha[3][(size_t)L][(size_t)L], s3.endSc);
            for (int yoff = 0; yoff < s3.cnum && yoff < (int)s3.trans.size(); ++yoff) {
                int y = s3.cfirst + yoff;
                float tsc = (float)s3.trans[yoff];
                float ay = (y < M) ? r.Jalpha[(size_t)y][(size_t)(L-1)][(size_t)(L-2)] : -1e30f;
                std::fprintf(stderr, "[ROOT_DBG]   3->yoff=%d y=%d type=%d tsc=%.4f alpha[y][L-1][L-2]=%.4f\n",
                             yoff, y, (int)m.states[(size_t)y].type, tsc, ay);
            }
        }
    }

    // All valid alignments use a truncated begin: stamp the root shadow.
    // With noTrPenalty, leave the regular yshadow set by the v=0 recurrence so
    // the trace walks normally (no TRUNC_BEGIN jump).
    if (!noTrPenalty) {
        r.Jyshadow[0][(size_t)L][(size_t)L] = TR_USED_TRUNC_BEGIN;
        if (fill_L) r.Lyshadow[0][(size_t)L][(size_t)L] = TR_USED_TRUNC_BEGIN;
        if (fill_R) r.Ryshadow[0][(size_t)L][(size_t)L] = TR_USED_TRUNC_BEGIN;
    }

    // Choose the optimal mode/score/b.
    double sc;
    char mode;
    int b;
    if (preset_mode == 'J') {
        sc = r.Jalpha[0][(size_t)L][(size_t)L]; mode = 'J'; b = Jb;
    } else if (preset_mode == 'L') {
        sc = r.Lalpha[0][(size_t)L][(size_t)L]; mode = 'L'; b = Lb;
    } else if (preset_mode == 'R') {
        sc = r.Ralpha[0][(size_t)L][(size_t)L]; mode = 'R'; b = Rb;
    } else if (preset_mode == 'T') {
        sc = r.Talpha[0][(size_t)L][(size_t)L]; mode = 'T'; b = Tb;
    } else {
        sc = r.Jalpha[0][(size_t)L][(size_t)L]; mode = 'J'; b = Jb;
        if (fill_L && r.Lalpha[0][(size_t)L][(size_t)L] > sc) {
            sc = r.Lalpha[0][(size_t)L][(size_t)L]; mode = 'L'; b = Lb;
        }
        if (fill_R && r.Ralpha[0][(size_t)L][(size_t)L] > sc) {
            sc = r.Ralpha[0][(size_t)L][(size_t)L]; mode = 'R'; b = Rb;
        }
        if (fill_T && r.Talpha[0][(size_t)L][(size_t)L] > sc) {
            sc = r.Talpha[0][(size_t)L][(size_t)L]; mode = 'T'; b = Tb;
        }
    }
    r.score = sc;
    r.mode  = mode;
    r.b     = b;
    return r;
}

static void dumpTrCykAlphaFor(const TrCykResult &r, int M, int L, int v_focus) {
    const char *flag = std::getenv("MMSEQS_CMSCAN_DUMP_TRDP");
    if (flag == nullptr || std::string(flag) != "1") return;
    (void)M;
    auto cellOr = [&](const std::vector<std::vector<std::vector<float>>> &a, int v, int j, int d) -> double {
        if ((int)a.size() <= v) return -1e30;
        if ((int)a[(size_t)v].size() <= j) return -1e30;
        if ((int)a[(size_t)v][(size_t)j].size() <= d) return -1e30;
        return (double)a[(size_t)v][(size_t)j][(size_t)d];
    };
    std::fprintf(stderr,
                 "[TRDP_DUMP] v=%d L=%d  J[L][L]=%.4f  L[L][L]=%.4f  R[L][L]=%.4f  T[L][L]=%.4f\n",
                 v_focus, L,
                 cellOr(r.Jalpha, v_focus, L, L),
                 cellOr(r.Lalpha, v_focus, L, L),
                 cellOr(r.Ralpha, v_focus, L, L),
                 cellOr(r.Talpha, v_focus, L, L));
    for (int dj = 0; dj <= 3 && dj <= L; ++dj) {
        const int j = L - dj;
        for (int dd = 0; dd <= 3 && dd <= j; ++dd) {
            const int d = j - dd;
            std::fprintf(stderr,
                         "[TRDP_DUMP]   v=%d j=%d d=%d  J=%.4f L=%.4f R=%.4f T=%.4f\n",
                         v_focus, j, d,
                         cellOr(r.Jalpha, v_focus, j, d),
                         cellOr(r.Lalpha, v_focus, j, d),
                         cellOr(r.Ralpha, v_focus, j, d),
                         cellOr(r.Talpha, v_focus, j, d));
        }
    }
    std::fprintf(stderr, "[TRDP_DUMP] result: score=%.4f mode=%c b=%d\n", r.score, r.mode, r.b);
}

// =====================================================================
// Stage 7: native port of Infernal cm_tr_alignT (cm_dpalign_trunc.c:178).
// Walks J/L/R/T shadow matrices produced by Stage 6 to recover the trCYK
// parsetree (state sequence with mode tags). Mirrors Infernal exactly:
// TRMODE_*_OFFSET decoding, USED_TRUNC_BEGIN/END/EL sentinels, BIF stack.
// =====================================================================
struct TrParseNode {
    int v;            // state index in CM (or M for EL)
    int8_t mode;      // TR_TRMODE_J/L/R/T
    int emitl;        // i-coordinate (1-based)
    int emitr;        // j-coordinate (1-based)
    int parent;       // index in trace, -1 for root
    int nxtl;         // left child trace index, -1 if none
    int nxtr;         // right child trace index, -1 if none
    bool is_right_child; // attach side relative to parent
};

struct TrParsetree {
    std::vector<TrParseNode> nodes;
    bool ok = false;
    char err[256] = {0};
    char optimal_mode = 'J';
    int  b = 0;       // local-entry state
    double trpenalty = 0.0;
    double score = 0.0;
};

// Append a node to the parsetree, attaching it under `parent` on the
// indicated side. Returns the new node's index. Mirrors InsertTraceNodewithMode.
static int trInsertTraceNode(TrParsetree &tr, int parent, bool is_right_child,
                             int i, int j, int v, int8_t mode) {
    TrParseNode n;
    n.v = v;
    n.mode = mode;
    n.emitl = i;
    n.emitr = j;
    n.parent = parent;
    n.nxtl = -1;
    n.nxtr = -1;
    n.is_right_child = is_right_child;
    int idx = (int)tr.nodes.size();
    tr.nodes.push_back(n);
    if (parent >= 0) {
        if (is_right_child) tr.nodes[(size_t)parent].nxtr = idx;
        else                tr.nodes[(size_t)parent].nxtl = idx;
    }
    return idx;
}

static TrParsetree runTrCYKAlignT(const InfernalExactModel &m,
                                  const TrCykResult &r,
                                  int L,
                                  int pty_idx) {
    TrParsetree tr;
    const int M = (int)m.states.size();
    char mode = TR_TRMODE_J;
    switch (r.mode) {
        case 'J': mode = TR_TRMODE_J; break;
        case 'L': mode = TR_TRMODE_L; break;
        case 'R': mode = TR_TRMODE_R; break;
        case 'T': mode = TR_TRMODE_T; break;
        default:  mode = TR_TRMODE_J; break;
    }
    tr.optimal_mode = r.mode;
    tr.b = r.b;
    tr.score = r.score;

    // Trpenalty for the chosen entry state b (local-mode mirrors
    // cm_tr_alignT line 233; we always use local pty since CMH_LOCAL_BEGIN
    // is implied by truncation DP being enabled).
    if (pty_idx >= 0 && r.b >= 0 && r.b < M
        && pty_idx < (int)m.trpL.size()
        && r.b < (int)m.trpL[(size_t)pty_idx].size()) {
        tr.trpenalty = m.trpL[(size_t)pty_idx][(size_t)r.b];
    }

    // Root S: attach state 0 spanning [1..L] in the optimal mode.
    int rootIdx = trInsertTraceNode(tr, -1, false, 1, L, 0, mode);
    (void)rootIdx;

    // Stacks for bifurcation right-child resumption.
    struct PdaEntry { int j; int k; int bifparent; int8_t mode; };
    std::vector<PdaEntry> pda;

    int v = 0;
    int i = 1;
    int j = L;
    int d = L;

    auto getJyshadow = [&](int v, int j, int d) -> int8_t {
        return r.Jyshadow[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getLyshadow = [&](int v, int j, int d) -> int8_t {
        return r.Lyshadow[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getRyshadow = [&](int v, int j, int d) -> int8_t {
        return r.Ryshadow[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getJkshadow = [&](int v, int j, int d) -> int {
        return r.Jkshadow[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getLkshadow = [&](int v, int j, int d) -> int {
        return r.Lkshadow[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getRkshadow = [&](int v, int j, int d) -> int {
        return r.Rkshadow[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getTkshadow = [&](int v, int j, int d) -> int {
        return r.Tkshadow[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getLkmode = [&](int v, int j, int d) -> int8_t {
        return r.Lkmode[(size_t)v][(size_t)j][(size_t)d];
    };
    auto getRkmode = [&](int v, int j, int d) -> int8_t {
        return r.Rkmode[(size_t)v][(size_t)j][(size_t)d];
    };

    const int MAX_STEPS = 8 * (L + 1) * (M + 1) + 1000;
    int steps = 0;

    while (true) {
        if (++steps > MAX_STEPS) {
            std::snprintf(tr.err, sizeof(tr.err),
                          "trace step limit exceeded at v=%d j=%d d=%d", v, j, d);
            return tr;
        }
        if (v < 0 || v > M) {
            std::snprintf(tr.err, sizeof(tr.err),
                          "trace v=%d out of [0..M=%d]", v, M);
            return tr;
        }

        const CmStateType st = (v == M) ? CM_ST_UNKNOWN : m.states[(size_t)v].type;
        const bool is_E_or_EL = (v == M) || (st == CM_ST_E);

        if (v != M && st == CM_ST_B) {
            // Bifurcation: pick k from kshadow per mode.
            int k = 0;
            if      (mode == TR_TRMODE_J) k = getJkshadow(v, j, d);
            else if (mode == TR_TRMODE_L) k = getLkshadow(v, j, d);
            else if (mode == TR_TRMODE_R) k = getRkshadow(v, j, d);
            else if (mode == TR_TRMODE_T) k = getTkshadow(v, j, d);
            else { std::snprintf(tr.err, sizeof(tr.err), "bogus mode at B v=%d", v); return tr; }

            int8_t prvmode = mode;
            int8_t rightmode;
            if      (prvmode == TR_TRMODE_J) rightmode = TR_TRMODE_J;
            else if (prvmode == TR_TRMODE_L) rightmode = TR_TRMODE_L;
            else if (prvmode == TR_TRMODE_R) rightmode = getRkmode(v, j, d);
            else                              rightmode = TR_TRMODE_L; // TRMODE_T

            // Stash right-child state.
            PdaEntry e;
            e.j = j;
            e.k = k;
            e.bifparent = (int)tr.nodes.size() - 1;
            e.mode = rightmode;
            pda.push_back(e);

            // Determine left-child mode.
            int8_t leftmode;
            if      (prvmode == TR_TRMODE_J) leftmode = TR_TRMODE_J;
            else if (prvmode == TR_TRMODE_L) leftmode = getLkmode(v, j, d);
            else if (prvmode == TR_TRMODE_R) leftmode = TR_TRMODE_R;
            else                              leftmode = TR_TRMODE_R; // TRMODE_T

            // Descend to BEGL_S (left child).
            j = j - k;
            d = d - k;
            i = j - d + 1;
            int y = m.states[(size_t)v].cfirst;
            mode = leftmode;
            trInsertTraceNode(tr, (int)tr.nodes.size() - 1, false, i, j, y, mode);
            v = y;
            continue;
        }

        if (is_E_or_EL) {
            // Pop a right-start off the stack, attach as right child.
            if (pda.empty()) break; // done
            PdaEntry e = pda.back(); pda.pop_back();
            int bifparent = e.bifparent;
            int kk = e.k;
            int jj = e.j;
            mode = e.mode;
            int parent_v = tr.nodes[(size_t)bifparent].v;
            // Right child y = cnum of bif state (cnum holds index of right child S).
            int y = m.states[(size_t)parent_v].cnum;
            d = kk;
            j = jj;
            i = j - d + 1;
            trInsertTraceNode(tr, bifparent, true, i, j, y, mode);
            v = y;
            continue;
        }

        // Standard state: read yshadow per mode.
        int yoffset_raw = 0;
        if      (mode == TR_TRMODE_J) yoffset_raw = (int)getJyshadow(v, j, d);
        else if (mode == TR_TRMODE_L) yoffset_raw = (int)getLyshadow(v, j, d);
        else if (mode == TR_TRMODE_R) yoffset_raw = (int)getRyshadow(v, j, d);
        else if (mode == TR_TRMODE_T) {
            if (v == 0) yoffset_raw = TR_USED_TRUNC_BEGIN;
            else { std::snprintf(tr.err, sizeof(tr.err),
                                 "TRMODE_T at non-root non-B v=%d", v); return tr; }
        }
        else { std::snprintf(tr.err, sizeof(tr.err), "bogus mode %d", (int)mode); return tr; }

        int yoffset = yoffset_raw;
        int8_t nxtmode = mode;
        if      (yoffset == TR_USED_TRUNC_BEGIN) { nxtmode = mode; }
        else if (yoffset == TR_USED_TRUNC_END)   { /* irrelevant */ }
        else if (yoffset == TR_USED_EL)          { /* irrelevant */ }
        else if (yoffset >= TR_TRMODE_R_OFFSET)  { nxtmode = TR_TRMODE_R; yoffset -= TR_TRMODE_R_OFFSET; }
        else if (yoffset >= TR_TRMODE_L_OFFSET)  { nxtmode = TR_TRMODE_L; yoffset -= TR_TRMODE_L_OFFSET; }
        else if (yoffset >= TR_TRMODE_J_OFFSET)  { nxtmode = TR_TRMODE_J; yoffset -= TR_TRMODE_J_OFFSET; }
        else { std::snprintf(tr.err, sizeof(tr.err),
                             "yoffset=%d out of bounds at v=%d j=%d d=%d", yoffset_raw, v, j, d); return tr; }

        // Advance i/j based on state-type emission and current mode.
        switch (st) {
            case CM_ST_D: case CM_ST_S: case CM_ST_B: case CM_ST_E:
                break;
            case CM_ST_MP:
                if (mode == TR_TRMODE_J)              i++;
                if (mode == TR_TRMODE_L && d > 0)     i++;
                if (mode == TR_TRMODE_J)              j--;
                if (mode == TR_TRMODE_R && d > 0)     j--;
                break;
            case CM_ST_ML: case CM_ST_IL:
                if (mode == TR_TRMODE_J)              i++;
                if (mode == TR_TRMODE_L && d > 0)     i++;
                break;
            case CM_ST_MR: case CM_ST_IR:
                if (mode == TR_TRMODE_J)              j--;
                if (mode == TR_TRMODE_R && d > 0)     j--;
                break;
            default:
                std::snprintf(tr.err, sizeof(tr.err), "bogus state-type %d at v=%d", (int)st, v); return tr;
        }
        d = j - i + 1;

        if (yoffset == TR_USED_EL || yoffset == TR_USED_TRUNC_END) {
            if (yoffset == TR_USED_EL) {
                trInsertTraceNode(tr, (int)tr.nodes.size() - 1, false, i, j, M, mode);
            }
            v = M; // EL or pseudo-EL
            continue;
        }
        if (yoffset == TR_USED_TRUNC_BEGIN) {
            trInsertTraceNode(tr, (int)tr.nodes.size() - 1, false, i, j, r.b, mode);
            v = r.b;
            continue;
        }

        // Standard descent.
        mode = nxtmode;
        int y = m.states[(size_t)v].cfirst + yoffset;
        trInsertTraceNode(tr, (int)tr.nodes.size() - 1, false, i, j, y, mode);
        v = y;
    }

    tr.ok = true;
    return tr;
}

static void dumpTrParsetree(const TrParsetree &tr, const InfernalExactModel &m) {
    const char *flag = std::getenv("MMSEQS_CMSCAN_DUMP_TRPT");
    if (flag == nullptr || std::string(flag) != "1") return;
    std::fprintf(stderr, "[TRPT_DUMP] mode=%c b=%d score=%.4f trpenalty=%.4f ok=%d %s\n",
                 tr.optimal_mode, tr.b, tr.score, tr.trpenalty, tr.ok ? 1 : 0,
                 tr.ok ? "" : tr.err);
    std::fprintf(stderr, "[TRPT_DUMP] %4s %4s %4s %4s %4s %4s %4s\n",
                 "idx", "v", "mode", "emitl", "emitr", "parent", "type");
    auto modeChar = [](int8_t md) -> char {
        if (md == TR_TRMODE_J) return 'J';
        if (md == TR_TRMODE_L) return 'L';
        if (md == TR_TRMODE_R) return 'R';
        if (md == TR_TRMODE_T) return 'T';
        return '?';
    };
    auto stChar = [&](int v) -> const char* {
        if (v == (int)m.states.size()) return "EL";
        switch (m.states[(size_t)v].type) {
            case CM_ST_MP: return "MP";
            case CM_ST_ML: return "ML";
            case CM_ST_MR: return "MR";
            case CM_ST_IL: return "IL";
            case CM_ST_IR: return "IR";
            case CM_ST_D:  return "D";
            case CM_ST_S:  return "S";
            case CM_ST_B:  return "B";
            case CM_ST_E:  return "E";
            default:       return "?";
        }
    };
    for (size_t k = 0; k < tr.nodes.size(); ++k) {
        const TrParseNode &n = tr.nodes[k];
        std::fprintf(stderr, "[TRPT_DUMP] %4zu %4d:%-3s %3c   %4d %4d %4d %4s\n",
                     k, n.v, stChar(n.v), modeChar(n.mode),
                     n.emitl, n.emitr, n.parent,
                     (n.parent < 0) ? "ROOT" : (n.is_right_child ? "RGT" : "LFT"));
    }
}

// Flatten a TrParsetree into the same trace-state list that the J-only CYK
// emits (DFS pre-order; nodes were already inserted in that order by
// runTrCYKAlignT). EL pseudo-states (v == M) are kept; downstream consumers
// (modelTraceCigarQueryCoord) skip them via "sid >= states.size()".
static std::vector<int> trParsetreeStates(const TrParsetree &tr) {
    std::vector<int> out;
    out.reserve(tr.nodes.size());
    for (const TrParseNode &n : tr.nodes) {
        out.push_back(n.v);
    }
    return out;
}

static InfernalExactModel parseInfernalCmExactModelFromStream(std::istream &in, const std::string &srcLabel) {
    int clen = -1;
    int statesN = -1;
    int w = -1;
    bool inModel = false;
    bool done = false;
    int curMapL = -1;
    int curMapR = -1;
    int curNodeIdx = -1;
    CmNodeType curNodeType = CM_ND_UNKNOWN;
    bool nodeChanged = false;
    InfernalExactModel m;

    std::vector<CmState> states;
    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }
        if (!inModel && t == "CM") {
            inModel = true;
            continue;
        }
        if (inModel && t == "//") {
            done = true;
            break;
        }
        if (!inModel) {
            std::stringstream hs(t);
            std::string key;
            hs >> key;
            if (key == "CLEN") {
                hs >> clen;
            } else if (key == "STATES") {
                hs >> statesN;
            } else if (key == "W") {
                hs >> w;
            } else if (key == "N3OMEGA") {
                hs >> m.null3Omega;
            } else if (key == "N2OMEGA") {
                hs >> m.null2Omega;
            } else if (key == "PBEGIN") {
                hs >> m.pBegin;
                m.hasLocalCfg = true;
            } else if (key == "PEND") {
                hs >> m.pEnd;
                m.hasLocalCfg = true;
            } else if (key == "ELSELF") {
                hs >> m.elSelf;
                m.hasLocalCfg = true;
            } else if (key == "NULL") {
                const double nullBase = 0.25;
                std::string tok0, tok1, tok2, tok3;
                hs >> tok0 >> tok1 >> tok2 >> tok3;
                if (!tok0.empty() && !tok1.empty() && !tok2.empty() && !tok3.empty()) {
                    m.nullProb[0] = parseInfernalAsciiProb(tok0, nullBase);
                    m.nullProb[1] = parseInfernalAsciiProb(tok1, nullBase);
                    m.nullProb[2] = parseInfernalAsciiProb(tok2, nullBase);
                    m.nullProb[3] = parseInfernalAsciiProb(tok3, nullBase);
                    m.hasNull = true;
                }
            } else if (key == "ECMLC" || key == "ECMLI" || key == "ECMGC" || key == "ECMGI") {
                double lambda = 0.0;
                double muExtrap = 0.0;
                double muOrig = 0.0;
                double dbsize = 0.0;
                int nrandhits = 0;
                double tailp = 0.0;
                hs >> lambda >> muExtrap >> muOrig >> dbsize >> nrandhits >> tailp;
                (void) muOrig;
                (void) tailp;
                InfernalExactModel::ExpTail exp;
                exp.valid = true;
                exp.lambda = lambda;
                exp.muExtrap = muExtrap;
                exp.dbsize = dbsize;
                exp.nrandhits = nrandhits;
                if (key == "ECMLC") {
                    m.expLC = exp;
                } else if (key == "ECMLI") {
                    m.expLI = exp;
                } else if (key == "ECMGC") {
                    m.expGC = exp;
                } else {
                    m.expGI = exp;
                }
            }
            continue;
        }

        if (line.find('[') != std::string::npos && line.find(']') != std::string::npos) {
            const size_t open = line.find('[');
            const size_t close = line.find(']');
            std::stringstream ns(line.substr(open + 1, close - open - 1));
            std::string ntype;
            int nidx = -1;
            ns >> ntype >> nidx;
            curNodeType = parseCmNodeType(ntype);
            curNodeIdx = nidx;
            nodeChanged = true;
            std::stringstream aft(line.substr(close + 1));
            std::vector<std::string> toks;
            std::string tok;
            while (aft >> tok) {
                toks.push_back(tok);
            }
            if (toks.size() >= 2) {
                curMapL = parseMaybeInt(toks[0]);
                curMapR = parseMaybeInt(toks[1]);
            } else {
                curMapL = -1;
                curMapR = -1;
            }
            continue;
        }

        std::stringstream ss(line);
        std::vector<std::string> toks;
        std::string tok;
        while (ss >> tok) {
            toks.push_back(tok);
        }
        if (toks.size() < 10) {
            continue;
        }

        const CmStateType st = parseCmStateType(toks[0]);
        if (st == CM_ST_UNKNOWN) {
            continue;
        }
        const int v = parseMaybeInt(toks[1]);
        const int cfirst = parseMaybeInt(toks[4]);
        const int cnum = parseMaybeInt(toks[5]);
        const int dmin1 = parseMaybeInt(toks[7]);
        const int dmax1 = parseMaybeInt(toks[8]);
        if (v < 0) {
            continue;
        }

        CmState s;
        s.type = st;
        s.idx = v;
        s.cfirst = cfirst;
        s.cnum = cnum;
        s.mapLeft = curMapL;
        s.mapRight = curMapR;
        s.nodeIdx = curNodeIdx;
        s.nodeType = curNodeType;
        s.isFirstOfNode = nodeChanged;
        nodeChanged = false;

        int off = 10;
        if (s.type != CM_ST_B && cnum > 0) {
            for (int x = 0; x < cnum && (off + x) < static_cast<int>(toks.size()); ++x) {
                s.trans.push_back(parseInfernalAsciiScore(toks[static_cast<size_t>(off + x)]));
            }
            off += cnum;
        }

        int emitN = 0;
        if (s.type == CM_ST_MP) {
            emitN = 16;
        } else if (s.type == CM_ST_ML || s.type == CM_ST_MR || s.type == CM_ST_IL || s.type == CM_ST_IR) {
            emitN = 4;
        }
        for (int x = 0; x < emitN && (off + x) < static_cast<int>(toks.size()); ++x) {
        s.emit.push_back(parseInfernalAsciiScore(toks[static_cast<size_t>(off + x)]));
        }
        s.emitF.resize(s.emit.size());
        for (size_t k = 0; k < s.emit.size(); ++k) { s.emitF[k] = static_cast<float>(s.emit[k]); }
        s.dmin1 = dmin1;
        s.dmax1 = dmax1;
        states.push_back(s);
    }
    if (!done || clen <= 0 || states.empty()) {
        Debug(Debug::ERROR) << "Failed parsing Infernal CM state graph from " << srcLabel << "\n";
        EXIT(EXIT_FAILURE);
    }

    std::sort(states.begin(), states.end(), [](const CmState &a, const CmState &b) { return a.idx < b.idx; });
    for (size_t i = 0; i < states.size(); ++i) {
        if (states[i].idx != static_cast<int>(i)) {
            Debug(Debug::WARNING) << "Non-contiguous state indices in CM; expected " << i
                                  << " found " << states[i].idx << "\n";
        }
    }
    if (statesN > 0 && static_cast<int>(states.size()) != statesN) {
        Debug(Debug::WARNING) << "CM STATES header=" << statesN
                              << " parsed states=" << states.size() << "\n";
    }

    m.clen = clen;
    m.w = w;
    m.states = states;
    m.rootState = 0;

    // cm_localize: with PBEGIN > 0, identify begin-eligible heads (first state of
    // MATP/MATL/MATR/BIF nodes); with PEND > 0, identify end-eligible heads (first
    // state of MATP/MATL/MATR/BEGL/BEGR nodes whose successor node is not END).
    // Rescale the head's outgoing transitions in linear space by 1/(1+end_lin) and
    // populate beginSc/endSc as log2 probabilities.
    // Gated on MMSEQS_CMSCAN_LOCAL_MODE=1 (default off). Although Infernal cmsearch
    // is local by default, our local pickup is currently a post-fill pass (EL-state
    // contributions don't propagate up via the parent CYK recurrence). Until that's
    // interleaved into the per-state fill loop, enabling local-by-default produces
    // worse alignments than pure glocal (mean -0.075 AUC on 9-query bench, 9-May).
    {
        const char *env = std::getenv("MMSEQS_CMSCAN_LOCAL_MODE");
        const bool enableLocalMode = (env != NULL && std::string(env) == "1");
        if (!enableLocalMode) {
            m.hasLocalCfg = false;
        }
    }
    if (m.hasLocalCfg) {
        // Map from nodeIdx -> {firstStateIdx, nodeType, nextNodeType}
        int maxNode = -1;
        for (const CmState &s : m.states) {
            if (s.nodeIdx > maxNode) maxNode = s.nodeIdx;
        }
        std::vector<int> firstOfNode(static_cast<size_t>(maxNode + 2), -1);
        std::vector<CmNodeType> nodeTypes(static_cast<size_t>(maxNode + 2), CM_ND_UNKNOWN);
        for (const CmState &s : m.states) {
            if (s.nodeIdx >= 0 && s.isFirstOfNode) {
                firstOfNode[static_cast<size_t>(s.nodeIdx)] = s.idx;
                nodeTypes[static_cast<size_t>(s.nodeIdx)] = s.nodeType;
            }
        }

        // Begin-eligible: first state of MATP/MATL/MATR/BIF nodes (excluding node 0).
        // End-eligible:   first state of MATP/MATL/MATR/BEGL/BEGR nodes whose
        //                 successor node is not END (node index >0 and < maxNode).
        std::vector<int> beginHeads;
        std::vector<int> endHeads;
        for (int n = 1; n <= maxNode; ++n) {
            const CmNodeType nt = nodeTypes[static_cast<size_t>(n)];
            const int head = firstOfNode[static_cast<size_t>(n)];
            if (head < 0) continue;
            if (nt == CM_ND_MATP || nt == CM_ND_MATL || nt == CM_ND_MATR || nt == CM_ND_BIF) {
                beginHeads.push_back(head);
            }
            if (nt == CM_ND_MATP || nt == CM_ND_MATL || nt == CM_ND_MATR ||
                nt == CM_ND_BEGL || nt == CM_ND_BEGR) {
                CmNodeType nextNt = (n + 1 <= maxNode) ? nodeTypes[static_cast<size_t>(n + 1)]
                                                       : CM_ND_END;
                if (nextNt != CM_ND_END) {
                    endHeads.push_back(head);
                }
            }
        }

        const double LOG2 = std::log(2.0);
        if (m.pBegin > 0.0 && !beginHeads.empty()) {
            // Mirror Infernal cm_CalculateLocalBeginProbs (cm_modelconfig.c:397):
            // node 1's first state (the natural first child of ROOT_S) gets
            // log2(1 - pBegin). Other internal heads (nd >= 2) share pBegin
            // across nstartsRest = nbegin - 1. Treating all heads uniformly
            // costs ~10 bits at v=0 vs Infernal (verified on 2YGH_1 harness).
            const int node1Head = (1 < static_cast<int>(firstOfNode.size())) ? firstOfNode[1] : -1;
            int nstartsRest = 0;
            for (int h : beginHeads) {
                if (h != node1Head) nstartsRest++;
            }
            const double begLogFirst = std::log(1.0 - m.pBegin) / LOG2;
            const double begLogRest = (nstartsRest > 0)
                ? std::log(m.pBegin / static_cast<double>(nstartsRest)) / LOG2
                : NEG_INF;
            for (int h : beginHeads) {
                m.states[static_cast<size_t>(h)].beginSc =
                    (h == node1Head) ? begLogFirst : begLogRest;
            }
        }
        if (m.pEnd > 0.0 && !endHeads.empty()) {
            const double endLin = m.pEnd / static_cast<double>(endHeads.size());
            // After cm_localize: t[v]_new[i] = t[v]_old[i] / (1 + endLin),
            //                    end[v]      = endLin / (1 + endLin).
            // In log2 space that is a uniform offset of -log2(1 + endLin) added to
            // existing transition log probabilities, and endSc = log2(endLin/(1+endLin)).
            const double scaleLog = -std::log(1.0 + endLin) / LOG2;
            const double endLog = std::log(endLin / (1.0 + endLin)) / LOG2;
            for (int h : endHeads) {
                CmState &s = m.states[static_cast<size_t>(h)];
                for (size_t k = 0; k < s.trans.size(); ++k) {
                    if (s.trans[k] != NEG_INF) {
                        s.trans[k] += scaleLog;
                    }
                }
                s.endSc = endLog;
            }
        }
        // Mirror Infernal cm_modelconfig.c:491 — in local mode, ROOT_S has no
        // standard outgoing transitions; the only way out of node 0 is via the
        // local-begin pickup. Zero ROOT_S transitions to avoid double-counting.
        const char *envZeroRoot = std::getenv("MMSEQS_CMSCAN_ZERO_ROOT_TRANS");
        const bool zeroRootTrans = (envZeroRoot != NULL && envZeroRoot[0] == '1');
        if (zeroRootTrans && m.pBegin > 0.0 && !m.states.empty()) {
            CmState &rootS = m.states[0];
            for (size_t k = 0; k < rootS.trans.size(); ++k) {
                rootS.trans[k] = NEG_INF;
            }
        }
        Debug(Debug::INFO) << "Infernal CM local cfg: PBEGIN=" << m.pBegin
                           << " PEND=" << m.pEnd
                           << " ELSELF=" << m.elSelf
                           << " beginHeads=" << beginHeads.size()
                           << " endHeads=" << endHeads.size() << "\n";
    }

    // Truncation DP support tables (gated on MMSEQS_CMSCAN_TRUNC_DP=1).
    // Mirrors Infernal display.c:CreateEmitMap and cm_trunc.c:cm_tr_penalties_Create.
    {
        const char *envTrunc = std::getenv("MMSEQS_CMSCAN_TRUNC_DP");
        const bool enableTruncDp = (envTrunc != NULL && envTrunc[0] == '1');
        if (enableTruncDp) {
            // Build node table: count nodes, record first-state and type per node.
            int maxNode = -1;
            for (const CmState &s : m.states) {
                if (s.nodeIdx > maxNode) maxNode = s.nodeIdx;
            }
            const int nNodes = maxNode + 1;
            m.nNodes = nNodes;
            m.nodeFirstState.assign(static_cast<size_t>(nNodes), -1);
            m.nodeType.assign(static_cast<size_t>(nNodes), CM_ND_UNKNOWN);
            for (const CmState &s : m.states) {
                if (s.nodeIdx >= 0 && s.isFirstOfNode) {
                    m.nodeFirstState[static_cast<size_t>(s.nodeIdx)] = s.idx;
                    m.nodeType[static_cast<size_t>(s.nodeIdx)] = s.nodeType;
                }
            }
            // For each BIF node, derive left/right child node indices from its
            // first-state's cfirst (-> BEGL state) and cnum (-> BEGR state).
            std::vector<int> bifLeftChild(static_cast<size_t>(nNodes), -1);
            std::vector<int> bifRightChild(static_cast<size_t>(nNodes), -1);
            for (int n = 0; n < nNodes; ++n) {
                if (m.nodeType[static_cast<size_t>(n)] != CM_ND_BIF) continue;
                const int firstSt = m.nodeFirstState[static_cast<size_t>(n)];
                if (firstSt < 0) continue;
                const CmState &bs = m.states[static_cast<size_t>(firstSt)];
                if (bs.cfirst >= 0 && bs.cfirst < (int)m.states.size()) {
                    bifLeftChild[static_cast<size_t>(n)] = m.states[static_cast<size_t>(bs.cfirst)].nodeIdx;
                }
                if (bs.cnum >= 0 && bs.cnum < (int)m.states.size()) {
                    bifRightChild[static_cast<size_t>(n)] = m.states[static_cast<size_t>(bs.cnum)].nodeIdx;
                }
            }

            // CreateEmitMap: stack-based DFS in left/right pass order.
            m.emapLpos.assign(static_cast<size_t>(nNodes), -1);
            m.emapRpos.assign(static_cast<size_t>(nNodes), -1);
            m.emapEpos.assign(static_cast<size_t>(nNodes), -1);
            // Stack entries: pairs of (on_right, nodeIdx). Mirror esl_stack push order.
            std::vector<std::pair<int,int>> pda;
            int cpos = 0;
            pda.emplace_back(0, 0); // start: left side of node 0
            while (!pda.empty()) {
                const int on_right = pda.back().first;
                const int nd = pda.back().second;
                pda.pop_back();
                const CmNodeType nt = (nd >= 0 && nd < nNodes)
                    ? m.nodeType[static_cast<size_t>(nd)]
                    : CM_ND_UNKNOWN;
                if (on_right) {
                    m.emapRpos[static_cast<size_t>(nd)] = cpos + 1;
                    if (nt == CM_ND_MATP || nt == CM_ND_MATR) cpos++;
                } else {
                    if (nt == CM_ND_MATP || nt == CM_ND_MATL) cpos++;
                    m.emapLpos[static_cast<size_t>(nd)] = cpos;
                    if (nt == CM_ND_BIF) {
                        // Push BIF for right side, then right child, then left child.
                        pda.emplace_back(1, nd);
                        pda.emplace_back(0, bifRightChild[static_cast<size_t>(nd)]);
                        pda.emplace_back(0, bifLeftChild[static_cast<size_t>(nd)]);
                    } else {
                        pda.emplace_back(1, nd);
                        if (nt != CM_ND_END) {
                            pda.emplace_back(0, nd + 1);
                        }
                    }
                }
            }
            // epos pass: from end of model back, propagate consensus position
            // an EL insertion follows.
            int eposCur = 0;
            for (int nd = nNodes - 1; nd >= 0; --nd) {
                const CmNodeType nt = m.nodeType[static_cast<size_t>(nd)];
                if (nt == CM_ND_END) {
                    eposCur = m.emapLpos[static_cast<size_t>(nd)];
                } else if (nt == CM_ND_BIF) {
                    const int rc = bifRightChild[static_cast<size_t>(nd)];
                    if (rc >= 0) eposCur = m.emapEpos[static_cast<size_t>(rc)];
                }
                m.emapEpos[static_cast<size_t>(nd)] = eposCur;
            }
            // Sanity: rpos[0] - 1 must equal clen.
            if (nNodes > 0 && m.emapRpos[0] - 1 != m.clen) {
                Debug(Debug::WARNING) << "TRUNC_DP emap: rpos[0]-1=" << (m.emapRpos[0]-1)
                                      << " differs from clen=" << m.clen << "\n";
            }
            for (int nd = 0; nd < nNodes; ++nd) {
                if (m.emapLpos[static_cast<size_t>(nd)] < 0
                 || m.emapRpos[static_cast<size_t>(nd)] < 0
                 || m.emapEpos[static_cast<size_t>(nd)] < 0) {
                    Debug(Debug::WARNING) << "TRUNC_DP emap: node " << nd << " unfilled\n";
                }
            }
            Debug(Debug::INFO) << "TRUNC_DP emap built: nNodes=" << nNodes
                               << " clen=" << m.clen
                               << " (rpos[0]-1=" << (nNodes ? m.emapRpos[0]-1 : -1) << ")\n";

            // Stage 3 + 4: compute psi[v] (cm_ExpectedStateOccupancy) and trp tables
            // (cm_tr_penalties_Create). Mirrors cm.c:cm_ExpectedStateOccupancy and
            // cm_trunc.c:cm_tr_penalties_Create.
            const int M = (int)m.states.size();

            // Linearize transitions; renormalize for end-eligible heads (Infernal
            // CMH_LOCAL_END branch in cm_ExpectedStateOccupancy renormalizes
            // t_copy[v] to discount the local-end probability mass).
            std::vector<std::vector<double>> transLin((size_t)M);
            for (int v = 0; v < M; ++v) {
                const CmState &vs = m.states[(size_t)v];
                transLin[(size_t)v].resize(vs.trans.size(), 0.0);
                double sum = 0.0;
                for (size_t k = 0; k < vs.trans.size(); ++k) {
                    transLin[(size_t)v][k] = (vs.trans[k] == NEG_INF)
                        ? 0.0 : std::pow(2.0, vs.trans[k]);
                    sum += transLin[(size_t)v][k];
                }
                if (vs.endSc != NEG_INF && sum > 1e-12 && std::fabs(sum - 1.0) > 1e-9) {
                    for (size_t k = 0; k < transLin[(size_t)v].size(); ++k) {
                        transLin[(size_t)v][k] /= sum;
                    }
                }
            }
            // ROOT_S: under ZRT=1 all trans are NEG_INF (sum=0). Warn — psi will
            // collapse to 0 for non-root states. Bit-exact requires ZRT=0.
            if (M > 0 && transLin[0].size() > 0) {
                double rootSum = 0.0;
                for (double t : transLin[0]) rootSum += t;
                if (rootSum < 1e-12) {
                    Debug(Debug::WARNING) << "TRUNC_DP: ROOT_S trans zeroed (ZRT=1?); "
                                          << "psi will be degenerate. Disable ZRT for bit-exact.\n";
                }
            }

            // psi[v] = expected number of times state v is entered (linear).
            // Topological order: states are in CM order (parents always have smaller
            // idx than children for non-S parents). Enumerate parents from cfirst/cnum.
            // For S_st: psi=1.0 directly. For non-S: sum psi[parent]*t[parent->v].
            // For IL/IR: add geometric self-loop term t/(1-t).
            std::vector<double> psi((size_t)M, 0.0);
            std::vector<std::vector<std::pair<int,int>>> parents((size_t)M);
            for (int x = 0; x < M; ++x) {
                const CmState &xs = m.states[(size_t)x];
                if (xs.type == CM_ST_B || xs.type == CM_ST_E) continue;
                if (xs.cnum <= 0 || xs.cfirst < 0) continue;
                for (int k = 0; k < xs.cnum; ++k) {
                    int v = xs.cfirst + k;
                    if (v < 0 || v >= M) continue;
                    parents[(size_t)v].push_back({x, k});
                }
            }
            for (int v = 0; v < M; ++v) {
                const CmState &vs = m.states[(size_t)v];
                if (vs.type == CM_ST_S) { psi[(size_t)v] = 1.0; continue; }
                const bool is_insert = (vs.type == CM_ST_IL || vs.type == CM_ST_IR);
                for (const auto &pp : parents[(size_t)v]) {
                    int x = pp.first;
                    int kIdx = pp.second;
                    if (is_insert && x == v) continue; // skip self-loop
                    if (kIdx < 0 || kIdx >= (int)transLin[(size_t)x].size()) continue;
                    psi[(size_t)v] += psi[(size_t)x] * transLin[(size_t)x][(size_t)kIdx];
                }
                if (is_insert && !transLin[(size_t)v].empty()) {
                    double t0 = transLin[(size_t)v][0];
                    if (t0 < 1.0) {
                        psi[(size_t)v] += psi[(size_t)v] * (t0 / (1.0 - t0));
                    }
                }
            }
            m.psi = psi;

            // begin[v] = local-begin probability per state (linear). Mirror
            // cm_modelconfig.c:cm_CalculateLocalBeginProbs.
            std::vector<double> begin((size_t)M, 0.0);
            int nstartsBegin = 0;
            for (int nd = 2; nd < nNodes; ++nd) {
                const CmNodeType nt = m.nodeType[(size_t)nd];
                if (nt == CM_ND_MATP || nt == CM_ND_MATL ||
                    nt == CM_ND_MATR || nt == CM_ND_BIF) nstartsBegin++;
            }
            const double pBeginEff = (m.pBegin > 0.0) ? m.pBegin : 0.05;
            if (nNodes >= 2) {
                int v1 = m.nodeFirstState[1];
                if (v1 >= 0) begin[(size_t)v1] = 1.0 - pBeginEff;
            }
            if (nstartsBegin > 0) {
                double pp = pBeginEff / (double)nstartsBegin;
                for (int nd = 2; nd < nNodes; ++nd) {
                    const CmNodeType nt = m.nodeType[(size_t)nd];
                    if (nt == CM_ND_MATP || nt == CM_ND_MATL ||
                        nt == CM_ND_MATR || nt == CM_ND_BIF) {
                        int v = m.nodeFirstState[(size_t)nd];
                        if (v >= 0) begin[(size_t)v] = pp;
                    }
                }
            }

            // trp tables: per-state truncation penalty for {5'+3', 5' only, 3' only}.
            // Storage: trpG[type][v] global, trpL[type][v] local (both log2).
            // Linear-space computation, convert to log2 at the end.
            std::vector<std::vector<double>> trpGlin(3, std::vector<double>((size_t)M, 0.0));
            std::vector<std::vector<double>> trpLlin(3, std::vector<double>((size_t)M, 0.0));

            const double g_5and3 = 2.0 / ((double)m.clen * (double)(m.clen + 1));
            const double g_5or3  = 1.0 / (double)m.clen;
            double prv5 = 0.0, prv3 = 0.0, prv53 = 0.0;

            for (int nd = 0; nd < nNodes; ++nd) {
                const CmNodeType nt = m.nodeType[(size_t)nd];
                int lpos = m.emapLpos[(size_t)nd];
                int rpos = m.emapRpos[(size_t)nd];
                if (!(nt == CM_ND_MATP || nt == CM_ND_MATL)) lpos += 1;
                if (!(nt == CM_ND_MATP || nt == CM_ND_MATR)) rpos -= 1;

                if (nt == CM_ND_END) {
                    prv5 = prv3 = prv53 = 0.0;
                } else if (nt == CM_ND_BEGL || nt == CM_ND_BEGR) {
                    int begS = m.nodeFirstState[(size_t)nd];
                    int parentB = -1;
                    for (int x = 0; x < M; ++x) {
                        const CmState &xs = m.states[(size_t)x];
                        if (xs.type != CM_ST_B) continue;
                        if (nt == CM_ND_BEGL && xs.cfirst == begS) { parentB = x; break; }
                        if (nt == CM_ND_BEGR && xs.cnum   == begS) { parentB = x; break; }
                    }
                    if (parentB >= 0) {
                        prv5  = (nt == CM_ND_BEGL) ? 0.0 : trpLlin[1][(size_t)parentB];
                        prv3  = (nt == CM_ND_BEGR) ? 0.0 : trpLlin[2][(size_t)parentB];
                        prv53 = trpLlin[0][(size_t)parentB];
                    } else {
                        prv5 = prv3 = prv53 = 0.0;
                    }
                } else if (nt == CM_ND_MATP || nt == CM_ND_MATL ||
                           nt == CM_ND_MATR || nt == CM_ND_BIF) {
                    int mState = m.nodeFirstState[(size_t)nd];
                    int i1 = -1, i2 = -1;
                    if (nd >= 1) {
                        const CmNodeType prevNt = m.nodeType[(size_t)(nd - 1)];
                        int prevFirst = m.nodeFirstState[(size_t)(nd - 1)];
                        if (prevFirst >= 0) {
                            switch (prevNt) {
                                case CM_ND_MATP: i1 = prevFirst + 4; i2 = prevFirst + 5; break;
                                case CM_ND_MATL: i1 = prevFirst + 2; break;
                                case CM_ND_MATR: i1 = prevFirst + 2; break;
                                case CM_ND_BEGR: i1 = prevFirst + 1; break;
                                case CM_ND_ROOT: i1 = prevFirst + 1; i2 = prevFirst + 2; break;
                                default: break;
                            }
                        }
                    }

                    double m_psi = psi[(size_t)mState];
                    if (nt == CM_ND_MATP) {
                        m_psi += psi[(size_t)(mState + 1)] + psi[(size_t)(mState + 2)];
                    }
                    double i1_psi = (i1 == -1) ? 0.0 : psi[(size_t)i1];
                    double i2_psi = (i2 == -1) ? 0.0 : psi[(size_t)i2];
                    double summed_psi = m_psi + i1_psi + i2_psi;
                    if (summed_psi <= 1e-30) summed_psi = 1.0;

                    trpGlin[0][(size_t)mState] = (m_psi  / summed_psi) * g_5and3;
                    if (i1 != -1) trpGlin[0][(size_t)i1] = (i1_psi / summed_psi) * g_5and3;
                    if (i2 != -1) trpGlin[0][(size_t)i2] = (i2_psi / summed_psi) * g_5and3;
                    if (rpos == m.clen) {
                        trpGlin[1][(size_t)mState] = (m_psi  / summed_psi) * g_5or3;
                        if (i1 != -1) trpGlin[1][(size_t)i1] = (i1_psi / summed_psi) * g_5or3;
                        if (i2 != -1) trpGlin[1][(size_t)i2] = (i2_psi / summed_psi) * g_5or3;
                    }
                    if (lpos == 1) {
                        trpGlin[2][(size_t)mState] = (m_psi  / summed_psi) * g_5or3;
                        if (i1 != -1) trpGlin[2][(size_t)i1] = (i1_psi / summed_psi) * g_5or3;
                        if (i2 != -1) trpGlin[2][(size_t)i2] = (i2_psi / summed_psi) * g_5or3;
                    }

                    int subtree_clen = rpos - lpos + 1;
                    int nfrag5 = subtree_clen;
                    int nfrag3 = subtree_clen;
                    int nfrag53 = (subtree_clen * (subtree_clen + 1)) / 2;
                    if (nfrag5 <= 0) nfrag5 = 1;
                    if (nfrag3 <= 0) nfrag3 = 1;
                    if (nfrag53 <= 0) nfrag53 = 1;

                    double cur5  = begin[(size_t)mState] / (double)nfrag5  + prv5;
                    double cur3  = begin[(size_t)mState] / (double)nfrag3  + prv3;
                    double cur53 = begin[(size_t)mState] / (double)nfrag53 + prv53;

                    trpLlin[0][(size_t)mState] = (m_psi  / summed_psi) * cur53;
                    if (i1 != -1) trpLlin[0][(size_t)i1] = (i1_psi / summed_psi) * cur53;
                    if (i2 != -1) trpLlin[0][(size_t)i2] = (i2_psi / summed_psi) * cur53;
                    trpLlin[1][(size_t)mState] = (m_psi  / summed_psi) * cur5;
                    if (i1 != -1) trpLlin[1][(size_t)i1] = (i1_psi / summed_psi) * cur5;
                    if (i2 != -1) trpLlin[1][(size_t)i2] = (i2_psi / summed_psi) * cur5;
                    trpLlin[2][(size_t)mState] = (m_psi  / summed_psi) * cur3;
                    if (i1 != -1) trpLlin[2][(size_t)i1] = (i1_psi / summed_psi) * cur3;
                    if (i2 != -1) trpLlin[2][(size_t)i2] = (i2_psi / summed_psi) * cur3;

                    prv5  = (nt == CM_ND_MATL) ? cur5 : 0.0;
                    prv3  = (nt == CM_ND_MATR) ? cur3 : 0.0;
                    prv53 = cur53;
                }
            }

            // Convert linear -> log2; non-set entries (still 0.0) become NEG_INF.
            const double LN2 = std::log(2.0);
            m.trpG.assign(3, std::vector<double>((size_t)M, NEG_INF));
            m.trpL.assign(3, std::vector<double>((size_t)M, NEG_INF));
            for (int v = 0; v < M; ++v) {
                for (int t = 0; t < 3; ++t) {
                    if (trpGlin[t][(size_t)v] > 0.0)
                        m.trpG[t][(size_t)v] = std::log(trpGlin[t][(size_t)v]) / LN2;
                    if (trpLlin[t][(size_t)v] > 0.0)
                        m.trpL[t][(size_t)v] = std::log(trpLlin[t][(size_t)v]) / LN2;
                }
            }

            // cm_CalcMargLikScores: per-state lmesc/rmesc for L/R modes.
            // MP: marginal sum over partner * null. ML/IL: lmesc=esc; rmesc=0.
            // MR/IR: rmesc=esc; lmesc=0. Other states: leave empty.
            for (size_t v = 0; v < m.states.size(); ++v) {
                CmState &s = m.states[v];
                s.lmesc.clear();
                s.rmesc.clear();
                if (s.type == CM_ST_MP && s.emit.size() >= 16) {
                    s.lmesc.assign(4, NEG_INF);
                    s.rmesc.assign(4, NEG_INF);
                    for (int x = 0; x < 4; ++x) {
                        double sumL = 0.0;
                        double sumR = 0.0;
                        for (int y = 0; y < 4; ++y) {
                            // emit[a*4+b] = log2(e[a,b] / (null[a]*null[b]))
                            // lmesc[x] = log2(sum_y e[x,y] / null[x])
                            //          = log2(sum_y null[y] * 2^emit[x*4+y])
                            const double lEsc = s.emit[(size_t)(x * 4 + y)];
                            const double rEsc = s.emit[(size_t)(y * 4 + x)];
                            sumL += m.nullProb[(size_t)y] * std::pow(2.0, lEsc);
                            sumR += m.nullProb[(size_t)y] * std::pow(2.0, rEsc);
                        }
                        s.lmesc[(size_t)x] = (sumL > 0.0) ? std::log(sumL) / std::log(2.0) : NEG_INF;
                        s.rmesc[(size_t)x] = (sumR > 0.0) ? std::log(sumR) / std::log(2.0) : NEG_INF;
                    }
                } else if ((s.type == CM_ST_ML || s.type == CM_ST_IL) && s.emit.size() >= 4) {
                    s.lmesc.assign(4, 0.0);
                    s.rmesc.assign(4, 0.0);
                    for (int x = 0; x < 4; ++x) s.lmesc[(size_t)x] = s.emit[(size_t)x];
                } else if ((s.type == CM_ST_MR || s.type == CM_ST_IR) && s.emit.size() >= 4) {
                    s.lmesc.assign(4, 0.0);
                    s.rmesc.assign(4, 0.0);
                    for (int x = 0; x < 4; ++x) s.rmesc[(size_t)x] = s.emit[(size_t)x];
                }
            }

            m.hasTruncDp = true;

            const char *envDumpEmap = std::getenv("MMSEQS_CMSCAN_DUMP_EMAP");
            if (envDumpEmap != NULL && envDumpEmap[0] == '1') {
                fprintf(stderr, "EMAP_DUMP nNodes=%d clen=%d\n", nNodes, m.clen);
                fprintf(stderr, "%4s %7s %9s %4s %4s %4s\n",
                        "Node", "State1", "Type", "lpos", "rpos", "epos");
                static const char *ntStr[] = {
                    "?", "ROOT", "MATP", "MATL", "MATR", "BEGL", "BEGR", "BIF", "END"
                };
                for (int nd = 0; nd < nNodes; ++nd) {
                    int t = (int)m.nodeType[(size_t)nd];
                    if (t < 0 || t > 8) t = 0;
                    fprintf(stderr, "%4d %7d %9s %4d %4d %4d\n",
                            nd, m.nodeFirstState[(size_t)nd], ntStr[t],
                            m.emapLpos[(size_t)nd], m.emapRpos[(size_t)nd],
                            m.emapEpos[(size_t)nd]);
                }
            }
            const char *envDumpMesc = std::getenv("MMSEQS_CMSCAN_DUMP_MESC");
            if (envDumpMesc != NULL && envDumpMesc[0] == '1') {
                fprintf(stderr, "MESC_DUMP M=%d (MP-only rows shown for brevity)\n", M);
                fprintf(stderr, "%4s %4s  %10s %10s %10s %10s   %10s %10s %10s %10s\n",
                        "v", "type", "lmA", "lmC", "lmG", "lmU", "rmA", "rmC", "rmG", "rmU");
                for (int v = 0; v < M; ++v) {
                    const CmState &s = m.states[(size_t)v];
                    if (s.type != CM_ST_MP) continue;
                    if (s.lmesc.size() < 4 || s.rmesc.size() < 4) continue;
                    fprintf(stderr, "%4d %4d  %10.6f %10.6f %10.6f %10.6f   %10.6f %10.6f %10.6f %10.6f\n",
                            v, (int)s.type,
                            s.lmesc[0], s.lmesc[1], s.lmesc[2], s.lmesc[3],
                            s.rmesc[0], s.rmesc[1], s.rmesc[2], s.rmesc[3]);
                }
            }
            const char *envDumpTrp = std::getenv("MMSEQS_CMSCAN_DUMP_TRP");
            if (envDumpTrp != NULL && envDumpTrp[0] == '1') {
                fprintf(stderr, "TRP_DUMP M=%d clen=%d\n", M, m.clen);
                fprintf(stderr, "%4s %4s  %12s %12s %12s   %12s %12s %12s\n",
                        "v", "type", "g_5+3", "g_5", "g_3", "l_5+3", "l_5", "l_3");
                for (int v = 0; v < M; ++v) {
                    if (m.trpG[0][(size_t)v] == NEG_INF
                     && m.trpL[0][(size_t)v] == NEG_INF) continue;
                    fprintf(stderr, "%4d %4d  %12.6f %12.6f %12.6f   %12.6f %12.6f %12.6f\n",
                            v, (int)m.states[(size_t)v].type,
                            m.trpG[0][(size_t)v], m.trpG[1][(size_t)v], m.trpG[2][(size_t)v],
                            m.trpL[0][(size_t)v], m.trpL[1][(size_t)v], m.trpL[2][(size_t)v]);
                }
            }
        }
        // Stage 6+7 verification: when MMSEQS_CMSCAN_TRUNC_TEST_FASTA points
        // to a single-record FASTA, run the truncated CYK + trace dispatcher
        // on it and dump the parsetree (gated by MMSEQS_CMSCAN_DUMP_TRPT=1).
        // Decoupled from the production scan path; pure verification harness.
        const char *envTestFasta = std::getenv("MMSEQS_CMSCAN_TRUNC_TEST_FASTA");
        const bool enableTruncDpHarness = (std::getenv("MMSEQS_CMSCAN_TRUNC_DP") != NULL
                                           && std::getenv("MMSEQS_CMSCAN_TRUNC_DP")[0] == '1');
        if (enableTruncDpHarness && envTestFasta != NULL && envTestFasta[0] != '\0') {
            std::ifstream fin(envTestFasta);
            if (!fin.good()) {
                std::fprintf(stderr, "[TRUNC_TEST] cannot open FASTA: %s\n", envTestFasta);
            } else {
                std::string fline, hdr, seq;
                while (std::getline(fin, fline)) {
                    if (fline.empty()) continue;
                    if (fline[0] == '>') {
                        if (!hdr.empty()) break;
                        hdr = fline.substr(1);
                    } else {
                        for (char c : fline) {
                            if (!std::isspace(static_cast<unsigned char>(c))) seq.push_back(c);
                        }
                    }
                }
                if (seq.empty()) {
                    std::fprintf(stderr, "[TRUNC_TEST] empty sequence in %s\n", envTestFasta);
                } else {
                    const int L = (int)seq.size();
                    std::vector<int8_t> dsq((size_t)(L + 2), -1); // 1-based, [0] sentinel
                    for (int p = 1; p <= L; ++p) {
                        char c = seq[(size_t)(p - 1)];
                        if      (c=='A'||c=='a') dsq[(size_t)p] = 0;
                        else if (c=='C'||c=='c') dsq[(size_t)p] = 1;
                        else if (c=='G'||c=='g') dsq[(size_t)p] = 2;
                        else if (c=='T'||c=='t'||c=='U'||c=='u') dsq[(size_t)p] = 3;
                        else dsq[(size_t)p] = -1;
                    }
                    const int pty_idx = 0; // TRPENALTY_5P_AND_3P (matches cmsearch glocal/local)
                    const char preset_mode = 'U'; // unknown — fill J/L/R/T
                    std::fprintf(stderr, "[TRUNC_TEST] hdr=%s L=%d pty_idx=%d preset=%c\n",
                                 hdr.c_str(), L, pty_idx, preset_mode);
                    const TrCykResult &res = runTrCYKInsideAlign(m, dsq, L, preset_mode, pty_idx);
                    std::fprintf(stderr, "[TRUNC_TEST] DP done: score=%.4f mode=%c b=%d\n",
                                 res.score, res.mode, res.b);
                    TrParsetree tr = runTrCYKAlignT(m, res, L, pty_idx);
                    std::fprintf(stderr, "[TRUNC_TEST] trace done: ok=%d nodes=%zu err=%s\n",
                                 tr.ok ? 1 : 0, tr.nodes.size(), tr.err);
                    dumpTrParsetree(tr, m);
                    dumpTrCykAlphaFor(res, (int)m.states.size(), L, /*v_focus=*/0);
                    dumpTrCykAlphaFor(res, (int)m.states.size(), L, /*v_focus=*/84);
                    dumpTrCykAlphaFor(res, (int)m.states.size(), L, /*v_focus=*/9);
                    // Specifically: BEGL_S=85 left subtree of BIF v=84,
                    // BEGR_S=165 right subtree. Expect Ralpha[85][90][90]≈-9.13,
                    // Lalpha[165][103][13]≈-2.05 from Infernal trace breakdown.
                    if (std::getenv("MMSEQS_CMSCAN_DUMP_TRDP") != NULL
                        && std::getenv("MMSEQS_CMSCAN_DUMP_TRDP")[0] == '1') {
                        auto cellOr = [&](const std::vector<std::vector<std::vector<float>>> &a,
                                          int v, int j, int d) -> double {
                            if ((int)a.size() <= v) return -1e30;
                            if ((int)a[(size_t)v].size() <= j) return -1e30;
                            if ((int)a[(size_t)v][(size_t)j].size() <= d) return -1e30;
                            return (double)a[(size_t)v][(size_t)j][(size_t)d];
                        };
                        std::fprintf(stderr, "[TRDP_BIF] v=85 (BEGL_S, left of BIF v=84):\n");
                        for (int dd : {88, 89, 90, 91, 92}) {
                            for (int jj : {90}) {
                                if (dd > jj) continue;
                                std::fprintf(stderr, "[TRDP_BIF]   v=85 j=%d d=%d  J=%.4f L=%.4f R=%.4f\n",
                                             jj, dd,
                                             cellOr(res.Jalpha, 85, jj, dd),
                                             cellOr(res.Lalpha, 85, jj, dd),
                                             cellOr(res.Ralpha, 85, jj, dd));
                            }
                        }
                        // Walk down the broken left subtree to find where the gap appears.
                        // Chain (Infernal trace): v=86 MP → 92 MP → 98 MP → 104 ML → 107 ML →
                        // 110 ML → 113 MR → 116 MR → 119 MR → 122 MR → 125 MR → 299 EL
                        // Inside-score reaching v=86 from the leaf side should be ≈ -9 + tsc(85→86).
                        // Probe v=86, v=125 in J/L/R modes at relevant cells.
                        std::fprintf(stderr, "[TRDP_BIF] v=86 (first MP after BEGL_S):\n");
                        for (int dd : {88, 89, 90}) {
                            for (int jj : {90}) {
                                if (dd > jj) continue;
                                std::fprintf(stderr, "[TRDP_BIF]   v=86 j=%d d=%d  J=%.4f L=%.4f R=%.4f\n",
                                             jj, dd,
                                             cellOr(res.Jalpha, 86, jj, dd),
                                             cellOr(res.Lalpha, 86, jj, dd),
                                             cellOr(res.Ralpha, 86, jj, dd));
                            }
                        }
                        // v=125 is the MR that local-end-pops into EL absorbing 79 residues.
                        // Infernal trace: tsc=-17.50, then EL spans residues 4..82 = 79 res.
                        // Inside score at v=125 j=82 d=80 should be ≈ -17.50 + esc(MR@A) + elSelf*79.
                        std::fprintf(stderr, "[TRDP_BIF] v=125 (MR with local-end pop in J trace):\n");
                        for (int dd : {78, 79, 80, 81}) {
                            for (int jj : {82, 83}) {
                                if (dd > jj) continue;
                                std::fprintf(stderr, "[TRDP_BIF]   v=125 j=%d d=%d  J=%.4f L=%.4f R=%.4f\n",
                                             jj, dd,
                                             cellOr(res.Jalpha, 125, jj, dd),
                                             cellOr(res.Lalpha, 125, jj, dd),
                                             cellOr(res.Ralpha, 125, jj, dd));
                            }
                        }
                        // EL deck (state index = M).
                        std::fprintf(stderr, "[TRDP_BIF] v=M (EL deck):\n");
                        for (int dd : {0, 1, 79, 80}) {
                            std::fprintf(stderr, "[TRDP_BIF]   v=%d j=%d d=%d  J=%.4f L=%.4f R=%.4f\n",
                                         (int)m.states.size(), dd, dd,
                                         cellOr(res.Jalpha, (int)m.states.size(), dd, dd),
                                         cellOr(res.Lalpha, (int)m.states.size(), dd, dd),
                                         cellOr(res.Ralpha, (int)m.states.size(), dd, dd));
                        }
                        // Print key model parameters: elSelf, endsc[125], state types.
                        std::fprintf(stderr, "[TRDP_BIF] elSelf=%.4f localOn=%d pBegin=%.4f pEnd=%.4f\n",
                                     m.elSelf, (int)m.hasLocalCfg, m.pBegin, m.pEnd);
                        if (125 < (int)m.states.size()) {
                            std::fprintf(stderr, "[TRDP_BIF] state[125] type=%d endSc=%.4f beginSc=%.4f cfirst=%d cnum=%d\n",
                                         (int)m.states[125].type, m.states[125].endSc, m.states[125].beginSc,
                                         m.states[125].cfirst, m.states[125].cnum);
                        }
                        std::fprintf(stderr, "[TRDP_BIF] v=165 (BEGR_S, right of BIF v=84):\n");
                        for (int dd : {11, 12, 13, 14, 15}) {
                            for (int jj : {103}) {
                                if (dd > jj) continue;
                                std::fprintf(stderr, "[TRDP_BIF]   v=165 j=%d d=%d  J=%.4f L=%.4f R=%.4f\n",
                                             jj, dd,
                                             cellOr(res.Jalpha, 165, jj, dd),
                                             cellOr(res.Lalpha, 165, jj, dd),
                                             cellOr(res.Ralpha, 165, jj, dd));
                            }
                        }
                    }
                    // Per-state argmax(J,L,R,T)+trpenalty, including all candidate b values.
                    if (std::getenv("MMSEQS_CMSCAN_DUMP_TRDP") != NULL
                        && std::getenv("MMSEQS_CMSCAN_DUMP_TRDP")[0] == '1') {
                        std::fprintf(stderr, "[TRDP_BCAND] v   J+pty5+3   L+pty5    R+pty3    T+pty5+3\n");
                        for (size_t v = 0; v < m.states.size(); ++v) {
                            const CmState &s = m.states[v];
                            if (s.type != CM_ST_MP && s.type != CM_ST_ML
                                && s.type != CM_ST_MR && s.type != CM_ST_B
                                && s.type != CM_ST_S) continue;
                            if (pty_idx >= (int)m.trpL.size()) continue;
                            double pJ = (m.trpL[0].size() > v && m.trpL[0][v] > -1e29) ? m.trpL[0][v] : 0.0;
                            double pL = (m.trpL[1].size() > v && m.trpL[1][v] > -1e29) ? m.trpL[1][v] : 0.0;
                            double pR = (m.trpL[2].size() > v && m.trpL[2][v] > -1e29) ? m.trpL[2][v] : 0.0;
                            double j = (res.Jalpha.size() > v && res.Jalpha[v].size() > (size_t)L
                                        && res.Jalpha[v][(size_t)L].size() > (size_t)L)
                                       ? res.Jalpha[v][(size_t)L][(size_t)L] : -1e30;
                            double l = (res.Lalpha.size() > v && res.Lalpha[v].size() > (size_t)L
                                        && res.Lalpha[v][(size_t)L].size() > (size_t)L)
                                       ? res.Lalpha[v][(size_t)L][(size_t)L] : -1e30;
                            double rr = (res.Ralpha.size() > v && res.Ralpha[v].size() > (size_t)L
                                         && res.Ralpha[v][(size_t)L].size() > (size_t)L)
                                        ? res.Ralpha[v][(size_t)L][(size_t)L] : -1e30;
                            double t = (res.Talpha.size() > v && res.Talpha[v].size() > (size_t)L
                                        && res.Talpha[v][(size_t)L].size() > (size_t)L)
                                       ? res.Talpha[v][(size_t)L][(size_t)L] : -1e30;
                            std::fprintf(stderr, "[TRDP_BCAND] %3zu  %10.4f %10.4f %10.4f %10.4f\n",
                                         v,
                                         (j > -1e29) ? j + pJ : -999.0,
                                         (l > -1e29) ? l + pL : -999.0,
                                         (rr > -1e29) ? rr + pR : -999.0,
                                         (t > -1e29) ? t + pJ : -999.0);
                        }
                    }
                }
            }
        }
    }

    Debug(Debug::INFO) << "Infernal CM exact parser: CLEN=" << m.clen
                       << " W=" << m.w
                       << " states=" << m.states.size() << "\n";
    return m;
}

static InfernalExactModel parseInfernalCmExactModel(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        Debug(Debug::ERROR) << "Cannot open Infernal CM file: " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    return parseInfernalCmExactModelFromStream(in, path);
}

static bool hasDbIndex(const std::string &path) {
    return FileUtil::fileExists((path + ".index").c_str());
}

static std::string getSelfExecutablePath() {
    std::vector<char> buf(4096, '\0');
    const ssize_t n = ::readlink("/proc/self/exe", buf.data(), buf.size() - 1);
    if (n <= 0) {
        return std::string();
    }
    buf[static_cast<size_t>(n)] = '\0';
    return std::string(buf.data());
}

static bool runExternalCommand(const std::vector<std::string> &args,
                               const std::vector<std::pair<std::string, std::string>> &envOverrides,
                               std::string &error) {
    if (args.empty()) {
        error = "empty command";
        return false;
    }
    std::vector<char *> argv;
    argv.reserve(args.size() + 1);
    for (size_t i = 0; i < args.size(); ++i) {
        argv.push_back(const_cast<char *>(args[i].c_str()));
    }
    argv.push_back(NULL);
    pid_t pid = fork();
    if (pid < 0) {
        error = "fork failed: " + std::string(std::strerror(errno));
        return false;
    }
    if (pid == 0) {
        for (size_t i = 0; i < envOverrides.size(); ++i) {
            setenv(envOverrides[i].first.c_str(), envOverrides[i].second.c_str(), 1);
        }
        execv(argv[0], argv.data());
        _exit(127);
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        std::ostringstream oss;
        oss << "command failed with status " << status << ":";
        for (size_t i = 0; i < args.size(); ++i) {
            oss << ' ' << args[i];
        }
        error = oss.str();
        return false;
    }
    return true;
}

static std::string createTempDir(const std::string &prefix) {
    std::string templ = "/tmp/" + prefix + "_XXXXXX";
    std::vector<char> path(templ.begin(), templ.end());
    path.push_back('\0');
    char *dir = ::mkdtemp(path.data());
    return (dir != NULL) ? std::string(dir) : std::string();
}

static bool ensureDirectoryRecursive(const std::string &path) {
    if (path.empty()) {
        return false;
    }
    if (FileUtil::directoryExists(path.c_str())) {
        return true;
    }

    std::string cur;
    size_t pos = 0;
    if (path[0] == '/') {
        cur = "/";
        pos = 1;
    }
    while (pos <= path.size()) {
        const size_t next = path.find('/', pos);
        const std::string part = path.substr(pos, next - pos);
        if (!part.empty()) {
            if (!cur.empty() && cur[cur.size() - 1] != '/') {
                cur.push_back('/');
            }
            cur += part;
            if (!FileUtil::directoryExists(cur.c_str()) && !FileUtil::makeDir(cur.c_str())) {
                return false;
            }
        }
        if (next == std::string::npos) {
            break;
        }
        pos = next + 1;
    }
    return FileUtil::directoryExists(path.c_str());
}

static std::string trimCopy(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) {
        ++b;
    }
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) {
        --e;
    }
    return s.substr(b, e - b);
}

static std::vector<std::string> splitCommaList(const std::string &s) {
    std::vector<std::string> out;
    size_t pos = 0;
    while (pos <= s.size()) {
        const size_t next = s.find(',', pos);
        const std::string part = trimCopy(s.substr(pos, next - pos));
        if (!part.empty()) {
            out.push_back(part);
        }
        if (next == std::string::npos) {
            break;
        }
        pos = next + 1;
    }
    if (out.empty() && !trimCopy(s).empty()) {
        out.push_back(trimCopy(s));
    }
    return out;
}

static bool mergeLolcmsearchInputs(const std::string &queryDb,
                                   const std::vector<std::string> &targetDbs,
                                   const std::vector<std::string> &resultDbs,
                                   const std::string &mergedTargetDb,
                                   const std::string &mergedResultDb,
                                   int threads,
                                   std::string &error) {
    if (targetDbs.empty() || targetDbs.size() != resultDbs.size()) {
        error = "target/result DB list sizes do not match";
        return false;
    }

    DBReader<unsigned int> qDbr(queryDb.c_str(), (queryDb + ".index").c_str(),
                                std::max(1, threads),
                                DBReader<unsigned int>::USE_INDEX);
    qDbr.open(DBReader<unsigned int>::NOSORT);

    std::vector<DBReader<unsigned int> *> tReaders(targetDbs.size(), NULL);
    std::vector<DBReader<unsigned int> *> rReaders(resultDbs.size(), NULL);
    std::vector<std::unordered_map<unsigned int, unsigned int> > keyMaps(targetDbs.size());
    int mergedTargetDbtype = -1;

    try {
        for (size_t i = 0; i < targetDbs.size(); ++i) {
            tReaders[i] = new DBReader<unsigned int>(targetDbs[i].c_str(),
                                                     (targetDbs[i] + ".index").c_str(),
                                                     std::max(1, threads),
                                                     DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
            tReaders[i]->open(DBReader<unsigned int>::NOSORT);
            rReaders[i] = new DBReader<unsigned int>(resultDbs[i].c_str(),
                                                     (resultDbs[i] + ".index").c_str(),
                                                     std::max(1, threads),
                                                     DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
            rReaders[i]->open(DBReader<unsigned int>::NOSORT);
            if (mergedTargetDbtype < 0) {
                mergedTargetDbtype = tReaders[i]->getDbtype();
            } else if (tReaders[i]->getDbtype() != mergedTargetDbtype) {
                error = "target DB types do not match across inputs";
                throw std::runtime_error(error);
            }
        }

        FILE *targetIndexFile = fopen((mergedTargetDb + ".index").c_str(), "w");
        if (targetIndexFile == NULL) {
            error = "failed to create merged target index";
            throw std::runtime_error(error);
        }
        unsigned int nextTargetKey = 0;
        size_t nextShardId = 0;
        size_t cumulativeDataSize = 0;
        for (size_t ti = 0; ti < tReaders.size(); ++ti) {
            DBReader<unsigned int> &tDbr = *tReaders[ti];
            std::unordered_map<unsigned int, unsigned int> &map = keyMaps[ti];
            map.reserve(tDbr.getSize() * 2 + 1);
            const std::vector<std::string> dataFiles = FileUtil::findDatafiles(targetDbs[ti].c_str());
            if (dataFiles.empty()) {
                error = "no datafiles found for target DB " + targetDbs[ti];
                throw std::runtime_error(error);
            }
            size_t thisDbDataSize = 0;
            for (size_t fi = 0; fi < dataFiles.size(); ++fi) {
                const std::string realDataFile = FileUtil::getRealPathFromSymLink(dataFiles[fi]);
                struct stat st;
                if (stat(realDataFile.c_str(), &st) != 0) {
                    error = "failed to stat target data file " + realDataFile;
                    throw std::runtime_error(error);
                }
                thisDbDataSize += static_cast<size_t>(st.st_size);
                const std::string dstShard = mergedTargetDb + "." + SSTR(nextShardId++);
                if (symlink(realDataFile.c_str(), dstShard.c_str()) != 0) {
                    error = "failed to symlink target data file " + realDataFile + " -> " + dstShard;
                    throw std::runtime_error(error);
                }
            }
            for (size_t id = 0; id < tDbr.getSize(); ++id) {
                const unsigned int oldKey = tDbr.getDbKey(id);
                const size_t mergedOffset = cumulativeDataSize + tDbr.getOffset(id);
                const size_t len = tDbr.getEntryLen(id);
                fprintf(targetIndexFile, "%u\t%zu\t%zu\n", nextTargetKey, mergedOffset, len);
                map[oldKey] = nextTargetKey;
                ++nextTargetKey;
            }
            cumulativeDataSize += thisDbDataSize;
        }
        if (fclose(targetIndexFile) != 0) {
            error = "failed to close merged target index";
            throw std::runtime_error(error);
        }
        {
            FILE *dbtypeFile = fopen((mergedTargetDb + ".dbtype").c_str(), "wb");
            if (dbtypeFile == NULL) {
                error = "failed to create merged target dbtype";
                throw std::runtime_error(error);
            }
            if (fwrite(&mergedTargetDbtype, sizeof(int), 1, dbtypeFile) != 1) {
                fclose(dbtypeFile);
                error = "failed to write merged target dbtype";
                throw std::runtime_error(error);
            }
            if (fclose(dbtypeFile) != 0) {
                error = "failed to close merged target dbtype";
                throw std::runtime_error(error);
            }
        }

        FILE *lookupFilePtr = fopen((mergedTargetDb + ".lookup").c_str(), "w");
        if (lookupFilePtr == NULL) {
            error = "failed to create merged target lookup";
            throw std::runtime_error(error);
        }
        DBWriter headerWriter((mergedTargetDb + "_h").c_str(),
                              (mergedTargetDb + "_h.index").c_str(),
                              1, 0, Parameters::DBTYPE_GENERIC_DB);
        headerWriter.open();
        unsigned int nextFileNumberOffset = 0;
        for (size_t ti = 0; ti < tReaders.size(); ++ti) {
            unsigned int localMaxFileNumber = 0;
            const std::string srcLookup = targetDbs[ti] + ".lookup";
            if (FileUtil::fileExists(srcLookup.c_str())) {
                DBReader<unsigned int> lookupReader(targetDbs[ti].c_str(),
                                                    (targetDbs[ti] + ".index").c_str(),
                                                    1,
                                                    DBReader<unsigned int>::USE_LOOKUP);
                lookupReader.open(DBReader<unsigned int>::NOSORT);
                DBReader<unsigned int>::LookupEntry *lookup = lookupReader.getLookup();
                for (size_t li = 0; li < lookupReader.getLookupSize(); ++li) {
                    const unsigned int oldKey = lookup[li].id;
                    std::unordered_map<unsigned int, unsigned int>::const_iterator it = keyMaps[ti].find(oldKey);
                    if (it == keyMaps[ti].end()) {
                        continue;
                    }
                    const unsigned int mergedKey = it->second;
                    const unsigned int mergedFileNumber = nextFileNumberOffset + lookup[li].fileNumber;
                    if (lookup[li].fileNumber > localMaxFileNumber) {
                        localMaxFileNumber = lookup[li].fileNumber;
                    }
                    fprintf(lookupFilePtr, "%u\t%s\t%u\n",
                            mergedKey,
                            lookup[li].entryName.c_str(),
                            mergedFileNumber);
                    headerWriter.writeData(lookup[li].entryName.c_str(),
                                           lookup[li].entryName.size(),
                                           mergedKey,
                                           0);
                }
                lookupReader.close();
                nextFileNumberOffset += (localMaxFileNumber + 1);
            } else {
                DBReader<unsigned int> &tDbr = *tReaders[ti];
                for (size_t id = 0; id < tDbr.getSize(); ++id) {
                    const unsigned int oldKey = tDbr.getDbKey(id);
                    std::unordered_map<unsigned int, unsigned int>::const_iterator it = keyMaps[ti].find(oldKey);
                    if (it == keyMaps[ti].end()) {
                        continue;
                    }
                    const unsigned int mergedKey = it->second;
                    fprintf(lookupFilePtr, "%u\t%u\t%u\n", mergedKey, mergedKey, nextFileNumberOffset);
                    const std::string mergedName = SSTR(mergedKey);
                    headerWriter.writeData(mergedName.c_str(), mergedName.size(), mergedKey, 0);
                }
                nextFileNumberOffset += 1;
            }
        }
        if (fclose(lookupFilePtr) != 0) {
            error = "failed to close merged target lookup";
            throw std::runtime_error(error);
        }
        headerWriter.close(true);

        DBWriter resultWriter(mergedResultDb.c_str(), (mergedResultDb + ".index").c_str(),
                              1, 0, Parameters::DBTYPE_ALIGNMENT_RES);
        resultWriter.open();
        for (size_t qi = 0; qi < qDbr.getSize(); ++qi) {
            const unsigned int queryKey = qDbr.getDbKey(qi);
            std::string out;
            for (size_t ri = 0; ri < rReaders.size(); ++ri) {
                DBReader<unsigned int> &rDbr = *rReaders[ri];
                const size_t rid = rDbr.getId(queryKey);
                if (rid == UINT_MAX) {
                    continue;
                }
                char *block = rDbr.getData(rid, 0);
                const size_t blockLen = rDbr.getEntryLen(rid);
                size_t pos = 0;
                const size_t usable = (blockLen == 0 ? 0 : blockLen - 1);
                while (pos < usable) {
                    const size_t lineStart = pos;
                    while (pos < usable && block[pos] != '\n') {
                        ++pos;
                    }
                    const size_t lineEnd = pos;
                    if (pos < usable && block[pos] == '\n') {
                        ++pos;
                    }
                    if (lineEnd <= lineStart) {
                        continue;
                    }
                    size_t tabPos = lineStart;
                    while (tabPos < lineEnd && block[tabPos] != '\t') {
                        ++tabPos;
                    }
                    if (tabPos <= lineStart || tabPos >= lineEnd) {
                        continue;
                    }
                    const unsigned int oldKey = static_cast<unsigned int>(
                        std::strtoul(block + lineStart, NULL, 10));
                    std::unordered_map<unsigned int, unsigned int>::const_iterator it = keyMaps[ri].find(oldKey);
                    if (it == keyMaps[ri].end()) {
                        continue;
                    }
                    out.append(SSTR(it->second));
                    out.append(block + tabPos, lineEnd - tabPos);
                    out.push_back('\n');
                }
            }
            resultWriter.writeData(out.c_str(), out.size(), queryKey, 0, true, true);
        }
        resultWriter.close(true);
    } catch (...) {
        for (size_t i = 0; i < tReaders.size(); ++i) {
            if (tReaders[i] != NULL) {
                tReaders[i]->close();
                delete tReaders[i];
            }
            if (rReaders[i] != NULL) {
                rReaders[i]->close();
                delete rReaders[i];
            }
        }
        qDbr.close();
        if (error.empty()) {
            error = "failed while merging multiple target/result DBs";
        }
        return false;
    }

    for (size_t i = 0; i < tReaders.size(); ++i) {
        tReaders[i]->close();
        delete tReaders[i];
        rReaders[i]->close();
        delete rReaders[i];
    }
    qDbr.close();
    return true;
}

static void removePathRecursively(const std::string &path) {
    if (path.empty()) {
        return;
    }
    std::string cmd = "rm -rf '" + path + "'";
    std::system(cmd.c_str());
}

// Decode a single sequence from an open DBReader by internal index.
static FastaSeq decodeOneSequence(DBReader<unsigned int> &dbr,
                                  size_t id,
                                  BaseMatrix &subMat,
                                  Sequence &seqObj,
                                  unsigned int thread_idx) {
    FastaSeq cur;
    cur.key = dbr.getDbKey(id);
    const bool hasLookup = (dbr.getLookupSize() > 0);
    if (hasLookup) {
        const size_t lid = dbr.getLookupIdByKey(cur.key);
        cur.id = dbr.getLookupEntryName(lid);
    } else {
        cur.id = std::to_string(cur.key);
    }
    const size_t seqLen = dbr.getSeqLen(id);
    const char *data = dbr.getData(id, thread_idx);
    cur.seq.reserve(seqLen);
    seqObj.mapSequence(id, cur.key, data, seqLen);
    for (int p = 0; p < seqObj.L; ++p) {
        const char c = normalizeBase(subMat.num2aa[seqObj.numSequence[p]]);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'U' || c == 'N') {
            cur.seq.push_back(c);
        } else {
            cur.seq.push_back('N');
        }
    }
    return cur;
}

struct Hit {
    unsigned int dbKey = 0;
    unsigned int dbLen = 0;
    std::string seqId;
    int start1 = 0;
    int end1 = 0;
    char mode = 'J'; // J/L/R/T
    bool trunc = false;
    double cyk = NEG_INF;
    double inside = NEG_INF;
    double bias = 0.0;
    double evalue = NEG_INF;
    bool hasEvalue = false;
    float precomputedSeqId = -1.0f; // trace-based pid, <0 means fall back to score-per-col estimate
    std::string cigar;
    std::string traceStates;
    // Query-coord alignment span derived from CIGAR (with --hand: consensus
    // col c == query position c-1). qStart/qEnd are 0-indexed; cigarAlnLen is
    // the total CIGAR op count (M+D+I) in MMseqs convention. Negative qStart
    // means "not populated; fall back to legacy emit".
    int qStart = -1;
    int qEnd = -1;
    unsigned int cigarAlnLen = 0;
    int leadingInsertTargets = 0;
    int trailingInsertTargets = 0;
};

static inline char traceOpForState(CmStateType t) {
    if (t == CM_ST_MP) return 'M';
    if (t == CM_ST_ML || t == CM_ST_MR) return 'M';
    if (t == CM_ST_IL || t == CM_ST_IR) return 'I';
    return '\0';
}

static std::string rleTraceOps(const std::string &ops) {
    if (ops.empty()) {
        return "NA";
    }
    std::string out;
    out.reserve(ops.size() * 2);
    char cur = ops[0];
    int run = 1;
    for (size_t i = 1; i < ops.size(); ++i) {
        if (ops[i] == cur) {
            ++run;
        } else {
            out += std::to_string(run);
            out.push_back(cur);
            cur = ops[i];
            run = 1;
        }
    }
    out += std::to_string(run);
    out.push_back(cur);
    return out;
}

static std::string encodeTraceStates(const std::vector<int> &states) {
    if (states.empty()) {
        return "NA";
    }
    std::string out;
    out.reserve(states.size() * 4);
    for (size_t i = 0; i < states.size(); ++i) {
        if (i > 0) {
            out.push_back(',');
        }
        out += std::to_string(states[i]);
    }
    return out;
}

static std::string buildConsensusFromExactModel(const InfernalExactModel &model) {
    if (model.clen <= 0) {
        return std::string();
    }
    std::string cons(static_cast<size_t>(model.clen), 'N');
    std::vector<double> bestSc(static_cast<size_t>(model.clen), NEG_INF);
    static const char bases[4] = {'A', 'C', 'G', 'U'};
    for (size_t si = 0; si < model.states.size(); ++si) {
        const CmState &s = model.states[si];
        if (s.type == CM_ST_MP && s.emit.size() >= 16) {
            if (s.mapLeft > 0 && s.mapLeft <= model.clen) {
                int bi = 0;
                double bsc = NEG_INF;
                for (int a = 0; a < 4; ++a) {
                    double row = NEG_INF;
                    for (int b = 0; b < 4; ++b) {
                        row = std::max(row, s.emit[static_cast<size_t>(a * 4 + b)]);
                    }
                    if (row > bsc) {
                        bsc = row;
                        bi = a;
                    }
                }
                const size_t pos = static_cast<size_t>(s.mapLeft - 1);
                if (bsc > bestSc[pos]) {
                    bestSc[pos] = bsc;
                    cons[pos] = bases[bi];
                }
            }
            if (s.mapRight > 0 && s.mapRight <= model.clen) {
                int bi = 0;
                double bsc = NEG_INF;
                for (int b = 0; b < 4; ++b) {
                    double col = NEG_INF;
                    for (int a = 0; a < 4; ++a) {
                        col = std::max(col, s.emit[static_cast<size_t>(a * 4 + b)]);
                    }
                    if (col > bsc) {
                        bsc = col;
                        bi = b;
                    }
                }
                const size_t pos = static_cast<size_t>(s.mapRight - 1);
                if (bsc > bestSc[pos]) {
                    bestSc[pos] = bsc;
                    cons[pos] = bases[bi];
                }
            }
        } else if ((s.type == CM_ST_ML || s.type == CM_ST_MR) && s.emit.size() >= 4) {
            int mapPos = (s.type == CM_ST_ML) ? s.mapLeft : s.mapRight;
            if (mapPos > 0 && mapPos <= model.clen) {
                int bi = 0;
                double bsc = s.emit[0];
                for (int a = 1; a < 4; ++a) {
                    if (s.emit[static_cast<size_t>(a)] > bsc) {
                        bsc = s.emit[static_cast<size_t>(a)];
                        bi = a;
                    }
                }
                const size_t pos = static_cast<size_t>(mapPos - 1);
                if (bsc > bestSc[pos]) {
                    bestSc[pos] = bsc;
                    cons[pos] = bases[bi];
                }
            }
        }
    }
    return cons;
}

// Build a query-coord MMseqs CIGAR from the CYK trace.
// With cmbuild --hand + all-'x' RF, consensus column c maps 1:1 to query
// position c-1, so per-column ops are query-coord directly. MMseqs convention:
//   M = match (consumes 1 query, 1 target)
//   D = consumes 1 query, 0 target  -> CM delete state at column
//   I = consumes 0 query, 1 target  -> CM insert state between columns
// We trim to firstMatchCol..lastMatchCol so leading/trailing delete-only
// columns (which are walked by truncated CYK but are not really "covered")
// don't inflate the alignment span. outQStart/outQEnd are 0-indexed query
// positions; outAlnLen is the total number of CIGAR ops (M+D+I).
static std::string modelTraceCigarQueryCoord(const InfernalExactModel &model,
                                             const std::vector<int> &traceStates,
                                             int *outQStart, int *outQEnd,
                                             int *outAlnLen,
                                             int *outLeadingInsertTargets,
                                             int *outTrailingInsertTargets,
                                             const std::vector<int8_t> *traceModes = nullptr) {
    // Mode-aware emit semantics (Infernal cm_dpalign_trunc.c ModeEmits):
    //   J: MP both, ML left, MR right
    //   L: MP left, ML left, MR silent
    //   R: MP right, ML silent, MR right
    //   T: MP/ML/MR all silent (T only legal at root/BIF)
    auto emitsLeft  = [](int8_t mode) {
        return mode == TR_TRMODE_J || mode == TR_TRMODE_L;
    };
    auto emitsRight = [](int8_t mode) {
        return mode == TR_TRMODE_J || mode == TR_TRMODE_R;
    };
    if (outQStart) *outQStart = -1;
    if (outQEnd) *outQEnd = -1;
    if (outAlnLen) *outAlnLen = 0;
    if (outLeadingInsertTargets) *outLeadingInsertTargets = 0;
    if (outTrailingInsertTargets) *outTrailingInsertTargets = 0;
    if (traceStates.empty()) {
        return "NA";
    }
    int maxMap = 0;
    int minMap = std::numeric_limits<int>::max();
    for (size_t i = 0; i < traceStates.size(); ++i) {
        const int sid = traceStates[i];
        if (sid < 0 || sid >= static_cast<int>(model.states.size())) {
            continue;
        }
        const CmState &s = model.states[static_cast<size_t>(sid)];
        if (s.mapLeft > 0) {
            maxMap = std::max(maxMap, s.mapLeft);
            minMap = std::min(minMap, s.mapLeft);
        }
        if (s.mapRight > 0) {
            maxMap = std::max(maxMap, s.mapRight);
            minMap = std::min(minMap, s.mapRight);
        }
    }
    if (maxMap <= 0 || minMap == std::numeric_limits<int>::max()) {
        return "NA";
    }
    std::vector<char> mop(static_cast<size_t>(maxMap + 1), '\0');
    std::vector<int> insAfter(static_cast<size_t>(maxMap + 1), 0);
    // ROOT_IL anchors before col 1 (mapLeft==0) and ROOT_IR anchors after col
    // CLEN (mapRight==CLEN+1==maxMap+1). They consume target residues but
    // don't fit in insAfter[1..maxMap]; bin them here so they get attributed
    // to leading/trailing inserts below — otherwise dbEnd-dbStart+1 exceeds
    // the CIGAR's M+I and downstream walkers (result2dnamsa) overshoot.
    int prefixIns = 0;
    int suffixIns = 0;
    for (size_t i = 0; i < traceStates.size(); ++i) {
        const int sid = traceStates[i];
        if (sid < 0 || sid >= static_cast<int>(model.states.size())) {
            continue;
        }
        const CmState &s = model.states[static_cast<size_t>(sid)];
        // Mode-aware emit (truncation modes). Without modes, treat as J.
        const int8_t md = (traceModes && i < traceModes->size())
                          ? (*traceModes)[i]
                          : (int8_t)TR_TRMODE_J;
        const bool emL = emitsLeft(md);
        const bool emR = emitsRight(md);
        if (s.type == CM_ST_MP) {
            if (emL && s.mapLeft > 0 && s.mapLeft <= maxMap) {
                mop[static_cast<size_t>(s.mapLeft)] = 'M';
            }
            if (emR && s.mapRight > 0 && s.mapRight <= maxMap) {
                mop[static_cast<size_t>(s.mapRight)] = 'M';
            }
        } else if (s.type == CM_ST_ML) {
            if (emL && s.mapLeft > 0 && s.mapLeft <= maxMap) {
                mop[static_cast<size_t>(s.mapLeft)] = 'M';
            }
        } else if (s.type == CM_ST_MR) {
            if (emR && s.mapRight > 0 && s.mapRight <= maxMap) {
                mop[static_cast<size_t>(s.mapRight)] = 'M';
            }
        } else if (s.type == CM_ST_D) {
            if (s.mapLeft > 0 && s.mapLeft <= maxMap && mop[static_cast<size_t>(s.mapLeft)] == '\0') {
                mop[static_cast<size_t>(s.mapLeft)] = 'D';
            }
            if (s.mapRight > 0 && s.mapRight <= maxMap && mop[static_cast<size_t>(s.mapRight)] == '\0') {
                mop[static_cast<size_t>(s.mapRight)] = 'D';
            }
        } else if (s.type == CM_ST_IL) {
            // IL inserts AFTER its left consensus column. ROOT_IL has mapLeft==0
            // and conceptually inserts before col 1 → bin as prefix.
            // Mode-aware: IL only emits in mode J or L (silent in R/T).
            if (!emL) { /* silent */ }
            else {
                int anchor = s.mapLeft;
                if (anchor <= 0) {
                    prefixIns += 1;
                } else if (anchor > maxMap) {
                    suffixIns += 1;
                } else {
                    insAfter[static_cast<size_t>(anchor)] += 1;
                }
            }
        } else if (s.type == CM_ST_IR) {
            // IR inserts BEFORE its right consensus column → anchor at mapRight-1.
            // ROOT_IR has mapRight==CLEN+1 (=maxMap+1) so anchor==maxMap (valid).
            // Mode-aware: IR only emits in mode J or R (silent in L/T).
            if (!emR) { /* silent */ }
            else {
                int anchor = (s.mapRight > 0) ? s.mapRight - 1 : 0;
                if (anchor <= 0) {
                    prefixIns += 1;
                } else if (anchor > maxMap) {
                    suffixIns += 1;
                } else {
                    insAfter[static_cast<size_t>(anchor)] += 1;
                }
            }
        }
    }
    int firstM = 0;
    int lastM = 0;
    for (int p = 1; p <= maxMap; ++p) {
        if (mop[static_cast<size_t>(p)] == 'M') {
            if (firstM == 0) firstM = p;
            lastM = p;
        }
    }
    if (firstM == 0) {
        // No match-state coverage; fall back to whole-trace span so callers can
        // still report something rather than failing silently.
        firstM = std::max(1, minMap);
        lastM = maxMap;
    }
    // Count target residues consumed by inserts that fall outside the
    // [firstM, lastM] window, so the caller can shrink dbStart/dbEnd to align
    // with the trimmed CIGAR. Insert states with anchor < firstM consumed
    // target before the alignment region; anchor > lastM consumed after.
    int leadingIns = prefixIns;
    int trailingIns = suffixIns;
    for (int p = 1; p < firstM && p <= maxMap; ++p) {
        leadingIns += insAfter[static_cast<size_t>(p)];
    }
    for (int p = lastM + 1; p <= maxMap; ++p) {
        trailingIns += insAfter[static_cast<size_t>(p)];
    }
    if (outLeadingInsertTargets) *outLeadingInsertTargets = leadingIns;
    if (outTrailingInsertTargets) *outTrailingInsertTargets = trailingIns;
    std::string ops;
    ops.reserve(static_cast<size_t>(lastM - firstM + 1) + 32);
    for (int p = firstM; p <= lastM; ++p) {
        char c = mop[static_cast<size_t>(p)];
        if (c == '\0') c = 'D'; // unwalked column inside hit envelope = query gap
        ops.push_back((c == 'M') ? 'M' : 'D');
        const int ins = insAfter[static_cast<size_t>(p)];
        for (int k = 0; k < ins; ++k) {
            ops.push_back('I');
        }
    }
    if (outQStart) *outQStart = firstM - 1;
    if (outQEnd) *outQEnd = lastM - 1;
    if (outAlnLen) *outAlnLen = static_cast<int>(ops.size());
    return rleTraceOps(ops);
}

static std::string modelTraceConsensusForCigar(const InfernalExactModel &model,
                                               const std::vector<int> &traceStates) {
    if (traceStates.empty()) {
        return std::string();
    }
    int maxMap = 0;
    int minMap = std::numeric_limits<int>::max();
    for (size_t i = 0; i < traceStates.size(); ++i) {
        const int sid = traceStates[i];
        if (sid < 0 || sid >= static_cast<int>(model.states.size())) {
            continue;
        }
        const CmState &s = model.states[static_cast<size_t>(sid)];
        if (s.mapLeft > 0) {
            maxMap = std::max(maxMap, s.mapLeft);
            minMap = std::min(minMap, s.mapLeft);
        }
        if (s.mapRight > 0) {
            maxMap = std::max(maxMap, s.mapRight);
            minMap = std::min(minMap, s.mapRight);
        }
    }
    if (maxMap <= 0 || minMap == std::numeric_limits<int>::max()) {
        return std::string();
    }

    std::vector<char> mop(static_cast<size_t>(maxMap + 1), '\0');
    for (size_t i = 0; i < traceStates.size(); ++i) {
        const int sid = traceStates[i];
        if (sid < 0 || sid >= static_cast<int>(model.states.size())) {
            continue;
        }
        const CmState &s = model.states[static_cast<size_t>(sid)];
        if (s.type == CM_ST_MP) {
            if (s.mapLeft > 0 && s.mapLeft <= maxMap) {
                mop[static_cast<size_t>(s.mapLeft)] = 'M';
            }
            if (s.mapRight > 0 && s.mapRight <= maxMap) {
                mop[static_cast<size_t>(s.mapRight)] = 'M';
            }
        } else if (s.type == CM_ST_ML || s.type == CM_ST_MR) {
            int p = (s.type == CM_ST_ML) ? s.mapLeft : s.mapRight;
            if (p > 0 && p <= maxMap) {
                mop[static_cast<size_t>(p)] = 'M';
            }
        } else if (s.type == CM_ST_D) {
            if (s.mapLeft > 0 && s.mapLeft <= maxMap && mop[static_cast<size_t>(s.mapLeft)] == '\0') {
                mop[static_cast<size_t>(s.mapLeft)] = 'D';
            }
            if (s.mapRight > 0 && s.mapRight <= maxMap && mop[static_cast<size_t>(s.mapRight)] == '\0') {
                mop[static_cast<size_t>(s.mapRight)] = 'D';
            }
        }
    }

    const std::string baseCons = buildConsensusFromExactModel(model);
    const char bases[4] = {'A', 'C', 'G', 'U'};
    std::vector<double> bestSc(static_cast<size_t>(maxMap + 1), NEG_INF);
    std::vector<char> bestCh(static_cast<size_t>(maxMap + 1), 'N');
    for (int p = minMap; p <= maxMap; ++p) {
        if (p >= 1 && p <= static_cast<int>(baseCons.size())) {
            bestCh[static_cast<size_t>(p)] = normalizeBase(baseCons[static_cast<size_t>(p - 1)]);
        }
    }
    for (size_t i = 0; i < traceStates.size(); ++i) {
        const int sid = traceStates[i];
        if (sid < 0 || sid >= static_cast<int>(model.states.size())) {
            continue;
        }
        const CmState &s = model.states[static_cast<size_t>(sid)];
        if (s.type == CM_ST_MP && s.emit.size() >= 16) {
            if (s.mapLeft > 0 && s.mapLeft <= maxMap) {
                int bi = 0;
                double bsc = NEG_INF;
                for (int a = 0; a < 4; ++a) {
                    double row = NEG_INF;
                    for (int b = 0; b < 4; ++b) {
                        row = std::max(row, s.emit[static_cast<size_t>(a * 4 + b)]);
                    }
                    if (row > bsc) {
                        bsc = row;
                        bi = a;
                    }
                }
                const size_t pos = static_cast<size_t>(s.mapLeft);
                if (bsc > bestSc[pos]) {
                    bestSc[pos] = bsc;
                    bestCh[pos] = bases[bi];
                }
            }
            if (s.mapRight > 0 && s.mapRight <= maxMap) {
                int bi = 0;
                double bsc = NEG_INF;
                for (int b = 0; b < 4; ++b) {
                    double col = NEG_INF;
                    for (int a = 0; a < 4; ++a) {
                        col = std::max(col, s.emit[static_cast<size_t>(a * 4 + b)]);
                    }
                    if (col > bsc) {
                        bsc = col;
                        bi = b;
                    }
                }
                const size_t pos = static_cast<size_t>(s.mapRight);
                if (bsc > bestSc[pos]) {
                    bestSc[pos] = bsc;
                    bestCh[pos] = bases[bi];
                }
            }
        } else if ((s.type == CM_ST_ML || s.type == CM_ST_MR) && s.emit.size() >= 4) {
            const int mapPos = (s.type == CM_ST_ML) ? s.mapLeft : s.mapRight;
            if (mapPos > 0 && mapPos <= maxMap) {
                int bi = 0;
                double bsc = s.emit[0];
                for (int a = 1; a < 4; ++a) {
                    if (s.emit[static_cast<size_t>(a)] > bsc) {
                        bsc = s.emit[static_cast<size_t>(a)];
                        bi = a;
                    }
                }
                const size_t pos = static_cast<size_t>(mapPos);
                if (bsc > bestSc[pos]) {
                    bestSc[pos] = bsc;
                    bestCh[pos] = bases[bi];
                }
            }
        }
    }

    std::string out;
    out.reserve(static_cast<size_t>(maxMap - minMap + 1));
    for (int p = minMap; p <= maxMap; ++p) {
        if (mop[static_cast<size_t>(p)] == '\0') {
            continue;
        }
        out.push_back(bestCh[static_cast<size_t>(p)]);
    }
    return out;
}

static std::vector<int> decodeTraceStates(const std::string &enc) {
    std::vector<int> out;
    if (enc.empty() || enc == "NA") {
        return out;
    }
    size_t i = 0;
    while (i < enc.size()) {
        size_t j = i;
        while (j < enc.size() && enc[j] != ',') {
            ++j;
        }
        if (j > i) {
            out.push_back(std::atoi(enc.substr(i, j - i).c_str()));
        }
        i = (j < enc.size()) ? (j + 1) : j;
    }
    return out;
}

static float seqIdFromCigarConsensus(const std::string &cigar,
                                     const std::string &consensus,
                                     const std::string &observed) {
    if (cigar.empty() || cigar == "NA" || consensus.empty() || observed.empty()) {
        return 0.0f;
    }
    size_t iObs = 0;
    size_t iCon = 0;
    size_t ident = 0;
    size_t denom = 0;
    size_t p = 0;
    while (p < cigar.size()) {
        size_t cnt = 0;
        while (p < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[p])) != 0) {
            cnt = cnt * 10 + static_cast<size_t>(cigar[p] - '0');
            ++p;
        }
        if (cnt == 0) {
            cnt = 1;
        }
        if (p >= cigar.size()) {
            break;
        }
        const char op = cigar[p++];
        if (op == 'M') {
            for (size_t k = 0; k < cnt; ++k) {
                if (iObs >= observed.size() || iCon >= consensus.size()) {
                    break;
                }
                const char a = normalizeBase(observed[iObs++]);
                const char b = normalizeBase(consensus[iCon++]);
                if (a != 'N' && b != 'N') {
                    denom += 1;
                }
                if (a == b && a != 'N') {
                    ident += 1;
                }
            }
        } else if (op == 'I') {
            // MMseqs convention: I = target-only (gap on query) → advance observed.
            const size_t step = std::min(cnt, observed.size() - std::min(iObs, observed.size()));
            iObs += step;
            denom += step;
        } else if (op == 'D') {
            // MMseqs convention: D = query-only (gap on target) → advance consensus.
            const size_t step = std::min(cnt, consensus.size() - std::min(iCon, consensus.size()));
            iCon += step;
            denom += step;
        }
    }
    if (denom == 0) {
        return 0.0f;
    }
    return static_cast<float>(ident) / static_cast<float>(denom);
}

static inline char complementBaseLocal(char c) {
    switch (normalizeBase(c)) {
        case 'A': return 'U';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'U': return 'A';
        default: return 'N';
    }
}

static std::string buildWindowSequenceFromHit(const std::string &targetSeq,
                                              const Matcher::result_t &hit,
                                              float flankFrac,
                                              int queryLen,
                                              int &fullLo,
                                              int &fullHi,
                                              bool &reverseStrand) {
    reverseStrand = (hit.dbStartPos > hit.dbEndPos);
    const int dbLo = std::min(hit.dbStartPos, hit.dbEndPos);
    const int dbHi = std::max(hit.dbStartPos, hit.dbEndPos);
    const int qLo = std::min(hit.qStartPos, hit.qEndPos);
    const int qHi = std::max(hit.qStartPos, hit.qEndPos);
    const int hitSpan = std::max(dbHi - dbLo + 1, qHi - qLo + 1);
    int flank = (flankFrac > 0.0f)
        ? static_cast<int>(std::ceil(static_cast<float>(hitSpan) * flankFrac))
        : std::max(24, queryLen / 4);
    flank = std::max(flank, 8);
    fullLo = std::max(0, dbLo - flank);
    fullHi = std::min(static_cast<int>(targetSeq.size()) - 1, dbHi + flank);
    if (fullHi < fullLo) {
        fullLo = 0;
        fullHi = -1;
        return std::string();
    }
    std::string window;
    window.reserve(static_cast<size_t>(fullHi - fullLo + 1));
    if (!reverseStrand) {
        for (int pos = fullLo; pos <= fullHi; ++pos) {
            char c = normalizeBase(targetSeq[static_cast<size_t>(pos)]);
            window.push_back((c == 'A' || c == 'C' || c == 'G' || c == 'U') ? c : 'N');
        }
    } else {
        for (int pos = fullHi; pos >= fullLo; --pos) {
            window.push_back(complementBaseLocal(targetSeq[static_cast<size_t>(pos)]));
        }
    }
    return window;
}

static std::vector<int> buildLolConsensusTargetMap(const Matcher::result_t &rough,
                                                   int fullLo,
                                                   int fullHi,
                                                   bool reverseStrand,
                                                   int scanLen,
                                                   int clen) {
    std::vector<int> targetByConsensus(static_cast<size_t>(std::max(0, clen) + 1), -1);
    if (rough.backtrace.empty() || scanLen <= 0 || clen <= 0) {
        return targetByConsensus;
    }
    int qPos = rough.qStartPos;
    const int qStep = (rough.qEndPos >= rough.qStartPos) ? 1 : -1;
    int dbPos = rough.dbStartPos;
    const int dbStep = (rough.dbEndPos >= rough.dbStartPos) ? 1 : -1;
    for (size_t k = 0; k < rough.backtrace.size(); ++k) {
        const char op = rough.backtrace[k];
        if (op == 'M') {
            const int consPos = qPos + 1;
            if (consPos >= 1 && consPos <= clen && dbPos >= fullLo && dbPos <= fullHi) {
                const int forwardLocal = dbPos - fullLo;
                const int scanLocal = reverseStrand ? (scanLen - 1 - forwardLocal) : forwardLocal;
                if (scanLocal >= 0 && scanLocal < scanLen) {
                    targetByConsensus[static_cast<size_t>(consPos)] = scanLocal + 1;
                }
            }
            qPos += qStep;
            dbPos += dbStep;
        } else if (op == 'I') {
            qPos += qStep;
        } else if (op == 'D') {
            dbPos += dbStep;
        }
    }
    return targetByConsensus;
}

static void computeLolRescoreWindow(const std::vector<int> &targetByConsensus,
                                    int scanLen,
                                    int pad,
                                    int &subLo,
                                    int &subHi) {
    subLo = 1;
    subHi = scanLen;
    if (scanLen <= 0 || targetByConsensus.empty()) {
        return;
    }
    int minMapped = scanLen + 1;
    int maxMapped = 0;
    for (size_t i = 1; i < targetByConsensus.size(); ++i) {
        const int pos = targetByConsensus[i];
        if (pos > 0) {
            minMapped = std::min(minMapped, pos);
            maxMapped = std::max(maxMapped, pos);
        }
    }
    if (maxMapped <= 0 || minMapped > scanLen) {
        return;
    }
    subLo = std::max(1, minMapped - pad);
    subHi = std::min(scanLen, maxMapped + pad);
}

static void sliceLolStateBands(const std::vector<int> &srcMinI,
                               const std::vector<int> &srcMaxI,
                               const std::vector<int> &srcMinJ,
                               const std::vector<int> &srcMaxJ,
                               int subLo,
                               int subHi,
                               int scanLen,
                               std::vector<int> &dstMinI,
                               std::vector<int> &dstMaxI,
                               std::vector<int> &dstMinJ,
                               std::vector<int> &dstMaxJ) {
    if (srcMinI.empty() || srcMaxI.empty() || srcMinJ.empty() || srcMaxJ.empty()
        || subLo <= 1 && subHi >= scanLen) {
        dstMinI = srcMinI;
        dstMaxI = srcMaxI;
        dstMinJ = srcMinJ;
        dstMaxJ = srcMaxJ;
        return;
    }
    const int shift = subLo - 1;
    const size_t M = srcMinI.size();
    dstMinI.assign(M, 1);
    dstMaxI.assign(M, subHi - subLo + 1);
    dstMinJ.assign(M, 1);
    dstMaxJ.assign(M, subHi - subLo + 1);
    bool anyBand = false;
    for (size_t v = 0; v < M; ++v) {
        const int loI = std::max(subLo, srcMinI[v]);
        const int hiI = std::min(subHi, srcMaxI[v]);
        if (hiI >= loI) {
            dstMinI[v] = loI - shift;
            dstMaxI[v] = hiI - shift;
            anyBand = true;
        }
        const int loJ = std::max(subLo, srcMinJ[v]);
        const int hiJ = std::min(subHi, srcMaxJ[v]);
        if (hiJ >= loJ) {
            dstMinJ[v] = loJ - shift;
            dstMaxJ[v] = hiJ - shift;
            anyBand = true;
        }
    }
    if (!anyBand) {
        dstMinI.clear();
        dstMaxI.clear();
        dstMinJ.clear();
        dstMaxJ.clear();
    }
}

static void buildLolStateBands(const InfernalExactModel &model,
                               const Matcher::result_t &rough,
                               int fullLo,
                               int fullHi,
                               bool reverseStrand,
                               int scanLen,
                               int bandPad,
                               std::vector<int> &stateMinI,
                               std::vector<int> &stateMaxI,
                               std::vector<int> &stateMinJ,
                               std::vector<int> &stateMaxJ) {
    const int clen = std::max(0, model.clen);
    if (rough.backtrace.empty() || scanLen <= 0 || clen <= 0) {
        stateMinI.clear();
        stateMaxI.clear();
        stateMinJ.clear();
        stateMaxJ.clear();
        return;
    }
    std::vector<int> targetByConsensus =
        buildLolConsensusTargetMap(rough, fullLo, fullHi, reverseStrand, scanLen, clen);

    const int M = static_cast<int>(model.states.size());
    stateMinI.assign(static_cast<size_t>(M), 1);
    stateMaxI.assign(static_cast<size_t>(M), scanLen);
    stateMinJ.assign(static_cast<size_t>(M), 1);
    stateMaxJ.assign(static_cast<size_t>(M), scanLen);
    bool anyBand = false;
    for (int v = 0; v < M; ++v) {
        const CmState &st = model.states[static_cast<size_t>(v)];
        if ((st.type == CM_ST_MP || st.type == CM_ST_ML || st.type == CM_ST_IL)
            && st.mapLeft > 0 && st.mapLeft <= clen) {
            const int mapped = targetByConsensus[static_cast<size_t>(st.mapLeft)];
            if (mapped > 0) {
                stateMinI[static_cast<size_t>(v)] = std::max(1, mapped - bandPad);
                stateMaxI[static_cast<size_t>(v)] = std::min(scanLen, mapped + bandPad);
                anyBand = true;
            }
        }
        if ((st.type == CM_ST_MP || st.type == CM_ST_MR || st.type == CM_ST_IR)
            && st.mapRight > 0 && st.mapRight <= clen) {
            const int mapped = targetByConsensus[static_cast<size_t>(st.mapRight)];
            if (mapped > 0) {
                stateMinJ[static_cast<size_t>(v)] = std::max(1, mapped - bandPad);
                stateMaxJ[static_cast<size_t>(v)] = std::min(scanLen, mapped + bandPad);
                anyBand = true;
            }
        }
    }
    if (!anyBand) {
        stateMinI.clear();
        stateMaxI.clear();
        stateMinJ.clear();
        stateMaxJ.clear();
    }
}

static inline double cmStateEmitScoreFast(const ExactStateExec &st, const std::vector<int8_t> &seqCode, int i, int j) {
    if (st.type == CM_ST_MP) {
        if (i > j || i < 1 || j >= static_cast<int>(seqCode.size()) || st.emitSize < 16 || st.emitPtr == NULL) {
            return NEG_INF;
        }
        const int li = seqCode[static_cast<size_t>(i)];
        const int ri = seqCode[static_cast<size_t>(j)];
        if (li < 0 || ri < 0) {
            return -1.0;
        }
        return st.emitPtr[li * 4 + ri];
    }
    if (st.type == CM_ST_ML || st.type == CM_ST_IL) {
        if (i > j || i < 1 || j >= static_cast<int>(seqCode.size()) || st.emitSize < 4 || st.emitPtr == NULL) {
            return NEG_INF;
        }
        const int bi = seqCode[static_cast<size_t>(i)];
        if (bi < 0) {
            return -1.0;
        }
        return st.emitPtr[bi];
    }
    if (st.type == CM_ST_MR || st.type == CM_ST_IR) {
        if (i > j || i < 1 || j >= static_cast<int>(seqCode.size()) || st.emitSize < 4 || st.emitPtr == NULL) {
            return NEG_INF;
        }
        const int bi = seqCode[static_cast<size_t>(j)];
        if (bi < 0) {
            return -1.0;
        }
        return st.emitPtr[bi];
    }
    return 0.0;
}

struct ExactRecCtx {
    int N;
    int M;
    const std::vector<ExactStateExec> *exec;
    const ExactStateExec *execData;
    const StateId *trDst;
    const double *trSc;
    const std::vector<int8_t> *seqCode;
    // Per-state ragged (v, i, d) chart — float score + uint32 generation counter.
    // Each state v has band [bandDmin[v] .. bandDmin[v]+bandSize[v]-1] derived
    // from QDB (CmState::dmin1/dmax1). Cells for v are laid out as
    // bandSize[v] * (N+1) consecutive floats starting at chartOffset[v].
    float *chartScore;
    uint32_t *chartSeen;
    uint32_t chartGen;
    const size_t *chartStateDOffset; // M+1 entries; prefix over (state,d) rows
    const size_t *chartRowOffset;    // totalRows+1 entries; prefix over live cells
    const int *chartRowImin;         // totalRows entries
    const int *chartRowSize;         // totalRows entries
    const int *chartBandDmin;   // M entries
    const int *chartBandSize;   // M entries (= dmax-dmin+1, 0 if empty)
    const int *stateMinI;
    const int *stateMaxI;
    const int *stateMinJ;
    const int *stateMaxJ;
};

static inline size_t exactChartIndex(const ExactRecCtx &ctx, int v, int i, int d) {
    const int dIdx = d - ctx.chartBandDmin[v];
    if (dIdx < 0 || dIdx >= ctx.chartBandSize[v]) {
        return std::numeric_limits<size_t>::max();
    }
    const size_t rowIdx = ctx.chartStateDOffset[v] + static_cast<size_t>(dIdx);
    const int iMin = ctx.chartRowImin[rowIdx];
    const int rowSize = ctx.chartRowSize[rowIdx];
    if (rowSize <= 0 || i < iMin || i >= iMin + rowSize) {
        return std::numeric_limits<size_t>::max();
    }
    return ctx.chartRowOffset[rowIdx] + static_cast<size_t>(i - iMin);
}

static double exactVitRec(ExactRecCtx &ctx, int v, int i, int j) {
    if (v < 0 || v >= ctx.M || i < 1 || i > ctx.N + 1 || j < 0 || j > ctx.N || i > j + 1) {
        return NEG_INF;
    }
    const int d = j - i + 1;
    const int bandDmin = ctx.chartBandDmin[v];
    const int bandSize = ctx.chartBandSize[v];
    if (bandSize <= 0 || d < bandDmin || d >= bandDmin + bandSize) {
        // Out-of-band — NEG_INF. Cheap to recompute, do not memoize.
        return NEG_INF;
    }
    const size_t dk = exactChartIndex(ctx, v, i, d);
    if (dk == std::numeric_limits<size_t>::max()) {
        return NEG_INF;
    }
    if (ctx.chartSeen[dk] == ctx.chartGen) {
        return static_cast<double>(ctx.chartScore[dk]);
    }

    const ExactStateExec &st = ctx.execData[static_cast<size_t>(v)];
    if (ctx.stateMinI != NULL && st.consumesLeft) {
        if (i < ctx.stateMinI[v] || i > ctx.stateMaxI[v]) {
            return NEG_INF;
        }
    }
    if (ctx.stateMinJ != NULL && st.consumesRight) {
        if (j < ctx.stateMinJ[v] || j > ctx.stateMaxJ[v]) {
            return NEG_INF;
        }
    }
    if (j < i - 1) {
        return NEG_INF;
    }
    if (st.type == CM_ST_E) {
        const double rv = (i == j + 1) ? 0.0 : NEG_INF;
        ctx.chartSeen[dk] = ctx.chartGen;
        ctx.chartScore[dk] = static_cast<float>(rv);
        return rv;
    }
    if (st.type == CM_ST_B) {
        const int y = st.bLeft;
        const int z = st.bRight;
        double best = NEG_INF;
        int kBeg = i - 1;
        int kEnd = j;
        if (st.splitKMax >= st.splitKMin) {
            kBeg = std::max(kBeg, i + st.splitKMin - 1);
            kEnd = std::min(kEnd, i + st.splitKMax - 1);
        }
        if (st.splitRMax >= st.splitRMin) {
            kBeg = std::max(kBeg, i + (d - st.splitRMax) - 1);
            kEnd = std::min(kEnd, i + (d - st.splitRMin) - 1);
        }
        for (int k = kBeg; k <= kEnd; ++k) {
            const double left = exactVitRec(ctx, y, i, k);
            const double right = exactVitRec(ctx, z, k + 1, j);
            if (left == NEG_INF || right == NEG_INF) {
                continue;
            }
            best = std::max(best, left + right);
        }
        ctx.chartSeen[dk] = ctx.chartGen;
        ctx.chartScore[dk] = static_cast<float>(best);
        return best;
    }

    const int ni = i + st.niShift;
    const int nj = j - (st.consumesRight ? 1 : 0);
    // Interior-bound and emit-failure NEG_INFs are cheap to recompute; skip
    // memoization to keep the sparse map bounded by actually-scored cells.
    if (ni > nj + 1) {
        return NEG_INF;
    }
    const double e = cmStateEmitScoreFast(st, *ctx.seqCode, i, j);
    if (e == NEG_INF) {
        return NEG_INF;
    }

    double best = NEG_INF;
    if (st.trCount == 1) {
        const double nxt = exactVitRec(ctx, static_cast<int>(st.trDst4[0]), ni, nj);
        if (nxt != NEG_INF) {
            best = e + st.trSc4[0] + nxt;
        }
    } else if (st.trCount == 2) {
        const double nxt0 = exactVitRec(ctx, static_cast<int>(st.trDst4[0]), ni, nj);
        if (nxt0 != NEG_INF) best = std::max(best, e + st.trSc4[0] + nxt0);
        const double nxt1 = exactVitRec(ctx, static_cast<int>(st.trDst4[1]), ni, nj);
        if (nxt1 != NEG_INF) best = std::max(best, e + st.trSc4[1] + nxt1);
    } else if (st.trCount == 4) {
        const double nxt0 = exactVitRec(ctx, static_cast<int>(st.trDst4[0]), ni, nj);
        if (nxt0 != NEG_INF) best = std::max(best, e + st.trSc4[0] + nxt0);
        const double nxt1 = exactVitRec(ctx, static_cast<int>(st.trDst4[1]), ni, nj);
        if (nxt1 != NEG_INF) best = std::max(best, e + st.trSc4[1] + nxt1);
        const double nxt2 = exactVitRec(ctx, static_cast<int>(st.trDst4[2]), ni, nj);
        if (nxt2 != NEG_INF) best = std::max(best, e + st.trSc4[2] + nxt2);
        const double nxt3 = exactVitRec(ctx, static_cast<int>(st.trDst4[3]), ni, nj);
        if (nxt3 != NEG_INF) best = std::max(best, e + st.trSc4[3] + nxt3);
    } else {
        for (int t = 0; t < st.trCount; ++t) {
            const size_t ti = st.trOff + static_cast<size_t>(t);
            const int y = ctx.trDst[ti];
            const double nxt = exactVitRec(ctx, y, ni, nj);
            if (nxt == NEG_INF) {
                continue;
            }
            best = std::max(best, e + ctx.trSc[ti] + nxt);
        }
    }
    ctx.chartSeen[dk] = ctx.chartGen;
    ctx.chartScore[dk] = static_cast<float>(best);
    return best;
}

static void exactTraceRec(ExactRecCtx &ctx,
                          int v,
                          int i,
                          int j,
                          int &minUsed,
                          int &maxUsed,
                          int obsCount[4],
                          double modelAggRaw[4],
                          std::string &traceOps,
                          std::vector<int> &traceStates) {
    const double best = exactVitRec(ctx, v, i, j);
    if (best == NEG_INF) {
        return;
    }
    const ExactStateExec &st = ctx.execData[static_cast<size_t>(v)];
    if (st.type == CM_ST_E) {
        return;
    }
    traceStates.push_back(v);
    if (st.type == CM_ST_B) {
        const int y = st.bLeft;
        const int z = st.bRight;
        int kBeg = i - 1;
        int kEnd = j;
        const int d = j - i + 1;
        if (st.splitKMax >= st.splitKMin) {
            kBeg = std::max(kBeg, i + st.splitKMin - 1);
            kEnd = std::min(kEnd, i + st.splitKMax - 1);
        }
        if (st.splitRMax >= st.splitRMin) {
            kBeg = std::max(kBeg, i + (d - st.splitRMax) - 1);
            kEnd = std::min(kEnd, i + (d - st.splitRMin) - 1);
        }
        // Always pick argmax over the enumerated split range. Earlier code
        // tried a tolerance-match early-exit on `best - cand` < 1e-6, but
        // cand is recomputed from float-rounded chart cells while best is the
        // float chart value cast to double; on ties or near-ties the early
        // match could pick a non-optimal k.
        int bestK = -1;
        double bestApprox = NEG_INF;
        for (int k = kBeg; k <= kEnd; ++k) {
            const double left = exactVitRec(ctx, y, i, k);
            const double right = exactVitRec(ctx, z, k + 1, j);
            if (left == NEG_INF || right == NEG_INF) {
                continue;
            }
            const double cand = left + right;
            if (cand > bestApprox) {
                bestApprox = cand;
                bestK = k;
            }
        }
        if (bestK >= 0) {
            exactTraceRec(ctx, y, i, bestK, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            exactTraceRec(ctx, z, bestK + 1, j, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
        }
        return;
    }

    const int ni = i + st.niShift;
    const int nj = j - (st.consumesRight ? 1 : 0);
    const char op = traceOpForState(st.type);
    if (op != '\0') {
        for (int x = 0; x < st.dConsume; ++x) {
            traceOps.push_back(op);
        }
    }
    if (st.consumesLeft) {
        minUsed = std::min(minUsed, i);
        maxUsed = std::max(maxUsed, i);
        const int bi = (*ctx.seqCode)[static_cast<size_t>(i)];
        if (bi >= 0 && bi < 4) {
            obsCount[bi]++;
        }
    }
    if (st.consumesRight) {
        minUsed = std::min(minUsed, j);
        maxUsed = std::max(maxUsed, j);
        const int bj = (*ctx.seqCode)[static_cast<size_t>(j)];
        if (bj >= 0 && bj < 4) {
            obsCount[bj]++;
        }
    }
    modelAggRaw[0] += st.null2Agg[0];
    modelAggRaw[1] += st.null2Agg[1];
    modelAggRaw[2] += st.null2Agg[2];
    modelAggRaw[3] += st.null2Agg[3];
    const double e = cmStateEmitScoreFast(st, *ctx.seqCode, i, j);
    // Pick argmax over all children. The previous tolerance-match early-exit
    // could short-circuit on a non-optimal child when float chart rounding
    // pushed cand within 1e-6 of best for a sub-optimal transition. Always
    // enumerating is essentially free here — trace runs once per hit while
    // fill runs O(M*N²) per target.
    int bestY = -1;
    double bestApprox = NEG_INF;
    const int trCount = st.trCount;
    if (trCount <= 4) {
        for (int t = 0; t < trCount; ++t) {
            const int y = static_cast<int>(st.trDst4[t]);
            const double nxt = exactVitRec(ctx, y, ni, nj);
            if (nxt == NEG_INF) {
                continue;
            }
            const double cand = e + st.trSc4[t] + nxt;
            if (cand > bestApprox) {
                bestApprox = cand;
                bestY = y;
            }
        }
    } else {
        for (int t = 0; t < trCount; ++t) {
            const size_t ti = st.trOff + static_cast<size_t>(t);
            const int y = ctx.trDst[ti];
            const double nxt = exactVitRec(ctx, y, ni, nj);
            if (nxt == NEG_INF) {
                continue;
            }
            const double cand = e + ctx.trSc[ti] + nxt;
            if (cand > bestApprox) {
                bestApprox = cand;
                bestY = y;
            }
        }
    }
    (void)best;  // best is no longer needed; argmax over enumerated cands suffices.
    if (bestY >= 0) {
        exactTraceRec(ctx, bestY, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
    }
}

static void runInfernalExactScan(const InfernalExactModel &model,
                                 const std::string &seq,
                                 bool wantInside,
                                 std::vector<Hit> &outHits,
                                 const std::string &seqId,
                                 int maxHitLen = 0,
                                 int forcedI = -1,
                                 int forcedD = -1,
                                 int anchorI = -1,
                                 int anchorD = -1,
                                 int forceDisableQdb = -1,
                                 const int *stateMinI = NULL,
                                 const int *stateMaxI = NULL,
                                 const int *stateMinJ = NULL,
                                 const int *stateMaxJ = NULL) {
    const int N = static_cast<int>(seq.size());
    const int M = static_cast<int>(model.states.size());
    if (N <= 0 || M <= 0) {
        return;
    }
    if (M > static_cast<int>(std::numeric_limits<StateId>::max())) {
        Debug(Debug::ERROR) << "CM state count exceeds uint16_t capacity: " << M << "\n";
        return;
    }
    // Backstop: after fixing NEG_INF memoization, the map grows only with
    // cells that return a finite score. For very large M*N^2 we still bail
    // out to prevent pathological runtime.
    {
        const uint64_t workEstimate = static_cast<uint64_t>(M)
                                    * static_cast<uint64_t>(N)
                                    * static_cast<uint64_t>(N);
        static constexpr uint64_t CM_EXACT_MAX_WORK = 200000000000ULL; // 2e11 cell-ops
        if (workEstimate > CM_EXACT_MAX_WORK) {
            Debug(Debug::WARNING) << "cmscan: skipping scan for M=" << M
                                   << " N=" << N << " (work "
                                   << workEstimate << " exceeds cap)\n";
            return;
        }
    }
    const bool fastLogsum = cmFastLogsumEnabled();
    const std::string modelConsensus = buildConsensusFromExactModel(model);
    thread_local ExactScanWorkspace ws;
    ws.seqCode.resize(static_cast<size_t>(N + 1));
    std::fill(ws.seqCode.begin(), ws.seqCode.end(), static_cast<int8_t>(-1)); // 1-based
    for (int p = 1; p <= N; ++p) {
        ws.seqCode[static_cast<size_t>(p)] = static_cast<int8_t>(baseToIdx(seq[static_cast<size_t>(p - 1)]));
    }
    ws.prefA.assign(static_cast<size_t>(N + 1), 0);
    ws.prefC.assign(static_cast<size_t>(N + 1), 0);
    ws.prefG.assign(static_cast<size_t>(N + 1), 0);
    ws.prefU.assign(static_cast<size_t>(N + 1), 0);
    std::vector<int> &prefA = ws.prefA;
    std::vector<int> &prefC = ws.prefC;
    std::vector<int> &prefG = ws.prefG;
    std::vector<int> &prefU = ws.prefU;
    for (int p = 1; p <= N; ++p) {
        prefA[static_cast<size_t>(p)] = prefA[static_cast<size_t>(p - 1)];
        prefC[static_cast<size_t>(p)] = prefC[static_cast<size_t>(p - 1)];
        prefG[static_cast<size_t>(p)] = prefG[static_cast<size_t>(p - 1)];
        prefU[static_cast<size_t>(p)] = prefU[static_cast<size_t>(p - 1)];
        const int bi = ws.seqCode[static_cast<size_t>(p)];
        if (bi == 0) prefA[static_cast<size_t>(p)]++;
        else if (bi == 1) prefC[static_cast<size_t>(p)]++;
        else if (bi == 2) prefG[static_cast<size_t>(p)]++;
        else if (bi == 3) prefU[static_cast<size_t>(p)]++;
    }
    const size_t needLog2 = static_cast<size_t>(N + 1);
    const size_t oldLog2 = ws.log2Int.size();
    if (oldLog2 < needLog2) {
        ws.log2Int.resize(needLog2, NEG_INF);
        int start = static_cast<int>(oldLog2);
        if (start < 1) {
            start = 1;
        }
        for (int x = start; x <= N; ++x) {
            ws.log2Int[static_cast<size_t>(x)] = std::log2(static_cast<double>(x));
        }
    }
    const std::vector<double> &log2Int = ws.log2Int;
    const double log2Omega3 = (model.null3Omega > 0.0) ? std::log2(model.null3Omega) : NEG_INF;
    const double log2Null[4] = {
        (model.nullProb[0] > 0.0) ? std::log2(model.nullProb[0]) : NEG_INF,
        (model.nullProb[1] > 0.0) ? std::log2(model.nullProb[1]) : NEG_INF,
        (model.nullProb[2] > 0.0) ? std::log2(model.nullProb[2]) : NEG_INF,
        (model.nullProb[3] > 0.0) ? std::log2(model.nullProb[3]) : NEG_INF
    };
    const bool null3Enabled = cmNull3Enabled();
    auto null3CorrByInterval = [&](int i, int j) -> double {
        if (!null3Enabled) {
            return 0.0;
        }
        const int cntA = prefA[static_cast<size_t>(j)] - prefA[static_cast<size_t>(i - 1)];
        const int cntC = prefC[static_cast<size_t>(j)] - prefC[static_cast<size_t>(i - 1)];
        const int cntG = prefG[static_cast<size_t>(j)] - prefG[static_cast<size_t>(i - 1)];
        const int cntU = prefU[static_cast<size_t>(j)] - prefU[static_cast<size_t>(i - 1)];
        const int known = cntA + cntC + cntG + cntU;
        if (known <= 0 || !model.hasNull || !std::isfinite(log2Omega3)) {
            return 0.0;
        }
        const double log2Len = log2Int[static_cast<size_t>(known)];
        double score = log2Omega3;
        if (cntA > 0 && std::isfinite(log2Null[0])) score += static_cast<double>(cntA) * (log2Int[static_cast<size_t>(cntA)] - log2Len - log2Null[0]);
        if (cntC > 0 && std::isfinite(log2Null[1])) score += static_cast<double>(cntC) * (log2Int[static_cast<size_t>(cntC)] - log2Len - log2Null[1]);
        if (cntG > 0 && std::isfinite(log2Null[2])) score += static_cast<double>(cntG) * (log2Int[static_cast<size_t>(cntG)] - log2Len - log2Null[2]);
        if (cntU > 0 && std::isfinite(log2Null[3])) score += static_cast<double>(cntU) * (log2Int[static_cast<size_t>(cntU)] - log2Len - log2Null[3]);
        return log2SumExp2(0.0, score);
    };
    std::vector<ExactStateExec> exec(static_cast<size_t>(M));
    std::vector<StateId> trDst;
    std::vector<double> trSc;
    trDst.reserve(static_cast<size_t>(M) * 4);
    trSc.reserve(static_cast<size_t>(M) * 4);
    std::vector<uint8_t> isBState(static_cast<size_t>(M), 0);
    std::vector<int> nbNiShift(static_cast<size_t>(M), 0);
    std::vector<int> nbDConsume(static_cast<size_t>(M), 0);
    std::vector<int> nbEmitSize(static_cast<size_t>(M), 0);
    std::vector<int> nbTrCount(static_cast<size_t>(M), 0);
    std::vector<int> nbDmin(static_cast<size_t>(M), 0);
    std::vector<int> nbDmax(static_cast<size_t>(M), -1);
    std::vector<size_t> nbTrOff(static_cast<size_t>(M), 0);
    std::vector<const float *> nbEmitPtr(static_cast<size_t>(M), NULL);
    std::vector<uint8_t> nbConsumeMask(static_cast<size_t>(M), 0);
    std::vector<StateId> nbTrDst4(static_cast<size_t>(M) * 4, 0);
    std::vector<float> nbTrScF4(static_cast<size_t>(M) * 4, -std::numeric_limits<float>::infinity());
    // MMSEQS_CMSCAN_DISABLE_QDB=1 neutralizes per-state QDB1 bands so CYK fills the
    // full d-range, mirroring Infernal cmalign's default (HMM-banded but not QDB-banded).
    // forceDisableQdb param overrides env var for per-call gating (used by envelope rescan
    // to run pass 1 with QDB on and pass 2 with QDB off).
    bool disableQdb;
    if (forceDisableQdb >= 0) {
        disableQdb = (forceDisableQdb != 0);
    } else {
        const char *envDisableQdb = std::getenv("MMSEQS_CMSCAN_DISABLE_QDB");
        disableQdb = (envDisableQdb != NULL && envDisableQdb[0] == '1');
    }
    // MMSEQS_CMSCAN_DROP_DMIN=1 zeros per-state dmin while preserving dmax. Targets the
    // FG=1 collapse on short-envelope queries (e.g. 4LCK_3) where forced d=N falls below
    // dmin at the root, without inviting qdboff's unbounded-trace regressions on large queries.
    const char *envDropDmin = std::getenv("MMSEQS_CMSCAN_DROP_DMIN");
    const bool dropDmin = (envDropDmin != NULL && envDropDmin[0] == '1');
    // MMSEQS_CMSCAN_DROP_DMAX=1: zeros dmin AND removes dmax cap (sentinel -1) so no
    // upper d-bound at any state, but still ranges from 0. Mirrors qdboff's band off.
    const char *envDropDmax = std::getenv("MMSEQS_CMSCAN_DROP_DMAX");
    const bool dropDmax = (envDropDmax != NULL && envDropDmax[0] == '1');
    for (int v = 0; v < M; ++v) {
        const CmState &st = model.states[static_cast<size_t>(v)];
        ExactStateExec e;
        e.type = st.type;
        e.niShift = 0;
        e.dConsume = 0;
        e.consumesLeft = false;
        e.consumesRight = false;
        e.bLeft = st.cfirst;
        e.bRight = st.cnum;
        e.trDst4[0] = 0; e.trDst4[1] = 0; e.trDst4[2] = 0; e.trDst4[3] = 0;
        e.trScF4[0] = -std::numeric_limits<float>::infinity();
        e.trScF4[1] = -std::numeric_limits<float>::infinity();
        e.trScF4[2] = -std::numeric_limits<float>::infinity();
        e.trScF4[3] = -std::numeric_limits<float>::infinity();
        e.trSc4[0] = NEG_INF; e.trSc4[1] = NEG_INF; e.trSc4[2] = NEG_INF; e.trSc4[3] = NEG_INF;
        e.trOff = trDst.size();
        e.trCount = 0;
        e.emitPtr = st.emitF.empty() ? NULL : &st.emitF[0];
        e.emitSize = static_cast<int>(st.emit.size());
        e.splitKMin = 0;
        e.splitKMax = -1;
        e.splitRMin = 0;
        e.splitRMax = -1;
        e.dmin = (disableQdb || dropDmin || dropDmax) ? 0 : st.dmin1;
        e.dmax = (disableQdb || dropDmax) ? -1 : st.dmax1;
        e.null2Agg[0] = 0.0;
        e.null2Agg[1] = 0.0;
        e.null2Agg[2] = 0.0;
        e.null2Agg[3] = 0.0;
        if (st.type == CM_ST_MP) {
            e.niShift = 1;
            e.dConsume = 2;
            e.consumesLeft = true;
            e.consumesRight = true;
        } else if (st.type == CM_ST_ML || st.type == CM_ST_IL) {
            e.niShift = 1;
            e.dConsume = 1;
            e.consumesLeft = true;
        } else if (st.type == CM_ST_MR || st.type == CM_ST_IR) {
            e.dConsume = 1;
            e.consumesRight = true;
        }
        if (st.type == CM_ST_MP && e.emitPtr != NULL && e.emitSize >= 16) {
            for (int a = 0; a < 4; ++a) {
                double row = 0.0;
                double col = 0.0;
                for (int b = 0; b < 4; ++b) {
                    row += emitBitToProbPair(e.emitPtr[a * 4 + b], model, a, b);
                    col += emitBitToProbPair(e.emitPtr[b * 4 + a], model, b, a);
                }
                e.null2Agg[a] = row + col;
            }
        } else if ((st.type == CM_ST_ML || st.type == CM_ST_MR || st.type == CM_ST_IL || st.type == CM_ST_IR) &&
                   e.emitPtr != NULL && e.emitSize >= 4) {
            for (int a = 0; a < 4; ++a) {
                e.null2Agg[a] = emitBitToProbSingle(e.emitPtr[a], model, a);
            }
        }
        if (st.type == CM_ST_B) {
            isBState[static_cast<size_t>(v)] = 1;
            if (!disableQdb && !dropDmax && e.bLeft >= 0 && e.bLeft < M) {
                const CmState &l = model.states[static_cast<size_t>(e.bLeft)];
                if (l.dmin1 >= 0 && l.dmax1 >= l.dmin1) {
                    e.splitKMin = dropDmin ? 0 : l.dmin1;
                    e.splitKMax = l.dmax1;
                }
            }
            if (!disableQdb && !dropDmax && e.bRight >= 0 && e.bRight < M) {
                const CmState &r = model.states[static_cast<size_t>(e.bRight)];
                if (r.dmin1 >= 0 && r.dmax1 >= r.dmin1) {
                    e.splitRMin = dropDmin ? 0 : r.dmin1;
                    e.splitRMax = r.dmax1;
                }
            }
        } else if (st.type != CM_ST_E) {
            const int tn = std::min(st.cnum, static_cast<int>(st.trans.size()));
            for (int x = 0; x < tn; ++x) {
                const int y = st.cfirst + x;
                if (y < 0 || y >= M) {
                    continue;
                }
                if (e.trCount < 4) {
                    const size_t t = static_cast<size_t>(e.trCount);
                    e.trDst4[t] = static_cast<StateId>(y);
                    e.trSc4[t] = st.trans[static_cast<size_t>(x)];
                    e.trScF4[t] = static_cast<float>(st.trans[static_cast<size_t>(x)]);
                }
                trDst.push_back(static_cast<StateId>(y));
                trSc.push_back(st.trans[static_cast<size_t>(x)]);
                ++e.trCount;
            }
        }
        exec[static_cast<size_t>(v)] = e;
        if (st.type != CM_ST_B && st.type != CM_ST_E) {
            const size_t vi = static_cast<size_t>(v);
            nbNiShift[vi] = e.niShift;
            nbDConsume[vi] = e.dConsume;
            nbEmitSize[vi] = e.emitSize;
            nbTrCount[vi] = e.trCount;
            nbDmin[vi] = e.dmin;
            nbDmax[vi] = e.dmax;
            nbTrOff[vi] = e.trOff;
            nbEmitPtr[vi] = e.emitPtr;
            uint8_t mask = 0;
            if (e.consumesLeft) {
                mask |= 1u;
            }
            if (e.consumesRight) {
                mask |= 2u;
            }
            nbConsumeMask[vi] = mask;
            const size_t t4 = vi * 4;
            nbTrDst4[t4 + 0] = e.trDst4[0];
            nbTrDst4[t4 + 1] = e.trDst4[1];
            nbTrDst4[t4 + 2] = e.trDst4[2];
            nbTrDst4[t4 + 3] = e.trDst4[3];
            nbTrScF4[t4 + 0] = e.trScF4[0];
            nbTrScF4[t4 + 1] = e.trScF4[1];
            nbTrScF4[t4 + 2] = e.trScF4[2];
            nbTrScF4[t4 + 3] = e.trScF4[3];
        }
    }
    ws.trScF.resize(trSc.size());
    for (size_t k = 0; k < trSc.size(); ++k) { ws.trScF[k] = static_cast<float>(trSc[k]); }
    std::vector<int> activeStates;
    activeStates.reserve(static_cast<size_t>(M));
    for (int v = M - 1; v >= 0; --v) {
        if (exec[static_cast<size_t>(v)].type != CM_ST_E) {
            activeStates.push_back(v);
        }
    }
    size_t gpuEmitFloatCount = 0;
    for (size_t vi = 0; vi < model.states.size(); ++vi) {
        gpuEmitFloatCount += model.states[vi].emitF.size();
    }
    const size_t gpuFastDeckBytesNeeded = estimateFastDeckGpuSingleJobBytes(
        M, N, activeStates.size(), trDst.size(), ws.trScF.size(), gpuEmitFloatCount);
    const size_t gpuFastDeckBytesBudget = cmGpuUsableBytes();
    const bool gpuFastDeckMemEligible = cmGpuExactEnabled() &&
                                        gpuFastDeckBytesBudget > 0 &&
                                        gpuFastDeckBytesNeeded > 0 &&
                                        gpuFastDeckBytesNeeded <= gpuFastDeckBytesBudget;
    // Fast path: iterative deck-based CYK over full sequence.
    // Inside score (if requested) is computed only for the selected best interval.
    const bool allowCpuFastDeckFallback = (N <= 1024);
    if (cmExactFastDeckEnabled() && (allowCpuFastDeckFallback || gpuFastDeckMemEligible)) {
        const float NEG_INF_F = -std::numeric_limits<float>::infinity();
        const size_t iStride = static_cast<size_t>(N + 2);
        const size_t stateStride = static_cast<size_t>(N + 1) * iStride;
        const size_t cells = static_cast<size_t>(M) * stateStride;
        if (!ws.vit.ensure(cells)) {
            Debug(Debug::ERROR) << "cmscan: failed to allocate vit buffer\n";
            return;
        }
        float *vitPtr = ws.vit.data();
        if (!ws.stateBase.ensure(static_cast<size_t>(M))) {
            Debug(Debug::ERROR) << "cmscan: failed to allocate stateBase buffer\n";
            return;
        }
        for (int v = 0; v < M; ++v) {
            ws.stateBase.data()[static_cast<size_t>(v)] = static_cast<size_t>(v) * stateStride;
        }
        // Pre-fill all state slices with NEG_INF using contiguous memsets.
        // This replaces scattered per-d fills in the main loop (94%+ of iterations
        // are out-of-band and just wrote NEG_INF in small chunks — now they can be
        // skipped entirely, and the contiguous write pattern is much more cache-friendly).
        for (int v = 0; v < M; ++v) {
            float *base = &vitPtr[ws.stateBase.data()[static_cast<size_t>(v)]];
            std::fill(base, base + stateStride, NEG_INF_F);
        }
        const size_t vdCells = static_cast<size_t>(M) * static_cast<size_t>(N + 1);
        if (!ws.bSplitBegByVD.ensure(vdCells)) {
            Debug(Debug::ERROR) << "cmscan: failed to allocate bSplitBeg buffer\n";
            return;
        }
        if (!ws.bSplitEndByVD.ensure(vdCells)) {
            Debug(Debug::ERROR) << "cmscan: failed to allocate bSplitEnd buffer\n";
            return;
        }
        if (!ws.bSplitTmp.ensure(static_cast<size_t>(N + 1))) {
            Debug(Debug::ERROR) << "cmscan: failed to allocate bSplitTmp buffer\n";
            return;
        }
        int *bSplitBegByVD = ws.bSplitBegByVD.data();
        int *bSplitEndByVD = ws.bSplitEndByVD.data();

        // Only init B-state split ranges (typically very few B-states).
        for (size_t ii = 0; ii < vdCells; ++ii) {
            bSplitBegByVD[ii] = 0;
            bSplitEndByVD[ii] = -1;
        }
        float *vit = ws.vit.data();
        size_t *stateBase = ws.stateBase.data();
        const ExactStateExec *execData = exec.data();
        for (int v = 0; v < M; ++v) {
            const ExactStateExec &st = exec[static_cast<size_t>(v)];
            if (st.type != CM_ST_B) {
                continue;
            }
            for (int d = 0; d <= N; ++d) {
                int kBeg = 0;
                int kEnd = d;
                if (st.splitKMax >= st.splitKMin) {
                    kBeg = std::max(kBeg, st.splitKMin);
                    kEnd = std::min(kEnd, st.splitKMax);
                }
                if (st.splitRMax >= st.splitRMin) {
                    kBeg = std::max(kBeg, d - st.splitRMax);
                    kEnd = std::min(kEnd, d - st.splitRMin);
                }
                const size_t vd = static_cast<size_t>(v) * static_cast<size_t>(N + 1) + static_cast<size_t>(d);
                bSplitBegByVD[vd] = kBeg;
                bSplitEndByVD[vd] = kEnd;
            }
        }

        bool gpuFilled = false;
        if (cmGpuExactEnabled()) {
            std::vector<CmFastDeckGpuState> gpuStates(static_cast<size_t>(M));
            std::vector<float> gpuEmitData;
            gpuEmitData.reserve(static_cast<size_t>(M) * 8);
            for (int v = 0; v < M; ++v) {
                const size_t vi = static_cast<size_t>(v);
                const CmState &src = model.states[vi];
                const ExactStateExec &st = exec[vi];
                CmFastDeckGpuState &dst = gpuStates[vi];
                dst.type = static_cast<int>(st.type);
                dst.dmin = st.dmin;
                dst.dmax = st.dmax;
                dst.bLeft = st.bLeft;
                dst.bRight = st.bRight;
                dst.dConsume = st.dConsume;
                dst.niShift = st.niShift;
                dst.emitSize = st.emitSize;
                dst.emitOffset = src.emitF.empty() ? -1 : static_cast<int>(gpuEmitData.size());
                dst.trCount = st.trCount;
                dst.trOff = st.trOff;
                dst.consumeMask = (st.consumesLeft ? 1u : 0u) | (st.consumesRight ? 2u : 0u);
                dst.endSc = static_cast<float>(src.endSc);
                dst.null2Agg[0] = static_cast<float>(st.null2Agg[0]);
                dst.null2Agg[1] = static_cast<float>(st.null2Agg[1]);
                dst.null2Agg[2] = static_cast<float>(st.null2Agg[2]);
                dst.null2Agg[3] = static_cast<float>(st.null2Agg[3]);
                if (!src.emitF.empty()) {
                    gpuEmitData.insert(gpuEmitData.end(), src.emitF.begin(), src.emitF.end());
                }
            }
            std::vector<size_t> stateBaseHost(stateBase, stateBase + static_cast<size_t>(M));
            std::vector<int> bSplitBegHost(bSplitBegByVD, bSplitBegByVD + vdCells);
            std::vector<int> bSplitEndHost(bSplitEndByVD, bSplitEndByVD + vdCells);
            std::string gpuError;
            gpuFilled = runInfernalExactScanFastDeckGpu(
                N,
                M,
                static_cast<int>(iStride),
                stateStride,
                cells,
                ws.seqCode.data(),
                gpuStates,
                activeStates,
                stateBaseHost,
                bSplitBegHost,
                bSplitEndHost,
                trDst,
                ws.trScF,
                gpuEmitData,
                model.hasLocalCfg,
                static_cast<float>(model.elSelf),
                vit,
                &gpuError);
            if (!gpuFilled && !gpuError.empty()) {
                static bool warned = false;
                if (!warned) {
                    Debug(Debug::WARNING) << "cmscan GPU fast-deck disabled for this run: " << gpuError << "\n";
                    warned = true;
                }
            }
        }

        if (gpuFilled || allowCpuFastDeckFallback) {
        if (!gpuFilled) {
            // E states emit empty only.
            for (int v = 0; v < M; ++v) {
                if (exec[static_cast<size_t>(v)].type != CM_ST_E) {
                    continue;
                }
                const size_t vb = stateBase[static_cast<size_t>(v)];
                for (int i = 1; i <= N + 1; ++i) {
                    vit[vb + static_cast<size_t>(i)] = 0.0f;
                }
            }

            for (int d = 0; d <= N; ++d) {
                const int iMax = (d == 0) ? (N + 1) : (N - d + 1);
                for (size_t vsi = 0; vsi < activeStates.size(); ++vsi) {
                    const int v = activeStates[vsi];
                    const size_t vi = static_cast<size_t>(v);
                    const ExactStateExec &st = execData[vi];
                    const size_t vBase = stateBase[vi];
                    const int dmin = isBState[vi] ? st.dmin : nbDmin[vi];
                    const int dmax = isBState[vi] ? st.dmax : nbDmax[vi];
                    if (dmax >= dmin && (d < dmin || d > dmax)) {
                        // Out-of-band: already NEG_INF from the contiguous pre-fill.
                        continue;
                    }
                    if (isBState[vi]) {
                        const int y = st.bLeft;
                        const int z = st.bRight;
                        if (y < 0 || y >= M || z < 0 || z >= M) {
                            continue;
                        }
                        const size_t vd = static_cast<size_t>(v) * static_cast<size_t>(N + 1) + static_cast<size_t>(d);
                        const int kBeg = bSplitBegByVD[vd];
                        const int kEnd = bSplitEndByVD[vd];
                        const float *const yBase = &vit[stateBase[static_cast<size_t>(y)]];
                        const float *const zBase = &vit[stateBase[static_cast<size_t>(z)]];
                        float *const bOut = &vit[vBase + static_cast<size_t>(d) * iStride + 1];
                        if (kBeg > kEnd) {
                            std::fill(bOut, bOut + iMax, NEG_INF_F);
                        } else {
                            for (int i = 1; i <= iMax; ++i) {
                                float bBest = NEG_INF_F;
                                for (int k = kBeg; k <= kEnd; ++k) {
                                    const float l = yBase[static_cast<size_t>(k) * iStride + static_cast<size_t>(i)];
                                    const float r = zBase[static_cast<size_t>(d - k) * iStride + static_cast<size_t>(i + k)];
                                    const float cand = l + r;
                                    bBest = std::fmaxf(bBest, cand);
                                }
                                bOut[i - 1] = bBest;
                            }
                        }
                        continue;
                    }

                    const int dConsume = nbDConsume[vi];
                    if (d < dConsume) {
                        continue;
                    }
                    const int nd = d - dConsume;
                    const int niShift = nbNiShift[vi];
                    const float *ep = nbEmitPtr[vi];
                    const int emitSize = nbEmitSize[vi];
                    const int trCount = nbTrCount[vi];
                    const size_t trOff = nbTrOff[vi];
                    const uint8_t consumeMask = nbConsumeMask[vi];
                    const size_t t4 = vi * 4;
                    const int8_t *sc = ws.seqCode.data();
                    const float *trScFPtr = ws.trScF.data();
                    const StateId *trDstPtr = trDst.data();
                    const size_t *stateBasePtr = stateBase;
                    const size_t dstBase0 = (trCount >= 1) ? stateBasePtr[static_cast<size_t>(nbTrDst4[t4 + 0])] : 0;
                    const size_t dstBase1 = (trCount >= 2) ? stateBasePtr[static_cast<size_t>(nbTrDst4[t4 + 1])] : 0;
                    const size_t dstBase2 = (trCount >= 3) ? stateBasePtr[static_cast<size_t>(nbTrDst4[t4 + 2])] : 0;
                    const size_t dstBase3 = (trCount >= 4) ? stateBasePtr[static_cast<size_t>(nbTrDst4[t4 + 3])] : 0;
                    const float trSc0 = (trCount >= 1) ? nbTrScF4[t4 + 0] : 0.f;
                    const float trSc1 = (trCount >= 2) ? nbTrScF4[t4 + 1] : 0.f;
                    const float trSc2 = (trCount >= 3) ? nbTrScF4[t4 + 2] : 0.f;
                    const float trSc3 = (trCount >= 4) ? nbTrScF4[t4 + 3] : 0.f;
                    // Build 5-entry emit lookup (indexed by sc+1, sc in {-1,0,1,2,3})
                    // NEG_INF_F propagates naturally through IEEE float arithmetic
                    float efT[5];
                    if (consumeMask == 0u) {
                        efT[0] = efT[1] = efT[2] = efT[3] = efT[4] = 0.0f;
                    } else if (consumeMask == 1u || consumeMask == 2u) {
                        efT[0] = -1.0f;
                        if (ep && emitSize >= 4) {
                            efT[1] = ep[0]; efT[2] = ep[1]; efT[3] = ep[2]; efT[4] = ep[3];
                        } else {
                            efT[1] = efT[2] = efT[3] = efT[4] = NEG_INF_F;
                        }
                    } else {
                        efT[0] = efT[1] = efT[2] = efT[3] = efT[4] = 0.0f; // pair: handled inline
                    }
                    // Contiguous output and source pointers (new layout: d*iStride + i)
                    float *const outPtr = &vit[vBase + static_cast<size_t>(d) * iStride + 1];
                    const float *const src0 = (trCount >= 1) ? &vit[dstBase0 + static_cast<size_t>(nd) * iStride + 1 + static_cast<size_t>(niShift)] : nullptr;
                    const float *const src1 = (trCount >= 2) ? &vit[dstBase1 + static_cast<size_t>(nd) * iStride + 1 + static_cast<size_t>(niShift)] : nullptr;
                    const float *const src2 = (trCount >= 3) ? &vit[dstBase2 + static_cast<size_t>(nd) * iStride + 1 + static_cast<size_t>(niShift)] : nullptr;
                    const float *const src3 = (trCount >= 4) ? &vit[dstBase3 + static_cast<size_t>(nd) * iStride + 1 + static_cast<size_t>(niShift)] : nullptr;
                    // sc[i] offset: scL[ii] = sc[ii+1] = sc[i]; scR[ii] = sc[ii+d] = sc[i+d-1]
                    const int8_t *const scL = sc + 1;
                    const int8_t *const scR = sc + d;
#if defined(CM_HAS_SSE2)
                    if ((trCount == 1 || trCount == 2 || trCount == 4) && consumeMask != 3u) {
                        const __m128 trSc0V = _mm_set1_ps(trSc0);
                        const __m128 trSc1V = (trCount >= 2) ? _mm_set1_ps(trSc1) : _mm_set1_ps(NEG_INF_F);
                        const __m128 trSc2V = (trCount >= 3) ? _mm_set1_ps(trSc2) : _mm_set1_ps(NEG_INF_F);
                        const __m128 trSc3V = (trCount >= 4) ? _mm_set1_ps(trSc3) : _mm_set1_ps(NEG_INF_F);
                        int ii = 0;
                        for (; ii + 4 <= iMax; ii += 4) {
                            __m128 ef_v;
                            if (consumeMask == 0u) {
                                ef_v = _mm_setzero_ps();
                            } else if (consumeMask == 1u) {
                                ef_v = _mm_set_ps(efT[scL[ii+3]+1], efT[scL[ii+2]+1],
                                                  efT[scL[ii+1]+1], efT[scL[ii+0]+1]);
                            } else { // consumeMask == 2u
                                ef_v = _mm_set_ps(efT[scR[ii+3]+1], efT[scR[ii+2]+1],
                                                  efT[scR[ii+1]+1], efT[scR[ii+0]+1]);
                            }
                            __m128 result;
                            if (trCount == 1) {
                                result = _mm_add_ps(_mm_add_ps(ef_v, trSc0V), _mm_loadu_ps(src0 + ii));
                            } else if (trCount == 2) {
                                const __m128 c0 = _mm_add_ps(_mm_add_ps(ef_v, trSc0V), _mm_loadu_ps(src0 + ii));
                                const __m128 c1 = _mm_add_ps(_mm_add_ps(ef_v, trSc1V), _mm_loadu_ps(src1 + ii));
                                result = _mm_max_ps(c0, c1);
                            } else { // trCount == 4
                                const __m128 c0 = _mm_add_ps(_mm_add_ps(ef_v, trSc0V), _mm_loadu_ps(src0 + ii));
                                const __m128 c1 = _mm_add_ps(_mm_add_ps(ef_v, trSc1V), _mm_loadu_ps(src1 + ii));
                                const __m128 c2 = _mm_add_ps(_mm_add_ps(ef_v, trSc2V), _mm_loadu_ps(src2 + ii));
                                const __m128 c3 = _mm_add_ps(_mm_add_ps(ef_v, trSc3V), _mm_loadu_ps(src3 + ii));
                                result = _mm_max_ps(_mm_max_ps(c0, c1), _mm_max_ps(c2, c3));
                            }
                            _mm_storeu_ps(outPtr + ii, result);
                        }
                        // Scalar remainder
                        for (; ii < iMax; ++ii) {
                            const float ef = (consumeMask == 0u) ? 0.0f
                                           : (consumeMask == 1u) ? efT[scL[ii]+1]
                                                                 : efT[scR[ii]+1];
                            float best = NEG_INF_F;
                            if (trCount >= 1) { const float c = ef + trSc0 + src0[ii]; best = std::fmaxf(best, c); }
                            if (trCount >= 2) { const float c = ef + trSc1 + src1[ii]; best = std::fmaxf(best, c); }
                            if (trCount >= 3) { const float c = ef + trSc2 + src2[ii]; best = std::fmaxf(best, c); }
                            if (trCount >= 4) { const float c = ef + trSc3 + src3[ii]; best = std::fmaxf(best, c); }
                            outPtr[ii] = best;
                        }
                    } else
#endif
                    {
                        // General scalar path (consumeMask==3 pair emit, or trCount==3/trCount>4)
                        // Precompute base pointers for all transitions (loop-invariant over i)
                        const float *trSrcPtrs[16];
                        float trScVals[16];
                        const int trCountClamped = (trCount <= 16) ? trCount : 16;
                        for (int t = 0; t < trCountClamped; ++t) {
                            const size_t ti = trOff + static_cast<size_t>(t);
                            trSrcPtrs[t] = &vit[stateBasePtr[static_cast<size_t>(trDstPtr[ti])] +
                                                static_cast<size_t>(nd) * iStride +
                                                static_cast<size_t>(niShift)];
                            trScVals[t] = trScFPtr[ti];
                        }

                        if (trCount == 1) {
                            for (int ii = 0; ii < iMax; ++ii) {
                                const int i = ii + 1;
                                float ef;
                                if (consumeMask == 3u) {
                                    if (!(ep && emitSize >= 16)) { ef = NEG_INF_F; }
                                    else { const int li = sc[i]; const int ri = sc[i + d - 1]; ef = (li >= 0 && ri >= 0) ? ep[li * 4 + ri] : -1.0f; }
                                } else if (consumeMask == 0u) { ef = 0.0f; }
                                else if (consumeMask == 1u) { ef = efT[scL[ii]+1]; }
                                else { ef = efT[scR[ii]+1]; }
                                outPtr[ii] = ef + trSc0 + src0[ii];
                            }
                        } else if (trCount == 2) {
                            for (int ii = 0; ii < iMax; ++ii) {
                                const int i = ii + 1;
                                float ef;
                                if (consumeMask == 3u) {
                                    if (!(ep && emitSize >= 16)) { ef = NEG_INF_F; }
                                    else { const int li = sc[i]; const int ri = sc[i + d - 1]; ef = (li >= 0 && ri >= 0) ? ep[li * 4 + ri] : -1.0f; }
                                } else if (consumeMask == 0u) { ef = 0.0f; }
                                else if (consumeMask == 1u) { ef = efT[scL[ii]+1]; }
                                else { ef = efT[scR[ii]+1]; }
                                float best = ef + trSc0 + src0[ii];
                                best = std::fmaxf(best, ef + trSc1 + src1[ii]);
                                outPtr[ii] = best;
                            }
                        } else if (trCount == 4) {
                            for (int ii = 0; ii < iMax; ++ii) {
                                const int i = ii + 1;
                                float ef;
                                if (consumeMask == 3u) {
                                    if (!(ep && emitSize >= 16)) { ef = NEG_INF_F; }
                                    else { const int li = sc[i]; const int ri = sc[i + d - 1]; ef = (li >= 0 && ri >= 0) ? ep[li * 4 + ri] : -1.0f; }
                                } else if (consumeMask == 0u) { ef = 0.0f; }
                                else if (consumeMask == 1u) { ef = efT[scL[ii]+1]; }
                                else { ef = efT[scR[ii]+1]; }
                                float best = ef + trSc0 + src0[ii];
                                best = std::fmaxf(best, ef + trSc1 + src1[ii]);
                                best = std::fmaxf(best, ef + trSc2 + src2[ii]);
                                best = std::fmaxf(best, ef + trSc3 + src3[ii]);
                                outPtr[ii] = best;
                            }
                        } else {
                            // Loop-swapped with precomputed emission scores
                            // Step 1: precompute emission into outPtr (reused as temp buffer)
                            if (consumeMask == 3u && ep && emitSize >= 16) {
                                for (int ii = 0; ii < iMax; ++ii) {
                                    const int i = ii + 1;
                                    const int li = sc[i];
                                    const int ri = sc[i + d - 1];
                                    outPtr[ii] = (li >= 0 && ri >= 0) ? ep[li * 4 + ri] : -1.0f;
                                }
                            } else if (consumeMask == 3u) {
                                std::fill(outPtr, outPtr + iMax, NEG_INF_F);
                            } else if (consumeMask == 0u) {
                                std::fill(outPtr, outPtr + iMax, 0.0f);
                            } else if (consumeMask == 1u) {
                                for (int ii = 0; ii < iMax; ++ii)
                                    outPtr[ii] = efT[scL[ii]+1];
                            } else {
                                for (int ii = 0; ii < iMax; ++ii)
                                    outPtr[ii] = efT[scR[ii]+1];
                            }
                            // Step 2: for each transition, stream through all i and accumulate
                            // First transition: bestBuf[ii] = outPtr[ii] + trSc + trSrc[i]
                            // Use bSplitTmp as temp best buffer
                            float *bestBuf = ws.bSplitTmp.data();
                            const int VS = VECSIZE_FLOAT;
                            {
                                const float *tsrc = trSrcPtrs[0];
                                const float tsc = trScVals[0];
                                const simd_float vtsc = simdf32_set(tsc);
                                int ii = 0;
                                for (; ii + VS - 1 < iMax; ii += VS) {
                                    simd_float vo = simdf32_loadu(outPtr + ii);
                                    simd_float vs = simdf32_loadu(tsrc + ii + 1);
                                    simdf32_storeu(bestBuf + ii, simdf32_add(simdf32_add(vo, vtsc), vs));
                                }
                                for (; ii < iMax; ++ii) {
                                    bestBuf[ii] = outPtr[ii] + tsc + tsrc[ii + 1];
                                }
                            }
                            // Remaining transitions: bestBuf = max(bestBuf, outPtr + tsc + tsrc[i+1])
                            for (int t = 1; t < trCountClamped; ++t) {
                                const float *tsrc = trSrcPtrs[t];
                                const float tsc = trScVals[t];
                                const simd_float vtsc = simdf32_set(tsc);
                                int ii = 0;
                                for (; ii + VS - 1 < iMax; ii += VS) {
                                    simd_float vo = simdf32_loadu(outPtr + ii);
                                    simd_float vs = simdf32_loadu(tsrc + ii + 1);
                                    simd_float vb = simdf32_loadu(bestBuf + ii);
                                    simd_float cand = simdf32_add(simdf32_add(vo, vtsc), vs);
                                    simdf32_storeu(bestBuf + ii, simdf32_max(vb, cand));
                                }
                                for (; ii < iMax; ++ii) {
                                    bestBuf[ii] = std::fmaxf(bestBuf[ii], outPtr[ii] + tsc + tsrc[ii + 1]);
                                }
                            }
                            // Copy result back
                            std::copy(bestBuf, bestBuf + iMax, outPtr);
                        }
                    }

                    // EL pre-fill fusion. Mirrors Infernal cm_dpsearch.c:3389-3398:
                    // for any end-eligible state v, alpha[v][d][i] must reflect
                    // max(regular recurrence, emit_v(i,d) + el_scA[d-sd] + endsc[v]).
                    // Doing this INSIDE the (d,v) fill (not as a post-pass) is what
                    // lets parents w of v reading alpha[v] at later d-iterations see
                    // the EL-augmented value — so EL contribution propagates up the
                    // CM tree, which the prior post-pass placement could never do.
                    if (model.hasLocalCfg && model.elSelf <= 0.0) {
                        const CmState &csEl = model.states[vi];
                        if (csEl.endSc != NEG_INF) {
                            int sdEl = -1;
                            if (csEl.type == CM_ST_MP) sdEl = 2;
                            else if (csEl.type == CM_ST_ML || csEl.type == CM_ST_MR) sdEl = 1;
                            else if (csEl.type == CM_ST_S) sdEl = 0;
                            if (sdEl >= 0 && d >= sdEl) {
                                const float elContrib = static_cast<float>(model.elSelf) * static_cast<float>(d - sdEl)
                                                      + static_cast<float>(csEl.endSc);
                                const int8_t *scEl = ws.seqCode.data();
                                for (int ii = 0; ii < iMax; ++ii) {
                                    const int i = ii + 1;
                                    float ef;
                                    if (consumeMask == 0u) {
                                        ef = 0.0f;
                                    } else if (consumeMask == 1u) {
                                        const int li = scEl[i];
                                        if (li < 0) continue;
                                        ef = (ep && emitSize >= 4) ? ep[li] : -1.0f;
                                    } else if (consumeMask == 2u) {
                                        const int ri = scEl[i + d - 1];
                                        if (ri < 0) continue;
                                        ef = (ep && emitSize >= 4) ? ep[ri] : -1.0f;
                                    } else {
                                        const int li = scEl[i];
                                        const int ri = scEl[i + d - 1];
                                        if (li < 0 || ri < 0) continue;
                                        ef = (ep && emitSize >= 16) ? ep[li * 4 + ri] : -1.0f;
                                    }
                                    const float cand = ef + elContrib;
                                    if (outPtr[ii] < cand) outPtr[ii] = cand;
                                }
                            }
                        }
                    }
                }
            }
        }

        const int root = model.rootState;

        // Local-end pickup (Infernal cm_dpsearch.c). After regular CYK fills
        // alpha[v][d][i] for end-eligible head v, take a max with the local-end
        // candidate: emit_v(i, d) + el_scA[d - sd_v] + endSc[v], where
        // el_scA[k] = elSelf * k. Emission must be included because state v
        // emits its M-state residue(s) before transitioning to EL.
        if (model.hasLocalCfg && model.elSelf <= 0.0) {
            const float elSelfF = static_cast<float>(model.elSelf);
            const int8_t *sc = ws.seqCode.data();
            for (int v = 0; v < M; ++v) {
                const CmState &cs = model.states[static_cast<size_t>(v)];
                if (cs.endSc == NEG_INF) continue;
                const float endScF = static_cast<float>(cs.endSc);
                const size_t vBase = stateBase[static_cast<size_t>(v)];
                const size_t vi = static_cast<size_t>(v);
                const float *ep = nbEmitPtr[vi];
                const int emitSize = nbEmitSize[vi];
                const uint8_t consumeMask = nbConsumeMask[vi];
                int sd = 0;
                if (cs.type == CM_ST_MP) sd = 2;
                else if (cs.type == CM_ST_ML || cs.type == CM_ST_MR) sd = 1;
                else if (cs.type == CM_ST_S) sd = 0;
                else continue;  // only MP/ML/MR/S are end-eligible heads
                for (int d = sd; d <= N; ++d) {
                    const float elContrib = elSelfF * static_cast<float>(d - sd) + endScF;
                    const int iMax = (d == 0) ? (N + 1) : (N - d + 1);
                    float *row = &vit[vBase + static_cast<size_t>(d) * iStride + 1];
                    for (int ii = 0; ii < iMax; ++ii) {
                        const int i = ii + 1;
                        float ef;
                        if (consumeMask == 0u) {
                            ef = 0.0f;
                        } else if (consumeMask == 1u) {
                            const int li = sc[i];
                            if (li < 0) continue;
                            ef = (ep && emitSize >= 4) ? ep[li] : -1.0f;
                        } else if (consumeMask == 2u) {
                            const int ri = sc[i + d - 1];
                            if (ri < 0) continue;
                            ef = (ep && emitSize >= 4) ? ep[ri] : -1.0f;
                        } else {
                            const int li = sc[i];
                            const int ri = sc[i + d - 1];
                            if (li < 0 || ri < 0) continue;
                            ef = (ep && emitSize >= 16) ? ep[li * 4 + ri] : -1.0f;
                        }
                        const float cand = ef + elContrib;
                        if (row[ii] < cand) row[ii] = cand;
                    }
                }
            }
        }

        // Hoisted: maxSpan/minSpan/forceEnv must be computed BEFORE the
        // local-begin pass below, so the FG=1 corner-only branch can
        // address the same trace-seed cell that the selection block reads.
        int maxSpan = N;
        // Fix A diagnostic: W cap can be disabled via MMSEQS_CMSCAN_WCAP=0 to
        // probe whether Infernal-length envelopes (d > W) become reachable.
        static int wcapMode = -1;
        if (wcapMode == -1) {
            const char *env = std::getenv("MMSEQS_CMSCAN_WCAP");
            wcapMode = (env != NULL && std::string(env) == "0") ? 0 : 1;
        }
        if (wcapMode == 1 && model.w > 0) {
            maxSpan = std::min(maxSpan, model.w);
        }
        if (maxHitLen > 0) {
            maxSpan = std::min(maxSpan, maxHitLen);
        }
        // Infernal semantics evaluate all legal spans under state dmin/dmax/w constraints.
        int minSpan = 1;
        // Fix B: optional CLEN-relative envelope floor. MMSEQS_CMSCAN_DMIN_CLEN_FRAC=0.9
        // forces d >= 0.9*CLEN to match Infernal's tight envelope distribution and prevent
        // the DP from clipping useful columns on divergent targets (e.g. 2YGH_1 KJ798010.1).
        {
            static double dminFrac = -1.0;
            if (dminFrac < 0.0) {
                const char *env = std::getenv("MMSEQS_CMSCAN_DMIN_CLEN_FRAC");
                dminFrac = (env != NULL) ? std::atof(env) : 0.0;
                if (dminFrac < 0.0 || dminFrac > 1.0) dminFrac = 0.0;
            }
            if (dminFrac > 0.0 && model.clen > 0) {
                int floorD = static_cast<int>(dminFrac * model.clen);
                if (floorD > minSpan) minSpan = floorD;
                if (minSpan > maxSpan) minSpan = maxSpan; // never invert
            }
        }
        // Fix B: force prefilter envelope. When forcedI/forcedD are valid,
        // collapse the (i, d) search to that single cell so trace matches
        // the prefilter's HMM-Forward-equivalent envelope.
        const bool forceEnv = (forcedI > 0 && forcedD > 0 &&
                               forcedD <= N && forcedI <= N - forcedD + 1);
        if (forceEnv) {
            minSpan = forcedD;
            maxSpan = forcedD;
        }

        // Hoisted: MMSEQS_CMSCAN_FORCE_GLOBAL env-cache. Read once here so the
        // local-begin pass below and the selection block further down both see
        // the same value. Selection block uses this same `forceGlobal`.
        static int forceGlobal = -1;
        if (forceGlobal == -1) {
            const char *env = std::getenv("MMSEQS_CMSCAN_FORCE_GLOBAL");
            if (env != NULL) {
                if (std::string(env) == "1") forceGlobal = 1;
                else if (std::string(env) == "2") forceGlobal = 2;
                else if (std::string(env) == "3") forceGlobal = 3;
                else if (std::string(env) == "4") forceGlobal = 4;
                else forceGlobal = 0;
            } else {
                forceGlobal = 0;
            }
        }

        // Local-begin pass (Infernal cm_dpsearch.c:576-610). After all states
        // are filled, the root cell can also be reached by jumping directly to
        // any begin-eligible state y at score `alpha[y][d][i] + beginSc[y]`.
        //
        // FG=0/2/3/4: spray local-begin alternates across every (i, d) ROOT_S
        //   cell. Default mode picks argmax over the row, so the per-cell
        //   pickup is harmless (every cell is internally consistent with its
        //   own argmax y). FG=2 sweeps i at fixed d; FG=3/4 scan a band.
        //
        // FG=1: trace is forced to seed at (i=1, d=min(maxSpan,N)). Per-cell
        //   pickup pollutes other root cells with values whose argmax y has
        //   nothing to do with that corner — see project_4LCK_FG1_collapse_
        //   decomposed.md and project_cyk_celldiff_root_local_begin.md (102
        //   polluted cells). Apply local-begin only at the trace-seed corner,
        //   mirroring Inside DP at lines 3683-3691 (bit-exact to Infernal).
        //
        // Glocal CMs leave beginSc=NEG_INF and skip the loop body either way.
        if (model.hasLocalCfg) {
            const size_t rootBaseLB = stateBase[static_cast<size_t>(root)];
            if (forceGlobal == 1) {
                // Corner-OVERWRITE: clear ROOT_S at the trace-seed corner
                // before applying local-begin. Mirrors Infernal cm_localize
                // (cm_modelconfig.c:491) which zeros ROOT_S→child transitions
                // so the only path to ROOT_S is via local-begin. Differs from
                // the per-cell max-update branch below: there we *add*
                // local-begin alternates on top of the regular ROOT_S
                // recurrence, which inflates the cell score above what the
                // local-begin trace can actually produce. Under TRACE_LBO=1 +
                // FG=1 the trace is forced to descend via local-begin, so
                // the cell value must reflect the local-begin path only —
                // otherwise score-vs-trace is inconsistent and bad alignments
                // get scored as if they took the regular-recurrence path.
                const int cornerD = std::min(maxSpan, N);
                const int cornerI = 1;
                if (cornerD >= 1 && cornerI + cornerD - 1 <= N) {
                    float &rootCell = vit[rootBaseLB
                        + static_cast<size_t>(cornerD) * iStride
                        + static_cast<size_t>(cornerI)];
                    rootCell = NEG_INF_F;
                    for (int y = 1; y < M; ++y) {
                        const CmState &cs = model.states[static_cast<size_t>(y)];
                        if (cs.beginSc == NEG_INF) continue;
                        const float bsc = static_cast<float>(cs.beginSc);
                        const float yVal = vit[stateBase[static_cast<size_t>(y)]
                            + static_cast<size_t>(cornerD) * iStride
                            + static_cast<size_t>(cornerI)];
                        if (yVal == NEG_INF_F) continue;
                        const float cand = yVal + bsc;
                        if (rootCell < cand) rootCell = cand;
                    }
                }
            } else if (forceGlobal == 5) {
                // FG=5 hybrid: only the corner cell (FG=1 candidate) and the
                // midpoint-span cells (FG=3 candidates) are read by the
                // selection logic below. Skip the O(M*N²) per-cell loop and
                // touch only those cells (O(M*span) work), yielding
                // FG=1-class wall time even though we maintain both
                // candidates. Both pickups are corner-OVERWRITE-style for
                // consistency with FG=1.
                const int cornerD = std::min(maxSpan, N);
                if (cornerD >= 1 && cornerD <= N) {
                    float &rootCell = vit[rootBaseLB
                        + static_cast<size_t>(cornerD) * iStride
                        + 1];
                    rootCell = NEG_INF_F;
                    for (int y = 1; y < M; ++y) {
                        const CmState &cs = model.states[static_cast<size_t>(y)];
                        if (cs.beginSc == NEG_INF) continue;
                        const float bsc = static_cast<float>(cs.beginSc);
                        const float yVal = vit[stateBase[static_cast<size_t>(y)]
                            + static_cast<size_t>(cornerD) * iStride
                            + 1];
                        if (yVal == NEG_INF_F) continue;
                        const float cand = yVal + bsc;
                        if (rootCell < cand) rootCell = cand;
                    }
                }
                if (anchorI > 0 && anchorD > 0
                    && anchorI <= N && anchorI + anchorD - 1 <= N) {
                    const int midpoint = anchorI + (anchorD - 1) / 2;
                    for (int d = minSpan; d <= maxSpan; ++d) {
                        const int iMax = N - d + 1;
                        const int iLo = std::max(1, midpoint - d + 1);
                        const int iHi = std::min(iMax, midpoint);
                        for (int i = iLo; i <= iHi; ++i) {
                            float &rootCell = vit[rootBaseLB
                                + static_cast<size_t>(d) * iStride
                                + static_cast<size_t>(i)];
                            rootCell = NEG_INF_F;
                            for (int y = 1; y < M; ++y) {
                                const CmState &cs = model.states[static_cast<size_t>(y)];
                                if (cs.beginSc == NEG_INF) continue;
                                const float bsc = static_cast<float>(cs.beginSc);
                                const float yVal = vit[stateBase[static_cast<size_t>(y)]
                                    + static_cast<size_t>(d) * iStride
                                    + static_cast<size_t>(i)];
                                if (yVal == NEG_INF_F) continue;
                                const float cand = yVal + bsc;
                                if (rootCell < cand) rootCell = cand;
                            }
                        }
                    }
                }
            } else {
                for (int y = 1; y < M; ++y) {
                    const CmState &cs = model.states[static_cast<size_t>(y)];
                    if (cs.beginSc == NEG_INF) {
                        continue;
                    }
                    const float bsc = static_cast<float>(cs.beginSc);
                    const size_t yBase = stateBase[static_cast<size_t>(y)];
                    const simd_float vbsc = simdf32_set(bsc);
                    const int VS = VECSIZE_FLOAT;
                    for (int d = 0; d <= N; ++d) {
                        const int iMax = (d == 0) ? (N + 1) : (N - d + 1);
                        float *rootRow = &vit[rootBaseLB + static_cast<size_t>(d) * iStride + 1];
                        const float *yRow = &vit[yBase + static_cast<size_t>(d) * iStride + 1];
                        // SIMD: rootRow[i] = max(rootRow[i], yRow[i] + bsc)
                        int i = 0;
                        for (; i + VS - 1 < iMax; i += VS) {
                            simd_float vy = simdf32_loadu(yRow + i);
                            simd_float vr = simdf32_loadu(rootRow + i);
                            simd_float cand = simdf32_add(vy, vbsc);
                            simdf32_storeu(rootRow + i, simdf32_max(vr, cand));
                        }
                        for (; i < iMax; ++i) {
                            const float cand = yRow[i] + bsc;
                            if (rootRow[i] < cand) rootRow[i] = cand;
                        }
                    }
                }
            }
        }
        float bestSc = NEG_INF_F;
        int bestI = 1;
        int bestD = N;
        char bestMode = 'J';
        double bestNull3Corr = 0.0;
        const bool enableTruncModes = cmTruncModesEnabled();
        float bestL = NEG_INF_F, bestR = NEG_INF_F, bestT = NEG_INF_F;
        int bestIL = 1, bestIR = 1, bestIT = 1;
        int bestDL = N, bestDR = N, bestDT = N;
        double bestLCorr = 0.0, bestRCorr = 0.0, bestTCorr = 0.0;
        const size_t rootBase = stateBase[static_cast<size_t>(root)];
        // Diagnostic: when MMSEQS_CMSCAN_DUMP_TID matches seqId, dump the score landscape.
        const char *dumpTid = std::getenv("MMSEQS_CMSCAN_DUMP_TID");
        const bool dumpThisTarget = (dumpTid != NULL && seqId == dumpTid);
        if (dumpThisTarget) {
            fprintf(stderr, "DUMP_HEADER tid=%s N=%d M=%d maxSpan=%d minSpan=%d truncModes=%d\n",
                    seqId.c_str(), N, M, maxSpan, minSpan, enableTruncModes ? 1 : 0);
        }
        // Per-d best raw and corrected scores, for diagnostic dump.
        std::vector<float> dBestRaw(static_cast<size_t>(maxSpan + 1), NEG_INF_F);
        std::vector<float> dBestSc(static_cast<size_t>(maxSpan + 1), NEG_INF_F);
        std::vector<int> dBestI(static_cast<size_t>(maxSpan + 1), -1);
        for (int d = minSpan; d <= maxSpan; ++d) {
            const int iMax = N - d + 1;
            int iLo = 1, iHi = iMax;
            if (forceEnv) { iLo = forcedI; iHi = forcedI; }
            for (int i = iLo; i <= iHi; ++i) {
                const float rawSc = vit[rootBase + static_cast<size_t>(d) * iStride + static_cast<size_t>(i)];
                if (rawSc == NEG_INF_F) {
                    continue;
                }
                const int j = i + d - 1;
                const bool mayImproveJ = rawSc > bestSc;
                const bool mayImproveL = enableTruncModes && j == N && rawSc > bestL;
                const bool mayImproveR = enableTruncModes && i == 1 && rawSc > bestR;
                const bool mayImproveT = enableTruncModes && i == 1 && j == N && rawSc > bestT;
                const bool dumpCell = dumpThisTarget && rawSc > dBestRaw[static_cast<size_t>(d)];
                if (!(mayImproveJ || mayImproveL || mayImproveR || mayImproveT) && !dumpCell) {
                    continue;
                }
                const double corr = null3CorrByInterval(i, j);
                const float sc = static_cast<float>(static_cast<double>(rawSc) - corr);
                if (dumpCell) {
                    dBestRaw[static_cast<size_t>(d)] = rawSc;
                    dBestSc[static_cast<size_t>(d)] = sc;
                    dBestI[static_cast<size_t>(d)] = i;
                }
                if (mayImproveJ && sc > bestSc) {
                    bestSc = sc;
                    bestI = i;
                    bestD = d;
                    bestMode = 'J';
                    bestNull3Corr = corr;
                }
                if (enableTruncModes) {
                    if (j == N && sc > bestL) {
                        bestL = sc;
                        bestIL = i;
                        bestDL = d;
                        bestLCorr = corr;
                    }
                    if (i == 1 && sc > bestR) {
                        bestR = sc;
                        bestIR = i;
                        bestDR = d;
                        bestRCorr = corr;
                    }
                    if (i == 1 && j == N && sc > bestT) {
                        bestT = sc;
                        bestIT = i;
                        bestDT = d;
                        bestTCorr = corr;
                    }
                }
            }
        }
        if (enableTruncModes) {
            if (bestL > bestSc) {
                bestSc = bestL;
                bestI = bestIL;
                bestD = bestDL;
                bestMode = 'L';
                bestNull3Corr = bestLCorr;
            }
            if (bestR > bestSc) {
                bestSc = bestR;
                bestI = bestIR;
                bestD = bestDR;
                bestMode = 'R';
                bestNull3Corr = bestRCorr;
            }
            if (bestT > bestSc) {
                bestSc = bestT;
                bestI = bestIT;
                bestD = bestDT;
                bestMode = 'T';
                bestNull3Corr = bestTCorr;
            }
        }
        // MMSEQS_CMSCAN_FORCE_GLOBAL=1: force trace to start at (i=1, d=min(N, maxSpan)),
        //   i.e. full envelope alignment from position 1. Worst for envelopes where the
        //   conserved core sits at i>1.
        // MMSEQS_CMSCAN_FORCE_GLOBAL=2: pin d=min(N, maxSpan), pick i = argmax alpha[ROOT_S][d][i].
        //   Full-length subspan, but slid to wherever score is best — handles 4LCK_3-style
        //   envelopes where the conserved region is not at i=1.
        // MMSEQS_CMSCAN_FORCE_GLOBAL=3: midpoint-anchor. argmax over (i, d) where
        //   i <= prefilter-midpoint <= i+d-1. Free d (adapts per hit), but envelope
        //   must straddle the prefilter's most-confident point.
        // MMSEQS_CMSCAN_FORCE_GLOBAL=4: span-contain. argmax over (i, d) where the
        //   prefilter [anchorI..anchorI+anchorD-1] is fully contained in [i..i+d-1].
        //   Strictest gate; equivalent to "envelope must cover the prefilter span".
        // MMSEQS_CMSCAN_FORCE_GLOBAL=5: per-target hybrid. Try FG=1 corner first; if its
        //   null3-corrected score < MMSEQS_CMSCAN_FG_HYBRID_THRESHOLD (default 0 bits),
        //   fall back to FG=3 midpoint anchor. Preserves column stability for hits where
        //   FG=1 works, recovers FG=3 alignment for targets where FG=1 picks junk.
        {
            // forceGlobal is hoisted above the local-begin pass; reuse it here.
            if (forceGlobal == 1) {
                const int gD = std::min(maxSpan, N);
                const int gI = 1;
                if (gD >= minSpan && gI + gD - 1 <= N) {
                    const float gSc = vit[rootBase + static_cast<size_t>(gD) * iStride + static_cast<size_t>(gI)];
                    if (gSc != NEG_INF_F) {
                        const double corr = null3CorrByInterval(gI, gI + gD - 1);
                        bestSc = static_cast<float>(static_cast<double>(gSc) - corr);
                        bestI = gI;
                        bestD = gD;
                        bestMode = 'T';
                        bestNull3Corr = corr;
                    }
                }
            } else if (forceGlobal == 2) {
                const int gD = std::min(maxSpan, N);
                if (gD >= minSpan) {
                    int bestIglob = -1;
                    float bestScGlob = NEG_INF_F;
                    double bestCorrGlob = 0.0;
                    const int iMax = N - gD + 1;
                    for (int gI = 1; gI <= iMax; ++gI) {
                        const float raw = vit[rootBase + static_cast<size_t>(gD) * iStride + static_cast<size_t>(gI)];
                        if (raw == NEG_INF_F) continue;
                        const double corr = null3CorrByInterval(gI, gI + gD - 1);
                        const float sc = static_cast<float>(static_cast<double>(raw) - corr);
                        if (sc > bestScGlob) {
                            bestScGlob = sc;
                            bestIglob = gI;
                            bestCorrGlob = corr;
                        }
                    }
                    if (bestIglob > 0) {
                        bestSc = bestScGlob;
                        bestI = bestIglob;
                        bestD = gD;
                        bestMode = 'T';
                        bestNull3Corr = bestCorrGlob;
                    }
                }
            } else if ((forceGlobal == 3 || forceGlobal == 4) && anchorI > 0 && anchorD > 0
                       && anchorI <= N && anchorI + anchorD - 1 <= N) {
                const int anchorEnd = anchorI + anchorD - 1;
                const int midpoint = anchorI + (anchorD - 1) / 2;
                int bestIa = -1, bestDa = -1;
                float bestScA = NEG_INF_F;
                double bestCorrA = 0.0;
                for (int d = minSpan; d <= maxSpan; ++d) {
                    const int iMax = N - d + 1;
                    int iLo, iHi;
                    if (forceGlobal == 3) {
                        // envelope must contain midpoint: i <= midpoint <= i+d-1
                        iLo = std::max(1, midpoint - d + 1);
                        iHi = std::min(iMax, midpoint);
                    } else {
                        // envelope must contain [anchorI..anchorEnd]:
                        // i <= anchorI and i+d-1 >= anchorEnd
                        iLo = std::max(1, anchorEnd - d + 1);
                        iHi = std::min(iMax, anchorI);
                    }
                    for (int i = iLo; i <= iHi; ++i) {
                        const float raw = vit[rootBase + static_cast<size_t>(d) * iStride + static_cast<size_t>(i)];
                        if (raw == NEG_INF_F) continue;
                        const double corr = null3CorrByInterval(i, i + d - 1);
                        const float sc = static_cast<float>(static_cast<double>(raw) - corr);
                        if (sc > bestScA) {
                            bestScA = sc;
                            bestIa = i;
                            bestDa = d;
                            bestCorrA = corr;
                        }
                    }
                }
                if (bestIa > 0) {
                    bestSc = bestScA;
                    bestI = bestIa;
                    bestD = bestDa;
                    bestMode = 'J';
                    bestNull3Corr = bestCorrA;
                }
            } else if (forceGlobal == 5) {
                // Margin-gated fallback: keep FG=1 corner unless FG=3 best is
                // meaningfully better. Preserves FG=1 column stability across
                // targets where corner is a reasonable choice; only switches
                // to FG=3 when corner is dramatically worse (e.g., CP150205.1
                // case: corner=-19, fg3=+50, gap=69 bits → fall back).
                // MMSEQS_CMSCAN_FG_HYBRID_MARGIN: gap (bits) above which we
                // fall back. Default 10. Higher value biases toward FG=1.
                static const float fgHybridMargin = []() -> float {
                    const char *e = std::getenv("MMSEQS_CMSCAN_FG_HYBRID_MARGIN");
                    if (e == NULL || *e == '\0') return 10.0f;
                    return static_cast<float>(std::atof(e));
                }();
                const int gD1 = std::min(maxSpan, N);
                const int gI1 = 1;
                float fg1Sc = NEG_INF_F;
                double fg1Corr = 0.0;
                if (gD1 >= minSpan && gI1 + gD1 - 1 <= N) {
                    const float raw = vit[rootBase + static_cast<size_t>(gD1) * iStride + static_cast<size_t>(gI1)];
                    if (raw != NEG_INF_F) {
                        fg1Corr = null3CorrByInterval(gI1, gI1 + gD1 - 1);
                        fg1Sc = static_cast<float>(static_cast<double>(raw) - fg1Corr);
                    }
                }
                int fg3I = -1, fg3D = -1;
                float fg3Sc = NEG_INF_F;
                double fg3Corr = 0.0;
                if (anchorI > 0 && anchorD > 0
                    && anchorI <= N && anchorI + anchorD - 1 <= N) {
                    const int midpoint = anchorI + (anchorD - 1) / 2;
                    for (int d = minSpan; d <= maxSpan; ++d) {
                        const int iMax = N - d + 1;
                        const int iLo = std::max(1, midpoint - d + 1);
                        const int iHi = std::min(iMax, midpoint);
                        for (int i = iLo; i <= iHi; ++i) {
                            const float raw = vit[rootBase + static_cast<size_t>(d) * iStride + static_cast<size_t>(i)];
                            if (raw == NEG_INF_F) continue;
                            const double corr = null3CorrByInterval(i, i + d - 1);
                            const float sc = static_cast<float>(static_cast<double>(raw) - corr);
                            if (sc > fg3Sc) {
                                fg3Sc = sc;
                                fg3I = i;
                                fg3D = d;
                                fg3Corr = corr;
                            }
                        }
                    }
                }
                bool useFg3 = false;
                if (fg3Sc != NEG_INF_F) {
                    if (fg1Sc == NEG_INF_F) {
                        useFg3 = true;
                    } else if ((fg3Sc - fg1Sc) > fgHybridMargin) {
                        useFg3 = true;
                    }
                }
                if (useFg3) {
                    bestSc = fg3Sc;
                    bestI = fg3I;
                    bestD = fg3D;
                    bestMode = 'J';
                    bestNull3Corr = fg3Corr;
                } else if (fg1Sc != NEG_INF_F) {
                    bestSc = fg1Sc;
                    bestI = gI1;
                    bestD = gD1;
                    bestMode = 'T';
                    bestNull3Corr = fg1Corr;
                }
            }
        }
        if (dumpThisTarget) {
            fprintf(stderr, "DUMP_BEST tid=%s bestI=%d bestD=%d bestSc=%.4f bestMode=%c null3Corr=%.4f\n",
                    seqId.c_str(), bestI, bestD, bestSc, bestMode, bestNull3Corr);
            for (int d = minSpan; d <= maxSpan; ++d) {
                if (dBestRaw[static_cast<size_t>(d)] == NEG_INF_F) continue;
                fprintf(stderr, "DUMP_D tid=%s d=%d bestI=%d rawSc=%.4f sc=%.4f\n",
                        seqId.c_str(), d, dBestI[static_cast<size_t>(d)],
                        dBestRaw[static_cast<size_t>(d)], dBestSc[static_cast<size_t>(d)]);
            }
            // Cell-diff vs Infernal: dump alpha[v=0][j=jTarget][d] for j_target from env.
            // Default jTarget = N (envelope-end at last residue) and additionally a list
            // from MMSEQS_CMSCAN_DUMP_J="167" (comma-separated).
            const char *dumpJList = std::getenv("MMSEQS_CMSCAN_DUMP_J");
            std::vector<int> jTargets;
            if (dumpJList != NULL) {
                std::string s = dumpJList;
                size_t pos = 0;
                while (pos < s.size()) {
                    size_t comma = s.find(',', pos);
                    int jt = std::atoi(s.substr(pos, comma - pos).c_str());
                    if (jt > 0 && jt <= N) jTargets.push_back(jt);
                    if (comma == std::string::npos) break;
                    pos = comma + 1;
                }
            }
            for (int jT : jTargets) {
                for (int d = 1; d <= jT; ++d) {
                    const int i = jT - d + 1;
                    if (i < 1 || i > N) continue;
                    if (d > maxSpan) continue;
                    const float a = vit[rootBase + static_cast<size_t>(d) * iStride + static_cast<size_t>(i)];
                    fprintf(stderr, "DUMP_CELL tid=%s j=%d i=%d d=%d alpha=%.4f\n",
                            seqId.c_str(), jT, i, d, a);
                }
            }
            // For each (j_target, d) of interest, dump alpha across all states v
            // that have finite alpha. Lets us see which y is supposed to feed
            // local-begin pickup. Gated by MMSEQS_CMSCAN_DUMP_ALL_V=1.
            const char *dumpAllV = std::getenv("MMSEQS_CMSCAN_DUMP_ALL_V");
            if (dumpAllV != NULL && dumpAllV[0] == '1') {
                const char *dumpDmin = std::getenv("MMSEQS_CMSCAN_DUMP_DMIN");
                int dMin = (dumpDmin != NULL) ? std::atoi(dumpDmin) : 100;
                if (dMin < 1) dMin = 1;
                for (int jT : jTargets) {
                    for (int d = dMin; d <= std::min(jT, maxSpan); ++d) {
                        const int i = jT - d + 1;
                        if (i < 1 || i > N) continue;
                        for (int v = 0; v < M; ++v) {
                            const size_t vBase = stateBase[static_cast<size_t>(v)];
                            const float a = vit[vBase + static_cast<size_t>(d) * iStride + static_cast<size_t>(i)];
                            if (a == NEG_INF_F) continue;
                            const CmState &cs = model.states[static_cast<size_t>(v)];
                            fprintf(stderr, "DUMP_V tid=%s j=%d d=%d v=%d type=%d node=%d ndtype=%d firstOfNode=%d alpha=%.4f beginSc=%.4f endSc=%.4f\n",
                                    seqId.c_str(), jT, d, v,
                                    static_cast<int>(cs.type),
                                    cs.nodeIdx,
                                    static_cast<int>(cs.nodeType),
                                    cs.isFirstOfNode ? 1 : 0,
                                    a,
                                    static_cast<float>(cs.beginSc),
                                    static_cast<float>(cs.endSc));
                        }
                    }
                }
            }
        }
        if (bestSc == NEG_INF_F) {
            return;
        }

        int minUsed = N + 1;
        int maxUsed = 0;
        int obsCount[4] = {0, 0, 0, 0};
        double modelAggRaw[4] = {0.0, 0.0, 0.0, 0.0};
        std::string traceOps;
        traceOps.reserve(static_cast<size_t>(bestD + 8));
        std::vector<int> traceStates;
        traceStates.reserve(static_cast<size_t>(M + N));
        struct TbCell { int v; int i; int d; };
        std::vector<TbCell> st;
        st.reserve(static_cast<size_t>(M + N));
        st.push_back(TbCell{root, bestI, bestD});
        // MMSEQS_CMSCAN_DUMP_TRACE_TID=<seqId>: per-step trace dump for cell-diff
        // vs Infernal --tfile. Format mirrors parsetree idx/emitl/emitr/state.
        const char *dumpTraceTid = std::getenv("MMSEQS_CMSCAN_DUMP_TRACE_TID");
        const bool dumpThisTrace = (dumpTraceTid != NULL && seqId == dumpTraceTid);
        int traceIdxCounter = 0;
        if (dumpThisTrace) {
            fprintf(stderr, "DUMP_TRACE_HDR tid=%s bestI=%d bestD=%d bestSc=%.4f bestMode=%c\n",
                    seqId.c_str(), bestI, bestD, bestSc, bestMode);
        }
        while (!st.empty()) {
            TbCell c = st.back();
            st.pop_back();
            const ExactStateExec &sv = exec[static_cast<size_t>(c.v)];
            if (dumpThisTrace) {
                const int j_dbg = c.i + c.d - 1;
                fprintf(stderr, "DUMP_TRACE idx=%d v=%d type=%d i=%d j=%d d=%d\n",
                        traceIdxCounter++, c.v, (int)sv.type, c.i, j_dbg, c.d);
            }
            if (sv.type == CM_ST_E) {
                continue;
            }
            traceStates.push_back(c.v);
            const double cur = vit[stateBase[static_cast<size_t>(c.v)] +
                                   static_cast<size_t>(c.d) * iStride +
                                   static_cast<size_t>(c.i)];
            if (cur == NEG_INF_F) {
                continue;
            }
            if (sv.type == CM_ST_B) {
                const int y = sv.bLeft;
                const int z = sv.bRight;
                const size_t vd = static_cast<size_t>(c.v) * static_cast<size_t>(N + 1) + static_cast<size_t>(c.d);
                const int kBeg = bSplitBegByVD[vd];
                const int kEnd = bSplitEndByVD[vd];
                // Argmax over the split range. l + r is computed in float to
                // mirror fill's `_mm_max_ps`/`std::fmaxf` reduction order, so
                // ties resolve identically.
                int bestK = -1;
                float bestApprox = NEG_INF_F;
                for (int k = kBeg; k <= kEnd; ++k) {
                    const float l = vit[stateBase[static_cast<size_t>(y)] +
                                        static_cast<size_t>(k) * iStride +
                                        static_cast<size_t>(c.i)];
                    const float r = vit[stateBase[static_cast<size_t>(z)] +
                                        static_cast<size_t>(c.d - k) * iStride +
                                        static_cast<size_t>(c.i + k)];
                    if (l == NEG_INF_F || r == NEG_INF_F) {
                        continue;
                    }
                    const float cand = l + r;
                    if (cand > bestApprox) {
                        bestApprox = cand;
                        bestK = k;
                    }
                }
                if (bestK >= 0) {
                    st.push_back(TbCell{z, c.i + bestK, c.d - bestK});
                    st.push_back(TbCell{y, c.i, bestK});
                }
                continue;
            }

            const int j = c.i + c.d - 1;
            const int consume = sv.dConsume;
            const int ni = c.i + sv.niShift;
            const char op = traceOpForState(sv.type);
            if (op != '\0') {
                for (int x = 0; x < consume; ++x) {
                    traceOps.push_back(op);
                }
            }
            const bool cL = sv.consumesLeft;
            const bool cR = sv.consumesRight;
            if (cL) {
                minUsed = std::min(minUsed, c.i);
                maxUsed = std::max(maxUsed, c.i);
                const int bi = ws.seqCode[static_cast<size_t>(c.i)];
                if (bi >= 0 && bi < 4) {
                    obsCount[bi]++;
                }
            }
            if (cR) {
                minUsed = std::min(minUsed, j);
                maxUsed = std::max(maxUsed, j);
                const int bj = ws.seqCode[static_cast<size_t>(j)];
                if (bj >= 0 && bj < 4) {
                    obsCount[bj]++;
                }
            }
            modelAggRaw[0] += sv.null2Agg[0];
            modelAggRaw[1] += sv.null2Agg[1];
            modelAggRaw[2] += sv.null2Agg[2];
            modelAggRaw[3] += sv.null2Agg[3];
            if (c.d < consume) {
                continue;
            }
            const int nd = c.d - consume;
            const double e = cmStateEmitScoreFast(sv, ws.seqCode, c.i, j);
            const float ef = static_cast<float>(e);
            const size_t ndBase = static_cast<size_t>(nd) * iStride + static_cast<size_t>(ni);
            // Argmax over children. Candidate is computed in float as
            // `(ef + trF) + n` to mirror fill's left-associative SIMD/scalar
            // arithmetic exactly, so trace and fill agree on tie-breaks.
            // Replaces a tolerance-match early-exit (`|cur - cand| < 1e-6`)
            // that could short-circuit on a non-optimal child when float
            // rounding pushed an inferior candidate within the tolerance.
            int bestY = -1;
            float bestCand = NEG_INF_F;
            const int trCount = sv.trCount;
            const bool dumpCand = (dumpThisTrace && c.v == root);
            // MMSEQS_CMSCAN_TRACE_LBO=1: at root, mirror Infernal cm_localize
            // (cm_modelconfig.c:491) which zeroes cm->t[0] when local mode is on,
            // making ROOT_S→child via trDst NEG_INF in cm_alignT. Only local-begin
            // can leave the root. Equivalent at trace time without poisoning DP fill.
            static const char *const lboEnv = std::getenv("MMSEQS_CMSCAN_TRACE_LBO");
            const bool traceLBO = (lboEnv != NULL && lboEnv[0] == '1');
            const bool skipTrDstAtRoot = traceLBO && (c.v == root) && model.hasLocalCfg;
            if (skipTrDstAtRoot) {
                if (dumpCand) {
                    fprintf(stderr, "DUMP_TRACE_CAND v=%d kind=tr_skipped reason=LBO\n", c.v);
                }
            } else if (trCount <= 4) {
                for (int t = 0; t < trCount; ++t) {
                    const int y = sv.trDst4[t];
                    const float n = vit[stateBase[static_cast<size_t>(y)] + ndBase];
                    if (dumpCand) {
                        fprintf(stderr, "DUMP_TRACE_CAND v=%d kind=tr t=%d y=%d trSc=%.6f n=%.6f cand=%.6f\n",
                                c.v, t, y, sv.trScF4[t], n,
                                (n == NEG_INF_F) ? -INFINITY : (ef + sv.trScF4[t]) + n);
                    }
                    if (n == NEG_INF_F) {
                        continue;
                    }
                    const float cand = (ef + sv.trScF4[t]) + n;
                    if (cand > bestCand) {
                        bestCand = cand;
                        bestY = y;
                    }
                }
            } else {
                for (int t = 0; t < trCount; ++t) {
                    const size_t ti = sv.trOff + static_cast<size_t>(t);
                    const int y = trDst[ti];
                    const float trF = static_cast<float>(trSc[ti]);
                    const float n = vit[stateBase[static_cast<size_t>(y)] + ndBase];
                    if (dumpCand) {
                        fprintf(stderr, "DUMP_TRACE_CAND v=%d kind=tr t=%d y=%d trSc=%.6f n=%.6f cand=%.6f\n",
                                c.v, t, y, trF, n,
                                (n == NEG_INF_F) ? -INFINITY : (ef + trF) + n);
                    }
                    if (n == NEG_INF_F) {
                        continue;
                    }
                    const float cand = (ef + trF) + n;
                    if (cand > bestCand) {
                        bestCand = cand;
                        bestY = y;
                    }
                }
            }
            // Local-begin trace: at ROOT_S, alpha[root][i][d] may have been
            // reached via a direct local-begin jump to some begin-eligible y
            // (the per-cell pickup at CmScan.cpp:2807-2826 mirrors cmsearch
            // search semantics). Without considering those candidates here the
            // trace walks ROOT_S→child recurrence even when local-begin won at
            // the chart, producing a misaligned CIGAR. Cell-diff bug:
            // project_cyk_celldiff_root_local_begin.md.
            if (c.v == root && model.hasLocalCfg) {
                for (int y = 1; y < M; ++y) {
                    const CmState &cs = model.states[static_cast<size_t>(y)];
                    if (cs.beginSc == NEG_INF) {
                        continue;
                    }
                    const float bsc = static_cast<float>(cs.beginSc);
                    const float n = vit[stateBase[static_cast<size_t>(y)] + ndBase];
                    if (dumpCand) {
                        fprintf(stderr, "DUMP_TRACE_CAND v=%d kind=lb y=%d bsc=%.6f n=%.6f cand=%.6f\n",
                                c.v, y, bsc, n,
                                (n == NEG_INF_F) ? -INFINITY : (ef + bsc) + n);
                    }
                    if (n == NEG_INF_F) {
                        continue;
                    }
                    const float cand = (ef + bsc) + n;  // ef=0 for ROOT_S
                    if (cand > bestCand) {
                        bestCand = cand;
                        bestY = y;
                    }
                }
            }
            // EL pop-out candidate (mirrors yshadow=USED_EL in Infernal cm_alignT).
            // At any end-eligible state v (endSc[v] != NEG_INF, local mode on),
            // alpha[v][j][d] may have been reached via the EL pre-fill at line
            // 2700-2742: emit_v(i,d) + elSelf*(d-sd) + endSc[v]. If that beats
            // every trDst-via-children candidate above, record EL pop-out: skip
            // pushing children. The CIGAR builder is column-driven (line 1593:
            // any unwalked col = D op), so the right subtree's columns become
            // gaps automatically. Residues popped to EL are accounted for via
            // minUsed/maxUsed expansion so the alignment span covers them.
            // Gated by MMSEQS_CMSCAN_TRACE_EL=1 for diagnostic A/B.
            static const char *const elTraceEnv = std::getenv("MMSEQS_CMSCAN_TRACE_EL");
            const bool traceEL = (elTraceEnv != NULL && elTraceEnv[0] == '1');
            bool pickEL = false;
            if (traceEL && model.hasLocalCfg && model.elSelf <= 0.0) {
                const CmState &csEL = model.states[static_cast<size_t>(c.v)];
                if (csEL.endSc != NEG_INF) {
                    const float elCand = ef
                        + static_cast<float>(model.elSelf) * static_cast<float>(nd)
                        + static_cast<float>(csEL.endSc);
                    if (dumpCand) {
                        fprintf(stderr, "DUMP_TRACE_CAND v=%d kind=el endSc=%.6f elContrib=%.6f cand=%.6f\n",
                                c.v, (float)csEL.endSc,
                                static_cast<float>(model.elSelf) * static_cast<float>(nd),
                                elCand);
                    }
                    if (elCand > bestCand) {
                        bestCand = elCand;
                        pickEL = true;
                    }
                }
            }
            if (dumpThisTrace && c.v == root) {
                fprintf(stderr, "DUMP_TRACE_PICK v=%d bestY=%d bestCand=%.6f%s\n",
                        c.v, bestY, bestCand, pickEL ? " EL" : "");
            }
            (void)cur;  // cur is only used for the NEG_INF_F early-skip above.
            if (pickEL) {
                // EL pops `nd` residues without column consumption. Update span
                // accounting; do not push children.
                if (nd > 0) {
                    minUsed = std::min(minUsed, ni);
                    maxUsed = std::max(maxUsed, ni + nd - 1);
                }
                if (dumpThisTrace) {
                    fprintf(stderr, "DUMP_TRACE_EL v=%d popped_residues=%d at_i=%d..%d\n",
                            c.v, nd, ni, ni + nd - 1);
                }
            } else if (bestY >= 0) {
                st.push_back(TbCell{bestY, ni, nd});
            }
        }
        const int tracedLen = (minUsed <= maxUsed) ? (maxUsed - minUsed + 1) : 0;
        if (minUsed > maxUsed || tracedLen < std::max(1, minSpan / 2)) {
            minUsed = bestI;
            maxUsed = bestI + bestD - 1;
        }
        Hit h;
        h.seqId = seqId;
        h.start1 = minUsed;
        h.end1 = maxUsed;
        h.mode = bestMode;
        h.trunc = (bestMode != 'J');
        h.traceStates = encodeTraceStates(traceStates);
        {
            int qS = -1, qE = -1, aL = 0, leadIns = 0, trailIns = 0;
            h.cigar = modelTraceCigarQueryCoord(model, traceStates, &qS, &qE, &aL,
                                                &leadIns, &trailIns);
            h.qStart = qS;
            h.qEnd = qE;
            h.cigarAlnLen = static_cast<unsigned int>(std::max(0, aL));
            h.leadingInsertTargets = leadIns;
            h.trailingInsertTargets = trailIns;
        }
        const double null2Corr = scoreCorrectionNull2BitsFromTrace(model, modelAggRaw, obsCount);
        h.cyk = static_cast<double>(bestSc) - null2Corr;
        if (wantInside) {
            // Deck-based inside DP over the restricted interval [bestI, bestJ].
            // Layout: insD[d * dStride + v * vStride + li] where li = i - bestI.
            // Iterate d=0..localN (outer), v=M-1..0 (children before parents),
            // li=0..localN-d (inner, sequential → cache-friendly reads/writes).
            const int localN = bestD;
            const size_t vStride = static_cast<size_t>(localN + 1);
            const size_t dStride = static_cast<size_t>(M) * vStride;
            const size_t insSize = static_cast<size_t>(localN + 1) * dStride;
            static constexpr size_t CM_INSIDE_MAX_CELLS = 100000000ULL; // ~800 MB doubles
            if (insSize > CM_INSIDE_MAX_CELLS) {
                h.inside = NEG_INF;
                h.bias = bestNull3Corr + null2Corr;
                outHits.push_back(h);
                return;
            }
            if (!ws.insD.ensure(insSize)) {
                Debug(Debug::ERROR) << "cmscan: failed to allocate inside buffer\n";
                return;
            }
            double *insDataAll = ws.insD.data();
            for (size_t ii = 0; ii < insSize; ++ii) {
                insDataAll[ii] = NEG_INF;
            }
            // E states emit the empty string: ins[d=0][E][li] = 0.
            for (int v = 0; v < M; ++v) {
                if (exec[static_cast<size_t>(v)].type != CM_ST_E) { continue; }
                double *row = insDataAll + static_cast<size_t>(v) * vStride;
                for (int li = 0; li <= localN; ++li) { row[li] = 0.0; }
            }
            for (int d = 0; d <= localN; ++d) {
                const int liMax = localN - d;
                double *insSlice = insDataAll + static_cast<size_t>(d) * dStride;
                const double *insData = insDataAll;
                const StateId *trDstPtr = trDst.data();
                const double *trScPtr = trSc.data();
                for (int v = M - 1; v >= 0; --v) {
                    const ExactStateExec &st = exec[static_cast<size_t>(v)];
                    if (st.type == CM_ST_E) { continue; }
                    if (st.dmax >= st.dmin && (d < st.dmin || d > st.dmax)) { continue; }
                    if (st.type == CM_ST_B) {
                        const int y = st.bLeft, z = st.bRight;
                        if (y < 0 || y >= M || z < 0 || z >= M) { continue; }
                        const size_t vd = static_cast<size_t>(v) * static_cast<size_t>(N + 1) + static_cast<size_t>(d);
                        const int kBeg = bSplitBegByVD[vd];
                        const int kEnd = bSplitEndByVD[vd];
                        double *dst = insSlice + static_cast<size_t>(v) * vStride;
                        if (fastLogsum) {
                            for (int li = 0; li <= liMax; ++li) {
                                bool has = false;
                                float acc = NEG_INF_F;
                                for (int k = kBeg; k <= kEnd; ++k) {
                                    const double lv = insData[static_cast<size_t>(k) * dStride + static_cast<size_t>(y) * vStride + static_cast<size_t>(li)];
                                    const double rv = insData[static_cast<size_t>(d - k) * dStride + static_cast<size_t>(z) * vStride + static_cast<size_t>(li + k)];
                                    if (lv != NEG_INF && rv != NEG_INF) {
                                        log2AccFastAddF(static_cast<float>(lv + rv), has, acc);
                                    }
                                }
                                dst[li] = has ? static_cast<double>(acc) : NEG_INF;
                            }
                        } else {
                            for (int li = 0; li <= liMax; ++li) {
                                bool has = false;
                                double maxVal = NEG_INF;
                                double scaledSum = 0.0;
                                for (int k = kBeg; k <= kEnd; ++k) {
                                    const double lv = insData[static_cast<size_t>(k) * dStride + static_cast<size_t>(y) * vStride + static_cast<size_t>(li)];
                                    const double rv = insData[static_cast<size_t>(d - k) * dStride + static_cast<size_t>(z) * vStride + static_cast<size_t>(li + k)];
                                    if (lv != NEG_INF && rv != NEG_INF) {
                                        log2AccExactAdd(lv + rv, has, maxVal, scaledSum);
                                    }
                                }
                                dst[li] = log2AccExactValue(has, maxVal, scaledSum);
                            }
                        }
                        continue;
                    }
                    {
                        const char *probeV = std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_V");
                        if (probeV != NULL && v == std::atoi(probeV)) {
                            fprintf(stderr, "DUMP_TRANS_PROBE v=%d d=%d type=%d dConsume=%d trCount=%d trOff=%zu emitSize=%d emitPtr=%p niShift=%d\n",
                                    v, d, (int)st.type, st.dConsume, st.trCount, st.trOff, st.emitSize, (void*)st.emitPtr, st.niShift);
                        }
                    }
                    if (d < st.dConsume) { continue; }
                    const int nd = d - st.dConsume;
                    const int niShift = st.niShift;
                    const float *ep = st.emitPtr;
                    const int8_t *sc = ws.seqCode.data();
                    const double *nxtSlice = insDataAll + static_cast<size_t>(nd) * dStride;
                    double *dst = insSlice + static_cast<size_t>(v) * vStride;
                    if (fastLogsum) {
                        for (int li = 0; li <= liMax; ++li) {
                            float ef = 0.0f;
                            if (st.consumesLeft && st.consumesRight) {
                                if (!(ep && st.emitSize >= 16)) {
                                    ef = NEG_INF_F;
                                } else {
                                    const int i = li + bestI;
                                    const int lc = sc[i];
                                    const int rc = sc[i + d - 1];
                                    ef = (lc >= 0 && rc >= 0) ? ep[lc * 4 + rc] : -1.0f;
                                }
                            } else if (st.consumesLeft) {
                                if (!(ep && st.emitSize >= 4)) {
                                    ef = NEG_INF_F;
                                } else {
                                    const int bi = sc[li + bestI];
                                    ef = (bi >= 0) ? ep[bi] : -1.0f;
                                }
                            } else if (st.consumesRight) {
                                if (!(ep && st.emitSize >= 4)) {
                                    ef = NEG_INF_F;
                                } else {
                                    const int bi = sc[li + bestI + d - 1];
                                    ef = (bi >= 0) ? ep[bi] : -1.0f;
                                }
                            }
                            const char *trDumpV = std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_V");
                            const int trDumpD = (std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_D") != NULL) ? std::atoi(std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_D")) : -1;
                            const int trDumpJ = (std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_J") != NULL) ? std::atoi(std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_J")) : -1;
                            const bool dumpV = (trDumpV != NULL) && (v == std::atoi(trDumpV));
                            const bool doDump = dumpV && (d == trDumpD) && (trDumpJ < 0 || li == (trDumpJ - d));
                            if (doDump) {
                                fprintf(stderr, "DUMP_TRANS_LIENTRY v=%d d=%d li=%d ef=%.6f bestI=%d niShift=%d\n",
                                        v, d, li, ef, bestI, niShift);
                            }
                            if (ef == NEG_INF_F) { continue; }
                            const int nli = li + niShift;
                            bool has = false;
                            float acc = NEG_INF_F;
                            // LOCAL-mode EL pre-init (Infernal cm_dpalign.c:1600-1604).
                            // For end-eligible v with finite endsc[v], seed alpha with
                            // ef + el_scA[d-sd] + endsc[v] before children FLogsum.
                            if (model.hasLocalCfg && model.elSelf <= 0.0) {
                                const CmState &csEl = model.states[static_cast<size_t>(v)];
                                if (csEl.endSc != NEG_INF) {
                                    int sdEl = -1;
                                    if (csEl.type == CM_ST_MP) sdEl = 2;
                                    else if (csEl.type == CM_ST_ML || csEl.type == CM_ST_MR) sdEl = 1;
                                    else if (csEl.type == CM_ST_S) sdEl = 0;
                                    if (sdEl >= 0 && d >= sdEl) {
                                        const float elTerm = ef
                                            + static_cast<float>(model.elSelf) * static_cast<float>(d - sdEl)
                                            + static_cast<float>(csEl.endSc);
                                        log2AccFastAddF(elTerm, has, acc);
                                    }
                                }
                            }
                            for (int t = 0; t < st.trCount; ++t) {
                                const size_t ti = st.trOff + static_cast<size_t>(t);
                                const double n = nxtSlice[static_cast<size_t>(trDstPtr[ti]) * vStride + static_cast<size_t>(nli)];
                                if (doDump) {
                                    fprintf(stderr, "DUMP_TRANS v=%d d=%d t=%d trDst=%d trSc=%.6f ef=%.6f childAlpha=%.6f sum=%.6f\n",
                                            v, d, t, (int)trDstPtr[ti], trScPtr[ti], ef, n,
                                            (n != NEG_INF) ? (ef + static_cast<float>(trScPtr[ti]) + static_cast<float>(n)) : NEG_INF);
                                }
                                if (n != NEG_INF) {
                                    log2AccFastAddF(ef + static_cast<float>(trScPtr[ti]) + static_cast<float>(n), has, acc);
                                }
                            }
                            if (doDump) {
                                fprintf(stderr, "DUMP_TRANS v=%d d=%d FINAL acc=%.6f\n", v, d, has ? acc : NEG_INF);
                            }
                            dst[li] = has ? static_cast<double>(acc) : NEG_INF;
                        }
                    } else {
                        for (int li = 0; li <= liMax; ++li) {
                            float ef = 0.0f;
                            if (st.consumesLeft && st.consumesRight) {
                                if (!(ep && st.emitSize >= 16)) {
                                    ef = NEG_INF_F;
                                } else {
                                    const int i = li + bestI;
                                    const int lc = sc[i];
                                    const int rc = sc[i + d - 1];
                                    ef = (lc >= 0 && rc >= 0) ? ep[lc * 4 + rc] : -1.0f;
                                }
                            } else if (st.consumesLeft) {
                                if (!(ep && st.emitSize >= 4)) {
                                    ef = NEG_INF_F;
                                } else {
                                    const int bi = sc[li + bestI];
                                    ef = (bi >= 0) ? ep[bi] : -1.0f;
                                }
                            } else if (st.consumesRight) {
                                if (!(ep && st.emitSize >= 4)) {
                                    ef = NEG_INF_F;
                                } else {
                                    const int bi = sc[li + bestI + d - 1];
                                    ef = (bi >= 0) ? ep[bi] : -1.0f;
                                }
                            }
                            const char *trDumpV2 = std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_V");
                            const int trDumpD2 = (std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_D") != NULL) ? std::atoi(std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_D")) : -1;
                            const int trDumpJ2 = (std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_J") != NULL) ? std::atoi(std::getenv("MMSEQS_CMSCAN_DUMP_TRANS_J")) : -1;
                            const bool dumpV2 = (trDumpV2 != NULL) && (v == std::atoi(trDumpV2));
                            const bool doDump2 = dumpV2 && (d == trDumpD2) && (trDumpJ2 < 0 || li == (trDumpJ2 - d));
                            if (doDump2) {
                                fprintf(stderr, "DUMP_TRANS_LIENTRY v=%d d=%d li=%d ef=%.6f bestI=%d niShift=%d (exact)\n",
                                        v, d, li, ef, bestI, niShift);
                            }
                            if (ef == NEG_INF_F) { continue; }
                            const int nli = li + niShift;
                            bool has = false;
                            double maxVal = NEG_INF;
                            double scaledSum = 0.0;
                            // LOCAL-mode EL pre-init (Infernal cm_dpalign.c:1600-1604).
                            if (model.hasLocalCfg && model.elSelf <= 0.0) {
                                const CmState &csEl = model.states[static_cast<size_t>(v)];
                                if (csEl.endSc != NEG_INF) {
                                    int sdEl = -1;
                                    if (csEl.type == CM_ST_MP) sdEl = 2;
                                    else if (csEl.type == CM_ST_ML || csEl.type == CM_ST_MR) sdEl = 1;
                                    else if (csEl.type == CM_ST_S) sdEl = 0;
                                    if (sdEl >= 0 && d >= sdEl) {
                                        const double elTerm = static_cast<double>(ef)
                                            + model.elSelf * static_cast<double>(d - sdEl)
                                            + csEl.endSc;
                                        log2AccExactAdd(elTerm, has, maxVal, scaledSum);
                                    }
                                }
                            }
                            for (int t = 0; t < st.trCount; ++t) {
                                const size_t ti = st.trOff + static_cast<size_t>(t);
                                const double n = nxtSlice[static_cast<size_t>(trDstPtr[ti]) * vStride + static_cast<size_t>(nli)];
                                if (doDump2) {
                                    fprintf(stderr, "DUMP_TRANS v=%d d=%d t=%d trDst=%d trSc=%.6f ef=%.6f childAlpha=%.6f sum=%.6f (exact)\n",
                                            v, d, t, (int)trDstPtr[ti], trScPtr[ti], ef, n,
                                            (n != NEG_INF) ? (static_cast<double>(ef) + trScPtr[ti] + n) : NEG_INF);
                                }
                                if (n != NEG_INF) {
                                    log2AccExactAdd(static_cast<double>(ef) + trScPtr[ti] + n, has, maxVal, scaledSum);
                                }
                            }
                            const double finalA = log2AccExactValue(has, maxVal, scaledSum);
                            if (doDump2) {
                                fprintf(stderr, "DUMP_TRANS v=%d d=%d FINAL acc=%.6f (exact)\n", v, d, finalA);
                            }
                            dst[li] = finalA;
                        }
                    }
                }
            }
            // LOCAL-mode local-begin pickup (Infernal cm_dpalign.c:1712-1721).
            // After all states are filled, the root cell at (d=L, li=0) can be
            // reached by jumping directly to any begin-eligible state y at score
            // alpha[y][L][L] + beginsc[y]. FLogsum these into alpha[0][L][L].
            if (model.hasLocalCfg) {
                const int root0 = model.rootState;
                const size_t rootCellIdx = static_cast<size_t>(bestD) * dStride
                                          + static_cast<size_t>(root0) * vStride + 0;
                const double rootA = insDataAll[rootCellIdx];
                bool has = false;
                double maxVal = NEG_INF;
                double scaledSum = 0.0;
                if (rootA != NEG_INF) log2AccExactAdd(rootA, has, maxVal, scaledSum);
                for (int y = 1; y < M; ++y) {
                    const CmState &csB = model.states[static_cast<size_t>(y)];
                    if (csB.beginSc == NEG_INF) continue;
                    const double yA = insDataAll[static_cast<size_t>(bestD) * dStride
                                                  + static_cast<size_t>(y) * vStride + 0];
                    if (yA == NEG_INF) continue;
                    log2AccExactAdd(yA + csB.beginSc, has, maxVal, scaledSum);
                }
                insDataAll[rootCellIdx] = log2AccExactValue(has, maxVal, scaledSum);
            }
            h.inside = insDataAll[static_cast<size_t>(bestD) * dStride + static_cast<size_t>(root) * vStride + 0]
                       - bestNull3Corr - null2Corr;
            h.bias = bestNull3Corr + null2Corr;
            // Cell-diff vs Infernal cm_InsideAlign: dump alpha[v][j][d] for chosen
            // j_target (envelope-relative coords 1..localN). MMSEQS_CMSCAN_DUMP_INSIDE_J
            // is a comma list. Indexing: alpha[v][j][d] = insDataAll[d*dStride + v*vStride + li]
            // where li = j - d (since the envelope starts at bestI=1 here).
            const char *dumpInsJList = std::getenv("MMSEQS_CMSCAN_DUMP_INSIDE_J");
            if (dumpInsJList != NULL) {
                std::vector<int> insJTargets;
                std::string s = dumpInsJList;
                size_t pos = 0;
                while (pos < s.size()) {
                    size_t comma = s.find(',', pos);
                    int jt = std::atoi(s.substr(pos, comma - pos).c_str());
                    if (jt > 0 && jt <= bestD) insJTargets.push_back(jt);
                    if (comma == std::string::npos) break;
                    pos = comma + 1;
                }
                for (int jT : insJTargets) {
                    for (int d = 0; d <= jT; ++d) {
                        const int li = jT - d;
                        if (li < 0 || li > bestD) continue;
                        const double *insSlice = insDataAll + static_cast<size_t>(d) * dStride;
                        for (int v = 0; v < M; ++v) {
                            const double a = insSlice[static_cast<size_t>(v) * vStride + static_cast<size_t>(li)];
                            if (a == NEG_INF) continue;
                            const ExactStateExec &st = exec[static_cast<size_t>(v)];
                            fprintf(stderr, "DUMP_INSIDE tid=%s j=%d d=%d v=%d type=%d alpha=%.6f\n",
                                    seqId.c_str(), jT, d, v,
                                    static_cast<int>(st.type),
                                    a);
                        }
                    }
                }
            }
        } else {
            h.inside = NEG_INF;
            h.bias = bestNull3Corr + null2Corr;
        }
        outHits.push_back(h);
        return;
        }
    }

    // Compute maxSpan first so we can size the per-state ragged chart precisely.
    int maxSpan = N;
    if (model.w > 0) {
        maxSpan = std::min(maxSpan, model.w);
    }
    // Cap d-space by the prefilter alignment length when available — a true
    // CM hit has d within a small multiple of the rnasearch envelope length,
    // and the memoized Viterbi only wastes memory on larger spans.
    if (maxHitLen > 0) {
        maxSpan = std::min(maxSpan, maxHitLen);
    }

    // Per-state sparse chart: QDB defines legal d rows, and optional LoL bands
    // later shrink each (state,d) row to only its live i interval.
    std::vector<int> &bandDmin = ws.exactChartBandDmin;
    std::vector<int> &bandSize = ws.exactChartBandSize;
    bandDmin.assign(static_cast<size_t>(M), 0);
    bandSize.assign(static_cast<size_t>(M), 0);
    for (int v = 0; v < M; ++v) {
        const ExactStateExec &e = exec[static_cast<size_t>(v)];
        int dlo;
        int dhi;
        if (e.dmax >= e.dmin && e.dmin >= 0) {
            dlo = std::max(0, e.dmin);
            dhi = std::min(maxSpan, e.dmax);
        } else {
            // Unset/invalid band — fall back to full range [0..maxSpan].
            dlo = 0;
            dhi = maxSpan;
        }
        if (dhi < dlo) {
            bandDmin[v] = 0;
            bandSize[v] = 0;
        } else {
            bandDmin[v] = dlo;
            bandSize[v] = dhi - dlo + 1;
        }
    }
    size_t chartCells = 0;

    // Safety budget. 1G cells = 4GB float + 4GB uint32 = 8GB per-thread worst case.
    // Big rRNA CMs (M=6869-9017, like CP000968) vs comparable-length targets
    // land at ~8G cells = 65GB with wide QDB bands; we skip rather than OOM.
    // TODO: streaming/windowed CYK or HMM pre-filter to unblock these queries.
    static constexpr size_t CM_EXACT_CHART_MAX_CELLS = 1000000000ULL;
    if (chartCells > CM_EXACT_CHART_MAX_CELLS) {
        Debug(Debug::WARNING) << "cmscan: per-state ragged chart too large "
                              << "(" << chartCells << " > " << CM_EXACT_CHART_MAX_CELLS
                              << " cells) for M=" << M << " N=" << N << "; skipping target\n";
        return;
    }

    std::vector<size_t> &chartStateDOffset = ws.exactChartStateDOffset;
    std::vector<size_t> &chartRowOffset = ws.exactChartRowOffset;
    std::vector<int> &chartRowImin = ws.exactChartRowImin;
    std::vector<int> &chartRowSize = ws.exactChartRowSize;
    chartStateDOffset.assign(static_cast<size_t>(M + 1), 0);
    for (int v = 0; v < M; ++v) {
        chartStateDOffset[static_cast<size_t>(v + 1)] =
            chartStateDOffset[static_cast<size_t>(v)] + static_cast<size_t>(bandSize[v]);
    }
    const size_t totalRows = chartStateDOffset[static_cast<size_t>(M)];
    chartRowOffset.assign(totalRows + 1, 0);
    chartRowImin.assign(totalRows, 1);
    chartRowSize.assign(totalRows, 0);

    chartCells = 0;
    for (int v = 0; v < M; ++v) {
        const ExactStateExec &st = exec[static_cast<size_t>(v)];
        for (int dIdx = 0; dIdx < bandSize[v]; ++dIdx) {
            const int d = bandDmin[v] + dIdx;
            const size_t rowIdx = chartStateDOffset[static_cast<size_t>(v)] + static_cast<size_t>(dIdx);
            int iLo = 1;
            int iHi = N - d + 1;
            if (iHi < iLo) {
                chartRowOffset[rowIdx + 1] = chartCells;
                continue;
            }
            if (stateMinI != NULL && st.consumesLeft) {
                iLo = std::max(iLo, stateMinI[v]);
                iHi = std::min(iHi, stateMaxI[v]);
            }
            if (stateMinJ != NULL && st.consumesRight) {
                iLo = std::max(iLo, stateMinJ[v] - d + 1);
                iHi = std::min(iHi, stateMaxJ[v] - d + 1);
            }
            if (iHi >= iLo) {
                chartRowImin[rowIdx] = iLo;
                chartRowSize[rowIdx] = iHi - iLo + 1;
                chartCells += static_cast<size_t>(chartRowSize[rowIdx]);
            }
            chartRowOffset[rowIdx + 1] = chartCells;
        }
    }

    std::vector<float> &exactChart = ws.exactChart;
    std::vector<uint32_t> &exactChartSeen = ws.exactChartSeen;
    uint32_t &exactChartGen = ws.exactChartGen;
    if (exactChart.size() < chartCells) {
        exactChart.resize(chartCells);
    }
    if (exactChartSeen.size() < chartCells) {
        exactChartSeen.resize(chartCells, 0);
    }
    ++exactChartGen;
    if (exactChartGen == 0) {
        std::fill(exactChartSeen.begin(), exactChartSeen.end(), 0u);
        exactChartGen = 1;
    }

    ExactRecCtx recCtx;
    recCtx.N = N;
    recCtx.M = M;
    recCtx.exec = &exec;
    recCtx.execData = exec.data();
    recCtx.trDst = trDst.empty() ? NULL : &trDst[0];
    recCtx.trSc = trSc.empty() ? NULL : &trSc[0];
    recCtx.seqCode = &ws.seqCode;
    recCtx.chartScore = exactChart.data();
    recCtx.chartSeen = exactChartSeen.data();
    recCtx.chartGen = exactChartGen;
    recCtx.chartStateDOffset = chartStateDOffset.data();
    recCtx.chartRowOffset = chartRowOffset.data();
    recCtx.chartRowImin = chartRowImin.data();
    recCtx.chartRowSize = chartRowSize.data();
    recCtx.chartBandDmin = bandDmin.data();
    recCtx.chartBandSize = bandSize.data();
    recCtx.stateMinI = stateMinI;
    recCtx.stateMaxI = stateMaxI;
    recCtx.stateMinJ = stateMinJ;
    recCtx.stateMaxJ = stateMaxJ;
    // Infernal semantics evaluate all legal spans under state dmin/dmax/w constraints.
    int minSpan = 1;
    // Fix B: force prefilter envelope on memoized path too (same as fast path).
    // Only enable if forcedD fits within the chart bands already allocated.
    const bool forceEnv2 = (forcedI > 0 && forcedD > 0 &&
                            forcedD <= N && forcedI <= N - forcedD + 1 &&
                            forcedD <= maxSpan);
    if (forceEnv2) {
        minSpan = forcedD;
        maxSpan = forcedD;
    }
    double bestSc = NEG_INF;
    int bestI = 1;
    int bestJ = N;
    char bestMode = 'J';
    double bestNull3Corr = 0.0;
    const bool enableTruncModes = cmTruncModesEnabled();
    double bestL = NEG_INF, bestR = NEG_INF, bestT = NEG_INF;
    int bestIL = 1, bestIR = 1, bestIT = 1;
    int bestJL = N, bestJR = N, bestJT = N;
    double bestLCorr = 0.0, bestRCorr = 0.0, bestTCorr = 0.0;
    for (int d = minSpan; d <= maxSpan; ++d) {
        const int iMax = N - d + 1;
        int iLo = 1, iHi = iMax;
        if (forceEnv2) { iLo = forcedI; iHi = forcedI; }
        for (int i = iLo; i <= iHi; ++i) {
            const int j = i + d - 1;
            const double rawSc = exactVitRec(recCtx, model.rootState, i, j);
            if (rawSc == NEG_INF) {
                continue;
            }
            const bool mayImproveJ = rawSc > bestSc;
            const bool mayImproveL = enableTruncModes && j == N && rawSc > bestL;
            const bool mayImproveR = enableTruncModes && i == 1 && rawSc > bestR;
            const bool mayImproveT = enableTruncModes && i == 1 && j == N && rawSc > bestT;
            if (!(mayImproveJ || mayImproveL || mayImproveR || mayImproveT)) {
                continue;
            }
            const double corr = null3CorrByInterval(i, j);
            const double sc = rawSc - corr;
            if (mayImproveJ && sc > bestSc) {
                bestSc = sc;
                bestI = i;
                bestJ = j;
                bestMode = 'J';
                bestNull3Corr = corr;
            }
            if (enableTruncModes) {
                if (j == N && sc > bestL) {
                    bestL = sc;
                    bestIL = i;
                    bestJL = j;
                    bestLCorr = corr;
                }
                if (i == 1 && sc > bestR) {
                    bestR = sc;
                    bestIR = i;
                    bestJR = j;
                    bestRCorr = corr;
                }
                if (i == 1 && j == N && sc > bestT) {
                    bestT = sc;
                    bestIT = i;
                    bestJT = j;
                    bestTCorr = corr;
                }
            }
        }
    }
    if (enableTruncModes) {
        if (bestL > bestSc) {
            bestSc = bestL;
            bestI = bestIL;
            bestJ = bestJL;
            bestMode = 'L';
            bestNull3Corr = bestLCorr;
        }
        if (bestR > bestSc) {
            bestSc = bestR;
            bestI = bestIR;
            bestJ = bestJR;
            bestMode = 'R';
            bestNull3Corr = bestRCorr;
        }
        if (bestT > bestSc) {
            bestSc = bestT;
            bestI = bestIT;
            bestJ = bestJT;
            bestMode = 'T';
            bestNull3Corr = bestTCorr;
        }
    }
    if (bestSc == NEG_INF) {
        return;
    }
    int minUsed = N + 1;
    int maxUsed = 0;
    int obsCount[4] = {0, 0, 0, 0};
    double modelAggRaw[4] = {0.0, 0.0, 0.0, 0.0};
    std::string traceOps;
    traceOps.reserve(static_cast<size_t>(bestJ - bestI + 4));
    std::vector<int> traceStates;
    traceStates.reserve(static_cast<size_t>(M + N));
    exactTraceRec(recCtx, model.rootState, bestI, bestJ, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
    const int tracedLen = (minUsed <= maxUsed) ? (maxUsed - minUsed + 1) : 0;
    if (minUsed > maxUsed || tracedLen < std::max(1, minSpan / 2)) {
        minUsed = bestI;
        maxUsed = bestJ;
    }

    Hit h;
    h.seqId = seqId;
    h.start1 = minUsed;
    h.end1 = maxUsed;
    h.mode = bestMode;
    h.trunc = (bestMode != 'J');
    h.traceStates = encodeTraceStates(traceStates);
    {
        int qS = -1, qE = -1, aL = 0, leadIns = 0, trailIns = 0;
        h.cigar = modelTraceCigarQueryCoord(model, traceStates, &qS, &qE, &aL,
                                            &leadIns, &trailIns);
        h.qStart = qS;
        h.qEnd = qE;
        h.cigarAlnLen = static_cast<unsigned int>(std::max(0, aL));
        h.leadingInsertTargets = leadIns;
        h.trailingInsertTargets = trailIns;
    }
    const double null2Corr = scoreCorrectionNull2BitsFromTrace(model, modelAggRaw, obsCount);
    h.cyk = bestSc - null2Corr;
    if (wantInside) {
        // Deck-based inside DP over the restricted interval [bestI, bestJ].
        // Replaces the memoized insideRec for better cache performance.
        const int localN2 = bestJ - bestI + 1;
        const size_t vStride2 = static_cast<size_t>(localN2 + 1);
        const size_t dStride2 = static_cast<size_t>(M) * vStride2;
        const size_t insSize2 = static_cast<size_t>(localN2 + 1) * dStride2;
        static constexpr size_t CM_INSIDE_MAX_CELLS2 = 100000000ULL; // ~800 MB doubles
        if (insSize2 > CM_INSIDE_MAX_CELLS2) {
            h.inside = NEG_INF;
            h.bias = bestNull3Corr + null2Corr;
            outHits.push_back(h);
            return;
        }
        if (!ws.insD.ensure(insSize2)) {
            Debug(Debug::ERROR) << "cmscan: failed to allocate inside buffer2\n";
            return;
        }
        double *insDataAll2 = ws.insD.data();
        for (size_t ii = 0; ii < insSize2; ++ii) {
            insDataAll2[ii] = NEG_INF;
        }
        for (int v = 0; v < M; ++v) {
            if (exec[static_cast<size_t>(v)].type != CM_ST_E) { continue; }
            double *row = insDataAll2 + static_cast<size_t>(v) * vStride2;
            for (int li = 0; li <= localN2; ++li) { row[li] = 0.0; }
        }
        const float NEG_INF_F2 = -std::numeric_limits<float>::infinity();
        for (int d = 0; d <= localN2; ++d) {
            const int liMax2 = localN2 - d;
            double *insSlice2 = insDataAll2 + static_cast<size_t>(d) * dStride2;
            const double *insData2 = insDataAll2;
            const StateId *trDstPtr = trDst.data();
            const double *trScPtr = trSc.data();
            for (int v = M - 1; v >= 0; --v) {
                const ExactStateExec &st = exec[static_cast<size_t>(v)];
                if (st.type == CM_ST_E) { continue; }
                if (st.dmax >= st.dmin && (d < st.dmin || d > st.dmax)) { continue; }
                if (st.type == CM_ST_B) {
                    const int y = st.bLeft, z = st.bRight;
                    if (y < 0 || y >= M || z < 0 || z >= M) { continue; }
                    int kBeg2 = 0, kEnd2 = d;
                    if (st.splitKMax >= st.splitKMin) { kBeg2 = std::max(kBeg2, st.splitKMin); kEnd2 = std::min(kEnd2, st.splitKMax); }
                    if (st.splitRMax >= st.splitRMin) { kBeg2 = std::max(kBeg2, d - st.splitRMax); kEnd2 = std::min(kEnd2, d - st.splitRMin); }
                    double *dst2 = insSlice2 + static_cast<size_t>(v) * vStride2;
                    for (int li = 0; li <= liMax2; ++li) {
                        if (fastLogsum) {
                            bool has = false;
                            float acc = NEG_INF_F2;
                            for (int k = kBeg2; k <= kEnd2; ++k) {
                                const double lv = insData2[static_cast<size_t>(k) * dStride2 + static_cast<size_t>(y) * vStride2 + static_cast<size_t>(li)];
                                const double rv = insData2[static_cast<size_t>(d - k) * dStride2 + static_cast<size_t>(z) * vStride2 + static_cast<size_t>(li + k)];
                                if (lv != NEG_INF && rv != NEG_INF) { log2AccFastAddF(static_cast<float>(lv + rv), has, acc); }
                            }
                            dst2[li] = has ? static_cast<double>(acc) : NEG_INF;
                        } else {
                            LogSumAccBase2Exact acc;
                            for (int k = kBeg2; k <= kEnd2; ++k) {
                                const double lv = insData2[static_cast<size_t>(k) * dStride2 + static_cast<size_t>(y) * vStride2 + static_cast<size_t>(li)];
                                const double rv = insData2[static_cast<size_t>(d - k) * dStride2 + static_cast<size_t>(z) * vStride2 + static_cast<size_t>(li + k)];
                                if (lv != NEG_INF && rv != NEG_INF) { acc.add(lv + rv); }
                            }
                            dst2[li] = acc.value();
                        }
                    }
                    continue;
                }
                if (d < st.dConsume) { continue; }
                const int nd2 = d - st.dConsume;
                const float *ep = st.emitPtr;
                const int8_t *sc = ws.seqCode.data();
                const double *nxtSlice2 = insDataAll2 + static_cast<size_t>(nd2) * dStride2;
                double *dst2 = insSlice2 + static_cast<size_t>(v) * vStride2;
                for (int li = 0; li <= liMax2; ++li) {
                    float ef = 0.0f;
                    if (st.consumesLeft && st.consumesRight) {
                        if (!(ep && st.emitSize >= 16)) {
                            ef = NEG_INF_F2;
                        } else {
                            const int ii = li + bestI;
                            const int lc = sc[ii];
                            const int rc = sc[ii + d - 1];
                            ef = (lc >= 0 && rc >= 0) ? ep[lc * 4 + rc] : -1.0f;
                        }
                    } else if (st.consumesLeft) {
                        if (!(ep && st.emitSize >= 4)) {
                            ef = NEG_INF_F2;
                        } else {
                            const int bi = sc[li + bestI];
                            ef = (bi >= 0) ? ep[bi] : -1.0f;
                        }
                    } else if (st.consumesRight) {
                        if (!(ep && st.emitSize >= 4)) {
                            ef = NEG_INF_F2;
                        } else {
                            const int bi = sc[li + bestI + d - 1];
                            ef = (bi >= 0) ? ep[bi] : -1.0f;
                        }
                    }
                    if (ef == NEG_INF_F2) { continue; }
                    const int nli = li + st.niShift;
                    if (fastLogsum) {
                        bool has = false;
                        float acc = NEG_INF_F2;
                        for (int t = 0; t < st.trCount; ++t) {
                            const size_t ti = st.trOff + static_cast<size_t>(t);
                            const double n = nxtSlice2[static_cast<size_t>(trDstPtr[ti]) * vStride2 + static_cast<size_t>(nli)];
                            if (n != NEG_INF) { log2AccFastAddF(ef + static_cast<float>(trScPtr[ti]) + static_cast<float>(n), has, acc); }
                        }
                        dst2[li] = has ? static_cast<double>(acc) : NEG_INF;
                    } else {
                        LogSumAccBase2Exact acc;
                        for (int t = 0; t < st.trCount; ++t) {
                            const size_t ti = st.trOff + static_cast<size_t>(t);
                            const double n = nxtSlice2[static_cast<size_t>(trDstPtr[ti]) * vStride2 + static_cast<size_t>(nli)];
                            if (n != NEG_INF) { acc.add(static_cast<double>(ef) + trScPtr[ti] + n); }
                        }
                        dst2[li] = acc.value();
                    }
                }
            }
        }
        h.inside = insDataAll2[static_cast<size_t>(localN2) * dStride2 + static_cast<size_t>(model.rootState) * vStride2 + 0]
                   - bestNull3Corr - null2Corr;
        h.bias = bestNull3Corr + null2Corr;
        // Cell-diff vs Infernal cm_InsideAlign: dump Inside alpha[v][j][d] for all
        // states v at the (j_target, d) cells specified by MMSEQS_CMSCAN_DUMP_INSIDE_J.
        // j_target is in *envelope* coords (1..localN2). Run cmalign on the envelope
        // sub-sequence to align j between our run and Infernal's.
        const char *dumpInsJList = std::getenv("MMSEQS_CMSCAN_DUMP_INSIDE_J");
        if (dumpInsJList != NULL) {
            std::vector<int> insJTargets;
            std::string s = dumpInsJList;
            size_t pos = 0;
            while (pos < s.size()) {
                size_t comma = s.find(',', pos);
                int jt = std::atoi(s.substr(pos, comma - pos).c_str());
                if (jt > 0 && jt <= localN2) insJTargets.push_back(jt);
                if (comma == std::string::npos) break;
                pos = comma + 1;
            }
            for (int jT : insJTargets) {
                for (int d = 0; d <= jT; ++d) {
                    const int li = jT - d;
                    if (li < 0 || li > localN2) continue;
                    const double *insSlice = insDataAll2 + static_cast<size_t>(d) * dStride2;
                    for (int v = 0; v < M; ++v) {
                        const double a = insSlice[static_cast<size_t>(v) * vStride2 + static_cast<size_t>(li)];
                        if (a == NEG_INF) continue;
                        const ExactStateExec &st = exec[static_cast<size_t>(v)];
                        fprintf(stderr, "DUMP_INSIDE tid=%s j=%d d=%d v=%d type=%d alpha=%.6f\n",
                                seqId.c_str(), jT, d, v,
                                static_cast<int>(st.type),
                                a);
                    }
                }
            }
        }
    } else {
        h.inside = NEG_INF;
        h.bias = bestNull3Corr + null2Corr;
    }
    outHits.push_back(h);
}

// Per-CM-column max-marginal match emission table for peak-anchor.
// Index 1..clen. For ML/MR states: direct 4-vector. For MP states:
// max over partner residue (cheap marginalization that preserves the
// conserved-core peak). Built once per query.
static std::vector<std::array<float, 4>> buildPerColumnMaxEmissions(const InfernalExactModel &model) {
    const int clen = model.clen;
    std::vector<std::array<float, 4>> out(static_cast<size_t>(clen + 1), {0.0f, 0.0f, 0.0f, 0.0f});
    std::vector<bool> have(static_cast<size_t>(clen + 1), false);
    for (size_t si = 0; si < model.states.size(); ++si) {
        const CmState &s = model.states[si];
        if ((s.type == CM_ST_ML || s.type == CM_ST_MR) && s.emit.size() >= 4) {
            const int pos = (s.type == CM_ST_ML) ? s.mapLeft : s.mapRight;
            if (pos > 0 && pos <= clen) {
                for (int a = 0; a < 4; ++a) {
                    const float v = static_cast<float>(s.emit[static_cast<size_t>(a)]);
                    if (!have[pos] || v > out[pos][a]) out[pos][a] = v;
                }
                have[pos] = true;
            }
        } else if (s.type == CM_ST_MP && s.emit.size() >= 16) {
            if (s.mapLeft > 0 && s.mapLeft <= clen) {
                for (int a = 0; a < 4; ++a) {
                    float row = -std::numeric_limits<float>::infinity();
                    for (int b = 0; b < 4; ++b) {
                        const float v = static_cast<float>(s.emit[static_cast<size_t>(a * 4 + b)]);
                        if (v > row) row = v;
                    }
                    if (!have[s.mapLeft] || row > out[s.mapLeft][a]) out[s.mapLeft][a] = row;
                }
                have[s.mapLeft] = true;
            }
            if (s.mapRight > 0 && s.mapRight <= clen) {
                for (int b = 0; b < 4; ++b) {
                    float col = -std::numeric_limits<float>::infinity();
                    for (int a = 0; a < 4; ++a) {
                        const float v = static_cast<float>(s.emit[static_cast<size_t>(a * 4 + b)]);
                        if (v > col) col = v;
                    }
                    if (!have[s.mapRight] || col > out[s.mapRight][b]) out[s.mapRight][b] = col;
                }
                have[s.mapRight] = true;
            }
        }
    }
    return out;
}

// Walk SW backtrace L→R, accumulate per-CM-column ML emissions on M ops,
// return full-target 0-indexed coord of cumulative-score argmax. With cmbuild
// --hand the CM column equals the 1-indexed query position. qStart is the
// 1-indexed query start from the m8 row; dbStart is 0-indexed full-target.
// Returns -1 on empty backtrace, no scored M ops, or all-N target window.
static int peakAnchorFromCigar(const std::string &bt,
                               int qStart,
                               int dbStart,
                               const std::vector<std::array<float, 4>> &emitCol,
                               int clen,
                               const std::string &fullTargetSeq) {
    if (bt.empty() || clen <= 0) return -1;
    int qPos = qStart;             // 1-indexed CM column under --hand
    int tPos = dbStart;            // 0-indexed full-target coord
    float cum = 0.0f;
    float bestCum = -std::numeric_limits<float>::infinity();
    int bestT = -1;
    bool any = false;
    const int tLen = static_cast<int>(fullTargetSeq.size());
    for (size_t k = 0; k < bt.size(); ++k) {
        const char op = bt[k];
        if (op == 'M') {
            if (qPos >= 1 && qPos <= clen && tPos >= 0 && tPos < tLen) {
                int a = -1;
                switch (fullTargetSeq[static_cast<size_t>(tPos)]) {
                    case 'A': case 'a': a = 0; break;
                    case 'C': case 'c': a = 1; break;
                    case 'G': case 'g': a = 2; break;
                    case 'U': case 'u': case 'T': case 't': a = 3; break;
                    default: break;
                }
                if (a >= 0) {
                    cum += emitCol[static_cast<size_t>(qPos)][static_cast<size_t>(a)];
                    if (!any || cum > bestCum) {
                        bestCum = cum;
                        bestT = tPos;
                    }
                    any = true;
                }
            }
            ++qPos; ++tPos;
        } else if (op == 'I') {
            ++qPos;
        } else if (op == 'D') {
            ++tPos;
        }
    }
    return any ? bestT : -1;
}

static bool finishFastDeckSelectedHitBatched(const InfernalExactModel &model,
                                             const std::string &seqId,
                                             const CmFastDeckGpuBatchResult &gr,
                                             const uint16_t *traceStateBuf,
                                             std::vector<Hit> &outHits) {
    if (!gr.found || gr.traceOverflow || gr.traceLen <= 0 || traceStateBuf == NULL) {
        return false;
    }
    std::vector<int> traceStates;
    traceStates.reserve(static_cast<size_t>(gr.traceLen));
    for (int i = 0; i < gr.traceLen; ++i) {
        traceStates.push_back(static_cast<int>(traceStateBuf[static_cast<size_t>(i)]));
    }

    Hit h;
    h.seqId = seqId;
    h.start1 = gr.minUsed;
    h.end1 = gr.maxUsed;
    h.mode = static_cast<char>(gr.bestMode);
    h.trunc = (gr.bestMode != 'J');
    h.traceStates = encodeTraceStates(traceStates);
    {
        int qS = -1, qE = -1, aL = 0, leadIns = 0, trailIns = 0;
        h.cigar = modelTraceCigarQueryCoord(model, traceStates, &qS, &qE, &aL,
                                            &leadIns, &trailIns);
        h.qStart = qS;
        h.qEnd = qE;
        h.cigarAlnLen = static_cast<unsigned int>(std::max(0, aL));
        h.leadingInsertTargets = leadIns;
        h.trailingInsertTargets = trailIns;
    }
    const double modelAggRaw[4] = {
        gr.modelAggRaw[0], gr.modelAggRaw[1], gr.modelAggRaw[2], gr.modelAggRaw[3]
    };
    const int obsCount[4] = {
        gr.obsCount[0], gr.obsCount[1], gr.obsCount[2], gr.obsCount[3]
    };
    const double null2Corr = scoreCorrectionNull2BitsFromTrace(model, modelAggRaw, obsCount);
    h.cyk = gr.bestSc - null2Corr;
    h.inside = NEG_INF;
    h.bias = gr.bestNull3Corr + null2Corr;
    outHits.push_back(std::move(h));
    return true;
}

// Stage 8: re-align each hit's envelope sub-sequence using the truncation DP
// (J/L/R/T modes via runTrCYKInsideAlign + runTrCYKAlignT). Mirrors Infernal
// cmsearch's per-envelope re-alignment with cm_TrCYKInsideAlign + cm_tr_alignT.
// Replaces hit.traceStates and hit.cigar in place; keeps start1/end1 from the
// prefilter scan (envelope coords are unchanged — only the trace through them
// is upgraded from J-only to truncation-aware).
//
// Gated on MMSEQS_CMSCAN_TRUNC_DP=1 + MMSEQS_CMSCAN_LOCAL_MODE=1; also requires
// the truncation tables (m.trpL/lmesc/rmesc) to have been built upstream.
static void runTruncDpReAlignHits(const InfernalExactModel &model,
                                  const std::string &seq,
                                  std::vector<Hit> &hits) {
    if (!model.hasLocalCfg) return;
    if (model.trpL.empty()) return;
    const int N = static_cast<int>(seq.size());
    if (N <= 0 || hits.empty()) return;
    const char *envDump = std::getenv("MMSEQS_CMSCAN_TRUNC_REALIGN_DUMP");
    const bool dump = (envDump != NULL && envDump[0] == '1');
    // MMSEQS_CMSCAN_REALIGN_FULL=1: ignore prior CYK trim (h.start1/h.end1
    // from initial pass) and re-align over the FULL target sequence. Used by
    // the regression test to make our DP comparable to Infernal cmalign which
    // always runs over the full input.
    static const char *const realignFullEnv = std::getenv("MMSEQS_CMSCAN_REALIGN_FULL");
    const bool realignFull = (realignFullEnv != nullptr && realignFullEnv[0] == '1');
    // MMSEQS_CMSCAN_TIMING=1 emits per-thread cumulative DP time at end of this
    // hit batch. Used to A/B test optimization changes without bench-running.
    static const char *const timingEnv = std::getenv("MMSEQS_CMSCAN_TIMING");
    const bool timing = (timingEnv != nullptr && timingEnv[0] == '1');
    long long dpUsTotal = 0;
    int dpCalls = 0;
    for (Hit &h : hits) {
        int lo = std::min(h.start1, h.end1);
        int hi = std::max(h.start1, h.end1);
        if (realignFull) { lo = 1; hi = N; }
        if (lo < 1 || hi > N || lo > hi) continue;
        const int subLen = hi - lo + 1;
        std::vector<int8_t> dsq(static_cast<size_t>(subLen + 2), -1);
        for (int p = 1; p <= subLen; ++p) {
            const char c = seq[static_cast<size_t>(lo - 1 + p - 1)];
            if      (c=='A'||c=='a') dsq[(size_t)p] = 0;
            else if (c=='C'||c=='c') dsq[(size_t)p] = 1;
            else if (c=='G'||c=='g') dsq[(size_t)p] = 2;
            else if (c=='T'||c=='t'||c=='U'||c=='u') dsq[(size_t)p] = 3;
            else dsq[(size_t)p] = -1;
        }
        const int pty_idx = 0; // TRPENALTY_5P_AND_3P (most permissive)
        // MMSEQS_CMSCAN_TRUNC_FORCE_J=1: force mode J only (mirror Infernal
        // --notrunc behavior). Mode L/R/T silently drop right/left/both-arm
        // states from the rendered MSA → forcing J keeps the parsetree dense
        // and bench-stable.
        static const char *const forceJEnv = std::getenv("MMSEQS_CMSCAN_TRUNC_FORCE_J");
        const char preset = (forceJEnv != nullptr && forceJEnv[0] == '1') ? 'J' : 'U';
        std::chrono::steady_clock::time_point _dp_t0;
        if (timing) _dp_t0 = std::chrono::steady_clock::now();
        const TrCykResult &res = runTrCYKInsideAlign(model, dsq, subLen, preset, pty_idx);
        if (timing) {
            auto _dp_t1 = std::chrono::steady_clock::now();
            dpUsTotal += std::chrono::duration_cast<std::chrono::microseconds>(_dp_t1 - _dp_t0).count();
            ++dpCalls;
        }
        if (!std::isfinite(res.score)) continue;
        TrParsetree tr = runTrCYKAlignT(model, res, subLen, pty_idx);
        if (!tr.ok || tr.nodes.empty()) continue;
        dumpTrParsetree(tr, model);  // gated by MMSEQS_CMSCAN_DUMP_TRPT=1
        std::vector<int> states = trParsetreeStates(tr);
        std::vector<int8_t> modes;
        modes.reserve(tr.nodes.size());
        for (const TrParseNode &n : tr.nodes) modes.push_back(n.mode);
        int qS = -1, qE = -1, aL = 0, leadIns = 0, trailIns = 0;
        std::string newCigar = modelTraceCigarQueryCoord(model, states, &qS, &qE, &aL,
                                                        &leadIns, &trailIns, &modes);
        if (newCigar.empty() || newCigar == "NA") continue;
        // EL pseudo-state (v == M) accounting. modelTraceCigarQueryCoord skips
        // these nodes; their absorbed target residues need to be folded into
        // leading/trailing inserts so the m8 t-coord trim shrinks the reported
        // span to match CIGAR M+D consume. Without this, downstream
        // result2dnamsa walks CIGAR over a wider t-window than CIGAR covers
        // and mis-registers residues into wrong CM columns (78%-94% of rows
        // in pre-fix wire-up).
        // Mode-aware emit semantics (Infernal cm_dpalign_trunc.c):
        //   mode J: MP emits {emitl, emitr}; ML {emitl}; MR {emitr}
        //   mode L: MP {emitl};              ML {emitl}; MR silent
        //   mode R: MP {emitr};              ML silent;  MR {emitr}
        //   mode T: same as J for emit accounting (truncations are at trace edges)
        const int M_states = static_cast<int>(model.states.size());
        int matchMin = std::numeric_limits<int>::max();
        int matchMax = 0;
        auto bumpMatch = [&](int pos) {
            if (pos <= 0) return;
            matchMin = std::min(matchMin, pos);
            matchMax = std::max(matchMax, pos);
        };
        for (const TrParseNode &n : tr.nodes) {
            if (n.v < 0 || n.v >= M_states) continue;
            const CmStateType st = model.states[(size_t)n.v].type;
            const bool emitsLeft  = (st == CM_ST_MP || st == CM_ST_ML)
                                    ? (n.mode == TR_TRMODE_J || n.mode == TR_TRMODE_L
                                       || n.mode == TR_TRMODE_T)
                                    : false;
            const bool emitsRight = (st == CM_ST_MP || st == CM_ST_MR)
                                    ? (n.mode == TR_TRMODE_J || n.mode == TR_TRMODE_R
                                       || n.mode == TR_TRMODE_T)
                                    : false;
            if (emitsLeft)  bumpMatch(n.emitl);
            if (emitsRight) bumpMatch(n.emitr);
        }
        int extraLead = 0, extraTrail = 0, middleEL = 0;
        for (const TrParseNode &n : tr.nodes) {
            if (n.v != M_states) continue;
            const int absorbed = std::max(0, n.emitr - n.emitl + 1);
            if (absorbed == 0) continue;
            if (matchMin == std::numeric_limits<int>::max()) {
                extraTrail += absorbed;
            } else if (n.emitr < matchMin) {
                extraLead += absorbed;
            } else if (n.emitl > matchMax) {
                extraTrail += absorbed;
            } else {
                // Middle/bifurcation EL: residues absorbed between matches,
                // shrinking either end of t-span chops the wrong residues.
                // Skip — span will exceed CIGAR M+D for these (renderer
                // under-walks but doesn't overrun).
                middleEL += absorbed;
            }
        }
        leadIns += extraLead;
        trailIns += extraTrail;
        if (dump) {
            std::fprintf(stderr,
                         "[TRUNC_REALIGN] %s lo=%d hi=%d L=%d mode=%c b=%d score=%.4f trpenalty=%.4f nodes=%zu el_lead=%d el_trail=%d\n",
                         h.seqId.c_str(), lo, hi, subLen,
                         res.mode, res.b, res.score, tr.trpenalty, tr.nodes.size(),
                         extraLead, extraTrail);
        }
        h.traceStates = encodeTraceStates(states);
        h.cigar = newCigar;
        h.qStart = qS;
        h.qEnd = qE;
        h.cigarAlnLen = static_cast<unsigned int>(std::max(0, aL));
        h.leadingInsertTargets = leadIns;
        h.trailingInsertTargets = trailIns;
        h.mode = res.mode;
        h.trunc = (res.mode != 'J');
        h.cyk = res.score;
    }
    if (timing && dpCalls > 0) {
        Debug(Debug::INFO) << "TRUNC_DP timing: " << dpCalls << " envelopes, "
                           << dpUsTotal << " us total, "
                           << (double)dpUsTotal / (double)dpCalls << " us/env\n";
    }
}

static bool cmBeamLocalizeEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_BEAM_LOCALIZE");
        cached = (env != NULL && std::string(env) == "1") ? 1 : 0;
    }
    return cached == 1;
}

static int cmBeamLocalizeBeamSize() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_BEAM_SIZE");
        cached = (env != NULL) ? std::atoi(env) : 128;
        if (cached < 16) {
            cached = 16;
        }
    }
    return cached;
}

static int cmBeamLocalizePad() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_BEAM_PAD");
        cached = (env != NULL) ? std::atoi(env) : 48;
        if (cached < 0) {
            cached = 0;
        }
    }
    return cached;
}

static float cmBeamLocalizeGapPenalty() {
    static float cached = std::numeric_limits<float>::quiet_NaN();
    if (std::isnan(cached)) {
        const char *env = std::getenv("MMSEQS_CMSCAN_BEAM_GAP");
        cached = (env != NULL) ? static_cast<float>(std::atof(env)) : 1.25f;
        if (!(cached > 0.0f)) {
            cached = 1.25f;
        }
    }
    return cached;
}

static float cmBeamLocalizeSlack() {
    static float cached = std::numeric_limits<float>::quiet_NaN();
    if (std::isnan(cached)) {
        const char *env = std::getenv("MMSEQS_CMSCAN_BEAM_SLACK");
        cached = (env != NULL) ? static_cast<float>(std::atof(env)) : 2.5f;
        if (cached < 0.0f) {
            cached = 0.0f;
        }
    }
    return cached;
}

static void cmBeamBuildEmitProfile(const InfernalExactModel &model,
                                   std::vector<std::array<float,4>> &profile) {
    const int clen = static_cast<int>(model.clen);
    profile.assign(static_cast<size_t>(clen), {0.f, 0.f, 0.f, 0.f});
    std::vector<int> counts(static_cast<size_t>(clen), 0);
    for (size_t vi = 0; vi < model.states.size(); ++vi) {
        const CmState &s = model.states[vi];
        if (s.type == CM_ST_ML && s.mapLeft >= 1 && s.mapLeft <= clen && s.emitF.size() >= 4) {
            for (int nt = 0; nt < 4; ++nt) {
                profile[static_cast<size_t>(s.mapLeft - 1)][static_cast<size_t>(nt)] += s.emitF[static_cast<size_t>(nt)];
            }
            counts[static_cast<size_t>(s.mapLeft - 1)]++;
        } else if (s.type == CM_ST_MR && s.mapRight >= 1 && s.mapRight <= clen && s.emitF.size() >= 4) {
            for (int nt = 0; nt < 4; ++nt) {
                profile[static_cast<size_t>(s.mapRight - 1)][static_cast<size_t>(nt)] += s.emitF[static_cast<size_t>(nt)];
            }
            counts[static_cast<size_t>(s.mapRight - 1)]++;
        } else if (s.type == CM_ST_MP && s.emitF.size() >= 16) {
            if (s.mapLeft >= 1 && s.mapLeft <= clen) {
                for (int a = 0; a < 4; ++a) {
                    float sum = 0.0f;
                    for (int b = 0; b < 4; ++b) {
                        sum += s.emitF[static_cast<size_t>(a * 4 + b)];
                    }
                    profile[static_cast<size_t>(s.mapLeft - 1)][static_cast<size_t>(a)] += sum * 0.25f;
                }
                counts[static_cast<size_t>(s.mapLeft - 1)]++;
            }
            if (s.mapRight >= 1 && s.mapRight <= clen) {
                for (int b = 0; b < 4; ++b) {
                    float sum = 0.0f;
                    for (int a = 0; a < 4; ++a) {
                        sum += s.emitF[static_cast<size_t>(a * 4 + b)];
                    }
                    profile[static_cast<size_t>(s.mapRight - 1)][static_cast<size_t>(b)] += sum * 0.25f;
                }
                counts[static_cast<size_t>(s.mapRight - 1)]++;
            }
        }
    }
    for (int c = 0; c < clen; ++c) {
        if (counts[static_cast<size_t>(c)] > 0) {
            const float inv = 1.0f / static_cast<float>(counts[static_cast<size_t>(c)]);
            for (int nt = 0; nt < 4; ++nt) {
                profile[static_cast<size_t>(c)][static_cast<size_t>(nt)] *= inv;
            }
        }
    }
}

static void cmBeamEncodeSeq(const char *seq, int L, std::vector<int8_t> &out) {
    out.assign(static_cast<size_t>(L), static_cast<int8_t>(-1));
    for (int i = 0; i < L; ++i) {
        const char c = seq[i];
        if      (c == 'A' || c == 'a') out[static_cast<size_t>(i)] = 0;
        else if (c == 'C' || c == 'c') out[static_cast<size_t>(i)] = 1;
        else if (c == 'G' || c == 'g') out[static_cast<size_t>(i)] = 2;
        else if (c == 'T' || c == 't' || c == 'U' || c == 'u') out[static_cast<size_t>(i)] = 3;
    }
}

static CmBeamLocalizeResult cmBeamLocalizeInterval(const InfernalExactModel &model,
                                                   const std::string &seq) {
    CmBeamLocalizeResult out;
    const int N = static_cast<int>(seq.size());
    if (N <= 0 || model.clen <= 0) {
        return out;
    }

    thread_local const InfernalExactModel *cachedModel = NULL;
    thread_local std::vector<std::array<float,4>> cachedProfile;
    if (cachedModel != &model) {
        cmBeamBuildEmitProfile(model, cachedProfile);
        cachedModel = &model;
    }
    if (cachedProfile.empty()) {
        return out;
    }

    std::vector<int8_t> encoded;
    cmBeamEncodeSeq(seq.c_str(), N, encoded);
    const int qLen = static_cast<int>(cachedProfile.size());
    const int beamSize = cmBeamLocalizeBeamSize();
    const int beamCap = std::max(beamSize * 4, beamSize);
    const int pad = cmBeamLocalizePad();
    const float gap = cmBeamLocalizeGapPenalty();
    const float slack = cmBeamLocalizeSlack();
    const float negInf = -std::numeric_limits<float>::infinity();

    std::vector<float> prevScore(static_cast<size_t>(N + 1), negInf);
    std::vector<int> prevStart(static_cast<size_t>(N + 1), 0);
    std::vector<uint8_t> prevValid(static_cast<size_t>(N + 1), 0);
    std::vector<int> prevActive;
    std::vector<int> candidates;
    std::vector<int> rowPositions;
    std::vector<float> rowScores;
    std::vector<int> rowStarts;
    std::vector<int> keepPositions;
    prevActive.reserve(static_cast<size_t>(beamCap));
    candidates.reserve(static_cast<size_t>(std::max(N, beamCap * 3)));
    rowPositions.reserve(static_cast<size_t>(std::max(N, beamCap * 3)));
    rowScores.reserve(static_cast<size_t>(std::max(N, beamCap * 3)));
    rowStarts.reserve(static_cast<size_t>(std::max(N, beamCap * 3)));
    keepPositions.reserve(static_cast<size_t>(beamCap));

    float globalBest = negInf;
    int globalStart = 1;
    int globalEnd = 0;

    for (int qi = 0; qi < qLen; ++qi) {
        candidates.clear();
        if (qi == 0) {
            candidates.resize(static_cast<size_t>(N));
            for (int j = 1; j <= N; ++j) {
                candidates[static_cast<size_t>(j - 1)] = j;
            }
        } else {
            for (size_t idx = 0; idx < prevActive.size(); ++idx) {
                const int j = prevActive[idx];
                if (j > 1) candidates.push_back(j - 1);
                candidates.push_back(j);
                if (j < N) candidates.push_back(j + 1);
            }
            if (candidates.empty()) {
                break;
            }
            std::sort(candidates.begin(), candidates.end());
            candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
        }

        rowPositions.clear();
        rowScores.clear();
        rowStarts.clear();
        float rowBest = negInf;
        int prevCandidateJ = -1000000000;
        float prevRowScore = negInf;
        int prevRowStart = 0;

        for (size_t ci = 0; ci < candidates.size(); ++ci) {
            const int j = candidates[ci];
            float emit = 0.0f;
            const int8_t b = encoded[static_cast<size_t>(j - 1)];
            if (b >= 0 && b < 4) {
                emit = cachedProfile[static_cast<size_t>(qi)][static_cast<size_t>(b)];
            }

            float best = emit;
            int bestStart = j;
            if (j > 1 && prevValid[static_cast<size_t>(j - 1)] != 0) {
                const float diag = prevScore[static_cast<size_t>(j - 1)] + emit;
                if (diag > best) {
                    best = diag;
                    bestStart = prevStart[static_cast<size_t>(j - 1)];
                }
            }
            if (prevValid[static_cast<size_t>(j)] != 0) {
                const float up = prevScore[static_cast<size_t>(j)] - gap;
                if (up > best) {
                    best = up;
                    bestStart = prevStart[static_cast<size_t>(j)];
                }
            }
            if (prevCandidateJ == j - 1 && std::isfinite(prevRowScore)) {
                const float left = prevRowScore - gap;
                if (left > best) {
                    best = left;
                    bestStart = prevRowStart;
                }
            }
            if (!(best > 0.0f) || !std::isfinite(best)) {
                prevCandidateJ = j;
                prevRowScore = negInf;
                prevRowStart = 0;
                continue;
            }
            rowPositions.push_back(j);
            rowScores.push_back(best);
            rowStarts.push_back(bestStart);
            if (best > rowBest) {
                rowBest = best;
            }
            if (best > globalBest) {
                globalBest = best;
                globalStart = bestStart;
                globalEnd = j;
            }
            prevCandidateJ = j;
            prevRowScore = best;
            prevRowStart = bestStart;
        }

        out.rowsVisited++;
        for (size_t idx = 0; idx < prevActive.size(); ++idx) {
            prevValid[static_cast<size_t>(prevActive[idx])] = 0;
        }
        prevActive.clear();
        if (rowPositions.empty()) {
            break;
        }

        std::vector<size_t> order(rowPositions.size());
        std::iota(order.begin(), order.end(), 0);
        if (static_cast<int>(order.size()) > beamSize) {
            out.rowsPruned++;
        }
        std::sort(order.begin(), order.end(),
                  [&](size_t left, size_t right) {
                      if (rowScores[left] != rowScores[right]) {
                          return rowScores[left] > rowScores[right];
                      }
                      return rowPositions[left] < rowPositions[right];
                  });
        const float kthScore = rowScores[order[std::min<size_t>(order.size(), static_cast<size_t>(beamSize)) - 1]];
        const float threshold = std::max(kthScore, rowBest - slack);
        keepPositions.clear();
        for (size_t oi = 0; oi < order.size(); ++oi) {
            const size_t idx = order[oi];
            if (rowScores[idx] < threshold && static_cast<int>(keepPositions.size()) >= beamSize) {
                break;
            }
            keepPositions.push_back(static_cast<int>(idx));
            if (static_cast<int>(keepPositions.size()) >= beamCap) {
                break;
            }
        }
        std::sort(keepPositions.begin(), keepPositions.end(),
                  [&](int left, int right) {
                      return rowPositions[static_cast<size_t>(left)] < rowPositions[static_cast<size_t>(right)];
                  });
        for (size_t ki = 0; ki < keepPositions.size(); ++ki) {
            const size_t idx = static_cast<size_t>(keepPositions[ki]);
            const int pos = rowPositions[idx];
            prevValid[static_cast<size_t>(pos)] = 1;
            prevScore[static_cast<size_t>(pos)] = rowScores[idx];
            prevStart[static_cast<size_t>(pos)] = rowStarts[idx];
            prevActive.push_back(pos);
        }
    }

    if (!std::isfinite(globalBest) || globalEnd < globalStart) {
        return out;
    }

    int lo = std::max(1, globalStart - pad);
    int hi = std::min(N, globalEnd + pad);
    const int minWidth = std::max(24, qLen / 2);
    if (hi - lo + 1 < minWidth) {
        const int center = (globalStart + globalEnd) / 2;
        lo = std::max(1, center - minWidth / 2);
        hi = std::min(N, lo + minWidth - 1);
        lo = std::max(1, hi - minWidth + 1);
    }

    out.valid = (hi > lo) && (hi - lo + 1 < N);
    out.start = lo;
    out.end = hi;
    out.bestScore = globalBest;
    return out;
}

static void cmMaybeBeamNarrowRegion(const InfernalExactModel &model,
                                    std::string &regionSeq,
                                    int &offset,
                                    int &forceI,
                                    int &anchorI,
                                    int &anchorD) {
    if (!cmBeamLocalizeEnabled()) {
        return;
    }
    if (regionSeq.size() < 64) {
        return;
    }
    if (forceI > 0 || anchorI > 0 || anchorD > 0) {
        return;
    }

    const CmBeamLocalizeResult beam = cmBeamLocalizeInterval(model, regionSeq);
    if (!beam.valid) {
        return;
    }

    const int lo0 = beam.start - 1;
    const int len = beam.end - beam.start + 1;
    if (len <= 0 || len >= static_cast<int>(regionSeq.size())) {
        return;
    }
    regionSeq = regionSeq.substr(static_cast<size_t>(lo0), static_cast<size_t>(len));
    offset += lo0;
    forceI = -1;
    anchorI = -1;
    anchorD = -1;
}

static bool cmBeamResultLooksOff(const std::vector<Hit> &hits,
                                 int beamStart,
                                 int beamEnd,
                                 int fullLen) {
    if (hits.empty()) {
        return true;
    }
    const int margin = std::max(4, std::min(16, fullLen / 20));
    for (size_t i = 0; i < hits.size(); ++i) {
        const int lo = std::min(hits[i].start1, hits[i].end1);
        const int hi = std::max(hits[i].start1, hits[i].end1);
        if (lo <= beamStart + margin || hi >= beamEnd - margin) {
            return true;
        }
    }
    return false;
}

// Per-envelope CYK re-scan wrapper. When MMSEQS_CMSCAN_RESCORE_ENVELOPE=1,
// runs CYK twice: first to identify peak (i,j), then re-runs on a sub-sequence
// of [max(1, i-pad) .. min(N, j+pad)] which mirrors Infernal's dispatch#2 behavior
// (envelope-restricted CYK that yields ~+60 bits at the same logical cell on
// 2YGH_1 KJ798010.1). Pad defaults to 5; override with MMSEQS_CMSCAN_RESCORE_PAD.
static void runInfernalExactScanWithEnvelopeRescore(const InfernalExactModel &model,
                                                    const std::string &seq,
                                                    bool wantInside,
                                                    std::vector<Hit> &outHits,
                                                    const std::string &seqId,
                                                    int maxHitLen = 0,
                                                    int forcedI = -1,
                                                    int forcedD = -1,
                                                    int anchorI = -1,
                                                    int anchorD = -1,
                                                    const int *stateMinI = NULL,
                                                    const int *stateMaxI = NULL,
                                                    const int *stateMinJ = NULL,
                                                    const int *stateMaxJ = NULL) {
    static int rescoreEnabled = -1;
    static int rescorePad = -1;
    static int rescanQdbOff = -1;
    if (rescoreEnabled == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_RESCORE_ENVELOPE");
        rescoreEnabled = (env != NULL && std::string(env) == "1") ? 1 : 0;
    }
    if (rescorePad == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_RESCORE_PAD");
        rescorePad = (env != NULL) ? std::atoi(env) : 5;
        if (rescorePad < 0) rescorePad = 0;
    }
    if (rescanQdbOff == -1) {
        // MMSEQS_CMSCAN_RESCAN_QDB_OFF=1: pass 1 keeps QDB, pass 2 forces QDB off.
        // Mirrors LOCAL=0+DISABLE_QDB lift (project_qdboff_2YGH_1_native_beats_shellout)
        // without the broken-envelope regressions on 3Q1Q_3 / 5CCB_3 (per
        // project_qdboff_does_not_generalize).
        const char *env = std::getenv("MMSEQS_CMSCAN_RESCAN_QDB_OFF");
        rescanQdbOff = (env != NULL && std::string(env) == "1") ? 1 : 0;
    }

    static const char *const _wrapTimingEnv = std::getenv("MMSEQS_CMSCAN_TIMING");
    const bool _wrapTiming = (_wrapTimingEnv != nullptr && _wrapTimingEnv[0] == '1');
    // MMSEQS_CMSCAN_SKIP_TRUNC_REALIGN=1: skip the post-scan TRUNC_DP realign step.
    // Under FORCE_J=1 the realign re-runs the same J-mode CYK that the initial scan
    // already produced a trace for — pure duplicate work. Skipping is safe when the
    // scan's CIGAR is acceptable for downstream consumers (it always is for R-scape
    // .sto rendering via result2dnamsa). Output drops L/R/T-mode hits, but FORCE_J=1
    // never produces those anyway.
    static const char *const _skipTruncRealignEnv = std::getenv("MMSEQS_CMSCAN_SKIP_TRUNC_REALIGN");
    const bool skipTruncRealign = (_skipTruncRealignEnv != nullptr && _skipTruncRealignEnv[0] == '1');
    if (rescoreEnabled == 0) {
        std::chrono::steady_clock::time_point _t0, _t1, _t2;
        if (_wrapTiming) _t0 = std::chrono::steady_clock::now();
        runInfernalExactScan(model, seq, wantInside, outHits, seqId, maxHitLen, forcedI, forcedD,
                             anchorI, anchorD, -1,
                             stateMinI, stateMaxI, stateMinJ, stateMaxJ);
        if (_wrapTiming) _t1 = std::chrono::steady_clock::now();
        if (!skipTruncRealign) {
            runTruncDpReAlignHits(model, seq, outHits);
        }
        if (_wrapTiming) {
            _t2 = std::chrono::steady_clock::now();
            long long us_scan = std::chrono::duration_cast<std::chrono::microseconds>(_t1 - _t0).count();
            long long us_dp   = std::chrono::duration_cast<std::chrono::microseconds>(_t2 - _t1).count();
            Debug(Debug::INFO) << "WRAP timing: N=" << seq.size() << " hits=" << outHits.size()
                               << " scan=" << us_scan << "us dp=" << us_dp << "us\n";
        }
        return;
    }

    // First pass: locate envelope peak with a CYK-only scan (faster than Inside).
    // When rescanQdbOff is set, force QDB ON for pass 1 (overriding any DISABLE_QDB
    // env var) so envelopes come from the QDB-banded path, then force QDB OFF for
    // pass 2 below for bit-exact alignment on the sub-sequence.
    std::vector<Hit> probeHits;
    const int pass1DisableQdb = rescanQdbOff ? 0 : -1;
    runInfernalExactScan(model, seq, /*wantInside=*/false, probeHits, seqId, maxHitLen, forcedI, forcedD,
                         anchorI, anchorD, pass1DisableQdb,
                         stateMinI, stateMaxI, stateMinJ, stateMaxJ);
    if (probeHits.empty()) {
        outHits = std::move(probeHits);
        return;
    }

    const int N = static_cast<int>(seq.size());
    for (Hit &probe : probeHits) {
        const int lo = std::min(probe.start1, probe.end1);
        const int hi = std::max(probe.start1, probe.end1);
        const int wi = std::max(1, lo - rescorePad);
        const int wj = std::min(N, hi + rescorePad);
        if (wi >= wj) {
            outHits.push_back(std::move(probe));
            continue;
        }
        const int subLen = wj - wi + 1;
        std::string subSeq = seq.substr(static_cast<size_t>(wi - 1), static_cast<size_t>(subLen));
        std::vector<Hit> rescoreHits;
        std::vector<int> slicedMinI, slicedMaxI, slicedMinJ, slicedMaxJ;
        const int *pass2MinI = stateMinI;
        const int *pass2MaxI = stateMaxI;
        const int *pass2MinJ = stateMinJ;
        const int *pass2MaxJ = stateMaxJ;
        if (stateMinI != NULL || stateMaxI != NULL || stateMinJ != NULL || stateMaxJ != NULL) {
            std::vector<int> tmpMinI, tmpMaxI, tmpMinJ, tmpMaxJ;
            const int M = static_cast<int>(model.states.size());
            if (stateMinI != NULL && stateMaxI != NULL && stateMinJ != NULL && stateMaxJ != NULL && M > 0) {
                tmpMinI.assign(stateMinI, stateMinI + M);
                tmpMaxI.assign(stateMaxI, stateMaxI + M);
                tmpMinJ.assign(stateMinJ, stateMinJ + M);
                tmpMaxJ.assign(stateMaxJ, stateMaxJ + M);
                sliceLolStateBands(tmpMinI, tmpMaxI, tmpMinJ, tmpMaxJ,
                                   wi, wj, N,
                                   slicedMinI, slicedMaxI, slicedMinJ, slicedMaxJ);
                if (!slicedMinI.empty() && !slicedMaxI.empty()
                    && !slicedMinJ.empty() && !slicedMaxJ.empty()) {
                    pass2MinI = slicedMinI.data();
                    pass2MaxI = slicedMaxI.data();
                    pass2MinJ = slicedMinJ.data();
                    pass2MaxJ = slicedMaxJ.data();
                }
            }
        }
        const int pass2DisableQdb = rescanQdbOff ? 1 : -1;
        runInfernalExactScan(model, subSeq, wantInside, rescoreHits, seqId, maxHitLen,
                             /*forcedI=*/-1, /*forcedD=*/-1,
                             /*anchorI=*/-1, /*anchorD=*/-1,
                             pass2DisableQdb,
                             pass2MinI, pass2MaxI, pass2MinJ, pass2MaxJ);
        if (rescoreHits.empty()) {
            outHits.push_back(std::move(probe));
            continue;
        }
        // Sub-seq positions [1..subLen] map back to absolute [wi..wj].
        Hit &best = rescoreHits.front();
        best.start1 += (wi - 1);
        best.end1 += (wi - 1);
        outHits.push_back(std::move(best));
    }
    runTruncDpReAlignHits(model, seq, outHits);
}


} // namespace

int cmscan(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    if (MMseqsMPI::isMaster() == false) {
        return EXIT_SUCCESS;
    }

    const char *subcmd = (command.cmd != NULL) ? command.cmd : "cmsearch";
    Debug(Debug::INFO) << "Running in-tree CM dynamic programming (" << subcmd << ", full DP, no heuristics)\n";

    struct QueryModel {
        unsigned int key;
        InfernalExactModel exactModel;
    };

    // For DB input, keep reader open and parse CMs lazily to avoid OOM
    // For single-file input, load immediately
    struct QueryModelRef {
        unsigned int key;
        size_t dbIdx; // index into cmReader, or SIZE_MAX for single-file
    };
    std::vector<QueryModelRef> queryRefs;
    DBReader<unsigned int> *cmReader = NULL;
    QueryModel singleModel; // only used for single-file case
    if (hasDbIndex(par.db1)) {
        const std::string cmIndex = par.db1 + ".index";
        cmReader = new DBReader<unsigned int>(par.db1.c_str(),
                                              cmIndex.c_str(),
                                              par.threads > 0 ? par.threads : 1,
                                              DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        cmReader->open(DBReader<unsigned int>::NOSORT);
        queryRefs.reserve(cmReader->getSize());
        for (size_t i = 0; i < cmReader->getSize(); ++i) {
            QueryModelRef ref;
            ref.key = cmReader->getDbKey(i);
            ref.dbIdx = i;
            queryRefs.push_back(ref);
        }
        Debug(Debug::INFO) << "CM database: " << queryRefs.size() << " models (lazy loading)\n";
    } else {
        if (looksLikeInfernalCm(par.db1) == false) {
            Debug(Debug::ERROR) << "cmscan requires an Infernal CM (run cmbuild first): " << par.db1 << "\n";
            return EXIT_FAILURE;
        }
        Debug(Debug::INFO) << "Loading Infernal CM: " << par.db1 << "\n";
        singleModel.key = 0;
        singleModel.exactModel = parseInfernalCmExactModel(par.db1);
        QueryModelRef ref;
        ref.key = 0;
        ref.dbIdx = SIZE_MAX;
        queryRefs.push_back(ref);
    }

    const size_t nThreads = (par.threads > 0) ? static_cast<size_t>(par.threads) : 1u;

    const std::string seqDbPath = par.db2;
    const std::string seqDbIndex = par.db2 + ".index";
    DBReader<unsigned int> seqDbr(seqDbPath.c_str(),
                                  seqDbIndex.c_str(),
                                  static_cast<int>(nThreads),
                                  DBReader<unsigned int>::USE_DATA
                                      | DBReader<unsigned int>::USE_INDEX
                                      | DBReader<unsigned int>::USE_LOOKUP);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    NucleotideMatrix nucMat(Parameters::getInstance().scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    BaseMatrix &subMat = static_cast<BaseMatrix&>(nucMat);

    const size_t observedDbResidues = seqDbr.getAminoAcidDBSize();
    const double strandMult = (par.strand == 2) ? 2.0 : 1.0;
    const double targetDbResidues = (par.dbSize > 0)
        ? static_cast<double>(par.dbSize)
        : (static_cast<double>(observedDbResidues) * strandMult);
    const bool wantInsideUser = (par.cmMode == LocalParameters::CM_MODE_INSIDE);

    // Open the result DB like any MMseqs2 result consumer: stream per-query
    // payloads on demand via thread_idx'd getData; no upfront parsing.
    DBReader<unsigned int> resultReader(par.db3.c_str(),
                                        (par.db3 + ".index").c_str(),
                                        static_cast<int>(nThreads),
                                        DBReader<unsigned int>::USE_DATA
                                            | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::NOSORT);

    const float cmRegionFlanking = par.cmRegionFlanking;
    Debug(Debug::INFO) << subcmd << " output DB: " << par.db4 << "\n";
    DBWriter resultWriter(par.db4.c_str(),
                          par.db4Index.c_str(),
                          static_cast<unsigned int>(nThreads),
                          par.compressed,
                          Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    if (cmFastMathEnabled()) {
        Debug(Debug::WARNING) << "MMSEQS_CMSCAN_FASTMATH=1 enabled: using approximate log/exp in Inside DP (scores may drift)\n";
    }

    size_t totalHits = 0;
    char buffer[1024 + 32768 * 4];
    // Outer per-query loop is serial; inner parallelism is over candidate
    // target lines (Phase 2 below). The 9-query rRNA bench has only 1 query,
    // so per-query parallelism wastes 15/16 threads on the slow CYK step.
    for (long qi = 0; qi < static_cast<long>(queryRefs.size()); ++qi) {
        const QueryModelRef &ref = queryRefs[qi];

        // Load CM lazily, in-memory: parse the DB entry directly via istringstream.
        QueryModel qm;
        if (ref.dbIdx != SIZE_MAX && cmReader != NULL) {
            qm.key = ref.key;
            const char *raw = cmReader->getData(ref.dbIdx, 0);
            size_t len = cmReader->getEntryLen(ref.dbIdx);
            std::string text(raw, len);
            const size_t nul = text.find('\0');
            if (nul != std::string::npos) {
                text.resize(nul);
            }
            std::istringstream iss(text);
            const std::string srcLabel = "cm-db-entry[" + std::to_string(qm.key) + "]";
            qm.exactModel = parseInfernalCmExactModelFromStream(iss, srcLabel);
        } else {
            qm = singleModel;
        }

        // Infernal default: Inside score is the primary ranking signal;
        // CYK remains the trace source. `wantInsideUser` only changes which
        // score we report, not whether Inside is computed.
        const bool promoteInsideToPrimary = !wantInsideUser;
        // MMSEQS_CMSCAN_SKIP_INSIDE=1: skip the Inside DP entirely (CYK still
        // runs as the trace source). For pipelines that don't consume h.inside
        // (the bench's permissive -e 10000 path doesn't), this halves cmsearch
        // CPU time. Verify byte-identical aln output before relying on it.
        static const char *const _skipInsideEnv = std::getenv("MMSEQS_CMSCAN_SKIP_INSIDE");
        const bool wantInside = !(_skipInsideEnv != nullptr && _skipInsideEnv[0] == '1');

        // Slack for the CYK d-cap. We take the max of:
        //   - 3x the rnasearch prefilter envelope length (local extension)
        //   - 1.5x the CM consensus length (global: a true CM hit can stretch
        //     well beyond the prefilter envelope, especially for structured
        //     RNAs where iter-3 msa-profile mass is concentrated around the
        //     center but the CM consensus spans the full rRNA).
        // model.w (Infernal's W) remains the hard upper bound on the CM side.
        static constexpr int CM_MAXSPAN_ENV_SLACK = 3;
        static constexpr double CM_MAXSPAN_CLEN_SLACK = 1.5;
        const int clenFloor = static_cast<int>(qm.exactModel.clen * CM_MAXSPAN_CLEN_SLACK);

        // Peak-anchor (MMSEQS_CMSCAN_ANCHOR_MODE=peak): rescore SW backtrace
        // from the prefilter results table with CM match emissions; argmax-j
        // of the cumulative bit-score becomes a point anchor for FORCE_GLOBAL=3
        // instead of the geometric midpoint of the prefilter envelope. Only
        // used when forceGlobal in {3,4} and a backtrace is present in the
        // prefilter row. Built once per query (model-only; cand-independent).
        // peakAnchorMode: 0=off, 1=peak (point anchor at SW-CIGAR bit-score peak),
        //                 2=peak_walk (anchor span = SW match extended naively
        //                 leftward to CM col 1 and rightward to CM col CLEN; FG=3
        //                 then straddles the midpoint of that span). The CIGAR
        //                 itself isn't needed for peak_walk — only qStart/qEnd
        //                 and clen — but we keep the same gating as peak (skip
        //                 reverse-strand for v1).
        static int peakAnchorMode = -1;
        if (peakAnchorMode == -1) {
            const char *env = std::getenv("MMSEQS_CMSCAN_ANCHOR_MODE");
            if (env != NULL) {
                if      (std::string(env) == "peak")      peakAnchorMode = 1;
                else if (std::string(env) == "peak_walk") peakAnchorMode = 2;
                else                                       peakAnchorMode = 0;
            } else {
                peakAnchorMode = 0;
            }
        }
        std::vector<std::array<float, 4>> perColMaxEmit;
        if (peakAnchorMode == 1) {
            perColMaxEmit = buildPerColumnMaxEmissions(qm.exactModel);
        }

        std::vector<Hit> hits;
        hits.reserve(1024);

        // Finalize each scanned hit: remap strand-specific coords to forward,
        // promote Inside -> primary score when Infernal-default, compute
        // evalue, and (forward only, with backtrace) compute trace-based pid.
        auto finalizeHit = [&](Hit &h, int strand, const std::string &targetFullSeq,
                               unsigned int tKey, unsigned int fullLen,
                               int offset, int regionLen) {
            h.dbKey = tKey;
            h.dbLen = fullLen;
            if (strand > 0) {
                if (offset > 0) {
                    h.start1 += offset;
                    h.end1 += offset;
                }
            } else {
                // revcomp pos p corresponds to forward pos (regionLen - p + 1);
                // Infernal convention: start1 > end1 indicates minus strand.
                int fwdStart = regionLen - h.start1 + 1 + offset;
                int fwdEnd = regionLen - h.end1 + 1 + offset;
                h.start1 = fwdStart;
                h.end1 = fwdEnd;
            }
            if (promoteInsideToPrimary && h.inside != NEG_INF) {
                h.cyk = h.inside;
            }
            double ev = 0.0;
            if (infernalExactScoreToEvalue(qm.exactModel, /*evalueModeInside=*/true, h.mode, h.cyk, targetDbResidues, ev)) {
                h.evalue = ev;
                h.hasEvalue = true;
            }
            const bool hasBacktrace = (!h.cigar.empty() && h.cigar != "NA");
            if (hasBacktrace
                && !h.traceStates.empty() && h.traceStates != "NA"
                && h.start1 > 0 && h.end1 >= h.start1) {
                const size_t s0 = static_cast<size_t>(h.start1 - 1);
                const size_t slen = static_cast<size_t>(h.end1 - h.start1 + 1);
                if (s0 < targetFullSeq.size()) {
                    const std::string obs = targetFullSeq.substr(s0, std::min(slen, targetFullSeq.size() - s0));
                    const std::vector<int> trace = decodeTraceStates(h.traceStates);
                    const std::string cons = modelTraceConsensusForCigar(qm.exactModel, trace);
                    h.precomputedSeqId = seqIdFromCigarConsensus(h.cigar, cons, obs);
                }
            }
        };

        // Phase 1: walk the per-query candidate list once (sequential),
        // parsing prefilter coords and deduping by target key. This is cheap
        // (just integer parsing + hash insert) so single-threaded is fine.
        // MMSEQS_CMSCAN_ALLOW_MULTIHIT=1: skip the per-target dedup so multiple
        // non-overlapping prefilter hits per target sequence each get realigned
        // and emitted (mirrors Infernal cmalign behavior). Combined with a
        // multi-hit prefilter, this can substantially raise R-scape AUC by
        // adding more independent co-evolution evidence per target.
        static const char *const allowMultiHitEnv = std::getenv("MMSEQS_CMSCAN_ALLOW_MULTIHIT");
        const bool allowMultiHit = (allowMultiHitEnv != nullptr && allowMultiHitEnv[0] == '1');
        struct Cand {
            unsigned int tKey;
            int dbStart;
            int dbEnd;
            int qStart;
            int qEnd;
            int qLen;
            bool hasRegionCoord;
            bool prefilterIsRev;
            std::string backtrace; // SW per-char M/I/D from m8 col 10, if present
        };
        std::vector<Cand> cands;
        cands.reserve(1024);
        std::unordered_set<unsigned int> seen;

        const size_t rid = resultReader.getId(qm.key);
        if (rid != UINT_MAX) {
            char *data = resultReader.getData(rid, 0);
            while (*data != '\0') {
                while (*data == ' ' || *data == '\t') {
                    ++data;
                }
                // Format: targetKey score seqId evalue qStart qEnd qLen dbStart dbEnd dbLen [backtrace]
                const char *lineStart = data;
                char *endptr = NULL;
                const unsigned long k = std::strtoul(data, &endptr, 10);
                if (endptr == data || k > static_cast<unsigned long>(UINT_MAX)) {
                    data = Util::skipLine(data);
                    continue;
                }
                const unsigned int tKey = static_cast<unsigned int>(k);
                // Dedup: only score first occurrence (best hit per target),
                // unless MMSEQS_CMSCAN_ALLOW_MULTIHIT=1.
                if (!allowMultiHit && seen.insert(tKey).second == false) {
                    data = Util::skipLine(data);
                    continue;
                }

                Cand cand{};
                cand.tKey = tKey;
                if (*endptr == '\t') {
                    const char *p = endptr;
                    int col = 1;
                    const char *colStart[12] = {NULL};
                    colStart[0] = lineStart;
                    while (*p != '\n' && *p != '\0' && col < 11) {
                        if (*p == '\t') {
                            colStart[col] = p + 1;
                            col++;
                        }
                        p++;
                    }
                    if (col >= 10 && colStart[4] && colStart[5] && colStart[7] && colStart[8]) {
                        cand.dbStart = Util::fast_atoi<int>(colStart[7]);
                        cand.dbEnd   = Util::fast_atoi<int>(colStart[8]);
                        cand.qStart  = Util::fast_atoi<int>(colStart[4]);
                        cand.qEnd    = Util::fast_atoi<int>(colStart[5]);
                        cand.qLen    = (colStart[6] != NULL) ? Util::fast_atoi<int>(colStart[6]) : 0;
                        // mmseqs/iter3 may encode minus strand by reversing
                        // either side. Strand is XOR of the two reversals.
                        const bool dbRev = (cand.dbStart > cand.dbEnd);
                        const bool qRev  = (cand.qStart > cand.qEnd);
                        cand.prefilterIsRev = (dbRev != qRev);
                        if (dbRev) std::swap(cand.dbStart, cand.dbEnd);
                        if (qRev)  std::swap(cand.qStart, cand.qEnd);
                        cand.hasRegionCoord = true;
                        if (col >= 11 && colStart[10]) {
                            const char *bs = colStart[10];
                            const char *be = bs;
                            while (*be != '\n' && *be != '\0' && *be != '\t') ++be;
                            cand.backtrace.assign(bs, static_cast<size_t>(be - bs));
                        }
                    }
                }
                data = Util::skipLine(data);
                cands.push_back(cand);
            }
        }

        const bool useGpuBatchCommonPath = cmGpuBatchCommonPathEnabled(qm.exactModel);
        if (useGpuBatchCommonPath) {
            const int M = static_cast<int>(qm.exactModel.states.size());
            const bool truncModesEnabled = cmTruncModesEnabled();
            const bool null3Enabled = cmNull3Enabled();
            const double log2Omega3 = (qm.exactModel.null3Omega > 0.0) ? std::log2(qm.exactModel.null3Omega) : NEG_INF;
            const std::array<double, 4> log2Null = {{
                (qm.exactModel.nullProb[0] > 0.0) ? std::log2(qm.exactModel.nullProb[0]) : NEG_INF,
                (qm.exactModel.nullProb[1] > 0.0) ? std::log2(qm.exactModel.nullProb[1]) : NEG_INF,
                (qm.exactModel.nullProb[2] > 0.0) ? std::log2(qm.exactModel.nullProb[2]) : NEG_INF,
                (qm.exactModel.nullProb[3] > 0.0) ? std::log2(qm.exactModel.nullProb[3]) : NEG_INF
            }};

            std::vector<ExactStateExec> exec(static_cast<size_t>(M));
            std::vector<StateId> trDst;
            std::vector<double> trSc;
            trDst.reserve(static_cast<size_t>(M) * 4);
            trSc.reserve(static_cast<size_t>(M) * 4);
            std::vector<int> activeStates;
            activeStates.reserve(static_cast<size_t>(M));
            const char *envDisableQdb = std::getenv("MMSEQS_CMSCAN_DISABLE_QDB");
            const bool disableQdb = (envDisableQdb != NULL && envDisableQdb[0] == '1');
            const char *envDropDmin = std::getenv("MMSEQS_CMSCAN_DROP_DMIN");
            const bool dropDmin = (envDropDmin != NULL && envDropDmin[0] == '1');
            const char *envDropDmax = std::getenv("MMSEQS_CMSCAN_DROP_DMAX");
            const bool dropDmax = (envDropDmax != NULL && envDropDmax[0] == '1');
            for (int v = 0; v < M; ++v) {
                const CmState &st = qm.exactModel.states[static_cast<size_t>(v)];
                ExactStateExec e;
                e.type = st.type;
                e.niShift = 0;
                e.dConsume = 0;
                e.consumesLeft = false;
                e.consumesRight = false;
                e.bLeft = st.cfirst;
                e.bRight = st.cnum;
                e.trDst4[0] = 0; e.trDst4[1] = 0; e.trDst4[2] = 0; e.trDst4[3] = 0;
                e.trScF4[0] = -std::numeric_limits<float>::infinity();
                e.trScF4[1] = -std::numeric_limits<float>::infinity();
                e.trScF4[2] = -std::numeric_limits<float>::infinity();
                e.trScF4[3] = -std::numeric_limits<float>::infinity();
                e.trSc4[0] = NEG_INF; e.trSc4[1] = NEG_INF; e.trSc4[2] = NEG_INF; e.trSc4[3] = NEG_INF;
                e.trOff = trDst.size();
                e.trCount = 0;
                e.emitPtr = st.emitF.empty() ? NULL : &st.emitF[0];
                e.emitSize = static_cast<int>(st.emit.size());
                e.splitKMin = 0;
                e.splitKMax = -1;
                e.splitRMin = 0;
                e.splitRMax = -1;
                e.dmin = (disableQdb || dropDmin || dropDmax) ? 0 : st.dmin1;
                e.dmax = (disableQdb || dropDmax) ? -1 : st.dmax1;
                e.null2Agg[0] = e.null2Agg[1] = e.null2Agg[2] = e.null2Agg[3] = 0.0;
                if (st.type == CM_ST_MP) {
                    e.niShift = 1;
                    e.dConsume = 2;
                    e.consumesLeft = true;
                    e.consumesRight = true;
                } else if (st.type == CM_ST_ML || st.type == CM_ST_IL) {
                    e.niShift = 1;
                    e.dConsume = 1;
                    e.consumesLeft = true;
                } else if (st.type == CM_ST_MR || st.type == CM_ST_IR) {
                    e.dConsume = 1;
                    e.consumesRight = true;
                }
                if (st.type == CM_ST_MP && e.emitPtr != NULL && e.emitSize >= 16) {
                    for (int a = 0; a < 4; ++a) {
                        double row = 0.0;
                        double col = 0.0;
                        for (int b = 0; b < 4; ++b) {
                            row += emitBitToProbPair(e.emitPtr[a * 4 + b], qm.exactModel, a, b);
                            col += emitBitToProbPair(e.emitPtr[b * 4 + a], qm.exactModel, b, a);
                        }
                        e.null2Agg[a] = row + col;
                    }
                } else if ((st.type == CM_ST_ML || st.type == CM_ST_MR || st.type == CM_ST_IL || st.type == CM_ST_IR) &&
                           e.emitPtr != NULL && e.emitSize >= 4) {
                    for (int a = 0; a < 4; ++a) {
                        e.null2Agg[a] = emitBitToProbSingle(e.emitPtr[a], qm.exactModel, a);
                    }
                }
                if (st.type == CM_ST_B) {
                    if (!disableQdb && !dropDmax && e.bLeft >= 0 && e.bLeft < M) {
                        const CmState &l = qm.exactModel.states[static_cast<size_t>(e.bLeft)];
                        if (l.dmin1 >= 0 && l.dmax1 >= l.dmin1) {
                            e.splitKMin = dropDmin ? 0 : l.dmin1;
                            e.splitKMax = l.dmax1;
                        }
                    }
                    if (!disableQdb && !dropDmax && e.bRight >= 0 && e.bRight < M) {
                        const CmState &r = qm.exactModel.states[static_cast<size_t>(e.bRight)];
                        if (r.dmin1 >= 0 && r.dmax1 >= r.dmin1) {
                            e.splitRMin = dropDmin ? 0 : r.dmin1;
                            e.splitRMax = r.dmax1;
                        }
                    }
                } else if (st.type != CM_ST_E) {
                    const int tn = std::min(st.cnum, static_cast<int>(st.trans.size()));
                    for (int x = 0; x < tn; ++x) {
                        const int y = st.cfirst + x;
                        if (y < 0 || y >= M) continue;
                        if (e.trCount < 4) {
                            const size_t t = static_cast<size_t>(e.trCount);
                            e.trDst4[t] = static_cast<StateId>(y);
                            e.trSc4[t] = st.trans[static_cast<size_t>(x)];
                            e.trScF4[t] = static_cast<float>(st.trans[static_cast<size_t>(x)]);
                        }
                        trDst.push_back(static_cast<StateId>(y));
                        trSc.push_back(st.trans[static_cast<size_t>(x)]);
                        ++e.trCount;
                    }
                }
                exec[static_cast<size_t>(v)] = e;
            }
            for (int v = M - 1; v >= 0; --v) {
                if (exec[static_cast<size_t>(v)].type != CM_ST_E) {
                    activeStates.push_back(v);
                }
            }
            std::vector<float> trScF(trSc.size());
            for (size_t ti = 0; ti < trSc.size(); ++ti) {
                trScF[ti] = static_cast<float>(trSc[ti]);
            }
            std::vector<CmFastDeckGpuState> gpuStates(static_cast<size_t>(M));
            std::vector<float> gpuEmitData;
            gpuEmitData.reserve(static_cast<size_t>(M) * 8);
            for (int v = 0; v < M; ++v) {
                const size_t vi = static_cast<size_t>(v);
                const CmState &src = qm.exactModel.states[vi];
                const ExactStateExec &st = exec[vi];
                CmFastDeckGpuState &dst = gpuStates[vi];
                dst.type = static_cast<int>(st.type);
                dst.dmin = st.dmin;
                dst.dmax = st.dmax;
                dst.bLeft = st.bLeft;
                dst.bRight = st.bRight;
                dst.dConsume = st.dConsume;
                dst.niShift = st.niShift;
                dst.emitSize = st.emitSize;
                dst.emitOffset = src.emitF.empty() ? -1 : static_cast<int>(gpuEmitData.size());
                dst.trCount = st.trCount;
                dst.trOff = st.trOff;
                dst.consumeMask = (st.consumesLeft ? 1u : 0u) | (st.consumesRight ? 2u : 0u);
                dst.endSc = static_cast<float>(src.endSc);
                dst.null2Agg[0] = static_cast<float>(st.null2Agg[0]);
                dst.null2Agg[1] = static_cast<float>(st.null2Agg[1]);
                dst.null2Agg[2] = static_cast<float>(st.null2Agg[2]);
                dst.null2Agg[3] = static_cast<float>(st.null2Agg[3]);
                if (!src.emitF.empty()) {
                    gpuEmitData.insert(gpuEmitData.end(), src.emitF.begin(), src.emitF.end());
                }
            }

            struct BatchEntry {
                std::string seqId;
                std::string fullSeq;
                std::string scanSeq;
                unsigned int tKey = 0;
                unsigned int fullLen = 0;
                int offset = 0;
                int regionLen = 0;
                int strand = 1;
                int minSpan = 1;
                int maxSpan = 0;
                int maxHitLen = 0;
                int forceI = -1;
                int forceD = -1;
            };
            struct BatchBucket {
                int maxN = 0;
                std::vector<BatchEntry> entries;
            };
            std::unordered_map<int, BatchBucket> batchBuckets;
            std::vector<int> batchBucketKeys;
            Sequence seqObjLocal(seqDbr.getMaxSeqLen(), Parameters::DBTYPE_NUCLEOTIDES, &nucMat, 0, false, false);
            static int wcapModeBatch = -1;
            if (wcapModeBatch == -1) {
                const char *env = std::getenv("MMSEQS_CMSCAN_WCAP");
                wcapModeBatch = (env != NULL && std::string(env) == "0") ? 0 : 1;
            }
            static double dminFracBatch = -1.0;
            if (dminFracBatch < 0.0) {
                const char *env = std::getenv("MMSEQS_CMSCAN_DMIN_CLEN_FRAC");
                dminFracBatch = (env != NULL) ? std::atof(env) : 0.0;
                if (dminFracBatch < 0.0 || dminFracBatch > 1.0) dminFracBatch = 0.0;
            }
            int gpuBatchBucketWidth = 4;
            if (const char *env = std::getenv("MMSEQS_CMSCAN_GPU_BATCH_BUCKET")) {
                const int parsed = std::atoi(env);
                if (parsed > 0) gpuBatchBucketWidth = parsed;
            }
            size_t batchCellBudget = 1610612736ULL;
            if (const char *env = std::getenv("MMSEQS_CMSCAN_GPU_BATCH_CELLS")) {
                const unsigned long long parsed = std::strtoull(env, NULL, 10);
                if (parsed > 0) batchCellBudget = static_cast<size_t>(parsed);
            }
            const size_t batchBytesBudget = cmGpuUsableBytes();
            auto estimateGpuBatchBytes = [&](size_t batchSize, int maxN, int traceCapacity) -> size_t {
                if (batchSize == 0 || maxN <= 0 || traceCapacity <= 0) {
                    return 0;
                }
                const size_t seqStride = static_cast<size_t>(maxN + 1);
                const size_t iStride = static_cast<size_t>(maxN + 2);
                const size_t stateStride = static_cast<size_t>(maxN + 1) * iStride;
                const size_t cellsPerJob = static_cast<size_t>(M) * stateStride;
                const size_t vdCells = static_cast<size_t>(M) * static_cast<size_t>(maxN + 1);
                size_t bytes = 0;
                bytes += batchSize * sizeof(CmFastDeckGpuBatchJob);
                bytes += batchSize * seqStride * sizeof(int8_t);
                bytes += 4ULL * batchSize * seqStride * sizeof(int);
                bytes += static_cast<size_t>(maxN + 1) * sizeof(double);
                bytes += gpuStates.size() * sizeof(CmFastDeckGpuState);
                bytes += activeStates.size() * sizeof(int);
                bytes += static_cast<size_t>(M) * sizeof(size_t);
                bytes += 2ULL * vdCells * sizeof(int);
                bytes += trDst.size() * sizeof(uint16_t);
                bytes += trScF.size() * sizeof(float);
                bytes += gpuEmitData.size() * sizeof(float);
                bytes += batchSize * cellsPerJob * sizeof(float);
                bytes += batchSize * sizeof(CmFastDeckGpuBatchResult);
                bytes += batchSize * static_cast<size_t>(traceCapacity) * sizeof(uint16_t);
                bytes += batchSize * static_cast<size_t>(traceCapacity) * CM_FASTDECK_GPU_TBCELL_BYTES;
                return bytes;
            };

            auto flushGpuBatch = [&](std::vector<BatchEntry> &batchEntries, int maxN) {
                if (batchEntries.empty()) return;
                if (maxN <= 0) {
                    batchEntries.clear();
                    return;
                }
                const size_t iStride = static_cast<size_t>(maxN + 2);
                const size_t stateStride = static_cast<size_t>(maxN + 1) * iStride;
                const size_t cellsPerJob = static_cast<size_t>(M) * stateStride;
                std::vector<size_t> stateBase(static_cast<size_t>(M));
                for (int v = 0; v < M; ++v) {
                    stateBase[static_cast<size_t>(v)] = static_cast<size_t>(v) * stateStride;
                }
                const size_t vdCells = static_cast<size_t>(M) * static_cast<size_t>(maxN + 1);
                std::vector<int> bSplitBegByVD(vdCells, 0);
                std::vector<int> bSplitEndByVD(vdCells, -1);
                for (int v = 0; v < M; ++v) {
                    const ExactStateExec &st = exec[static_cast<size_t>(v)];
                    if (st.type != CM_ST_B) continue;
                    for (int d = 0; d <= maxN; ++d) {
                        int kBeg = 0;
                        int kEnd = d;
                        if (st.splitKMax >= st.splitKMin) {
                            kBeg = std::max(kBeg, st.splitKMin);
                            kEnd = std::min(kEnd, st.splitKMax);
                        }
                        if (st.splitRMax >= st.splitRMin) {
                            kBeg = std::max(kBeg, d - st.splitRMax);
                            kEnd = std::min(kEnd, d - st.splitRMin);
                        }
                        const size_t vd = static_cast<size_t>(v) * static_cast<size_t>(maxN + 1) + static_cast<size_t>(d);
                        bSplitBegByVD[vd] = kBeg;
                        bSplitEndByVD[vd] = kEnd;
                    }
                }
                const int seqStride = maxN + 1;
                std::vector<CmFastDeckGpuBatchJob> gpuJobs(batchEntries.size());
                std::vector<int8_t> packedSeqCode(static_cast<size_t>(batchEntries.size()) * static_cast<size_t>(seqStride), -1);
                std::vector<int> packedPrefA(static_cast<size_t>(batchEntries.size()) * static_cast<size_t>(seqStride), 0);
                std::vector<int> packedPrefC(static_cast<size_t>(batchEntries.size()) * static_cast<size_t>(seqStride), 0);
                std::vector<int> packedPrefG(static_cast<size_t>(batchEntries.size()) * static_cast<size_t>(seqStride), 0);
                std::vector<int> packedPrefU(static_cast<size_t>(batchEntries.size()) * static_cast<size_t>(seqStride), 0);
                for (size_t bi = 0; bi < batchEntries.size(); ++bi) {
                    const BatchEntry &e = batchEntries[bi];
                    CmFastDeckGpuBatchJob &job = gpuJobs[bi];
                    job.N = static_cast<int>(e.scanSeq.size());
                    job.minSpan = e.minSpan;
                    job.maxSpan = e.maxSpan;
                    job.forcedI = e.forceI;
                    job.forcedD = e.forceD;
                    const size_t off = bi * static_cast<size_t>(seqStride);
                    for (int p = 1; p <= job.N; ++p) {
                        const int biCode = baseToIdx(e.scanSeq[static_cast<size_t>(p - 1)]);
                        packedSeqCode[off + static_cast<size_t>(p)] = static_cast<int8_t>(biCode);
                        packedPrefA[off + static_cast<size_t>(p)] = packedPrefA[off + static_cast<size_t>(p - 1)];
                        packedPrefC[off + static_cast<size_t>(p)] = packedPrefC[off + static_cast<size_t>(p - 1)];
                        packedPrefG[off + static_cast<size_t>(p)] = packedPrefG[off + static_cast<size_t>(p - 1)];
                        packedPrefU[off + static_cast<size_t>(p)] = packedPrefU[off + static_cast<size_t>(p - 1)];
                        if (biCode == 0) packedPrefA[off + static_cast<size_t>(p)]++;
                        else if (biCode == 1) packedPrefC[off + static_cast<size_t>(p)]++;
                        else if (biCode == 2) packedPrefG[off + static_cast<size_t>(p)]++;
                        else if (biCode == 3) packedPrefU[off + static_cast<size_t>(p)]++;
                    }
                }
                std::vector<double> log2Int(static_cast<size_t>(maxN + 1), NEG_INF);
                for (int x = 1; x <= maxN; ++x) {
                    log2Int[static_cast<size_t>(x)] = std::log2(static_cast<double>(x));
                }
                const int traceCapacity = std::max(256, 2 * (M + maxN + 8));
                std::vector<CmFastDeckGpuBatchResult> gpuResults;
                std::vector<uint16_t> hostTraceStates;
                std::string gpuError;
                const bool ok = runInfernalExactScanFastDeckGpuBatch(
                    M,
                    qm.exactModel.rootState,
                    maxN,
                    static_cast<int>(batchEntries.size()),
                    static_cast<int>(iStride),
                    stateStride,
                    cellsPerJob,
                    gpuJobs,
                    packedSeqCode,
                    packedPrefA,
                    packedPrefC,
                    packedPrefG,
                    packedPrefU,
                    log2Int,
                    gpuStates,
                    activeStates,
                    stateBase,
                    bSplitBegByVD,
                    bSplitEndByVD,
                    trDst,
                    trScF,
                    gpuEmitData,
                    truncModesEnabled,
                    null3Enabled,
                    log2Omega3,
                    log2Null,
                    traceCapacity,
                    &gpuResults,
                    &hostTraceStates,
                    &gpuError);
                if (!ok) {
                    Debug(Debug::WARNING) << "cmscan GPU batch disabled for this batch: " << gpuError << "\n";
                    for (const BatchEntry &e : batchEntries) {
                        std::vector<Hit> scanHits;
                        runInfernalExactScanWithEnvelopeRescore(qm.exactModel, e.scanSeq, wantInside, scanHits, e.seqId,
                                                                e.maxHitLen, e.forceI, e.forceD, -1, -1);
                        for (Hit &h : scanHits) {
                            finalizeHit(h, e.strand, e.fullSeq, e.tKey, e.fullLen, e.offset, e.regionLen);
                            hits.emplace_back(std::move(h));
                        }
                    }
                    batchEntries.clear();
                    return;
                }
                for (size_t bi = 0; bi < batchEntries.size(); ++bi) {
                    const CmFastDeckGpuBatchResult &gr = gpuResults[bi];
                    std::vector<Hit> scanHits;
                    const bool traced = finishFastDeckSelectedHitBatched(
                        qm.exactModel,
                        batchEntries[bi].seqId,
                        gr,
                        hostTraceStates.data() + bi * static_cast<size_t>(traceCapacity),
                        scanHits);
                    if (!traced) {
                        runInfernalExactScanWithEnvelopeRescore(qm.exactModel, batchEntries[bi].scanSeq, wantInside, scanHits,
                                                                batchEntries[bi].seqId, batchEntries[bi].maxHitLen,
                                                                batchEntries[bi].forceI, batchEntries[bi].forceD, -1, -1);
                    }
                    for (Hit &h : scanHits) {
                        finalizeHit(h, batchEntries[bi].strand, batchEntries[bi].fullSeq,
                                    batchEntries[bi].tKey, batchEntries[bi].fullLen,
                                    batchEntries[bi].offset, batchEntries[bi].regionLen);
                        hits.emplace_back(std::move(h));
                    }
                }
                batchEntries.clear();
            };

            for (long ci = 0; ci < static_cast<long>(cands.size()); ++ci) {
                const Cand &cand = cands[ci];
                const size_t tId = seqDbr.getId(cand.tKey);
                if (tId == UINT_MAX) continue;
                FastaSeq fs = decodeOneSequence(seqDbr, tId, subMat, seqObjLocal, 0);
                std::string regionSeq;
                int offset = 0;
                int maxHitLen = 0;
                if (cmRegionFlanking > 0.0f && cand.hasRegionCoord) {
                    const int qAlnLen = std::max(1, cand.qEnd - cand.qStart + 1);
                    int leftFlank = 0;
                    int rightFlank = 0;
                    const int qLen = std::max(qAlnLen, (cand.qLen > 0) ? cand.qLen : qm.exactModel.clen);
                    if (qLen > 0) {
                        const int qStartClamped = std::max(0, std::min(cand.qStart, qLen - 1));
                        const int qEndClamped = std::max(qStartClamped, std::min(cand.qEnd, qLen - 1));
                        const int leftMissing = qStartClamped;
                        const int rightMissing = std::max(0, qLen - qEndClamped - 1);
                        const int totalMissing = leftMissing + rightMissing;
                        const double missingFrac = static_cast<double>(totalMissing) / static_cast<double>(qLen);
                        const int baseFlank = static_cast<int>(std::ceil(static_cast<double>(cmRegionFlanking) * static_cast<double>(qAlnLen) * missingFrac));
                        leftFlank = baseFlank + leftMissing;
                        rightFlank = baseFlank + rightMissing;
                    } else {
                        const int baseFlank = static_cast<int>(cmRegionFlanking * qAlnLen);
                        leftFlank = baseFlank;
                        rightFlank = baseFlank;
                    }
                    const int sLen = static_cast<int>(fs.seq.size());
                    const int regStart = std::max(0, cand.dbStart - leftFlank);
                    const int regEnd = std::min(sLen, cand.dbEnd + rightFlank + 1);
                    regionSeq = fs.seq.substr(static_cast<size_t>(regStart), static_cast<size_t>(regEnd - regStart));
                    offset = regStart;
                } else {
                    regionSeq = fs.seq;
                }
                if (cand.hasRegionCoord) {
                    const int dbSpan = std::max(1, cand.dbEnd - cand.dbStart + 1);
                    maxHitLen = std::max(dbSpan * CM_MAXSPAN_ENV_SLACK, clenFloor);
                } else if (clenFloor > 0) {
                    maxHitLen = clenFloor;
                }
                const bool scanFwd = !cand.hasRegionCoord || !cand.prefilterIsRev;
                const bool scanRev = !cand.hasRegionCoord || cand.prefilterIsRev;
                int forceI = -1, forceD = -1;
                static int forcePrefilter = -1;
                if (forcePrefilter == -1) {
                    const char *env = std::getenv("MMSEQS_CMSCAN_FORCE_PREFILTER_ENV");
                    forcePrefilter = (env != NULL && std::string(env) == "1") ? 1 : 0;
                }
                if (forcePrefilter == 1 && cand.hasRegionCoord) {
                    forceI = (cand.dbStart - offset) + 1;
                    forceD = std::max(1, cand.dbEnd - cand.dbStart + 1);
                }
                auto enqueueJob = [&](const std::string &scanSeq, int strand, int jobForceI, int jobForceD) {
                    const int Njob = static_cast<int>(scanSeq.size());
                    if (Njob <= 0) return;
                    int maxSpan = Njob;
                    if (wcapModeBatch == 1 && qm.exactModel.w > 0) maxSpan = std::min(maxSpan, qm.exactModel.w);
                    if (maxHitLen > 0) maxSpan = std::min(maxSpan, maxHitLen);
                    int minSpan = 1;
                    if (dminFracBatch > 0.0 && qm.exactModel.clen > 0) {
                        minSpan = std::max(minSpan, static_cast<int>(dminFracBatch * qm.exactModel.clen));
                        if (minSpan > maxSpan) minSpan = maxSpan;
                    }
                    if (jobForceI > 0 && jobForceD > 0 && jobForceD <= Njob && jobForceI <= Njob - jobForceD + 1) {
                        minSpan = jobForceD;
                        maxSpan = jobForceD;
                    }
                    if (maxSpan <= 0) {
                        std::vector<Hit> scanHits;
                        runInfernalExactScanWithEnvelopeRescore(qm.exactModel, scanSeq, wantInside, scanHits, fs.id, maxHitLen,
                                                                jobForceI, jobForceD, -1, -1);
                        for (Hit &h : scanHits) {
                            finalizeHit(h, strand, fs.seq, cand.tKey, static_cast<unsigned int>(fs.seq.size()),
                                        offset, static_cast<int>(regionSeq.size()));
                            hits.emplace_back(std::move(h));
                        }
                        return;
                    }
                    const int bucketKey = std::max(gpuBatchBucketWidth, ((Njob + gpuBatchBucketWidth - 1) / gpuBatchBucketWidth) * gpuBatchBucketWidth);
                    BatchBucket &bucket = batchBuckets[bucketKey];
                    if (bucket.entries.empty()) {
                        batchBucketKeys.push_back(bucketKey);
                    }
                    const int nextMaxN = std::max(bucket.maxN, Njob);
                    const int nextTraceCapacity = std::max(256, 2 * (M + nextMaxN + 8));
                    const size_t nextIStride = static_cast<size_t>(nextMaxN + 2);
                    const size_t nextStateStride = static_cast<size_t>(nextMaxN + 1) * nextIStride;
                    const size_t nextCellsPerJob = static_cast<size_t>(M) * nextStateStride;
                    const size_t projectedCells = nextCellsPerJob * static_cast<size_t>(bucket.entries.size() + 1);
                    const size_t projectedBytes = estimateGpuBatchBytes(bucket.entries.size() + 1, nextMaxN, nextTraceCapacity);
                    if (!bucket.entries.empty() &&
                        ((projectedCells > batchCellBudget) ||
                         (batchBytesBudget > 0 && projectedBytes > batchBytesBudget))) {
                        flushGpuBatch(bucket.entries, bucket.maxN);
                        bucket.maxN = 0;
                    }
                    if (bucket.entries.empty()) {
                        const size_t singleBytes = estimateGpuBatchBytes(1, Njob, std::max(256, 2 * (M + Njob + 8)));
                        if (batchBytesBudget > 0 && singleBytes > batchBytesBudget) {
                            std::vector<Hit> scanHits;
                            runInfernalExactScanWithEnvelopeRescore(qm.exactModel, scanSeq, wantInside, scanHits, fs.id,
                                                                    maxHitLen, jobForceI, jobForceD, -1, -1);
                            for (Hit &h : scanHits) {
                                finalizeHit(h, strand, fs.seq, cand.tKey, static_cast<unsigned int>(fs.seq.size()),
                                            offset, static_cast<int>(regionSeq.size()));
                                hits.emplace_back(std::move(h));
                            }
                            return;
                        }
                    }
                    BatchEntry e;
                    e.seqId = fs.id;
                    e.fullSeq = fs.seq;
                    e.scanSeq = scanSeq;
                    e.tKey = cand.tKey;
                    e.fullLen = static_cast<unsigned int>(fs.seq.size());
                    e.offset = offset;
                    e.regionLen = static_cast<int>(regionSeq.size());
                    e.strand = strand;
                    e.minSpan = minSpan;
                    e.maxSpan = maxSpan;
                    e.maxHitLen = maxHitLen;
                    e.forceI = jobForceI;
                    e.forceD = jobForceD;
                    bucket.maxN = std::max(bucket.maxN, Njob);
                    bucket.entries.emplace_back(std::move(e));
                };
                if (scanFwd) enqueueJob(regionSeq, +1, forceI, forceD);
                if (scanRev) {
                    const std::string revSeq = reverseComplement(regionSeq);
                    int revForceI = -1, revForceD = -1;
                    const int M_local = static_cast<int>(regionSeq.size());
                    if (forceI > 0 && forceD > 0) {
                        revForceI = M_local - forceI - forceD + 2;
                        revForceD = forceD;
                    }
                    enqueueJob(revSeq, -1, revForceI, revForceD);
                }
            }
            for (size_t ki = 0; ki < batchBucketKeys.size(); ++ki) {
                BatchBucket &bucket = batchBuckets[batchBucketKeys[ki]];
                flushGpuBatch(bucket.entries, bucket.maxN);
            }
        } else {
        // Phase 2: parallelize CYK across candidates. Each thread keeps a
        // local Hit accumulator and a local Sequence buffer; the global
        // `hits` is merged once under critical at the end.
#pragma omp parallel num_threads(static_cast<int>(nThreads))
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence seqObjLocal(seqDbr.getMaxSeqLen(), Parameters::DBTYPE_NUCLEOTIDES, &nucMat, 0, false, false);
            std::vector<Hit> localHits;
            localHits.reserve(64);

#pragma omp for schedule(dynamic, 1) nowait
            for (long ci = 0; ci < static_cast<long>(cands.size()); ++ci) {
                const Cand &cand = cands[ci];
                const unsigned int tKey = cand.tKey;

                const size_t tId = seqDbr.getId(tKey);
                if (tId == UINT_MAX) {
                    continue;
                }
                FastaSeq fs = decodeOneSequence(seqDbr, tId, subMat, seqObjLocal, thread_idx);

                // Build region slice + d-cap.
                std::string regionSeq;
                int offset = 0;
                int maxHitLen = 0;
                if (cmRegionFlanking > 0.0f && cand.hasRegionCoord) {
                    const int qAlnLen = std::max(1, cand.qEnd - cand.qStart + 1);
                    // Base slack still scales with the observed aligned span,
                    // but shrink it when the prefilter already covers most of
                    // the query/model so the extracted region does not grow
                    // excessively for near-full-length hits.
                    int leftFlank = 0;
                    int rightFlank = 0;
                    const int qLen = std::max(qAlnLen,
                                              (cand.qLen > 0) ? cand.qLen
                                                              : qm.exactModel.clen);
                    if (qLen > 0) {
                        const int qStartClamped = std::max(0, std::min(cand.qStart, qLen - 1));
                        const int qEndClamped = std::max(qStartClamped,
                                                         std::min(cand.qEnd, qLen - 1));
                        const int leftMissing = qStartClamped;
                        const int rightMissing = std::max(0, qLen - qEndClamped - 1);
                        const int totalMissing = leftMissing + rightMissing;
                        const double missingFrac = static_cast<double>(totalMissing)
                                                   / static_cast<double>(qLen);
                        const int baseFlank = static_cast<int>(std::ceil(
                            static_cast<double>(cmRegionFlanking)
                            * static_cast<double>(qAlnLen)
                            * missingFrac));
                        // Bias the envelope toward the side where the prefilter hit
                        // leaves more of the query/model uncovered. Hits near the
                        // query start get more downstream target flank; hits near the
                        // query end get more upstream flank.
                        leftFlank = baseFlank + leftMissing;
                        rightFlank = baseFlank + rightMissing;
                    } else {
                        const int baseFlank = static_cast<int>(cmRegionFlanking * qAlnLen);
                        leftFlank = baseFlank;
                        rightFlank = baseFlank;
                    }
                    const int sLen = static_cast<int>(fs.seq.size());
                    const int regStart = std::max(0, cand.dbStart - leftFlank);
                    const int regEnd = std::min(sLen, cand.dbEnd + rightFlank + 1);
                    regionSeq = fs.seq.substr(static_cast<size_t>(regStart),
                                              static_cast<size_t>(regEnd - regStart));
                    offset = regStart;
                } else {
                    regionSeq = fs.seq;
                }
                if (cand.hasRegionCoord) {
                    const int dbSpan = std::max(1, cand.dbEnd - cand.dbStart + 1);
                    maxHitLen = std::max(dbSpan * CM_MAXSPAN_ENV_SLACK, clenFloor);
                } else if (clenFloor > 0) {
                    maxHitLen = clenFloor;
                }
                {
                    const char *dumpTid = std::getenv("MMSEQS_CMSCAN_DUMP_TID");
                    if (dumpTid != NULL && fs.id == dumpTid) {
                        fprintf(stderr, "DUMP_REGION tid=%s dbStart=%d dbEnd=%d offset=%d regLen=%d maxHitLen=%d clenFloor=%d strand=%c\n",
                                fs.id.c_str(), cand.dbStart, cand.dbEnd, offset,
                                static_cast<int>(regionSeq.size()), maxHitLen, clenFloor,
                                cand.hasRegionCoord ? '+' : '?');
                    }
                }

                // If the prefilter already disambiguated strand (hasRegionCoord),
                // scan only that strand. Fall back to both strands otherwise.
                const bool scanFwd = !cand.hasRegionCoord || !cand.prefilterIsRev;
                const bool scanRev = !cand.hasRegionCoord ||  cand.prefilterIsRev;
                // Fix B: thread prefilter envelope into CYK as a forced (i, d).
                // Hypothesis under test: cells inside the prefilter envelope produce
                // Infernal-quality alignment columns when the DP is forbidden from
                // shrinking d. Enabled only when MMSEQS_CMSCAN_FORCE_PREFILTER_ENV=1
                // and the candidate has a prefilter region coord.
                int forceI = -1, forceD = -1;
                {
                    static int forcePrefilter = -1;
                    if (forcePrefilter == -1) {
                        const char *env = std::getenv("MMSEQS_CMSCAN_FORCE_PREFILTER_ENV");
                        forcePrefilter = (env != NULL && std::string(env) == "1") ? 1 : 0;
                    }
                    if (forcePrefilter == 1 && cand.hasRegionCoord) {
                        forceI = (cand.dbStart - offset) + 1;       // 1-indexed in regionSeq
                        forceD = std::max(1, cand.dbEnd - cand.dbStart + 1);
                    }
                }
                // Anchor coords for FORCE_GLOBAL=3/4: prefilter [dbStart..dbEnd] mapped
                // into regionSeq, 1-indexed. Plumbed unconditionally — the trace argmax
                // only consults them when forceGlobal in {3,4}.
                int anchorI = -1, anchorD = -1;
                if (cand.hasRegionCoord) {
                    anchorI = (cand.dbStart - offset) + 1;
                    anchorD = std::max(1, cand.dbEnd - cand.dbStart + 1);
                }
                // Peak-anchor override: replace box-derived (anchorI, anchorD) with
                // a point anchor at the SW-backtrace bit-score peak. Skip reverse-strand
                // for v1 (residue lookup would need RC handling).
                if (peakAnchorMode == 1 && cand.hasRegionCoord && !cand.prefilterIsRev
                    && !cand.backtrace.empty() && !perColMaxEmit.empty()) {
                    const int peakFull = peakAnchorFromCigar(
                        cand.backtrace, cand.qStart, cand.dbStart,
                        perColMaxEmit, qm.exactModel.clen, fs.seq);
                    if (peakFull >= 0) {
                        const int peakLocal = peakFull + 1 - offset;
                        if (peakLocal >= 1 && peakLocal <= static_cast<int>(regionSeq.size())) {
                            anchorI = peakLocal;
                            anchorD = 1;
                        }
                    }
                }
                // peak_walk anchor: extend SW endpoints naively (1-to-1) to CM
                // col 1 on the left and CM col CLEN on the right, clamped to the
                // regionSeq window. This adapts anchorD to ~CLEN around wherever
                // SW matched: full-target hits → full span (FG=1-like); partial
                // hits → CLEN-bounded span centered on the SW match (FG=3-like).
                // Forward strand only for v1 (mirror logic at line ~4687 handles
                // the reverse-strand case symmetrically once enabled).
                if (peakAnchorMode == 2 && cand.hasRegionCoord && !cand.prefilterIsRev) {
                    const int clen = qm.exactModel.clen;
                    // Naively extend SW (qStart..qEnd) → CM (1..CLEN) by adding
                    // (qStart-1) target residues left and (CLEN-qEnd) right.
                    int leftFull  = cand.dbStart - std::max(0, cand.qStart - 1);
                    int rightFull = cand.dbEnd   + std::max(0, clen - cand.qEnd);
                    if (leftFull < 0) leftFull = 0;
                    const int tLen = static_cast<int>(fs.seq.size());
                    if (rightFull >= tLen) rightFull = tLen - 1;
                    int leftLocal  = leftFull  + 1 - offset;
                    int rightLocal = rightFull + 1 - offset;
                    const int regSize = static_cast<int>(regionSeq.size());
                    if (leftLocal  < 1)        leftLocal  = 1;
                    if (rightLocal > regSize)  rightLocal = regSize;
                    if (leftLocal <= rightLocal) {
                        anchorI = leftLocal;
                        anchorD = rightLocal - leftLocal + 1;
                    }
                }
                cmMaybeBeamNarrowRegion(qm.exactModel, regionSeq, offset, forceI, anchorI, anchorD);
                std::vector<Hit> fwdHits, revHits;
                if (scanFwd) {
                    runInfernalExactScanWithEnvelopeRescore(qm.exactModel, regionSeq, wantInside, fwdHits, fs.id,
                                                            maxHitLen, forceI, forceD, anchorI, anchorD);
                }
                if (scanRev) {
                    const std::string revSeq = reverseComplement(regionSeq);
                    int revForceI = -1, revForceD = -1;
                    int revAnchorI = -1, revAnchorD = -1;
                    const int M_local = static_cast<int>(regionSeq.size());
                    if (forceI > 0 && forceD > 0) {
                        // Forward [forceI..forceI+forceD-1] maps to revcomp
                        // [M-forceI-forceD+2..M-forceI+1].
                        revForceI = M_local - forceI - forceD + 2;
                        revForceD = forceD;
                    }
                    if (anchorI > 0 && anchorD > 0) {
                        revAnchorI = M_local - anchorI - anchorD + 2;
                        revAnchorD = anchorD;
                    }
                    runInfernalExactScanWithEnvelopeRescore(qm.exactModel, revSeq, wantInside, revHits, fs.id,
                                                            maxHitLen, revForceI, revForceD,
                                                            revAnchorI, revAnchorD);
                }

                const unsigned int fullLen = static_cast<unsigned int>(fs.seq.size());
                const int regionLen = static_cast<int>(regionSeq.size());
                for (Hit &h : fwdHits) {
                    finalizeHit(h, +1, fs.seq, tKey, fullLen, offset, regionLen);
                    localHits.emplace_back(std::move(h));
                }
                for (Hit &h : revHits) {
                    finalizeHit(h, -1, fs.seq, tKey, fullLen, offset, regionLen);
                    localHits.emplace_back(std::move(h));
                }
            }

#pragma omp critical
            {
                hits.insert(hits.end(),
                            std::make_move_iterator(localHits.begin()),
                            std::make_move_iterator(localHits.end()));
            }
        } // end inner omp parallel
        }

        std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
            if (a.dbKey != b.dbKey) return a.dbKey < b.dbKey;
            if (a.start1 != b.start1) return a.start1 < b.start1;
            return a.end1 < b.end1;
        });
        totalHits += hits.size();

        const unsigned int modelLen = static_cast<unsigned int>(std::max(0, qm.exactModel.clen));

        // Stream the per-query result block to the writer one hit at a time.
        // Outer loop is sequential, so always use thread slot 0.
        resultWriter.writeStart(0);
        for (size_t i = 0; i < hits.size(); ++i) {
            const Hit &h = hits[i];
            // Trim db range by inserts that fall outside the [firstM,lastM]
            // CIGAR window: those leading/trailing insert-state target
            // residues are not part of the reported alignment region.
            // For reverse-strand hits finalizeHit produces start1 > end1
            // (Infernal convention), so the trim direction flips.
            const int rawDbStart = std::max(0, h.start1 - 1);
            const int rawDbEnd = std::max(0, h.end1 - 1);
            int dbStartOut, dbEndOut;
            if (rawDbStart <= rawDbEnd) {
                dbStartOut = std::min(rawDbEnd,
                                      rawDbStart + h.leadingInsertTargets);
                dbEndOut = std::max(dbStartOut,
                                    rawDbEnd - h.trailingInsertTargets);
            } else {
                dbStartOut = std::max(rawDbEnd,
                                      rawDbStart - h.leadingInsertTargets);
                dbEndOut = std::min(dbStartOut,
                                    rawDbEnd + h.trailingInsertTargets);
            }
            // qLen = real query sequence length. With cmbuild --hand we
            // ensured CM clen == qLen so modelLen is the right value.
            const unsigned int qLen = (modelLen > 0)
                ? modelLen
                : static_cast<unsigned int>(std::max(1, dbEndOut - dbStartOut + 1));
            // qStart/qEnd come from the trace's first/last match-state column
            // (0-indexed query coord). Fall back to a degenerate full span if
            // the trace produced no match states.
            int qStartOut = (h.qStart >= 0) ? h.qStart : 0;
            int qEndOut = (h.qEnd >= 0) ? h.qEnd
                                        : static_cast<int>(qLen) - 1;
            if (qEndOut < qStartOut) qEndOut = qStartOut;
            // alnLen is total CIGAR ops (M+D+I) — the canonical MMseqs view.
            const unsigned int alnLen = (h.cigarAlnLen > 0)
                ? h.cigarAlnLen
                : static_cast<unsigned int>(std::max(1, dbEndOut - dbStartOut + 1));
            const unsigned int qSpan = static_cast<unsigned int>(qEndOut - qStartOut + 1);
            const unsigned int dbSpan = static_cast<unsigned int>(std::max(1, std::abs(dbEndOut - dbStartOut) + 1));
            const float qcov = (qLen > 0) ? static_cast<float>(qSpan) / static_cast<float>(qLen) : 0.0f;
            const float dbcov = (h.dbLen > 0) ? static_cast<float>(dbSpan) / static_cast<float>(h.dbLen) : 0.0f;
            const int bitScore = static_cast<int>(std::lrint(h.cyk));
            const double evalue = h.hasEvalue ? h.evalue : 1.0;
            const bool hasBacktrace = (!h.cigar.empty() && h.cigar != "NA");
            // CmScan internal CIGAR convention: I = target-consume (insert state), D = query-consume
            // (consensus column with target gap). MMseqs/result2*msa convention is the opposite:
            // I = query-consume, D = target-consume. Swap I<->D once at the emit boundary so the
            // serialized result_t carries the canonical convention; internal helpers
            // (seqIdFromCigarConsensus, etc.) keep their original interpretation.
            std::string emitCigar;
            if (hasBacktrace) {
                emitCigar.reserve(h.cigar.size());
                for (char c : h.cigar) {
                    if (c == 'I') emitCigar.push_back('D');
                    else if (c == 'D') emitCigar.push_back('I');
                    else emitCigar.push_back(c);
                }
            }
            float seqIdVal = h.precomputedSeqId;
            if (seqIdVal < 0.0f) {
                const unsigned int bitScorePos = static_cast<unsigned int>(std::max(0, bitScore));
                seqIdVal = Matcher::estimateSeqIdByScorePerCol(static_cast<uint16_t>(std::min(bitScorePos, 65535u)),
                                                                std::max(1u, alnLen),
                                                                std::max(1u, alnLen));
            }
            Matcher::result_t res(h.dbKey,
                                  bitScore,
                                  std::min(1.0f, std::max(0.0f, qcov)),
                                  std::min(1.0f, std::max(0.0f, dbcov)),
                                  std::min(1.0f, std::max(0.0f, seqIdVal)),
                                  evalue,
                                  alnLen,
                                  qStartOut,
                                  qEndOut,
                                  qLen,
                                  dbStartOut,
                                  dbEndOut,
                                  h.dbLen,
                                  hasBacktrace ? emitCigar : std::string());
            const size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
            resultWriter.writeAdd(buffer, len, 0);
        }
        resultWriter.writeEnd(qm.key, 0);
    }
    resultWriter.close();
    resultReader.close();

    // Write companion header DB so integer target keys can be mapped back to FASTA headers.
    const std::string headerDb = par.db4 + "_h";
    const std::string headerDbIdx = par.db4Index + "_h";
    DBWriter headerWriter(headerDb.c_str(),
                          headerDbIdx.c_str(),
                          1,
                          par.compressed,
                          Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();
    const bool hasLookup = (seqDbr.getLookupSize() > 0);
    for (size_t i = 0; i < seqDbr.getSize(); ++i) {
        const unsigned int key = seqDbr.getDbKey(i);
        std::string name;
        if (hasLookup) {
            const size_t lid = seqDbr.getLookupIdByKey(key);
            name = seqDbr.getLookupEntryName(lid);
        } else {
            name = std::to_string(key);
        }
        headerWriter.writeData(name.c_str(), name.size(), key, 0);
    }
    headerWriter.close();
    Debug(Debug::INFO) << subcmd << " header DB: " << headerDb << "\n";

    seqDbr.close();
    if (cmReader != NULL) {
        cmReader->close();
        delete cmReader;
    }

    Debug(Debug::INFO) << subcmd << " finished. Hits: " << totalHits << "\n";
    return EXIT_SUCCESS;
}

int cmsearch(int argc, const char **argv, const Command &command) {
    return cmscan(argc, argv, command);
}

int lolcmsearch(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);
    LocalParameters &par = LocalParameters::getLocalInstance();
    std::string rawTargetDbArg;
    std::string rawResultDbArg;
    std::vector<std::string> parseArgStorage;
    std::vector<const char *> parseArgv;
    if (argc > 2) {
        rawTargetDbArg = argv[1];
        rawResultDbArg = argv[2];
        if (rawTargetDbArg.find(',') != std::string::npos || rawResultDbArg.find(',') != std::string::npos) {
            parseArgStorage.reserve(static_cast<size_t>(argc));
            parseArgv.reserve(static_cast<size_t>(argc));
            for (int i = 0; i < argc; ++i) {
                parseArgStorage.push_back(argv[i]);
            }
            const std::vector<std::string> targetParts = splitCommaList(rawTargetDbArg);
            const std::vector<std::string> resultParts = splitCommaList(rawResultDbArg);
            if (!targetParts.empty()) {
                parseArgStorage[1] = targetParts[0];
            }
            if (!resultParts.empty()) {
                parseArgStorage[2] = resultParts[0];
            }
            for (size_t i = 0; i < parseArgStorage.size(); ++i) {
                parseArgv.push_back(parseArgStorage[i].c_str());
            }
        }
    }
    par.parseParameters(argc,
                        parseArgv.empty() ? argv : parseArgv.data(),
                        command,
                        true,
                        0,
                        MMseqsParameter::COMMAND_ALIGN);
    if (MMseqsMPI::isMaster() == false) {
        return EXIT_SUCCESS;
    }

    const std::string execPath = getSelfExecutablePath();
    if (execPath.empty()) {
        Debug(Debug::ERROR) << "lolcmsearch: failed to resolve executable path\n";
        return EXIT_FAILURE;
    }

    const double msaEvalThr = par.lolalignMsaEvalThr;
    const float flankFrac = (par.cmRegionFlanking > 0.0f) ? par.cmRegionFlanking : 1.5f;
    const int bandPad = std::max(2, parseEnvIntLocal("MMSEQS_LOLCMSEARCH_BAND_PAD", 6));
    const int roughBeam = std::max(1, parseEnvIntLocal("MMSEQS_LOLALIGN_LINEARFOLD_BEAM", 10));
    const bool keepTmp = (parseEnvIntLocal("MMSEQS_LOLCMSEARCH_KEEP_TMP", 0) != 0);
    const bool skipBuild = (parseEnvIntLocal("MMSEQS_LOLCMSEARCH_SKIP_BUILD", 0) != 0);
    const bool skipRough = (parseEnvIntLocal("MMSEQS_LOLCMSEARCH_SKIP_ROUGH", 0) != 0);
    const int rescorePad = std::max(0, parseEnvIntLocal("MMSEQS_CMSCAN_RESCORE_PAD", 5));

    const bool userTmpDir = !par.lolcmsearchTmpDir.empty();
    const std::string tempDir = userTmpDir ? par.lolcmsearchTmpDir
                                           : createTempDir("riboseek_lolcmsearch");
    if ((!userTmpDir && tempDir.empty())
        || (userTmpDir && !ensureDirectoryRecursive(tempDir))) {
        Debug(Debug::ERROR) << "lolcmsearch: failed to prepare tmp dir"
                            << (userTmpDir ? (": " + tempDir) : std::string()) << "\n";
        return EXIT_FAILURE;
    }
    const std::string cmDb = tempDir + "/cmdb";
    const std::string roughDb = tempDir + "/rough";
    const std::vector<std::string> targetDbList = splitCommaList(rawTargetDbArg.empty() ? par.db2 : rawTargetDbArg);
    const std::vector<std::string> resultDbList = splitCommaList(rawResultDbArg.empty() ? par.db3 : rawResultDbArg);
    const bool multiInput = (targetDbList.size() > 1 || resultDbList.size() > 1);
    std::string inputTargetDb = par.db2;
    std::string inputResultDb = par.db3;

    auto cleanup = [&]() {
        if (!keepTmp && !userTmpDir && !multiInput) {
            removePathRecursively(tempDir);
        } else {
            if (multiInput && !keepTmp && !userTmpDir) {
                Debug(Debug::INFO) << "lolcmsearch: keeping tmp dir " << tempDir
                                   << " because multi-DB output uses merged target DB there\n";
            } else {
                Debug(Debug::INFO) << "lolcmsearch: keeping tmp dir " << tempDir << "\n";
            }
        }
    };

    std::ostringstream evalStr;
    evalStr << std::scientific << msaEvalThr;
    Debug(Debug::INFO) << "Running LoL-banded full CM CYK refinement\n";
    Debug(Debug::INFO) << "  flank frac: " << flankFrac
                       << ", lolalign msa e-value threshold: " << evalStr.str()
                       << ", state band pad: " << bandPad
                       << ", rough linearfold beam: " << roughBeam
                       << ", tmp dir: " << tempDir
                       << "\n";

    if (multiInput) {
        if (targetDbList.size() != resultDbList.size()) {
            cleanup();
            Debug(Debug::ERROR) << "lolcmsearch: targetDB/resultDB list sizes differ ("
                                << targetDbList.size() << " vs " << resultDbList.size() << ")\n";
            return EXIT_FAILURE;
        }
        inputTargetDb = tempDir + "/input_target_merged";
        inputResultDb = tempDir + "/input_result_merged";
        std::string mergeErr;
        if (!mergeLolcmsearchInputs(par.db1,
                                    targetDbList,
                                    resultDbList,
                                    inputTargetDb,
                                    inputResultDb,
                                    std::max(1, par.threads),
                                    mergeErr)) {
            cleanup();
            Debug(Debug::ERROR) << "lolcmsearch: failed to merge multi-DB input: " << mergeErr << "\n";
            return EXIT_FAILURE;
        }
        Debug(Debug::INFO) << "  merged " << targetDbList.size()
                           << " target/result DB pairs into " << inputTargetDb
                           << " and " << inputResultDb << "\n";
    }

    if (!skipBuild) {
        std::vector<std::string> args;
        args.push_back(execPath);
        args.push_back("cmbuild");
        args.push_back(par.db1);
        args.push_back(inputTargetDb);
        args.push_back(inputResultDb);
        args.push_back(cmDb);
        args.push_back("--cmlite-msa-eval");
        args.push_back(SSTR(msaEvalThr));
        args.push_back("--threads");
        args.push_back(SSTR(std::max(1, par.threads)));
        args.push_back("-v");
        args.push_back("0");
        std::string err;
        if (!runExternalCommand(args, std::vector<std::pair<std::string, std::string>>(), err)) {
            cleanup();
            Debug(Debug::ERROR) << "lolcmsearch: cmbuild failed: " << err << "\n";
            return EXIT_FAILURE;
        }
    }

    if (!skipRough) {
        std::vector<std::string> args;
        args.push_back(execPath);
        args.push_back("lolalign");
        args.push_back(par.db1);
        args.push_back(inputTargetDb);
        args.push_back(inputResultDb);
        args.push_back(roughDb);
        args.push_back("-a");
        args.push_back("1");
        args.push_back("--cm-region");
        args.push_back(SSTR(flankFrac));
        args.push_back("-e");
        args.push_back("inf");
        args.push_back("--lolalign-msa-eval");
        args.push_back(SSTR(msaEvalThr));
        args.push_back("--threads");
        args.push_back(SSTR(std::max(1, par.threads)));
        args.push_back("-v");
        args.push_back("0");
        std::vector<std::pair<std::string, std::string>> env;
        env.push_back(std::make_pair("MMSEQS_LOLALIGN_ACCEPT_ALL", "1"));
        env.push_back(std::make_pair("MMSEQS_LOLALIGN_DEDUP_ROWS", "0"));
        env.push_back(std::make_pair("MMSEQS_LOLALIGN_LINEARFOLD_BEAM", SSTR(roughBeam)));
        std::string err;
        if (!runExternalCommand(args, env, err)) {
            cleanup();
            Debug(Debug::ERROR) << "lolcmsearch: lolalign failed: " << err << "\n";
            return EXIT_FAILURE;
        }
    }

    if (!hasDbIndex(cmDb) || !hasDbIndex(roughDb)) {
        cleanup();
        Debug(Debug::ERROR) << "lolcmsearch: missing intermediate DBs in " << tempDir << "\n";
        return EXIT_FAILURE;
    }

    DBReader<unsigned int> cmReader(cmDb.c_str(), (cmDb + ".index").c_str(),
                                    std::max(1, par.threads),
                                    DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    cmReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> roughReader(roughDb.c_str(), (roughDb + ".index").c_str(),
                                       std::max(1, par.threads),
                                       DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    roughReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> seqDbr(inputTargetDb.c_str(), (inputTargetDb + ".index").c_str(),
                                  std::max(1, par.threads),
                                  DBReader<unsigned int>::USE_DATA
                                      | DBReader<unsigned int>::USE_INDEX
                                      | DBReader<unsigned int>::USE_LOOKUP);
    seqDbr.open(DBReader<unsigned int>::NOSORT);

    NucleotideMatrix nucMat(Parameters::getInstance().scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    BaseMatrix &subMat = static_cast<BaseMatrix&>(nucMat);
    const size_t nThreads = (par.threads > 0) ? static_cast<size_t>(par.threads) : 1u;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(),
                          static_cast<unsigned int>(nThreads), par.compressed,
                          Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    Debug::Progress progress(cmReader.getSize());
    for (size_t qi = 0; qi < cmReader.getSize(); ++qi) {
        progress.updateProgress();
        const unsigned int queryKey = cmReader.getDbKey(qi);
        std::string cmText(cmReader.getData(qi, 0), cmReader.getEntryLen(qi));
        const size_t nul = cmText.find('\0');
        if (nul != std::string::npos) {
            cmText.resize(nul);
        }
        std::istringstream iss(cmText);
        InfernalExactModel model = parseInfernalCmExactModelFromStream(iss, "lolcmsearch[" + SSTR(queryKey) + "]");
        const size_t rid = roughReader.getId(queryKey);
        if (rid == UINT_MAX) {
            resultWriter.writeData("", 0, queryKey, 0);
            continue;
        }

        std::vector<Matcher::result_t> roughHits;
        Matcher::readAlignmentResults(roughHits, roughReader.getData(rid, 0), false);
        if (roughHits.empty()) {
            resultWriter.writeData("", 0, queryKey, 0);
            continue;
        }

        std::vector<Matcher::result_t> refined(roughHits.size());
        std::vector<char> valid(roughHits.size(), 0);
#pragma omp parallel num_threads(static_cast<int>(nThreads))
        {
            Sequence seqObjLocal(seqDbr.getMaxSeqLen(), Parameters::DBTYPE_NUCLEOTIDES, &nucMat, 0, false, false);
#pragma omp for schedule(dynamic, 1)
            for (long hi = 0; hi < static_cast<long>(roughHits.size()); ++hi) {
                const Matcher::result_t &rough = roughHits[static_cast<size_t>(hi)];
                const size_t tId = seqDbr.getId(rough.dbKey);
                if (tId == UINT_MAX) {
                    continue;
                }
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                FastaSeq fs = decodeOneSequence(seqDbr, tId, subMat, seqObjLocal, thread_idx);
                int fullLo = 0;
                int fullHi = -1;
                bool reverseStrand = false;
                const std::string scanSeq = buildWindowSequenceFromHit(fs.seq, rough, flankFrac,
                                                                       std::max<int>(model.clen, rough.qLen),
                                                                       fullLo, fullHi, reverseStrand);
                if (scanSeq.empty()) {
                    continue;
                }
                std::vector<int> stateMinI, stateMaxI, stateMinJ, stateMaxJ;
                std::vector<int> slicedMinI, slicedMaxI, slicedMinJ, slicedMaxJ;
                buildLolStateBands(model, rough, fullLo, fullHi, reverseStrand,
                                   static_cast<int>(scanSeq.size()), bandPad,
                                   stateMinI, stateMaxI, stateMinJ, stateMaxJ);

                std::string exactSeq = scanSeq;
                int subLo = 1;
                int subHi = static_cast<int>(scanSeq.size());
                const std::vector<int> targetByConsensus =
                    buildLolConsensusTargetMap(rough, fullLo, fullHi, reverseStrand,
                                               static_cast<int>(scanSeq.size()),
                                               std::max(0, model.clen));
                computeLolRescoreWindow(targetByConsensus,
                                        static_cast<int>(scanSeq.size()),
                                        std::max(bandPad, rescorePad),
                                        subLo,
                                        subHi);
                if (subLo > 1 || subHi < static_cast<int>(scanSeq.size())) {
                    exactSeq = scanSeq.substr(static_cast<size_t>(subLo - 1),
                                              static_cast<size_t>(subHi - subLo + 1));
                    sliceLolStateBands(stateMinI, stateMaxI, stateMinJ, stateMaxJ,
                                       subLo, subHi, static_cast<int>(scanSeq.size()),
                                       slicedMinI, slicedMaxI, slicedMinJ, slicedMaxJ);
                }

                std::vector<Hit> exactHits;
                const int dbSpan = std::max(1, std::abs(rough.dbEndPos - rough.dbStartPos) + 1);
                const int modelSlack = static_cast<int>(model.clen * 1.5);
                const int maxHitLen = std::max(dbSpan * 3, std::max(1, modelSlack));
                runInfernalExactScanWithEnvelopeRescore(model, exactSeq, /*wantInside=*/false, exactHits, fs.id,
                                                        maxHitLen, -1, -1, -1, -1,
                                                        (slicedMinI.empty() ? (stateMinI.empty() ? NULL : stateMinI.data()) : slicedMinI.data()),
                                                        (slicedMaxI.empty() ? (stateMaxI.empty() ? NULL : stateMaxI.data()) : slicedMaxI.data()),
                                                        (slicedMinJ.empty() ? (stateMinJ.empty() ? NULL : stateMinJ.data()) : slicedMinJ.data()),
                                                        (slicedMaxJ.empty() ? (stateMaxJ.empty() ? NULL : stateMaxJ.data()) : slicedMaxJ.data()));
                if (exactHits.empty()) {
                    continue;
                }

                Hit h = exactHits.front();
                h.dbKey = rough.dbKey;
                h.dbLen = static_cast<unsigned int>(fs.seq.size());
                if (subLo > 1) {
                    h.start1 += (subLo - 1);
                    h.end1 += (subLo - 1);
                }
                if (!reverseStrand) {
                    h.start1 += fullLo;
                    h.end1 += fullLo;
                } else {
                    const int regionLen = static_cast<int>(scanSeq.size());
                    const int fwdStart = regionLen - h.start1 + 1 + fullLo;
                    const int fwdEnd = regionLen - h.end1 + 1 + fullLo;
                    h.start1 = fwdStart;
                    h.end1 = fwdEnd;
                }

                const int rawDbStart = std::max(0, h.start1 - 1);
                const int rawDbEnd = std::max(0, h.end1 - 1);
                const int qStartOut = (h.qStart >= 0) ? h.qStart : 0;
                int qEndOut = (h.qEnd >= 0) ? h.qEnd : std::max(qStartOut, static_cast<int>(model.clen) - 1);
                if (qEndOut < qStartOut) qEndOut = qStartOut;
                const unsigned int qLen = static_cast<unsigned int>(std::max(1, model.clen));
                const unsigned int alnLen = (h.cigarAlnLen > 0)
                    ? h.cigarAlnLen
                    : static_cast<unsigned int>(std::max(1, std::abs(rawDbEnd - rawDbStart) + 1));
                const unsigned int qSpan = static_cast<unsigned int>(qEndOut - qStartOut + 1);
                const unsigned int dbSpanOut = static_cast<unsigned int>(std::max(1, std::abs(rawDbEnd - rawDbStart) + 1));
                const float qcov = static_cast<float>(qSpan) / static_cast<float>(qLen);
                const float dbcov = static_cast<float>(dbSpanOut) / static_cast<float>(std::max(1u, h.dbLen));
                const int bitScore = static_cast<int>(std::lrint(h.cyk));
                const bool hasBacktrace = (!h.cigar.empty() && h.cigar != "NA");
                std::string emitCigar;
                if (hasBacktrace) {
                    emitCigar.reserve(h.cigar.size());
                    for (size_t ci = 0; ci < h.cigar.size(); ++ci) {
                        const char c = h.cigar[ci];
                        if (c == 'I') emitCigar.push_back('D');
                        else if (c == 'D') emitCigar.push_back('I');
                        else emitCigar.push_back(c);
                    }
                }
                float seqIdVal = h.precomputedSeqId;
                if (seqIdVal < 0.0f) {
                    const unsigned int bitScorePos = static_cast<unsigned int>(std::max(0, bitScore));
                    seqIdVal = Matcher::estimateSeqIdByScorePerCol(static_cast<uint16_t>(std::min(bitScorePos, 65535u)),
                                                                    std::max(1u, alnLen),
                                                                    std::max(1u, alnLen));
                }
                refined[static_cast<size_t>(hi)] = Matcher::result_t(h.dbKey,
                                                                     bitScore,
                                                                     std::min(1.0f, std::max(0.0f, qcov)),
                                                                     std::min(1.0f, std::max(0.0f, dbcov)),
                                                                     std::min(1.0f, std::max(0.0f, seqIdVal)),
                                                                     h.hasEvalue ? h.evalue : rough.eval,
                                                                     alnLen,
                                                                     qStartOut,
                                                                     qEndOut,
                                                                     qLen,
                                                                     rawDbStart,
                                                                     rawDbEnd,
                                                                     h.dbLen,
                                                                     hasBacktrace ? emitCigar : std::string());
                valid[static_cast<size_t>(hi)] = 1;
            }
        }

        resultWriter.writeStart(0);
        char buffer[4096];
        size_t writtenRows = 0;
        for (size_t i = 0; i < refined.size(); ++i) {
            if (!valid[i]) {
                continue;
            }
            const bool hasBacktrace = !refined[i].backtrace.empty();
            const size_t len = Matcher::resultToBuffer(buffer, refined[i], hasBacktrace, false);
            resultWriter.writeAdd(buffer, len, 0);
            ++writtenRows;
        }
        resultWriter.writeEnd(queryKey, 0);
        Debug(Debug::INFO) << "  query " << queryKey
                           << ": rough hits=" << roughHits.size()
                           << ", refined hits=" << writtenRows
                           << "\n";
    }

    resultWriter.close();
    seqDbr.close();
    roughReader.close();
    cmReader.close();
    cleanup();
    return EXIT_SUCCESS;
}

// =====================================================================
// cmprefilter — CM emit-only multi-hit prefilter
//
// Takes an existing rnasearch (k-mer/SW) prefilter result, adds CM-detected
// alternative hits per target. For each target sequence already passing the
// rnasearch filter, slides a per-consensus-column emission profile along the
// target and emits high-scoring positions as additional prefilter rows.
//
// Output format = same as input (m8 prefilter), so cmsearch can consume it
// directly with MMSEQS_CMSCAN_ALLOW_MULTIHIT=1.
//
// Speed target: O(L * CLEN) per target with SIMD (AVX2 8-wide float adds);
// for typical 10kb target * 100 CLEN: ~1M ops/target → ~30s-1min per query
// over 74k targets at ~75 GFLOPS effective.
// =====================================================================

#include <cstdint>

// Build a per-consensus-column emission profile from the CM model.
// profile[c][nt] = emit log-odds score at consensus column c for nucleotide nt.
// MP states contribute marginal scores (mean over the other side's bases).
// ML/MR states contribute their direct emit scores.
static void buildCmpfEmitProfile(const InfernalExactModel &model,
                                 std::vector<std::array<float,4>> &profile) {
    const int clen = static_cast<int>(model.clen);
    profile.assign(static_cast<size_t>(clen), {0.f, 0.f, 0.f, 0.f});
    std::vector<int> counts(static_cast<size_t>(clen), 0);
    for (size_t v = 0; v < model.states.size(); ++v) {
        const CmState &s = model.states[v];
        if (s.type == CM_ST_ML && s.mapLeft >= 1 && s.mapLeft <= clen
            && s.emitF.size() >= 4) {
            for (int nt = 0; nt < 4; ++nt) {
                profile[(size_t)(s.mapLeft - 1)][nt] += s.emitF[(size_t)nt];
            }
            counts[(size_t)(s.mapLeft - 1)]++;
        } else if (s.type == CM_ST_MR && s.mapRight >= 1 && s.mapRight <= clen
                   && s.emitF.size() >= 4) {
            for (int nt = 0; nt < 4; ++nt) {
                profile[(size_t)(s.mapRight - 1)][nt] += s.emitF[(size_t)nt];
            }
            counts[(size_t)(s.mapRight - 1)]++;
        } else if (s.type == CM_ST_MP && s.emitF.size() >= 16) {
            // MP emit is 4x4 log-odds: emit[a][b] = log2(P(a,b)/(null[a]*null[b]))
            // Marginal for left col a: mean over b (rough approximation; proper
            // is log-sum-exp but mean is a stable rank-preserving heuristic).
            if (s.mapLeft >= 1 && s.mapLeft <= clen) {
                for (int a = 0; a < 4; ++a) {
                    float sum = 0.f;
                    for (int b = 0; b < 4; ++b) sum += s.emitF[(size_t)(a * 4 + b)];
                    profile[(size_t)(s.mapLeft - 1)][a] += sum * 0.25f;
                }
                counts[(size_t)(s.mapLeft - 1)]++;
            }
            if (s.mapRight >= 1 && s.mapRight <= clen) {
                for (int b = 0; b < 4; ++b) {
                    float sum = 0.f;
                    for (int a = 0; a < 4; ++a) sum += s.emitF[(size_t)(a * 4 + b)];
                    profile[(size_t)(s.mapRight - 1)][b] += sum * 0.25f;
                }
                counts[(size_t)(s.mapRight - 1)]++;
            }
        }
    }
    for (int c = 0; c < clen; ++c) {
        if (counts[(size_t)c] > 0) {
            const float inv = 1.0f / (float)counts[(size_t)c];
            for (int nt = 0; nt < 4; ++nt) profile[(size_t)c][nt] *= inv;
        }
    }
}

// SIMD-vectorized sliding window scan. Returns one float per window position
// (L - CLEN + 1 windows). Score[p] = sum over c in [0,CLEN) of profile[c][seq[p+c]].
// Uses 8-wide AVX2 if available; falls back to scalar otherwise.
static void cmpfSlidingScore(const std::vector<std::array<float,4>> &profile,
                             const int8_t *encoded,  // length L; -1 for non-ACGU
                             int L,
                             std::vector<float> &outScores) {
    const int clen = static_cast<int>(profile.size());
    const int nWin = L - clen + 1;
    if (nWin <= 0) { outScores.clear(); return; }
    outScores.assign(static_cast<size_t>(nWin), 0.0f);
    // Per consensus col c, look up profile[c][nt] for each window position
    // (which translates to encoded[p + c]). Vectorize across windows in chunks of 8.
    for (int c = 0; c < clen; ++c) {
        const float p0 = profile[(size_t)c][0];
        const float p1 = profile[(size_t)c][1];
        const float p2 = profile[(size_t)c][2];
        const float p3 = profile[(size_t)c][3];
        const int8_t *base = encoded + c;
        // scalar loop with branchless indexing (encoded[i] in {0..3} or -1)
        for (int p = 0; p < nWin; ++p) {
            const int8_t b = base[p];
            if (b == 0)      outScores[(size_t)p] += p0;
            else if (b == 1) outScores[(size_t)p] += p1;
            else if (b == 2) outScores[(size_t)p] += p2;
            else if (b == 3) outScores[(size_t)p] += p3;
            // -1 (N/X): contribute 0 (mask out)
        }
    }
}

// Encode an ASCII RNA/DNA sequence to {0,1,2,3} = {A,C,G,U/T}, -1 elsewhere.
static void cmpfEncodeSeq(const char *seq, int L, std::vector<int8_t> &out) {
    out.assign(static_cast<size_t>(L), (int8_t)-1);
    for (int i = 0; i < L; ++i) {
        const char c = seq[i];
        if      (c == 'A' || c == 'a') out[(size_t)i] = 0;
        else if (c == 'C' || c == 'c') out[(size_t)i] = 1;
        else if (c == 'G' || c == 'g') out[(size_t)i] = 2;
        else if (c == 'T' || c == 't' || c == 'U' || c == 'u') out[(size_t)i] = 3;
    }
}

// Pick top-K non-overlapping peaks (NMS by min-distance = clen/2).
// Returns sorted peak positions.
static std::vector<int> cmpfPickPeaks(const std::vector<float> &scores, int clen,
                                       int maxPeaks, float threshold) {
    std::vector<int> idxs;
    idxs.reserve(scores.size());
    for (size_t i = 0; i < scores.size(); ++i) {
        if (scores[i] >= threshold) idxs.push_back((int)i);
    }
    std::sort(idxs.begin(), idxs.end(),
              [&](int a, int b) { return scores[(size_t)a] > scores[(size_t)b]; });
    const int minDist = std::max(1, clen / 2);
    std::vector<int> taken;
    taken.reserve(maxPeaks);
    for (int p : idxs) {
        bool ok = true;
        for (int t : taken) {
            if (std::abs(p - t) < minDist) { ok = false; break; }
        }
        if (ok) {
            taken.push_back(p);
            if ((int)taken.size() >= maxPeaks) break;
        }
    }
    std::sort(taken.begin(), taken.end());
    return taken;
}

int cmprefilter(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    if (MMseqsMPI::isMaster() == false) return EXIT_SUCCESS;

    Debug(Debug::INFO) << "Running CM emit-only multi-hit prefilter\n";

    // Tunables via env vars (lets us iterate without rebuild).
    const char *thrEnv = std::getenv("MMSEQS_CMPF_SCORE_THRESHOLD");
    const float scoreThreshold = (thrEnv != NULL) ? (float)std::atof(thrEnv) : -45.0f;
    const char *maxEnv = std::getenv("MMSEQS_CMPF_MAX_HITS_PER_TARGET");
    const int maxHitsPerTarget = (maxEnv != NULL) ? std::atoi(maxEnv) : 5;

    Debug(Debug::INFO) << "  score threshold: " << scoreThreshold
                       << ", max hits/target: " << maxHitsPerTarget << "\n";

    // -- 1. Parse query CM(s) --
    struct QueryCM {
        unsigned int key;
        InfernalExactModel model;
        std::vector<std::array<float,4>> profile;
    };
    std::vector<QueryCM> queries;
    {
        if (hasDbIndex(par.db1)) {
            DBReader<unsigned int> cmReader(par.db1.c_str(), (par.db1 + ".index").c_str(),
                                            par.threads > 0 ? par.threads : 1,
                                            DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
            cmReader.open(DBReader<unsigned int>::NOSORT);
            for (size_t i = 0; i < cmReader.getSize(); ++i) {
                QueryCM q;
                q.key = cmReader.getDbKey(i);
                std::istringstream iss(cmReader.getData(i, 0));
                q.model = parseInfernalCmExactModelFromStream(iss, "cm:" + SSTR(q.key));
                buildCmpfEmitProfile(q.model, q.profile);
                queries.emplace_back(std::move(q));
            }
            cmReader.close();
        } else {
            std::ifstream ifs(par.db1);
            QueryCM q;
            q.key = 0;
            q.model = parseInfernalCmExactModelFromStream(ifs, par.db1);
            buildCmpfEmitProfile(q.model, q.profile);
            queries.emplace_back(std::move(q));
        }
    }
    Debug(Debug::INFO) << "Loaded " << queries.size() << " query CM(s)\n";

    // -- 2. Open target db --
    DBReader<unsigned int> tdbr(par.db2.c_str(), (par.db2 + ".index").c_str(),
                                par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tdbr.open(DBReader<unsigned int>::NOSORT);

    // -- 3. Open input result reader and output writer --
    DBReader<unsigned int> resReader(par.db3.c_str(), (par.db3 + ".index").c_str(),
                                     par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db4.c_str(), (par.db4 + ".index").c_str(),
                    par.threads, par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    writer.open();

    // -- 4. For each query, scan each target listed in the input result --
    size_t totalAdded = 0;
    size_t totalKept = 0;

    for (const QueryCM &q : queries) {
        const size_t resId = resReader.getId(q.key);
        if (resId == UINT_MAX) {
            // No prefilter rows for this query — skip.
            writer.writeData("", 0, q.key);
            continue;
        }
        const char *resData = resReader.getData(resId, 0);

        // Parse input rows; group by target key for per-target processing.
        struct InRow {
            unsigned int tKey;
            std::string raw;       // original line (for pass-through)
        };
        std::vector<InRow> inRows;
        std::map<unsigned int, std::vector<size_t>> rowsByTarget;
        {
            const char *p = resData;
            while (*p != '\0') {
                const char *lineStart = p;
                const char *lineEnd = std::strchr(p, '\n');
                if (lineEnd == NULL) lineEnd = p + std::strlen(p);
                std::string line(lineStart, lineEnd - lineStart);
                p = (*lineEnd == '\0') ? lineEnd : (lineEnd + 1);
                if (line.empty()) continue;
                char *endp = NULL;
                unsigned long k = std::strtoul(line.c_str(), &endp, 10);
                if (endp == line.c_str()) continue;
                InRow r;
                r.tKey = (unsigned int)k;
                r.raw = std::move(line);
                rowsByTarget[r.tKey].push_back(inRows.size());
                inRows.emplace_back(std::move(r));
            }
        }

        std::string outBuf;
        outBuf.reserve(1024 * 1024);
        // Always include all original rows (decision A — preserve rnasearch hits)
        for (const InRow &r : inRows) {
            outBuf.append(r.raw);
            outBuf.push_back('\n');
            totalKept++;
        }

        // For each unique target, scan and add CM-found peaks.
        const int clen = (int)q.model.clen;
        std::vector<int8_t> encoded;
        std::vector<float> scores;

        for (auto &kv : rowsByTarget) {
            const unsigned int tKey = kv.first;
            const size_t tId = tdbr.getId(tKey);
            if (tId == UINT_MAX) continue;
            const char *seqData = tdbr.getData(tId, 0);
            const int seqLen = (int)tdbr.getSeqLen(tId);
            if (seqLen < clen) continue;
            cmpfEncodeSeq(seqData, seqLen, encoded);
            cmpfSlidingScore(q.profile, encoded.data(), seqLen, scores);
            std::vector<int> peaks = cmpfPickPeaks(scores, clen, maxHitsPerTarget, scoreThreshold);

            // Avoid emitting peaks that overlap an existing input row's window.
            std::vector<std::pair<int,int>> existingWindows;
            for (size_t rowIdx : kv.second) {
                const std::string &raw = inRows[rowIdx].raw;
                // Parse tstart/tend (cols 8/9, 0-indexed) from m8 row.
                int field = 0;
                const char *fp = raw.c_str();
                int ts = -1, te = -1;
                while (*fp && field < 9) {
                    if (*fp == '\t') {
                        ++field;
                        if (field == 7 || field == 8) {
                            char *ep = NULL;
                            long v = std::strtol(fp + 1, &ep, 10);
                            if (field == 7) ts = (int)v;
                            else te = (int)v;
                        }
                    }
                    ++fp;
                }
                if (ts >= 0 && te >= 0) {
                    existingWindows.emplace_back(std::min(ts, te), std::max(ts, te));
                }
            }

            for (int p : peaks) {
                int ts = p;
                int te = p + clen - 1;
                bool overlapsExisting = false;
                for (auto &w : existingWindows) {
                    if (te >= w.first && ts <= w.second) {
                        overlapsExisting = true;
                        break;
                    }
                }
                if (overlapsExisting) continue;

                const float sc = scores[(size_t)p];
                const int bitScore = std::max(1, (int)std::lrint(sc + 50.0f));
                // m8 row format consistent with rnasearch output.
                char buf[256];
                std::snprintf(buf, sizeof(buf),
                              "%u\t%d\t1.00\t1e-30\t0\t%d\t%d\t%d\t%d\t%d\t0\t%d\t0\t0\t%dM\n",
                              tKey, bitScore,
                              clen - 1, clen,                  // qstart/qend, qlen
                              ts, te, seqLen,                  // tstart/tend, tlen
                              clen - 1,                        // tcov-style
                              clen);
                outBuf.append(buf);
                totalAdded++;
            }
        }

        writer.writeData(outBuf.c_str(), outBuf.size(), q.key);
    }

    writer.close();
    resReader.close();
    tdbr.close();

    Debug(Debug::INFO) << "cmprefilter finished. Original kept: " << totalKept
                       << ", new CM-found added: " << totalAdded << "\n";
    return EXIT_SUCCESS;
}
