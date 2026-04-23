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

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <unistd.h>
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

struct UnaryRule {
    int lhs;
    char base;
    double logp;
};

struct BinaryRule {
    int lhs;
    int rhs1;
    int rhs2;
    double logp;
};

struct PairRule {
    int lhs;
    int inner; // -1 means epsilon inner (consumes only the two paired terminals)
    char leftBase;
    char rightBase;
    double logp;
};

struct EpsRule {
    int lhs;
    int rhs;
    double logp;
};

struct Grammar {
    std::map<std::string, int> ntToId;
    std::vector<std::string> idToNt;
    int start = -1;
    double calA = 1.0;
    double calB = 0.0;
    double calSplit = std::numeric_limits<double>::infinity();
    double calNegA = 1.0;
    double calNegB = 0.0;
    bool hasPiecewiseCal = false;
    int modelLen = -1;
    std::vector<UnaryRule> unaryRules;
    std::vector<BinaryRule> binaryRules;
    std::vector<PairRule> pairRules;
    std::vector<EpsRule> epsRules;
};

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
    ExpTail expLC;
    ExpTail expLI;
    ExpTail expGC;
    ExpTail expGI;
    std::vector<CmState> states;
    int rootState;
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
    // Per-state ragged (v, i, d) chart for Infernal-compat exactVitRec.
    // Each state v gets bandSize[v] * (N+1) cells, allocated contiguously
    // and indexed via cumulative chartOffset[v]. Total cells = chartOffset[M].
    // Bands derive from cmbuild's QDB output (CmState::dmin1/dmax1) and shrink
    // memory ~30x vs the previous flat M*N*maxSpan layout.
    std::vector<float> exactChart;
    std::vector<uint32_t> exactChartSeen;
    uint32_t exactChartGen = 1;
    std::vector<size_t> exactChartOffset;   // M+1 entries, exactChartOffset[M] = total cells
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

static bool cmRefineTraceEnabled() {
    static int cached = -1;
    if (cached == -1) {
        const char *env = std::getenv("MMSEQS_CMSCAN_REFINE");
        if (env != NULL) {
            cached = (std::string(env) == "0") ? 0 : 1;
        } else {
            // Skip refine traceback by default — no quality benefit, 2-3x cost.
            cached = 0;
        }
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

static inline double nativeCalibratedScoreToEvalue(double scoreBits,
                                                   double targetDbResidues) {
    // Native fallback calibration for uncalibrated grammars:
    // E ~= N * 2^(-S), where N is effective target DB residues and S is bit score.
    const double n = std::max(1.0, targetDbResidues);
    const double s = std::max(-1000.0, std::min(1000.0, scoreBits));
    const double ev = n * std::exp2(-s);
    return std::max(1e-300, std::min(1e300, ev));
}

static int getNtId(Grammar &g, const std::string &nt) {
    std::map<std::string, int>::iterator it = g.ntToId.find(nt);
    if (it != g.ntToId.end()) {
        return it->second;
    }
    int id = static_cast<int>(g.idToNt.size());
    g.ntToId[nt] = id;
    g.idToNt.push_back(nt);
    return id;
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

static InfernalExactModel parseInfernalCmExactModel(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        Debug(Debug::ERROR) << "Cannot open Infernal CM file: " << path << "\n";
        EXIT(EXIT_FAILURE);
    }

    int clen = -1;
    int statesN = -1;
    int w = -1;
    bool inModel = false;
    bool done = false;
    int curMapL = -1;
    int curMapR = -1;
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
            (void) nidx;
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
        Debug(Debug::ERROR) << "Failed parsing Infernal CM state graph from " << path << "\n";
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
    Debug(Debug::INFO) << "Infernal CM exact parser: CLEN=" << m.clen
                       << " W=" << m.w
                       << " states=" << m.states.size() << "\n";
    return m;
}

// CM grammar format (plain text):
// START <NT>
// UNARY <LHS> <BASE> <LOGP>
// BINARY <LHS> <RHS1> <RHS2> <LOGP>
// PAIR <LHS> <INNER|- > <LEFT_BASE> <RIGHT_BASE> <LOGP>
// EPS <LHS> <RHS> <LOGP>
static Grammar parseGrammar(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        Debug(Debug::ERROR) << "Cannot open CM grammar file: " << path << "\n";
        EXIT(EXIT_FAILURE);
    }

    Grammar g;
    std::string line;
    size_t lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            std::stringstream cs(line.substr(1));
            std::string k;
            cs >> k;
            if (k == "CAL_A") {
                double v = 1.0;
                if (cs >> v) {
                    g.calA = v;
                }
            } else if (k == "CAL_B") {
                double v = 0.0;
                if (cs >> v) {
                    g.calB = v;
                }
            } else if (k == "CAL_SPLIT") {
                double v = 0.0;
                if (cs >> v) {
                    g.calSplit = v;
                }
            } else if (k == "CAL_NEG_A") {
                double v = 1.0;
                if (cs >> v) {
                    g.calNegA = v;
                    g.hasPiecewiseCal = true;
                }
            } else if (k == "CAL_NEG_B") {
                double v = 0.0;
                if (cs >> v) {
                    g.calNegB = v;
                    g.hasPiecewiseCal = true;
                }
            } else if (k == "model_len") {
                int v = -1;
                if (cs >> v) {
                    g.modelLen = v;
                }
            }
            continue;
        }
        std::stringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "START") {
            std::string s;
            ss >> s;
            if (s.empty()) {
                Debug(Debug::ERROR) << "Invalid START at " << path << ":" << lineNo << "\n";
                EXIT(EXIT_FAILURE);
            }
            g.start = getNtId(g, s);
        } else if (tag == "UNARY") {
            std::string lhs;
            std::string b;
            double logp = 0.0;
            ss >> lhs >> b >> logp;
            if (lhs.empty() || b.size() != 1) {
                Debug(Debug::ERROR) << "Invalid UNARY at " << path << ":" << lineNo << "\n";
                EXIT(EXIT_FAILURE);
            }
            UnaryRule r;
            r.lhs = getNtId(g, lhs);
            r.base = normalizeBase(b[0]);
            r.logp = logp;
            g.unaryRules.push_back(r);
        } else if (tag == "BINARY") {
            std::string lhs, r1, r2;
            double logp = 0.0;
            ss >> lhs >> r1 >> r2 >> logp;
            if (lhs.empty() || r1.empty() || r2.empty()) {
                Debug(Debug::ERROR) << "Invalid BINARY at " << path << ":" << lineNo << "\n";
                EXIT(EXIT_FAILURE);
            }
            BinaryRule r;
            r.lhs = getNtId(g, lhs);
            r.rhs1 = getNtId(g, r1);
            r.rhs2 = getNtId(g, r2);
            r.logp = logp;
            g.binaryRules.push_back(r);
        } else if (tag == "PAIR") {
            std::string lhs, inner, lb, rb;
            double logp = 0.0;
            ss >> lhs >> inner >> lb >> rb >> logp;
            if (lhs.empty() || inner.empty() || lb.size() != 1 || rb.size() != 1) {
                Debug(Debug::ERROR) << "Invalid PAIR at " << path << ":" << lineNo << "\n";
                EXIT(EXIT_FAILURE);
            }
            PairRule r;
            r.lhs = getNtId(g, lhs);
            r.inner = (inner == "-") ? -1 : getNtId(g, inner);
            r.leftBase = normalizeBase(lb[0]);
            r.rightBase = normalizeBase(rb[0]);
            r.logp = logp;
            g.pairRules.push_back(r);
        } else if (tag == "EPS") {
            std::string lhs, rhs;
            double logp = 0.0;
            ss >> lhs >> rhs >> logp;
            if (lhs.empty() || rhs.empty()) {
                Debug(Debug::ERROR) << "Invalid EPS at " << path << ":" << lineNo << "\n";
                EXIT(EXIT_FAILURE);
            }
            EpsRule r;
            r.lhs = getNtId(g, lhs);
            r.rhs = getNtId(g, rhs);
            r.logp = logp;
            g.epsRules.push_back(r);
        } else {
            Debug(Debug::ERROR) << "Unknown grammar directive '" << tag << "' at " << path << ":" << lineNo << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    if (g.start < 0) {
        Debug(Debug::ERROR) << "Grammar has no START directive: " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    return g;
}

static std::vector<FastaSeq> readFasta(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        Debug(Debug::ERROR) << "Cannot open FASTA file: " << path << "\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<FastaSeq> seqs;
    FastaSeq cur;
    unsigned int nextKey = 1;
    std::set<unsigned int> usedKeys;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!cur.id.empty()) {
                seqs.push_back(cur);
                cur = FastaSeq();
            }
            cur.id = trim(line.substr(1));
            // If FASTA header is numeric, preserve it as target key; otherwise assign a stable local key.
            std::string idTok = cur.id;
            const size_t ws = idTok.find_first_of(" \t");
            if (ws != std::string::npos) {
                idTok = idTok.substr(0, ws);
            }
            char *endptr = NULL;
            const unsigned long parsed = std::strtoul(idTok.c_str(), &endptr, 10);
            if (endptr != idTok.c_str() && *endptr == '\0' && parsed <= static_cast<unsigned long>(UINT_MAX)) {
                cur.key = static_cast<unsigned int>(parsed);
            } else {
                cur.key = nextKey;
            }
            while (usedKeys.find(cur.key) != usedKeys.end()) {
                cur.key = nextKey;
                nextKey += 1;
            }
            usedKeys.insert(cur.key);
            nextKey += 1;
            continue;
        }
        std::string s = trim(line);
        for (size_t i = 0; i < s.size(); ++i) {
            char c = normalizeBase(s[i]);
            if (c == 'A' || c == 'C' || c == 'G' || c == 'U' || c == 'N') {
                cur.seq.push_back(c);
            } else {
                // Preserve sequence length for ambiguous IUPAC symbols by mapping to N.
                cur.seq.push_back('N');
            }
        }
    }
    if (!cur.id.empty()) {
        seqs.push_back(cur);
    }
    return seqs;
}

static bool hasDbIndex(const std::string &path) {
    return FileUtil::fileExists((path + ".index").c_str());
}

struct RegionCoord {
    int dbStart;  // 0-based start in target
    int dbEnd;    // 0-based end (inclusive) in target
    int qStart;   // 0-based start in query
    int qEnd;     // 0-based end (inclusive) in query
};

struct CandidateResult {
    std::unordered_map<unsigned int, std::unordered_set<unsigned int> > candidates;
    // (queryKey, targetKey) → region coordinate from prefilter alignment
    std::unordered_map<uint64_t, RegionCoord> regionCoords;

    static uint64_t coordKey(unsigned int qKey, unsigned int tKey) {
        return (static_cast<uint64_t>(qKey) << 32) | static_cast<uint64_t>(tKey);
    }
};

static CandidateResult
loadResultDbCandidates(const std::string &resultDbPath, int threads) {
    CandidateResult result;
    const std::string idx = resultDbPath + ".index";
    // Data may be in resultDbPath or resultDbPath.0 (split DB format)
    const std::string splitPath = resultDbPath + ".0";
    if ((!FileUtil::fileExists(resultDbPath.c_str()) && !FileUtil::fileExists(splitPath.c_str()))
        || !FileUtil::fileExists(idx.c_str())) {
        return result;
    }
    DBReader<unsigned int> rdb(resultDbPath.c_str(),
                               idx.c_str(),
                               threads > 0 ? threads : 1,
                               DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    rdb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    result.candidates.reserve(rdb.getSize() * 2 + 1);
    for (size_t i = 0; i < rdb.getSize(); ++i) {
        const unsigned int qKey = rdb.getDbKey(i);
        std::unordered_set<unsigned int> &cand = result.candidates[qKey];
        char *data = rdb.getData(i, 0);
        while (*data != '\0') {
            while (*data == ' ' || *data == '\t') {
                ++data;
            }
            // Try to parse as alignment result (tab-separated, ≥10 columns)
            // Format: targetKey score seqId evalue qStart qEnd qLen dbStart dbEnd dbLen [backtrace]
            const char *lineStart = data;
            char *endptr = NULL;
            const unsigned long k = std::strtoul(data, &endptr, 10);
            if (endptr != data && k <= static_cast<unsigned long>(UINT_MAX)) {
                const unsigned int tKey = static_cast<unsigned int>(k);
                cand.insert(tKey);
                // Check if this line has alignment fields (tab after first number)
                if (*endptr == '\t') {
                    // Count tabs to find dbStart (col 7) and dbEnd (col 8)
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
                    if (col >= 10 && colStart[4] != NULL && colStart[5] != NULL
                        && colStart[7] != NULL && colStart[8] != NULL) {
                        uint64_t ck = CandidateResult::coordKey(qKey, tKey);
                        // Only store the first (best-scoring) hit per (query, target)
                        if (result.regionCoords.find(ck) == result.regionCoords.end()) {
                            int dbStart = Util::fast_atoi<int>(colStart[7]);
                            int dbEnd = Util::fast_atoi<int>(colStart[8]);
                            if (dbStart > dbEnd) {
                                std::swap(dbStart, dbEnd);
                            }
                            int qStart = Util::fast_atoi<int>(colStart[4]);
                            int qEnd = Util::fast_atoi<int>(colStart[5]);
                            if (qStart > qEnd) {
                                std::swap(qStart, qEnd);
                            }
                            RegionCoord rc;
                            rc.dbStart = dbStart;
                            rc.dbEnd = dbEnd;
                            rc.qStart = qStart;
                            rc.qEnd = qEnd;
                            result.regionCoords[ck] = rc;
                        }
                    }
                }
            }
            data = Util::skipLine(data);
        }
    }
    rdb.close();
    return result;
}

// Decode a single sequence from an open DBReader by internal index.
static FastaSeq decodeOneSequence(DBReader<unsigned int> &dbr,
                                  size_t id,
                                  BaseMatrix &subMat,
                                  Sequence &seqObj) {
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
    const char *data = dbr.getData(id, 0);
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

static std::string writeTempTextFile(const std::string &text) {
    char templ[] = "/tmp/mmseqs_cm_XXXXXX";
    const int fd = mkstemp(templ);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot create temporary file for CM parsing\n";
        EXIT(EXIT_FAILURE);
    }
    size_t off = 0;
    while (off < text.size()) {
        const ssize_t n = ::write(fd, text.data() + off, text.size() - off);
        if (n <= 0) {
            ::close(fd);
            Debug(Debug::ERROR) << "Cannot write temporary CM parse file\n";
            EXIT(EXIT_FAILURE);
        }
        off += static_cast<size_t>(n);
    }
    ::close(fd);
    return std::string(templ);
}

static bool looksLikeInfernalCmText(const std::string &txt) {
    const size_t eol = txt.find('\n');
    const std::string first = trim(txt.substr(0, (eol == std::string::npos) ? txt.size() : eol));
    return first.find("INFERNAL1/") != std::string::npos;
}

struct Hit {
    unsigned int dbKey = 0;
    unsigned int dbLen = 0;
    size_t targetIdx = 0;
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
    std::string cigar;
    std::string traceStates;
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

static std::string modelTraceCigar(const InfernalExactModel &model, const std::vector<int> &traceStates) {
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
    std::vector<char> mop(static_cast<size_t>(maxMap + 1), '\0'); // raw map coordinate space
    std::vector<int> insAfter(static_cast<size_t>(maxMap + 1), 0);
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
        } else if (s.type == CM_ST_ML) {
            if (s.mapLeft > 0 && s.mapLeft <= maxMap) {
                mop[static_cast<size_t>(s.mapLeft)] = 'M';
            }
        } else if (s.type == CM_ST_MR) {
            if (s.mapRight > 0 && s.mapRight <= maxMap) {
                mop[static_cast<size_t>(s.mapRight)] = 'M';
            }
        } else if (s.type == CM_ST_D) {
            if (s.mapLeft > 0 && s.mapLeft <= maxMap && mop[static_cast<size_t>(s.mapLeft)] == '\0') {
                mop[static_cast<size_t>(s.mapLeft)] = 'D';
            }
            if (s.mapRight > 0 && s.mapRight <= maxMap && mop[static_cast<size_t>(s.mapRight)] == '\0') {
                mop[static_cast<size_t>(s.mapRight)] = 'D';
            }
        } else if (s.type == CM_ST_IL || s.type == CM_ST_IR) {
            int anchor = (s.mapLeft > 0) ? s.mapLeft : s.mapRight;
            if (anchor > 0 && anchor <= maxMap) {
                insAfter[static_cast<size_t>(anchor)] += 1;
            }
        }
    }
    std::string ops;
    ops.reserve(static_cast<size_t>((maxMap - minMap + 1) + 32));
    for (int p = minMap; p <= maxMap; ++p) {
        if (mop[static_cast<size_t>(p)] == '\0') {
            mop[static_cast<size_t>(p)] = 'D'; // unset consensus position = gap in target
        }
        const char o = (mop[static_cast<size_t>(p)] == 'D') ? 'I' : 'M';
        ops.push_back(o);
        const int ins = insAfter[static_cast<size_t>(p)];
        for (int k = 0; k < ins; ++k) {
            ops.push_back('D');
        }
    }
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
            const size_t step = std::min(cnt, consensus.size() - std::min(iCon, consensus.size()));
            iCon += step;
            denom += step;
        } else if (op == 'D') {
            const size_t step = std::min(cnt, observed.size() - std::min(iObs, observed.size()));
            iObs += step;
            denom += step;
        }
    }
    if (denom == 0) {
        return 0.0f;
    }
    return static_cast<float>(ident) / static_cast<float>(denom);
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
    const size_t *chartOffset;  // M+1 entries; chartOffset[M] = total cells
    const int *chartBandDmin;   // M entries
    const int *chartBandSize;   // M entries (= dmax-dmin+1, 0 if empty)
};

static inline size_t exactChartIndex(const ExactRecCtx &ctx, int v, int i, int d) {
    return ctx.chartOffset[v]
         + static_cast<size_t>(i) * static_cast<size_t>(ctx.chartBandSize[v])
         + static_cast<size_t>(d - ctx.chartBandDmin[v]);
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
    if (ctx.chartSeen[dk] == ctx.chartGen) {
        return static_cast<double>(ctx.chartScore[dk]);
    }

    const ExactStateExec &st = ctx.execData[static_cast<size_t>(v)];
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
            if (std::fabs(best - cand) < 1e-6) {
                exactTraceRec(ctx, y, i, k, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
                exactTraceRec(ctx, z, k + 1, j, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
                return;
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
    int bestY = -1;
    double bestApprox = NEG_INF;
    if (st.trCount == 1) {
        const int y = static_cast<int>(st.trDst4[0]);
        const double nxt = exactVitRec(ctx, y, ni, nj);
        if (nxt != NEG_INF) {
            bestApprox = e + st.trSc4[0] + nxt;
            bestY = y;
        }
        if (nxt != NEG_INF && std::fabs(best - (e + st.trSc4[0] + nxt)) < 1e-6) {
            exactTraceRec(ctx, y, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            return;
        }
    } else if (st.trCount == 2) {
        const int y0 = static_cast<int>(st.trDst4[0]);
        const double nxt0 = exactVitRec(ctx, y0, ni, nj);
        if (nxt0 != NEG_INF) {
            const double cand0 = e + st.trSc4[0] + nxt0;
            if (cand0 > bestApprox) { bestApprox = cand0; bestY = y0; }
        }
        if (nxt0 != NEG_INF && std::fabs(best - (e + st.trSc4[0] + nxt0)) < 1e-6) {
            exactTraceRec(ctx, y0, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            return;
        }
        const int y1 = static_cast<int>(st.trDst4[1]);
        const double nxt1 = exactVitRec(ctx, y1, ni, nj);
        if (nxt1 != NEG_INF) {
            const double cand1 = e + st.trSc4[1] + nxt1;
            if (cand1 > bestApprox) { bestApprox = cand1; bestY = y1; }
        }
        if (nxt1 != NEG_INF && std::fabs(best - (e + st.trSc4[1] + nxt1)) < 1e-6) {
            exactTraceRec(ctx, y1, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            return;
        }
    } else if (st.trCount == 4) {
        const int y0 = static_cast<int>(st.trDst4[0]);
        const double nxt0 = exactVitRec(ctx, y0, ni, nj);
        if (nxt0 != NEG_INF) {
            const double cand0 = e + st.trSc4[0] + nxt0;
            if (cand0 > bestApprox) { bestApprox = cand0; bestY = y0; }
        }
        if (nxt0 != NEG_INF && std::fabs(best - (e + st.trSc4[0] + nxt0)) < 1e-6) {
            exactTraceRec(ctx, y0, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            return;
        }
        const int y1 = static_cast<int>(st.trDst4[1]);
        const double nxt1 = exactVitRec(ctx, y1, ni, nj);
        if (nxt1 != NEG_INF) {
            const double cand1 = e + st.trSc4[1] + nxt1;
            if (cand1 > bestApprox) { bestApprox = cand1; bestY = y1; }
        }
        if (nxt1 != NEG_INF && std::fabs(best - (e + st.trSc4[1] + nxt1)) < 1e-6) {
            exactTraceRec(ctx, y1, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            return;
        }
        const int y2 = static_cast<int>(st.trDst4[2]);
        const double nxt2 = exactVitRec(ctx, y2, ni, nj);
        if (nxt2 != NEG_INF) {
            const double cand2 = e + st.trSc4[2] + nxt2;
            if (cand2 > bestApprox) { bestApprox = cand2; bestY = y2; }
        }
        if (nxt2 != NEG_INF && std::fabs(best - (e + st.trSc4[2] + nxt2)) < 1e-6) {
            exactTraceRec(ctx, y2, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            return;
        }
        const int y3 = static_cast<int>(st.trDst4[3]);
        const double nxt3 = exactVitRec(ctx, y3, ni, nj);
        if (nxt3 != NEG_INF) {
            const double cand3 = e + st.trSc4[3] + nxt3;
            if (cand3 > bestApprox) { bestApprox = cand3; bestY = y3; }
        }
        if (nxt3 != NEG_INF && std::fabs(best - (e + st.trSc4[3] + nxt3)) < 1e-6) {
            exactTraceRec(ctx, y3, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
            return;
        }
    } else {
        for (int t = 0; t < st.trCount; ++t) {
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
            if (std::fabs(best - cand) < 1e-6) {
                exactTraceRec(ctx, y, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
                return;
            }
        }
    }
    if (bestY >= 0) {
        exactTraceRec(ctx, bestY, ni, nj, minUsed, maxUsed, obsCount, modelAggRaw, traceOps, traceStates);
    }
}

struct RefineRecCtx {
    int N;
    int M;
    const std::vector<ExactStateExec> *exec;
    const StateId *trDst;
    const double *trSc;
    const std::vector<int8_t> *seqCode;
    const int *bSplitBegByVD;
    const int *bSplitEndByVD;
    bool useDense;
    std::vector<double> *memoDense;
    std::vector<uint32_t> *memoSeen;
    uint32_t memoGen;
    std::unordered_map<int64_t, double> *memoVit;
};

static inline int64_t refineRecKey(const RefineRecCtx &ctx, int v, int i, int j) {
    return (static_cast<int64_t>(v) * static_cast<int64_t>(ctx.N + 2) + i)
           * static_cast<int64_t>(ctx.N + 2) + j;
}

static double refineVitRec(RefineRecCtx &ctx, int v, int i, int j) {
    if (v < 0 || v >= ctx.M || i < 1 || i > ctx.N + 1 || j < 0 || j > ctx.N || i > j + 1) {
        return NEG_INF;
    }
    const int d = j - i + 1;
    const ExactStateExec &sv = (*ctx.exec)[static_cast<size_t>(v)];
    if (sv.dmax >= sv.dmin && (d < sv.dmin || d > sv.dmax)) {
        return NEG_INF;
    }
    const int64_t key = refineRecKey(ctx, v, i, j);
    if (ctx.useDense) {
        const size_t dk = static_cast<size_t>(key);
        if ((*ctx.memoSeen)[dk] == ctx.memoGen) {
            return (*ctx.memoDense)[dk];
        }
    } else {
        std::unordered_map<int64_t, double>::const_iterator it = ctx.memoVit->find(key);
        if (it != ctx.memoVit->end()) {
            return it->second;
        }
    }
    double rv = NEG_INF;
    if (sv.type == CM_ST_E) {
        rv = (i == j + 1) ? 0.0 : NEG_INF;
    } else if (sv.type == CM_ST_B) {
        const int y = sv.bLeft;
        const int z = sv.bRight;
        if (y >= 0 && y < ctx.M && z >= 0 && z < ctx.M) {
            const size_t vd = static_cast<size_t>(v) * static_cast<size_t>(ctx.N + 1) + static_cast<size_t>(d);
            const int kBeg = ctx.bSplitBegByVD[vd];
            const int kEnd = ctx.bSplitEndByVD[vd];
            for (int k = kBeg; k <= kEnd; ++k) {
                const double l = refineVitRec(ctx, y, i, i + k - 1);
                const double r = refineVitRec(ctx, z, i + k, j);
                if (l != NEG_INF && r != NEG_INF) {
                    rv = std::max(rv, l + r);
                }
            }
        }
    } else if (d >= sv.dConsume) {
        const double e = cmStateEmitScoreFast(sv, *ctx.seqCode, i, j);
        if (e != NEG_INF) {
            const int ni = i + sv.niShift;
            const int nj = j - (sv.dConsume - sv.niShift);
            if (sv.trCount == 1) {
                const double nxt = refineVitRec(ctx, static_cast<int>(sv.trDst4[0]), ni, nj);
                if (nxt != NEG_INF) {
                    rv = std::max(rv, e + sv.trSc4[0] + nxt);
                }
            } else if (sv.trCount == 2) {
                const double nxt0 = refineVitRec(ctx, static_cast<int>(sv.trDst4[0]), ni, nj);
                if (nxt0 != NEG_INF) {
                    rv = std::max(rv, e + sv.trSc4[0] + nxt0);
                }
                const double nxt1 = refineVitRec(ctx, static_cast<int>(sv.trDst4[1]), ni, nj);
                if (nxt1 != NEG_INF) {
                    rv = std::max(rv, e + sv.trSc4[1] + nxt1);
                }
            } else if (sv.trCount == 4) {
                const double nxt0 = refineVitRec(ctx, static_cast<int>(sv.trDst4[0]), ni, nj);
                if (nxt0 != NEG_INF) {
                    rv = std::max(rv, e + sv.trSc4[0] + nxt0);
                }
                const double nxt1 = refineVitRec(ctx, static_cast<int>(sv.trDst4[1]), ni, nj);
                if (nxt1 != NEG_INF) {
                    rv = std::max(rv, e + sv.trSc4[1] + nxt1);
                }
                const double nxt2 = refineVitRec(ctx, static_cast<int>(sv.trDst4[2]), ni, nj);
                if (nxt2 != NEG_INF) {
                    rv = std::max(rv, e + sv.trSc4[2] + nxt2);
                }
                const double nxt3 = refineVitRec(ctx, static_cast<int>(sv.trDst4[3]), ni, nj);
                if (nxt3 != NEG_INF) {
                    rv = std::max(rv, e + sv.trSc4[3] + nxt3);
                }
            } else {
                for (int t = 0; t < sv.trCount; ++t) {
                    const size_t ti = sv.trOff + static_cast<size_t>(t);
                    const double nxt = refineVitRec(ctx, ctx.trDst[ti], ni, nj);
                    if (nxt != NEG_INF) {
                        rv = std::max(rv, e + ctx.trSc[ti] + nxt);
                    }
                }
            }
        }
    }
    if (ctx.useDense) {
        const size_t dk = static_cast<size_t>(key);
        (*ctx.memoSeen)[dk] = ctx.memoGen;
        (*ctx.memoDense)[dk] = rv;
    } else {
        (*ctx.memoVit)[key] = rv;
    }
    return rv;
}

static void refineTraceRec(RefineRecCtx &ctx,
                           int v,
                           int i,
                           int j,
                           std::vector<int> &traceStatesRefined) {
    const double best = refineVitRec(ctx, v, i, j);
    if (best == NEG_INF) {
        return;
    }
    const ExactStateExec &sv = (*ctx.exec)[static_cast<size_t>(v)];
    if (sv.type == CM_ST_E) {
        return;
    }
    traceStatesRefined.push_back(v);
    if (sv.type == CM_ST_B) {
        const int d = j - i + 1;
        const int y = sv.bLeft;
        const int z = sv.bRight;
        const size_t vd = static_cast<size_t>(v) * static_cast<size_t>(ctx.N + 1) + static_cast<size_t>(d);
        const int kBeg = ctx.bSplitBegByVD[vd];
        const int kEnd = ctx.bSplitEndByVD[vd];
        int bestK = -1;
        double bestCand = NEG_INF;
        for (int k = kBeg; k <= kEnd; ++k) {
            const double l = refineVitRec(ctx, y, i, i + k - 1);
            const double r = refineVitRec(ctx, z, i + k, j);
            if (l == NEG_INF || r == NEG_INF) {
                continue;
            }
            const double cand = l + r;
            if (cand > bestCand) {
                bestCand = cand;
                bestK = k;
            }
            if (std::fabs(best - cand) < 1e-9) {
                refineTraceRec(ctx, y, i, i + k - 1, traceStatesRefined);
                refineTraceRec(ctx, z, i + k, j, traceStatesRefined);
                return;
            }
        }
        if (bestK >= 0) {
            refineTraceRec(ctx, y, i, i + bestK - 1, traceStatesRefined);
            refineTraceRec(ctx, z, i + bestK, j, traceStatesRefined);
        }
        return;
    }
    const int d = j - i + 1;
    if (d < sv.dConsume) {
        return;
    }
    const double e = cmStateEmitScoreFast(sv, *ctx.seqCode, i, j);
    if (e == NEG_INF) {
        return;
    }
    const int ni = i + sv.niShift;
    const int nj = j - (sv.dConsume - sv.niShift);
    int bestY = -1;
    double bestCand = NEG_INF;
    if (sv.trCount == 1) {
        const int y = static_cast<int>(sv.trDst4[0]);
        const double nxt = refineVitRec(ctx, y, ni, nj);
        if (nxt != NEG_INF) {
            const double cand = e + sv.trSc4[0] + nxt;
            bestCand = cand;
            bestY = y;
            if (std::fabs(best - cand) < 1e-9) {
                refineTraceRec(ctx, y, ni, nj, traceStatesRefined);
                return;
            }
        }
    } else if (sv.trCount == 2) {
        for (int t = 0; t < 2; ++t) {
            const int y = static_cast<int>(sv.trDst4[t]);
            const double nxt = refineVitRec(ctx, y, ni, nj);
            if (nxt == NEG_INF) {
                continue;
            }
            const double cand = e + sv.trSc4[t] + nxt;
            if (cand > bestCand) {
                bestCand = cand;
                bestY = y;
            }
            if (std::fabs(best - cand) < 1e-9) {
                refineTraceRec(ctx, y, ni, nj, traceStatesRefined);
                return;
            }
        }
    } else if (sv.trCount == 4) {
        for (int t = 0; t < 4; ++t) {
            const int y = static_cast<int>(sv.trDst4[t]);
            const double nxt = refineVitRec(ctx, y, ni, nj);
            if (nxt == NEG_INF) {
                continue;
            }
            const double cand = e + sv.trSc4[t] + nxt;
            if (cand > bestCand) {
                bestCand = cand;
                bestY = y;
            }
            if (std::fabs(best - cand) < 1e-9) {
                refineTraceRec(ctx, y, ni, nj, traceStatesRefined);
                return;
            }
        }
    } else {
        for (int t = 0; t < sv.trCount; ++t) {
            const size_t ti = sv.trOff + static_cast<size_t>(t);
            const int y = ctx.trDst[ti];
            const double nxt = refineVitRec(ctx, y, ni, nj);
            if (nxt == NEG_INF) {
                continue;
            }
            const double cand = e + ctx.trSc[ti] + nxt;
            if (cand > bestCand) {
                bestCand = cand;
                bestY = y;
            }
            if (std::fabs(best - cand) < 1e-9) {
                refineTraceRec(ctx, y, ni, nj, traceStatesRefined);
                return;
            }
        }
    }
    if (bestY >= 0) {
        refineTraceRec(ctx, bestY, ni, nj, traceStatesRefined);
    }
}

static void runInfernalExactScan(const InfernalExactModel &model,
                                 const std::string &seq,
                                 bool wantInside,
                                 std::vector<Hit> &outHits,
                                 const std::string &seqId,
                                 int maxHitLen = 0) {
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
    auto null3CorrByInterval = [&](int i, int j) -> double {
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
        e.dmin = st.dmin1;
        e.dmax = st.dmax1;
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
            if (e.bLeft >= 0 && e.bLeft < M) {
                const CmState &l = model.states[static_cast<size_t>(e.bLeft)];
                if (l.dmin1 >= 0 && l.dmax1 >= l.dmin1) {
                    e.splitKMin = l.dmin1;
                    e.splitKMax = l.dmax1;
                }
            }
            if (e.bRight >= 0 && e.bRight < M) {
                const CmState &r = model.states[static_cast<size_t>(e.bRight)];
                if (r.dmin1 >= 0 && r.dmax1 >= r.dmin1) {
                    e.splitRMin = r.dmin1;
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
    // Fast path: iterative deck-based CYK over full sequence.
    // Inside score (if requested) is computed only for the selected best interval.
    if (cmExactFastDeckEnabled() && N <= 1024) {
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
                        {
                            const float *tsrc = trSrcPtrs[0];
                            const float tsc = trScVals[0];
                            for (int ii = 0; ii < iMax; ++ii) {
                                bestBuf[ii] = outPtr[ii] + tsc + tsrc[ii + 1];
                            }
                        }
                        // Remaining transitions
                        for (int t = 1; t < trCountClamped; ++t) {
                            const float *tsrc = trSrcPtrs[t];
                            const float tsc = trScVals[t];
                            for (int ii = 0; ii < iMax; ++ii) {
                                bestBuf[ii] = std::fmaxf(bestBuf[ii], outPtr[ii] + tsc + tsrc[ii + 1]);
                            }
                        }
                        // Copy result back
                        std::copy(bestBuf, bestBuf + iMax, outPtr);
                    }
                }
            }
        }

        const int root = model.rootState;
        int maxSpan = N;
        if (model.w > 0) {
            maxSpan = std::min(maxSpan, model.w);
        }
        if (maxHitLen > 0) {
            maxSpan = std::min(maxSpan, maxHitLen);
        }
        // Infernal semantics evaluate all legal spans under state dmin/dmax/w constraints.
        int minSpan = 1;
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
        for (int d = minSpan; d <= maxSpan; ++d) {
            const int iMax = N - d + 1;
            for (int i = 1; i <= iMax; ++i) {
                const float rawSc = vit[rootBase + static_cast<size_t>(d) * iStride + static_cast<size_t>(i)];
                if (rawSc == NEG_INF_F) {
                    continue;
                }
                const int j = i + d - 1;
                const bool mayImproveJ = rawSc > bestSc;
                const bool mayImproveL = enableTruncModes && j == N && rawSc > bestL;
                const bool mayImproveR = enableTruncModes && i == 1 && rawSc > bestR;
                const bool mayImproveT = enableTruncModes && i == 1 && j == N && rawSc > bestT;
                if (!(mayImproveJ || mayImproveL || mayImproveR || mayImproveT)) {
                    continue;
                }
                const double corr = null3CorrByInterval(i, j);
                const float sc = static_cast<float>(static_cast<double>(rawSc) - corr);
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
        while (!st.empty()) {
            TbCell c = st.back();
            st.pop_back();
            const ExactStateExec &sv = exec[static_cast<size_t>(c.v)];
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
                bool found = false;
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
                    if (std::fabs(cur - (l + r)) < 1e-6) {
                        st.push_back(TbCell{z, c.i + k, c.d - k});
                        st.push_back(TbCell{y, c.i, k});
                        found = true;
                        break;
                    }
                }
                if (!found && bestK >= 0) {
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
            const size_t ndBase = static_cast<size_t>(nd) * iStride + static_cast<size_t>(ni);
            bool pushedNext = false;
            int bestY = -1;
            float bestCand = NEG_INF_F;
            if (sv.trCount == 1) {
                const int y = sv.trDst4[0];
                const float n = vit[stateBase[static_cast<size_t>(y)] + ndBase];
                if (n != NEG_INF_F) {
                    const float cand = static_cast<float>(e + sv.trSc4[0] + n);
                    bestCand = cand;
                    bestY = y;
                }
                if (n != NEG_INF_F && std::fabs(cur - (e + sv.trSc4[0] + n)) < 1e-6) {
                    st.push_back(TbCell{y, ni, nd});
                    pushedNext = true;
                }
            } else if (sv.trCount == 2) {
                const int y0 = sv.trDst4[0];
                const float n0 = vit[stateBase[static_cast<size_t>(y0)] + ndBase];
                if (n0 != NEG_INF_F) {
                    const float cand0 = static_cast<float>(e + sv.trSc4[0] + n0);
                    if (cand0 > bestCand) { bestCand = cand0; bestY = y0; }
                }
                if (n0 != NEG_INF_F && std::fabs(cur - (e + sv.trSc4[0] + n0)) < 1e-6) {
                    st.push_back(TbCell{y0, ni, nd});
                    pushedNext = true;
                } else {
                    const int y1 = sv.trDst4[1];
                    const float n1 = vit[stateBase[static_cast<size_t>(y1)] + ndBase];
                    if (n1 != NEG_INF_F) {
                        const float cand1 = static_cast<float>(e + sv.trSc4[1] + n1);
                        if (cand1 > bestCand) { bestCand = cand1; bestY = y1; }
                    }
                    if (n1 != NEG_INF_F && std::fabs(cur - (e + sv.trSc4[1] + n1)) < 1e-6) {
                        st.push_back(TbCell{y1, ni, nd});
                        pushedNext = true;
                    }
                }
            } else if (sv.trCount == 4) {
                const int y0 = sv.trDst4[0];
                const float n0 = vit[stateBase[static_cast<size_t>(y0)] + ndBase];
                if (n0 != NEG_INF_F) {
                    const float cand0 = static_cast<float>(e + sv.trSc4[0] + n0);
                    if (cand0 > bestCand) { bestCand = cand0; bestY = y0; }
                }
                if (n0 != NEG_INF_F && std::fabs(cur - (e + sv.trSc4[0] + n0)) < 1e-6) {
                    st.push_back(TbCell{y0, ni, nd});
                    pushedNext = true;
                } else {
                    const int y1 = sv.trDst4[1];
                    const float n1 = vit[stateBase[static_cast<size_t>(y1)] + ndBase];
                    if (n1 != NEG_INF_F) {
                        const float cand1 = static_cast<float>(e + sv.trSc4[1] + n1);
                        if (cand1 > bestCand) { bestCand = cand1; bestY = y1; }
                    }
                    if (n1 != NEG_INF_F && std::fabs(cur - (e + sv.trSc4[1] + n1)) < 1e-6) {
                        st.push_back(TbCell{y1, ni, nd});
                        pushedNext = true;
                    } else {
                        const int y2 = sv.trDst4[2];
                        const float n2 = vit[stateBase[static_cast<size_t>(y2)] + ndBase];
                        if (n2 != NEG_INF_F) {
                            const float cand2 = static_cast<float>(e + sv.trSc4[2] + n2);
                            if (cand2 > bestCand) { bestCand = cand2; bestY = y2; }
                        }
                        if (n2 != NEG_INF_F && std::fabs(cur - (e + sv.trSc4[2] + n2)) < 1e-6) {
                            st.push_back(TbCell{y2, ni, nd});
                            pushedNext = true;
                        } else {
                            const int y3 = sv.trDst4[3];
                            const float n3 = vit[stateBase[static_cast<size_t>(y3)] + ndBase];
                            if (n3 != NEG_INF_F) {
                                const float cand3 = static_cast<float>(e + sv.trSc4[3] + n3);
                                if (cand3 > bestCand) { bestCand = cand3; bestY = y3; }
                            }
                            if (n3 != NEG_INF_F && std::fabs(cur - (e + sv.trSc4[3] + n3)) < 1e-6) {
                                st.push_back(TbCell{y3, ni, nd});
                                pushedNext = true;
                            }
                        }
                    }
                }
            } else {
                for (int t = 0; t < sv.trCount; ++t) {
                    const size_t ti = sv.trOff + static_cast<size_t>(t);
                    const int y = trDst[ti];
                    const double tr = trSc[ti];
                    const float n = vit[stateBase[static_cast<size_t>(y)] + ndBase];
                    if (n == NEG_INF_F) {
                        continue;
                    }
                    const float cand = static_cast<float>(e + tr + n);
                    if (cand > bestCand) {
                        bestCand = cand;
                        bestY = y;
                    }
                    if (std::fabs(cur - (e + tr + n)) < 1e-6) {
                        st.push_back(TbCell{y, ni, nd});
                        pushedNext = true;
                        break;
                    }
                }
            }
            if (!pushedNext && bestY >= 0) {
                st.push_back(TbCell{bestY, ni, nd});
            }
        }
        std::vector<int> traceStatesRefined;
        if (cmRefineTraceEnabled()) {
            traceStatesRefined.reserve(traceStates.size() + 8);
            const size_t refineCells = static_cast<size_t>(M) * static_cast<size_t>(N + 2) * static_cast<size_t>(N + 2);
            const bool refineUseDense = (refineCells <= 12000000ULL); // ~96MB values + ~48MB seen
            std::vector<double> &memoRefDense = ws.memoVitDense;
            std::vector<uint32_t> &memoRefSeen = ws.memoVitSeen;
            uint32_t &memoRefGen = ws.memoVitGen;
            std::unordered_map<int64_t, double> &memoRefMap = ws.memoVitMap;
            uint32_t curRefGen = 0;
            if (refineUseDense) {
                if (memoRefDense.size() < refineCells) {
                    memoRefDense.resize(refineCells);
                }
                if (memoRefSeen.size() < refineCells) {
                    memoRefSeen.resize(refineCells, 0);
                }
                ++memoRefGen;
                if (memoRefGen == 0) {
                    std::fill(memoRefSeen.begin(), memoRefSeen.end(), 0u);
                    memoRefGen = 1;
                }
                curRefGen = memoRefGen;
                memoRefMap.clear();
            } else {
                memoRefDense.clear();
                memoRefSeen.clear();
                memoRefMap.clear();
                memoRefMap.reserve(static_cast<size_t>(M) * static_cast<size_t>(bestD + 8));
            }
            RefineRecCtx rctx;
            rctx.N = N;
            rctx.M = M;
            rctx.exec = &exec;
            rctx.trDst = trDst.empty() ? NULL : &trDst[0];
            rctx.trSc = trSc.empty() ? NULL : &trSc[0];
            rctx.seqCode = &ws.seqCode;
            rctx.bSplitBegByVD = bSplitBegByVD;
            rctx.bSplitEndByVD = bSplitEndByVD;
            rctx.useDense = refineUseDense;
            rctx.memoDense = &memoRefDense;
            rctx.memoSeen = &memoRefSeen;
            rctx.memoGen = curRefGen;
            rctx.memoVit = &memoRefMap;
            refineTraceRec(rctx, root, bestI, bestI + bestD - 1, traceStatesRefined);
        }
        const std::vector<int> &traceStatesForCigar = traceStatesRefined.empty() ? traceStates : traceStatesRefined;

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
        h.cigar = rleTraceOps(traceOps);
        h.traceStates = encodeTraceStates(traceStatesForCigar);
        h.cigar = modelTraceCigar(model, traceStatesForCigar);
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
                            if (ef == NEG_INF_F) { continue; }
                            const int nli = li + niShift;
                            bool has = false;
                            float acc = NEG_INF_F;
                            for (int t = 0; t < st.trCount; ++t) {
                                const size_t ti = st.trOff + static_cast<size_t>(t);
                                const double n = nxtSlice[static_cast<size_t>(trDstPtr[ti]) * vStride + static_cast<size_t>(nli)];
                                if (n != NEG_INF) {
                                    log2AccFastAddF(ef + static_cast<float>(trScPtr[ti]) + static_cast<float>(n), has, acc);
                                }
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
                            if (ef == NEG_INF_F) { continue; }
                            const int nli = li + niShift;
                            bool has = false;
                            double maxVal = NEG_INF;
                            double scaledSum = 0.0;
                            for (int t = 0; t < st.trCount; ++t) {
                                const size_t ti = st.trOff + static_cast<size_t>(t);
                                const double n = nxtSlice[static_cast<size_t>(trDstPtr[ti]) * vStride + static_cast<size_t>(nli)];
                                if (n != NEG_INF) {
                                    log2AccExactAdd(static_cast<double>(ef) + trScPtr[ti] + n, has, maxVal, scaledSum);
                                }
                            }
                            dst[li] = log2AccExactValue(has, maxVal, scaledSum);
                        }
                    }
                }
            }
            h.inside = insDataAll[static_cast<size_t>(bestD) * dStride + static_cast<size_t>(root) * vStride + 0]
                       - bestNull3Corr - null2Corr;
            h.bias = bestNull3Corr + null2Corr;
        } else {
            h.inside = NEG_INF;
            h.bias = bestNull3Corr + null2Corr;
        }
        outHits.push_back(h);
        return;
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

    // Per-state ragged chart: each state v gets bandSize[v]*(N+1) cells, with
    // bands derived from QDB (CmState::dmin1/dmax1 as copied into exec[v].dmin/dmax).
    // For big rRNA CMs (M ~ 9000) bands are typically tight (avg ~20) so the
    // ragged chart is ~30x smaller than a flat M*(N+1)*(maxSpan+1) chart.
    std::vector<size_t> &chartOffset = ws.exactChartOffset;
    std::vector<int> &bandDmin = ws.exactChartBandDmin;
    std::vector<int> &bandSize = ws.exactChartBandSize;
    chartOffset.assign(static_cast<size_t>(M) + 1, 0);
    bandDmin.assign(static_cast<size_t>(M), 0);
    bandSize.assign(static_cast<size_t>(M), 0);
    const size_t rowPlus1 = static_cast<size_t>(N + 1);
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
        chartOffset[v + 1] = chartOffset[v] +
            static_cast<size_t>(bandSize[v]) * rowPlus1;
    }
    size_t chartCells = chartOffset[M];

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
    recCtx.chartOffset = chartOffset.data();
    recCtx.chartBandDmin = bandDmin.data();
    recCtx.chartBandSize = bandSize.data();
    // Infernal semantics evaluate all legal spans under state dmin/dmax/w constraints.
    int minSpan = 1;
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
        for (int i = 1; i <= iMax; ++i) {
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
    h.cigar = rleTraceOps(traceOps);
    h.traceStates = encodeTraceStates(traceStates);
    h.cigar = modelTraceCigar(model, traceStates);
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
    } else {
        h.inside = NEG_INF;
        h.bias = bestNull3Corr + null2Corr;
    }
    outHits.push_back(h);
}

static void runDpScan(const Grammar &g,
                      const std::string &seq,
                      bool wantInside,
                      bool wantBacktrace,
                      std::vector<Hit> &outHits,
                      const std::string &seqId) {
    const int n = static_cast<int>(seq.size());
    const int m = static_cast<int>(g.idToNt.size());
    if (n == 0 || m == 0) {
        return;
    }

    // Banded CYK: cap span length at W to bound memory to O(m*n*W).
    // W mirrors Infernal's fallback "max emit length" (3*modelLen when known).
    int W = (g.modelLen > 0) ? std::min(n, 3 * g.modelLen) : n;
    if (W < 1) W = n;
    const int Wp1 = W + 1;
    const size_t cellCount = static_cast<size_t>(m) * static_cast<size_t>(n) * static_cast<size_t>(Wp1);
    std::vector<double> cyk(cellCount, NEG_INF);
    std::vector<double> inside;
    if (wantInside) {
        inside.assign(cellCount, NEG_INF);
    }

    const auto idx = [n, Wp1](int nt, int i, int j) -> size_t {
        return static_cast<size_t>(nt) * static_cast<size_t>(n) * static_cast<size_t>(Wp1)
             + static_cast<size_t>(i) * static_cast<size_t>(Wp1)
             + static_cast<size_t>(j - i + 1);
    };
    const auto scoreEq = [](double a, double b) -> bool {
        if (a == NEG_INF || b == NEG_INF) {
            return false;
        }
        const double d = std::fabs(a - b);
        const double s = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
        return d <= 1e-5 * s;
    };

    int mRoot = -1;
    int dRoot = -1;
    for (size_t nt = 0; nt < g.idToNt.size(); ++nt) {
        const std::string &name = g.idToNt[nt];
        if (name.find("M_1_") == 0) {
            mRoot = static_cast<int>(nt);
        } else if (name.find("D_1_") == 0) {
            dRoot = static_cast<int>(nt);
        }
    }

    const auto applyEpsilonClosure = [&](std::vector<double> &chart, int i, int j, bool insideMode) {
        if (g.epsRules.empty()) {
            return;
        }
        bool changed = true;
        int iter = 0;
        while (changed && iter < static_cast<int>(g.idToNt.size()) + 2) {
            changed = false;
            ++iter;
            for (size_t r = 0; r < g.epsRules.size(); ++r) {
                const EpsRule &er = g.epsRules[r];
                const size_t lhs = idx(er.lhs, i, j);
                const size_t rhs = idx(er.rhs, i, j);
                if (chart[rhs] == NEG_INF) {
                    continue;
                }
                const double cand = chart[rhs] + er.logp;
                if (!insideMode) {
                    if (cand > chart[lhs]) {
                        chart[lhs] = cand;
                        changed = true;
                    }
                } else {
                    const double old = chart[lhs];
                    chart[lhs] = logSumExpMaybeFast(chart[lhs], cand);
                    if (chart[lhs] != old) {
                        changed = true;
                    }
                }
            }
        }
    };

    for (int i = 0; i < n; ++i) {
        const char b = seq[i];
        for (size_t r = 0; r < g.unaryRules.size(); ++r) {
            const UnaryRule &ur = g.unaryRules[r];
            if (ur.base == b || b == 'N') {
                size_t at = idx(ur.lhs, i, i);
                cyk[at] = std::max(cyk[at], ur.logp);
                if (wantInside) {
                    inside[at] = logSumExpMaybeFast(inside[at], ur.logp);
                }
            }
        }
        applyEpsilonClosure(cyk, i, i, false);
        if (wantInside) {
            applyEpsilonClosure(inside, i, i, true);
        }
    }

    // Infernal-like local semantics: evaluate all spans up to the band.
    const int maxSpan = W;
    const int minSpan = 1;
    for (int len = 2; len <= maxSpan; ++len) {
        for (int i = 0; i + len - 1 < n; ++i) {
            const int j = i + len - 1;

            for (size_t br = 0; br < g.binaryRules.size(); ++br) {
                const BinaryRule &r = g.binaryRules[br];
                double best = NEG_INF;
                double sum = NEG_INF;
                for (int k = i; k < j; ++k) {
                    const double lCyk = cyk[idx(r.rhs1, i, k)];
                    const double rCyk = cyk[idx(r.rhs2, k + 1, j)];
                    if (lCyk != NEG_INF && rCyk != NEG_INF) {
                        best = std::max(best, lCyk + rCyk + r.logp);
                    }
                    if (wantInside) {
                        const double lIn = inside[idx(r.rhs1, i, k)];
                        const double rIn = inside[idx(r.rhs2, k + 1, j)];
                        if (lIn != NEG_INF && rIn != NEG_INF) {
                            sum = logSumExpMaybeFast(sum, lIn + rIn + r.logp);
                        }
                    }
                }
                size_t at = idx(r.lhs, i, j);
                cyk[at] = std::max(cyk[at], best);
                if (wantInside) {
                    inside[at] = logSumExpMaybeFast(inside[at], sum);
                }
            }

            for (size_t pr = 0; pr < g.pairRules.size(); ++pr) {
                const PairRule &r = g.pairRules[pr];
                const char lb = seq[i];
                const char rb = seq[j];
                if (!((lb == r.leftBase || lb == 'N') && (rb == r.rightBase || rb == 'N'))) {
                    continue;
                }
                double innerCyk = NEG_INF;
                double innerInside = NEG_INF;
                if (r.inner == -1) {
                    if (j == i + 1) {
                        innerCyk = 0.0;
                        innerInside = 0.0;
                    }
                } else if (i + 1 <= j - 1) {
                    innerCyk = cyk[idx(r.inner, i + 1, j - 1)];
                    if (wantInside) {
                        innerInside = inside[idx(r.inner, i + 1, j - 1)];
                    }
                }
                if (innerCyk != NEG_INF) {
                    const size_t at = idx(r.lhs, i, j);
                    cyk[at] = std::max(cyk[at], innerCyk + r.logp);
                    if (wantInside && innerInside != NEG_INF) {
                        inside[at] = logSumExpMaybeFast(inside[at], innerInside + r.logp);
                    }
                }
            }
            applyEpsilonClosure(cyk, i, j, false);
            if (wantInside) {
                applyEpsilonClosure(inside, i, j, true);
            }
        }
    }

    // Report best parse per sequence for native grammar mode.
    // This is closer to Infernal cmsearch top-hit behavior than emitting every finite window.
    double bestCyk = NEG_INF;
    double bestInside = NEG_INF;
    int bestI = -1;
    int bestJ = -1;
    for (int i = 0; i < n; ++i) {
        const int jMax = std::min(n - 1, i + maxSpan - 1);
        for (int j = i; j <= jMax; ++j) {
            const double cykSc = cyk[idx(g.start, i, j)];
            const double inSc = wantInside ? inside[idx(g.start, i, j)] : NEG_INF;
            if (cykSc == NEG_INF && inSc == NEG_INF) {
                continue;
            }
            const int len = j - i + 1;
            if (len < minSpan) {
                continue;
            }
            const bool strictlyBetter = (bestI < 0 || cykSc > bestCyk);
            if (strictlyBetter) {
                bestCyk = cykSc;
                bestInside = inSc;
                bestI = i;
                bestJ = j;
            }
        }
    }
    if (bestI >= 0) {
        if (cmDebugExtEnabled()) {
            double leftExt = NEG_INF;
            double rightExt = NEG_INF;
            const int bestLen = bestJ - bestI + 1;
            const bool extFits = (bestLen + 1) <= W;
            if (bestI > 0 && extFits) {
                leftExt = cyk[idx(g.start, bestI - 1, bestJ)];
            }
            if (bestJ + 1 < n && extFits) {
                rightExt = cyk[idx(g.start, bestI, bestJ + 1)];
            }
            const double mBest = (mRoot >= 0) ? cyk[idx(mRoot, bestI, bestJ)] : NEG_INF;
            const double dBest = (dRoot >= 0) ? cyk[idx(dRoot, bestI, bestJ)] : NEG_INF;
            const double mRight = (mRoot >= 0 && bestJ + 1 < n && extFits) ? cyk[idx(mRoot, bestI, bestJ + 1)] : NEG_INF;
            const double dRight = (dRoot >= 0 && bestJ + 1 < n && extFits) ? cyk[idx(dRoot, bestI, bestJ + 1)] : NEG_INF;
            const double mCoreR = (mRoot >= 0 && bestJ > bestI) ? cyk[idx(mRoot, bestI, bestJ - 1)] : NEG_INF;
            const double dCoreR = (dRoot >= 0 && bestJ > bestI) ? cyk[idx(dRoot, bestI, bestJ - 1)] : NEG_INF;
            const double mCoreL = (mRoot >= 0 && bestI + 1 <= bestJ) ? cyk[idx(mRoot, bestI + 1, bestJ)] : NEG_INF;
            const double dCoreL = (dRoot >= 0 && bestI + 1 <= bestJ) ? cyk[idx(dRoot, bestI + 1, bestJ)] : NEG_INF;
            Debug(Debug::INFO) << "cmsearch-ext-debug " << seqId
                               << " best=[" << (bestI + 1) << "," << (bestJ + 1) << "]"
                               << " score=" << bestCyk
                               << " leftExt=" << leftExt
                               << " rightExt=" << rightExt
                               << " dL=" << (leftExt == NEG_INF ? NEG_INF : (bestCyk - leftExt))
                               << " dR=" << (rightExt == NEG_INF ? NEG_INF : (bestCyk - rightExt))
                               << " mBest=" << mBest
                               << " dBest=" << dBest
                               << " mRight=" << mRight
                               << " dRight=" << dRight
                               << " mCoreR=" << mCoreR
                               << " dCoreR=" << dCoreR
                               << " mCoreL=" << mCoreL
                               << " dCoreL=" << dCoreL
                               << "\n";
        }
        Hit h;
        h.seqId = seqId;
        h.start1 = bestI + 1;
        h.end1 = bestJ + 1;
        h.cigar = "NA";
        h.traceStates = "NA";
        if (wantBacktrace) {
            std::string traceOps;
            traceOps.reserve(static_cast<size_t>(std::max(1, bestJ - bestI + 1)));
            std::unordered_set<uint64_t> epsStack;
            std::function<bool(int, int, int)> traceNt =
                [&](int nt, int i, int j) -> bool {
                    if (i < 0 || j < i || j >= n || nt < 0 || nt >= m) {
                        return false;
                    }
                    const double v = cyk[idx(nt, i, j)];
                    if (v == NEG_INF) {
                        return false;
                    }
                    if (i == j) {
                        const char b = seq[i];
                        for (size_t urIdx = 0; urIdx < g.unaryRules.size(); ++urIdx) {
                            const UnaryRule &ur = g.unaryRules[urIdx];
                            if (ur.lhs != nt) {
                                continue;
                            }
                            if (!(ur.base == b || b == 'N')) {
                                continue;
                            }
                            if (scoreEq(v, ur.logp)) {
                                traceOps.push_back('M');
                                return true;
                            }
                        }
                    }
                    for (size_t pr = 0; pr < g.pairRules.size(); ++pr) {
                        const PairRule &r = g.pairRules[pr];
                        if (r.lhs != nt) {
                            continue;
                        }
                        if (i >= j) {
                            continue;
                        }
                        const char lb = seq[i];
                        const char rb = seq[j];
                        if (!((lb == r.leftBase || lb == 'N') && (rb == r.rightBase || rb == 'N'))) {
                            continue;
                        }
                        if (r.inner == -1) {
                            if (j == i + 1 && scoreEq(v, r.logp)) {
                                traceOps.push_back('M');
                                traceOps.push_back('M');
                                return true;
                            }
                            continue;
                        }
                        if (i + 1 > j - 1) {
                            continue;
                        }
                        const double inV = cyk[idx(r.inner, i + 1, j - 1)];
                        if (inV == NEG_INF) {
                            continue;
                        }
                        if (!scoreEq(v, inV + r.logp)) {
                            continue;
                        }
                        const size_t save = traceOps.size();
                        traceOps.push_back('M');
                        if (traceNt(r.inner, i + 1, j - 1)) {
                            traceOps.push_back('M');
                            return true;
                        }
                        traceOps.resize(save);
                    }
                    for (size_t br = 0; br < g.binaryRules.size(); ++br) {
                        const BinaryRule &r = g.binaryRules[br];
                        if (r.lhs != nt) {
                            continue;
                        }
                        for (int k = i; k < j; ++k) {
                            const double lv = cyk[idx(r.rhs1, i, k)];
                            const double rv = cyk[idx(r.rhs2, k + 1, j)];
                            if (lv == NEG_INF || rv == NEG_INF) {
                                continue;
                            }
                            if (!scoreEq(v, lv + rv + r.logp)) {
                                continue;
                            }
                            const size_t save = traceOps.size();
                            if (traceNt(r.rhs1, i, k) && traceNt(r.rhs2, k + 1, j)) {
                                return true;
                            }
                            traceOps.resize(save);
                        }
                    }
                    const uint64_t stackKey = (static_cast<uint64_t>(nt) << 42)
                                            | (static_cast<uint64_t>(i) << 21)
                                            | static_cast<uint64_t>(j);
                    if (epsStack.find(stackKey) != epsStack.end()) {
                        return false;
                    }
                    epsStack.insert(stackKey);
                    for (size_t erIdx = 0; erIdx < g.epsRules.size(); ++erIdx) {
                        const EpsRule &er = g.epsRules[erIdx];
                        if (er.lhs != nt) {
                            continue;
                        }
                        const double rv = cyk[idx(er.rhs, i, j)];
                        if (rv == NEG_INF) {
                            continue;
                        }
                        if (!scoreEq(v, rv + er.logp)) {
                            continue;
                        }
                        const size_t save = traceOps.size();
                        traceOps.push_back('D');
                        if (traceNt(er.rhs, i, j)) {
                            epsStack.erase(stackKey);
                            return true;
                        }
                        traceOps.resize(save);
                    }
                    epsStack.erase(stackKey);
                    return false;
                };
            if (traceNt(g.start, bestI, bestJ)) {
                h.cigar = rleTraceOps(traceOps);
            } else {
                h.cigar = std::to_string(std::max(1, bestJ - bestI + 1)) + "M";
            }
        }
        const bool useNegCal = g.hasPiecewiseCal && bestCyk < g.calSplit;
        const double calA = useNegCal ? g.calNegA : g.calA;
        const double calB = useNegCal ? g.calNegB : g.calB;
        h.cyk = calA * bestCyk + calB;
        if (bestInside == NEG_INF) {
            h.inside = NEG_INF;
        } else {
            const bool useNegCalInside = g.hasPiecewiseCal && bestInside < g.calSplit;
            const double calAInside = useNegCalInside ? g.calNegA : g.calA;
            const double calBInside = useNegCalInside ? g.calNegB : g.calB;
            h.inside = calAInside * bestInside + calBInside;
        }
        outHits.push_back(h);
    }
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
        Grammar g;
        InfernalExactModel exactModel;
        bool useInfernalCompat;
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
                                              1,
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
        singleModel.key = 0;
        singleModel.useInfernalCompat = false;
        if (looksLikeInfernalCm(par.db1)) {
            Debug(Debug::INFO) << "Detected Infernal .cm file, loading internal exact state-graph model\n";
            singleModel.exactModel = parseInfernalCmExactModel(par.db1);
            singleModel.useInfernalCompat = true;
        } else {
            singleModel.g = parseGrammar(par.db1);
        }
        QueryModelRef ref;
        ref.key = 0;
        ref.dbIdx = SIZE_MAX;
        queryRefs.push_back(ref);
    }

    const std::string seqDbPath = par.db2;
    const std::string seqDbIndex = par.db2 + ".index";
    DBReader<unsigned int> seqDbr(seqDbPath.c_str(),
                                  seqDbIndex.c_str(),
                                  1,
                                  DBReader<unsigned int>::USE_DATA
                                      | DBReader<unsigned int>::USE_INDEX
                                      | DBReader<unsigned int>::USE_LOOKUP);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    NucleotideMatrix nucMat(Parameters::getInstance().scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    BaseMatrix &subMat = static_cast<BaseMatrix&>(nucMat);
    Sequence seqObj(seqDbr.getMaxSeqLen(), Parameters::DBTYPE_NUCLEOTIDES, &nucMat, 0, false, false);

    const size_t observedDbResidues = seqDbr.getAminoAcidDBSize();
    const double strandMult = (par.strand == 2) ? 2.0 : 1.0;
    const double targetDbResidues = (par.dbSize > 0)
        ? static_cast<double>(par.dbSize)
        : (static_cast<double>(observedDbResidues) * strandMult);
    const bool wantInsideUser = (par.cmMode == LocalParameters::CM_MODE_INSIDE);

    CandidateResult candidateResult = loadResultDbCandidates(par.db3, par.threads);
    std::unordered_map<unsigned int, std::unordered_set<unsigned int> > &candidatesByQuery = candidateResult.candidates;
    const float cmRegionFlanking = par.cmRegionFlanking;
    Debug(Debug::INFO) << subcmd << " output DB: " << par.db4 << "\n";
    DBWriter resultWriter(par.db4.c_str(),
                          par.db4Index.c_str(),
                          1,
                          par.compressed,
                          Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    const size_t nThreads = (par.threads > 0) ? static_cast<size_t>(par.threads) : 1u;
    char buffer[1024 + 32768 * 4];
    size_t totalHits = 0;
    for (size_t qi = 0; qi < queryRefs.size(); ++qi) {
        const QueryModelRef &ref = queryRefs[qi];

        // Load CM lazily: parse from DB on demand, or use preloaded single model
        QueryModel qm;
        if (ref.dbIdx != SIZE_MAX && cmReader != NULL) {
            qm.key = ref.key;
            qm.useInfernalCompat = false;
            const char *raw = cmReader->getData(ref.dbIdx, 0);
            size_t len = cmReader->getEntryLen(ref.dbIdx);
            std::string text(raw, len);
            const size_t nul = text.find('\0');
            if (nul != std::string::npos) {
                text.resize(nul);
            }
            const std::string tmpPath = writeTempTextFile(text);
            if (looksLikeInfernalCmText(text)) {
                qm.exactModel = parseInfernalCmExactModel(tmpPath);
                qm.useInfernalCompat = true;
            } else {
                qm.g = parseGrammar(tmpPath);
            }
            std::remove(tmpPath.c_str());
        } else {
            qm = singleModel;
        }

        std::unordered_map<unsigned int, std::unordered_set<unsigned int> >::const_iterator cq = candidatesByQuery.find(qm.key);
        if (cq == candidatesByQuery.end()) {
            resultWriter.writeData("", 0, qm.key, 0);
            continue;
        }

        // Decode only the candidate target sequences for this query
        std::vector<FastaSeq> targetSeqs;
        targetSeqs.reserve(cq->second.size());
        for (std::unordered_set<unsigned int>::const_iterator it = cq->second.begin(); it != cq->second.end(); ++it) {
            const size_t id = seqDbr.getId(*it);
            if (id == UINT_MAX) {
                continue;
            }
            targetSeqs.push_back(decodeOneSequence(seqDbr, id, subMat, seqObj));
        }
        if (targetSeqs.empty()) {
            resultWriter.writeData("", 0, qm.key, 0);
            continue;
        }

        std::vector<Hit> hits;
        hits.reserve(1024);
        const bool infernalDefaultPrimaryInside = qm.useInfernalCompat && !wantInsideUser;
        const bool wantInside = wantInsideUser || infernalDefaultPrimaryInside;
        if (wantInside && cmFastMathEnabled()) {
            Debug(Debug::WARNING) << "MMSEQS_CMSCAN_FASTMATH=1 enabled: using approximate log/exp in Inside DP (scores may drift)\n";
        }
        // Build per-target region offsets when --cm-region is active
        std::vector<int> regionOffsets(targetSeqs.size(), 0);
        std::vector<std::string> regionSeqs(targetSeqs.size());
        // Per-target maxHitLen: cap the CYK d-search by a multiple of the
        // rnasearch envelope length. A true CM hit has d within a small
        // multiple of the prefilter span, and the memoized Viterbi only
        // wastes memory on larger spans.
        std::vector<int> maxHitLens(targetSeqs.size(), 0);
        // Slack for the CYK d-cap. We take the max of:
        //   - 3x the rnasearch prefilter envelope length (local extension)
        //   - 1.5x the CM consensus length (global: a true CM hit can stretch
        //     well beyond the prefilter envelope, especially for structured
        //     RNAs where iter-3 msa-profile mass is concentrated around the
        //     center but the CM consensus spans the full rRNA).
        // model.w (Infernal's W) remains the hard upper bound on the CM side.
        static constexpr int CM_MAXSPAN_ENV_SLACK = 3;
        static constexpr double CM_MAXSPAN_CLEN_SLACK = 1.5;
        const int clenFloor = qm.useInfernalCompat
            ? static_cast<int>(qm.exactModel.clen * CM_MAXSPAN_CLEN_SLACK)
            : 0;
        for (size_t i = 0; i < targetSeqs.size(); ++i) {
            const std::string &fullSeq = targetSeqs[i].seq;
            uint64_t ck = CandidateResult::coordKey(qm.key, targetSeqs[i].key);
            std::unordered_map<uint64_t, RegionCoord>::const_iterator it = candidateResult.regionCoords.find(ck);
            if (cmRegionFlanking > 0.0f) {
                if (it != candidateResult.regionCoords.end()) {
                    const int hitStart = it->second.dbStart;
                    const int hitEnd = it->second.dbEnd;
                    const int qAlnLen = std::max(1, it->second.qEnd - it->second.qStart + 1);
                    const int flank = static_cast<int>(cmRegionFlanking * qAlnLen);
                    const int seqLen = static_cast<int>(fullSeq.size());
                    const int regStart = std::max(0, hitStart - flank);
                    const int regEnd = std::min(seqLen, hitEnd + flank + 1);
                    regionSeqs[i] = fullSeq.substr(static_cast<size_t>(regStart),
                                                   static_cast<size_t>(regEnd - regStart));
                    regionOffsets[i] = regStart;
                } else {
                    regionSeqs[i] = fullSeq;
                }
            } else {
                regionSeqs[i] = fullSeq;
            }
            if (it != candidateResult.regionCoords.end()) {
                const int dbSpan = std::max(1, it->second.dbEnd - it->second.dbStart + 1);
                maxHitLens[i] = std::max(dbSpan * CM_MAXSPAN_ENV_SLACK, clenFloor);
            } else if (clenFloor > 0) {
                maxHitLens[i] = clenFloor;
            }
        }

        // Scan both forward and reverse complement strands
        std::vector<std::string> revSeqs(targetSeqs.size());
        for (size_t i = 0; i < targetSeqs.size(); ++i) {
            revSeqs[i] = reverseComplement(regionSeqs[i]);
        }

        std::vector< std::vector<Hit> > perSeqHits(targetSeqs.size());
        std::vector< std::vector<Hit> > perSeqHitsRev(targetSeqs.size());
        if (nThreads <= 1 || targetSeqs.size() <= 1) {
            for (size_t i = 0; i < targetSeqs.size(); ++i) {
                if (qm.useInfernalCompat) {
                    runInfernalExactScan(qm.exactModel, regionSeqs[i], wantInside, perSeqHits[i], targetSeqs[i].id, maxHitLens[i]);
                    runInfernalExactScan(qm.exactModel, revSeqs[i], wantInside, perSeqHitsRev[i], targetSeqs[i].id, maxHitLens[i]);
                } else {
                    runDpScan(qm.g, regionSeqs[i], wantInside, par.addBacktrace, perSeqHits[i], targetSeqs[i].id);
                    runDpScan(qm.g, revSeqs[i], wantInside, par.addBacktrace, perSeqHitsRev[i], targetSeqs[i].id);
                }
            }
        } else {
            #pragma omp parallel for schedule(dynamic, 1) num_threads(static_cast<int>(nThreads))
            for (long i = 0; i < static_cast<long>(targetSeqs.size()); ++i) {
                const size_t ui = static_cast<size_t>(i);
                if (qm.useInfernalCompat) {
                    runInfernalExactScan(qm.exactModel, regionSeqs[ui], wantInside, perSeqHits[ui], targetSeqs[ui].id, maxHitLens[ui]);
                    runInfernalExactScan(qm.exactModel, revSeqs[ui], wantInside, perSeqHitsRev[ui], targetSeqs[ui].id, maxHitLens[ui]);
                } else {
                    runDpScan(qm.g, regionSeqs[ui], wantInside, par.addBacktrace, perSeqHits[ui], targetSeqs[ui].id);
                    runDpScan(qm.g, revSeqs[ui], wantInside, par.addBacktrace, perSeqHitsRev[ui], targetSeqs[ui].id);
                }
            }
        }
        for (size_t i = 0; i < targetSeqs.size(); ++i) {
            const int offset = regionOffsets[i];
            const int seqLen = static_cast<int>(regionSeqs[i].size());
            // Forward strand hits
            for (size_t j = 0; j < perSeqHits[i].size(); ++j) {
                perSeqHits[i][j].dbKey = targetSeqs[i].key;
                perSeqHits[i][j].dbLen = static_cast<unsigned int>(targetSeqs[i].seq.size());
                perSeqHits[i][j].targetIdx = i;
                if (offset > 0) {
                    perSeqHits[i][j].start1 += offset;
                    perSeqHits[i][j].end1 += offset;
                }
            }
            hits.insert(hits.end(), perSeqHits[i].begin(), perSeqHits[i].end());
            // Reverse strand hits: map coordinates back to forward strand
            // revcomp pos p corresponds to forward pos (seqLen - p + 1)
            for (size_t j = 0; j < perSeqHitsRev[i].size(); ++j) {
                Hit &h = perSeqHitsRev[i][j];
                h.dbKey = targetSeqs[i].key;
                h.dbLen = static_cast<unsigned int>(targetSeqs[i].seq.size());
                h.targetIdx = i;
                int fwdStart = seqLen - h.start1 + 1 + offset;
                int fwdEnd = seqLen - h.end1 + 1 + offset;
                // Use start > end to indicate minus strand (Infernal convention)
                h.start1 = fwdStart;
                h.end1 = fwdEnd;
            }
            hits.insert(hits.end(), perSeqHitsRev[i].begin(), perSeqHitsRev[i].end());
        }
        if (infernalDefaultPrimaryInside) {
            for (size_t i = 0; i < hits.size(); ++i) {
                if (hits[i].inside != NEG_INF) {
                    hits[i].cyk = hits[i].inside;
                }
            }
        }
        if (qm.useInfernalCompat) {
            const bool evalueModeInside = infernalDefaultPrimaryInside || wantInsideUser;
            for (size_t i = 0; i < hits.size(); ++i) {
                double ev = 0.0;
                if (infernalExactScoreToEvalue(qm.exactModel, evalueModeInside, hits[i].mode, hits[i].cyk, targetDbResidues, ev)) {
                    hits[i].evalue = ev;
                    hits[i].hasEvalue = true;
                }
            }
        } else {
            for (size_t i = 0; i < hits.size(); ++i) {
                hits[i].evalue = nativeCalibratedScoreToEvalue(hits[i].cyk, targetDbResidues);
                hits[i].hasEvalue = true;
            }
        }
        std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
            if (a.dbKey != b.dbKey) {
                return a.dbKey < b.dbKey;
            }
            if (a.start1 != b.start1) {
                return a.start1 < b.start1;
            }
            return a.end1 < b.end1;
        });
        totalHits += hits.size();

        const unsigned int modelLen = qm.useInfernalCompat
            ? static_cast<unsigned int>(std::max(0, qm.exactModel.clen))
            : static_cast<unsigned int>(std::max(0, qm.g.modelLen));
        std::string resultPayload;
        resultPayload.reserve(hits.size() * 96);
        for (size_t i = 0; i < hits.size(); ++i) {
            const int dbStart = std::max(0, hits[i].start1 - 1);
            const int dbEnd = std::max(0, hits[i].end1 - 1);
            const unsigned int alnLen = static_cast<unsigned int>(std::max(1, dbEnd - dbStart + 1));
            const unsigned int qLen = (modelLen > 0) ? modelLen : alnLen;
            const int qStart = 0;
            const int qEnd = static_cast<int>(std::min(qLen, alnLen)) - 1;
            const float qcov = (qLen > 0) ? static_cast<float>(alnLen) / static_cast<float>(qLen) : 0.0f;
            const float dbcov = (hits[i].dbLen > 0) ? static_cast<float>(alnLen) / static_cast<float>(hits[i].dbLen) : 0.0f;
            const int bitScore = static_cast<int>(std::lrint(hits[i].cyk));
            const double evalue = hits[i].hasEvalue ? hits[i].evalue : 1.0;
            const bool hasBacktrace = (hits[i].cigar.empty() == false && hits[i].cigar != "NA");
            float seqIdVal = 0.0f;
            bool haveTraceIdentity = false;
            if (qm.useInfernalCompat && hasBacktrace && hits[i].traceStates.empty() == false && hits[i].traceStates != "NA") {
                if (hits[i].targetIdx < targetSeqs.size() && hits[i].start1 > 0 && hits[i].end1 >= hits[i].start1) {
                    const FastaSeq &fs = targetSeqs[hits[i].targetIdx];
                    const size_t s0 = static_cast<size_t>(hits[i].start1 - 1);
                    const size_t slen = static_cast<size_t>(hits[i].end1 - hits[i].start1 + 1);
                    if (s0 < fs.seq.size()) {
                        const std::string obs = fs.seq.substr(s0, std::min(slen, fs.seq.size() - s0));
                        const std::vector<int> trace = decodeTraceStates(hits[i].traceStates);
                        const std::string cons = modelTraceConsensusForCigar(qm.exactModel, trace);
                        seqIdVal = seqIdFromCigarConsensus(hits[i].cigar, cons, obs);
                        haveTraceIdentity = true;
                    }
                }
            }
            if (!haveTraceIdentity) {
                const unsigned int bitScorePos = static_cast<unsigned int>(std::max(0, bitScore));
                seqIdVal = Matcher::estimateSeqIdByScorePerCol(static_cast<uint16_t>(std::min(bitScorePos, 65535u)),
                                                               std::max(1u, alnLen),
                                                               std::max(1u, alnLen));
            }
            Matcher::result_t res(hits[i].dbKey,
                                  bitScore,
                                  std::min(1.0f, std::max(0.0f, qcov)),
                                  std::min(1.0f, std::max(0.0f, dbcov)),
                                  std::min(1.0f, std::max(0.0f, seqIdVal)),
                                  evalue,
                                  alnLen,
                                  qStart,
                                  qEnd,
                                  qLen,
                                  dbStart,
                                  dbEnd,
                                  hits[i].dbLen,
                                  hasBacktrace ? hits[i].cigar : std::string());
            const bool addBacktrace = hasBacktrace;
            const size_t len = Matcher::resultToBuffer(buffer, res, addBacktrace, false);
            resultPayload.append(buffer, len);
        }
        resultWriter.writeData(resultPayload.c_str(), resultPayload.size(), qm.key, 0);
    }
    resultWriter.close();

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
