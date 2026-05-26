#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "MMseqsMPI.h"
#include "Util.h"
#include "RnaMatcher.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "commons/DinucleotideMapping.h"
#include "commons/RNAFoldBridge.h"
#ifdef HAVE_RNAFOLD_PREDICT
extern "C" {
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/partfunc/local.h>
}
#endif
#ifdef HAVE_INFERNAL_BRIDGE
#include "infernal/InfernalBridge.h"
#endif

#ifdef OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

static const float NEG_INF_F = -1e30f;

struct StructVec {
    float u;
    float l;
    float r;

    StructVec() : u(1.0f), l(0.0f), r(0.0f) {}
    StructVec(float uu, float ll, float rr) : u(uu), l(ll), r(rr) {}
};

struct SparsePairMatrix {
    int len = 0;
    std::vector<int> rowOffsets;
    std::vector<int> cols;
    std::vector<float> probs;
};

struct FoldResult {
    std::vector<StructVec> structVec;
    std::vector<int> partners;
    std::vector<float> pairMatrix;
    std::shared_ptr<const SparsePairMatrix> sparsePairMatrix;
    std::string backendUsed;
};

struct OrientedHit {
    int qStart = -1;
    int qEnd = -1;
    int dbStart = -1;
    int dbEnd = -1;
    bool reverseStrand = false;
    std::string backtrace;
};

struct WindowInfo {
    std::string seq;
    std::vector<StructVec> structVec;
    std::vector<int> partners;
    std::vector<float> pairMatrix;
    std::shared_ptr<const std::vector<float>> pairMatrixRef;
    std::shared_ptr<const SparsePairMatrix> sparsePairMatrix;
    int fullLo = 0;
    int fullHi = -1;
    bool reverseStrand = false;
    int localSeedStart = 0;
    int localSeedEnd = 0;
};

static std::shared_ptr<const SparsePairMatrix> makeSparsePairMatrixFromDense(int len,
                                                                              const std::vector<float> &dense) {
    if (len <= 0 || dense.size() != static_cast<size_t>(len) * static_cast<size_t>(len)) {
        return std::shared_ptr<const SparsePairMatrix>();
    }
    std::shared_ptr<SparsePairMatrix> matrix(new SparsePairMatrix());
    matrix->len = len;
    matrix->rowOffsets.assign(static_cast<size_t>(len + 1), 0);
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            if (i != j && dense[static_cast<size_t>(i) * static_cast<size_t>(len) + static_cast<size_t>(j)] > 0.0f) {
                matrix->rowOffsets[static_cast<size_t>(i + 1)]++;
            }
        }
    }
    for (int i = 1; i <= len; ++i) {
        matrix->rowOffsets[static_cast<size_t>(i)] += matrix->rowOffsets[static_cast<size_t>(i - 1)];
    }
    const int total = matrix->rowOffsets[static_cast<size_t>(len)];
    matrix->cols.assign(static_cast<size_t>(total), 0);
    matrix->probs.assign(static_cast<size_t>(total), 0.0f);
    std::vector<int> cursor = matrix->rowOffsets;
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            const float prob = dense[static_cast<size_t>(i) * static_cast<size_t>(len) + static_cast<size_t>(j)];
            if (i == j || prob <= 0.0f) {
                continue;
            }
            const int pos = cursor[static_cast<size_t>(i)]++;
            matrix->cols[static_cast<size_t>(pos)] = j;
            matrix->probs[static_cast<size_t>(pos)] = prob;
        }
    }
    return matrix;
}

static float sparsePairMatrixAt(const SparsePairMatrix &matrix, int i, int j) {
    if (i < 0 || j < 0 || i >= matrix.len || j >= matrix.len || matrix.rowOffsets.size() != static_cast<size_t>(matrix.len + 1)) {
        return 0.0f;
    }
    const int begin = matrix.rowOffsets[static_cast<size_t>(i)];
    const int end = matrix.rowOffsets[static_cast<size_t>(i + 1)];
    const std::vector<int>::const_iterator first = matrix.cols.begin() + begin;
    const std::vector<int>::const_iterator last = matrix.cols.begin() + end;
    const std::vector<int>::const_iterator it = std::lower_bound(first, last, j);
    if (it == last || *it != j) {
        return 0.0f;
    }
    return matrix.probs[static_cast<size_t>(it - matrix.cols.begin())];
}

static std::vector<float> sparsePairMatrixToDense(const SparsePairMatrix &matrix) {
    std::vector<float> dense(static_cast<size_t>(std::max(0, matrix.len * matrix.len)), 0.0f);
    if (matrix.len <= 0 || matrix.rowOffsets.size() != static_cast<size_t>(matrix.len + 1)) {
        return dense;
    }
    for (int row = 0; row < matrix.len; ++row) {
        for (int idx = matrix.rowOffsets[static_cast<size_t>(row)]; idx < matrix.rowOffsets[static_cast<size_t>(row + 1)]; ++idx) {
            dense[static_cast<size_t>(row * matrix.len + matrix.cols[static_cast<size_t>(idx)])] = matrix.probs[static_cast<size_t>(idx)];
        }
    }
    return dense;
}

static void clearWindowStructure(WindowInfo &window) {
    std::vector<StructVec>().swap(window.structVec);
    std::vector<int>().swap(window.partners);
    std::vector<float>().swap(window.pairMatrix);
    window.pairMatrixRef.reset();
    window.sparsePairMatrix.reset();
}

static void clearWindowMaterialized(WindowInfo &window) {
    clearWindowStructure(window);
    std::string().swap(window.seq);
}

static bool windowBoundsValid(const WindowInfo &window) {
    return window.fullHi >= window.fullLo && window.fullLo >= 0;
}

static const std::vector<float> *windowPairMatrixPtr(const WindowInfo &window) {
    if (!window.pairMatrix.empty()) {
        return &window.pairMatrix;
    }
    if (window.pairMatrixRef && !window.pairMatrixRef->empty()) {
        return window.pairMatrixRef.get();
    }
    return NULL;
}

static bool windowHasPairMatrix(const WindowInfo &window) {
    return windowPairMatrixPtr(window) != NULL
        || (window.sparsePairMatrix && !window.sparsePairMatrix->probs.empty());
}

static float windowPairMatrixAt(const WindowInfo &window, int i, int j) {
    const int targetLen = static_cast<int>(window.seq.size());
    if (i < 0 || j < 0 || i >= targetLen || j >= targetLen) {
        return 0.0f;
    }
    const std::vector<float> *pairMatrix = windowPairMatrixPtr(window);
    if (pairMatrix != NULL && pairMatrix->size() == static_cast<size_t>(targetLen * targetLen)) {
        return (*pairMatrix)[static_cast<size_t>(i * targetLen + j)];
    }
    if (window.sparsePairMatrix && window.sparsePairMatrix->len == targetLen) {
        return sparsePairMatrixAt(*window.sparsePairMatrix, i, j);
    }
    return 0.0f;
}

struct ProjectedAlignment {
    std::string row;
    std::vector<int> targetPosByQuery;
    std::vector<StructVec> structByQuery;
};

struct AlignmentResult {
    bool valid = false;
    float score = 0.0f;
    int qStart = -1;
    int qEnd = -1;
    int dbStart = -1;
    int dbEnd = -1;
    int matches = 0;
    std::string backtrace;
};

struct AlignmentSummary {
    int queryResidues = 0;
    int targetResidues = 0;
    int matchColumns = 0;
    float qcov = 0.0f;
    float dbcov = 0.0f;
    float seqId = 0.0f;
    float meanStructSim = 0.0f;
    float meanStemAnchor = 0.0f;
    float consistency = 0.0f;
};

struct CandidateHit {
    RnaMatcher::result_t raw;
    OrientedHit oriented;
    WindowInfo window;
};

struct ModelHit {
    size_t candidateIdx = 0;
    float seedWeight = 1.0f;
    float weight = 1.0f;
    AlignmentResult aln;
};

struct StructureSettings {
    std::string backend = "auto";
    int minLoop = 3;
    int linearfoldMinLength = 128;
    int linearfoldBeamSize = 100;
    int rnaPlfoldWindow = 0;
    int rnaPlfoldSpan = 0;
    bool keepPairMatrix = true;
    bool sparsePairMatrix = false;
    bool dotBracketOnly = false;
};

struct StemGuide {
    std::vector<int> partnerByQuery;
    std::vector<float> anchorPosByQuery;
    std::vector<float> supportByQuery;
};

struct StemPattern {
    int leftStart = -1;
    int rightOuter = -1;
    int stemLen = 0;
    int armLen = 0;
    int span = 0;
    float priority = 0.0f;
    float pairSupport = 0.0f;
    std::string leftArm;
    std::string rightArm;
    std::vector<std::array<float, 4>> leftProfile;
    std::vector<std::array<float, 4>> rightProfile;
};

struct StemMatch {
    int patternIdx = -1;
    int leftQuery = -1;
    int rightQuery = -1;
    int leftTarget = -1;
    int rightTarget = -1;
    int diagCenter = 0;
    float score = 0.0f;
    float support = 0.0f;
};

struct ChainState {
    float score = NEG_INF_F;
    int hitIdx = -1;
    int prevHitIdx = -1;
    int prevRank = -1;
};

struct AnchorBlock {
    int qStart = -1;
    int qEnd = -1;
    int tStart = -1;
    int tEnd = -1;
    float support = 0.0f;
};

struct ModuleBandPrior {
    bool valid = false;
    int expectedRightStart = -1;
    int leftDiagLo = 0;
    int leftDiagHi = 0;
    int rightDiagLo = 0;
    int rightDiagHi = 0;
    int spanLo = 0;
    int spanHi = 0;
    float leftDiagCenter = 0.0f;
    float rightDiagCenter = 0.0f;
    float spanCenter = 0.0f;
};

struct LolPosterior {
    int n = 0;
    int m = 0;
    float maxP = 0.0f;
    float macScore = 0.0f;
    std::vector<float> posterior;
};

struct LolAnchor {
    int q = -1;
    int t = -1;
    float support = 0.0f;
};

struct LolSectionBlock {
    bool stem = false;
    int qLo = -1;
    int qHi = -1;
    int tLo = -1;
    int tHi = -1;
    float support = 0.0f;
};

static ProjectedAlignment projectAlignment(int qLen,
                                           const WindowInfo &window,
                                           const AlignmentResult &aln);

static AlignmentResult stitchLolSectionConstrainedCYK(const std::string &querySeq,
                                                      const std::vector<std::array<float, 4>> &profile,
                                                      const std::vector<StructVec> &queryStruct,
                                                      const std::vector<int> &partners,
                                                      const WindowInfo &window,
                                                      const std::vector<StemPattern> &patterns,
                                                      const std::vector<StemMatch> &matches,
                                                      const std::vector<int> &chain,
                                                      float structScoreWeight);

static inline char normalizeBase(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (c == 'T') {
        c = 'U';
    }
    return c;
}

static inline char complementBase(char c) {
    c = normalizeBase(c);
    switch (c) {
        case 'A': return 'U';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'U': return 'A';
        default: return 'N';
    }
}

static inline bool isCanonicalPair(char a, char b) {
    a = normalizeBase(a);
    b = normalizeBase(b);
    return (a == 'A' && b == 'U')
        || (a == 'U' && b == 'A')
        || (a == 'C' && b == 'G')
        || (a == 'G' && b == 'C')
        || (a == 'G' && b == 'U')
        || (a == 'U' && b == 'G');
}

static inline int encodeBase(char c) {
    switch (normalizeBase(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'U': return 3;
        default: return -1;
    }
}

static inline float clampf(float x, float lo, float hi) {
    return std::max(lo, std::min(hi, x));
}

static float parseEnvFloat(const char *name, float defaultValue) {
    const char *env = std::getenv(name);
    return (env != NULL) ? static_cast<float>(std::atof(env)) : defaultValue;
}

static int parseEnvInt(const char *name, int defaultValue) {
    const char *env = std::getenv(name);
    return (env != NULL) ? std::atoi(env) : defaultValue;
}

static unsigned long long parseEnvUnsignedLongLong(const char *name, unsigned long long defaultValue) {
    const char *env = std::getenv(name);
    if (env == NULL || *env == '\0') {
        return defaultValue;
    }
    char *end = NULL;
    errno = 0;
    const unsigned long long parsed = std::strtoull(env, &end, 10);
    if (errno != 0 || end == env) {
        return defaultValue;
    }
    return parsed;
}

static unsigned long long densePairMatrixBytesForLength(size_t len) {
    return static_cast<unsigned long long>(len)
         * static_cast<unsigned long long>(len)
         * static_cast<unsigned long long>(sizeof(float));
}

static std::vector<int> dotBracketToPartners(const std::string &dot) {
    std::vector<int> partners(dot.size(), -1);
    std::vector<int> stack;
    stack.reserve(dot.size());
    for (size_t i = 0; i < dot.size(); ++i) {
        if (dot[i] == '(') {
            stack.push_back(static_cast<int>(i));
        } else if (dot[i] == ')' && !stack.empty()) {
            const int j = stack.back();
            stack.pop_back();
            partners[i] = j;
            partners[static_cast<size_t>(j)] = static_cast<int>(i);
        }
    }
    return partners;
}

static std::vector<StructVec> dotBracketToStructVec(const std::string &dot) {
    std::vector<StructVec> out(dot.size());
    for (size_t i = 0; i < dot.size(); ++i) {
        switch (dot[i]) {
            case '(':
                out[i] = StructVec(0.0f, 0.0f, 1.0f);
                break;
            case ')':
                out[i] = StructVec(0.0f, 1.0f, 0.0f);
                break;
            default:
                out[i] = StructVec(1.0f, 0.0f, 0.0f);
                break;
        }
    }
    return out;
}

static std::vector<StructVec> pairMatrixToStructVec(const std::vector<float> &pairMatrix,
                                                    const std::vector<int> &partners,
                                                    int len) {
    std::vector<StructVec> out(static_cast<size_t>(len), StructVec(1.0f, 0.0f, 0.0f));
    if (len <= 0 || pairMatrix.size() != static_cast<size_t>(len * len)) {
        return out;
    }
    for (int i = 0; i < len; ++i) {
        float leftMass = 0.0f;
        float rightMass = 0.0f;
        for (int j = 0; j < len; ++j) {
            const float p = pairMatrix[static_cast<size_t>(i * len + j)];
            if (p <= 0.0f) {
                continue;
            }
            if (j < i) {
                leftMass += p;
            } else if (j > i) {
                rightMass += p;
            }
        }
        float paired = leftMass + rightMass;
        if (paired > 1.0f && paired > 0.0f) {
            const float scale = 1.0f / paired;
            leftMass *= scale;
            rightMass *= scale;
            paired = 1.0f;
        }
        out[static_cast<size_t>(i)] = StructVec(1.0f - paired, leftMass, rightMass);
        if (partners.size() == static_cast<size_t>(len) && partners[static_cast<size_t>(i)] >= 0) {
            if (partners[static_cast<size_t>(i)] < i) {
                out[static_cast<size_t>(i)].l = std::max(out[static_cast<size_t>(i)].l, 0.5f);
            } else {
                out[static_cast<size_t>(i)].r = std::max(out[static_cast<size_t>(i)].r, 0.5f);
            }
            out[static_cast<size_t>(i)].u = std::min(out[static_cast<size_t>(i)].u, 0.5f);
        }
    }
    return out;
}

static std::vector<StructVec> sparsePairMatrixToStructVec(const SparsePairMatrix &pairMatrix,
                                                          const std::vector<int> &partners,
                                                          int len) {
    std::vector<StructVec> out(static_cast<size_t>(len), StructVec(1.0f, 0.0f, 0.0f));
    if (len <= 0 || pairMatrix.len != len || pairMatrix.rowOffsets.size() != static_cast<size_t>(len + 1)) {
        return out;
    }
    for (int i = 0; i < len; ++i) {
        float leftMass = 0.0f;
        float rightMass = 0.0f;
        for (int idx = pairMatrix.rowOffsets[static_cast<size_t>(i)]; idx < pairMatrix.rowOffsets[static_cast<size_t>(i + 1)]; ++idx) {
            const int j = pairMatrix.cols[static_cast<size_t>(idx)];
            const float prob = pairMatrix.probs[static_cast<size_t>(idx)];
            if (prob <= 0.0f) {
                continue;
            }
            if (j < i) {
                leftMass += prob;
            } else if (j > i) {
                rightMass += prob;
            }
        }
        float paired = leftMass + rightMass;
        if (paired > 1.0f && paired > 0.0f) {
            const float scale = 1.0f / paired;
            leftMass *= scale;
            rightMass *= scale;
            paired = 1.0f;
        }
        out[static_cast<size_t>(i)] = StructVec(1.0f - paired, leftMass, rightMass);
        if (partners.size() == static_cast<size_t>(len) && partners[static_cast<size_t>(i)] >= 0) {
            if (partners[static_cast<size_t>(i)] < i) {
                out[static_cast<size_t>(i)].l = std::max(out[static_cast<size_t>(i)].l, 0.5f);
            } else {
                out[static_cast<size_t>(i)].r = std::max(out[static_cast<size_t>(i)].r, 0.5f);
            }
            out[static_cast<size_t>(i)].u = std::min(out[static_cast<size_t>(i)].u, 0.5f);
        }
    }
    return out;
}

static float basePreferenceFromProfile(const std::array<float, 4> &column,
                                       int base) {
    if (base < 0) {
        return 0.0f;
    }
    float best = column[0];
    float worst = column[0];
    for (int b = 1; b < 4; ++b) {
        best = std::max(best, column[static_cast<size_t>(b)]);
        worst = std::min(worst, column[static_cast<size_t>(b)]);
    }
    if (best <= worst + 1e-6f) {
        return 0.5f;
    }
    return clampf((column[static_cast<size_t>(base)] - worst) / (best - worst), 0.0f, 1.0f);
}

static float scoreArmAt(const std::string &targetSeq,
                        int pos,
                        const std::vector<std::array<float, 4>> &profileCols,
                        float minMeanSupport) {
    if (pos < 0 || pos + static_cast<int>(profileCols.size()) > static_cast<int>(targetSeq.size())) {
        return -1.0f;
    }
    float total = 0.0f;
    for (size_t k = 0; k < profileCols.size(); ++k) {
        const int b = encodeBase(targetSeq[static_cast<size_t>(pos) + k]);
        if (b < 0) {
            return -1.0f;
        }
        total += basePreferenceFromProfile(profileCols[k], b);
    }
    const float mean = total / static_cast<float>(std::max<size_t>(1, profileCols.size()));
    return (mean >= minMeanSupport) ? mean : -1.0f;
}

static float scoreStemModuleAt(const StemPattern &pat,
                               const WindowInfo &window,
                               int leftPos,
                               int rightStart,
                               float leftScore,
                               float rightScore) {
    const int targetLen = static_cast<int>(window.seq.size());
    if (leftPos < 0 || rightStart < 0
        || leftPos + pat.armLen > targetLen
        || rightStart + pat.armLen > targetLen) {
        return -1.0f;
    }

    float pairStruct = 0.0f;
    float canonical = 0.0f;
    float profilePair = 0.0f;
    for (int k = 0; k < pat.armLen; ++k) {
        const int tl = leftPos + k;
        const int tr = rightStart + (pat.armLen - 1 - k);
        const char lb = window.seq[static_cast<size_t>(tl)];
        const char rb = window.seq[static_cast<size_t>(tr)];
        if (encodeBase(lb) < 0 || encodeBase(rb) < 0) {
            return -1.0f;
        }
        canonical += isCanonicalPair(lb, rb) ? 1.0f : 0.0f;

        const float pairProb = windowPairMatrixAt(window, tl, tr);
        if (pairProb > 0.0f) {
            pairStruct += pairProb;
        } else if (window.partners.size() == static_cast<size_t>(targetLen)
                   && window.partners[static_cast<size_t>(tl)] == tr) {
            pairStruct += 1.0f;
        }

        if (!pat.leftProfile.empty() && !pat.rightProfile.empty()) {
            const int lCode = encodeBase(lb);
            const int rCode = encodeBase(rb);
            profilePair += 0.5f * (basePreferenceFromProfile(pat.leftProfile[static_cast<size_t>(k)], lCode)
                                 + basePreferenceFromProfile(pat.rightProfile[static_cast<size_t>(pat.armLen - 1 - k)], rCode));
        }
    }

    const float denom = static_cast<float>(std::max(1, pat.armLen));
    canonical /= denom;
    pairStruct /= denom;
    profilePair /= denom;
    return 0.20f * (leftScore + rightScore)
         + 0.30f * canonical
         + 0.30f * pairStruct
         + 0.20f * profilePair;
}

static std::vector<float> buildPairSupportMatrix(int qLen,
                                                 const std::vector<ModelHit> &modelHits,
                                                 const std::vector<CandidateHit> &candidates) {
    std::vector<float> pairSupport(static_cast<size_t>(qLen * qLen), 0.0f);
    if (qLen <= 1 || modelHits.empty()) {
        return pairSupport;
    }

    float totalWeight = 0.0f;
    for (size_t h = 0; h < modelHits.size(); ++h) {
        const ModelHit &modelHit = modelHits[h];
        if (!modelHit.aln.valid) {
            continue;
        }
        const CandidateHit &cand = candidates[modelHit.candidateIdx];
        ProjectedAlignment proj = projectAlignment(qLen, cand.window, modelHit.aln);
        const int targetLen = static_cast<int>(cand.window.seq.size());
        if (targetLen <= 0) {
            continue;
        }
        totalWeight += modelHit.weight;
        for (int i = 0; i < qLen; ++i) {
            const int ti = proj.targetPosByQuery[static_cast<size_t>(i)];
            if (ti < 0) {
                continue;
            }
            for (int j = i + 1; j < qLen; ++j) {
                const int tj = proj.targetPosByQuery[static_cast<size_t>(j)];
                if (tj < 0) {
                    continue;
                }
                float support = 0.0f;
                support = windowPairMatrixAt(cand.window, ti, tj);
                if (support <= 0.0f
                    && cand.window.partners.size() == static_cast<size_t>(targetLen)
                    && cand.window.partners[static_cast<size_t>(ti)] == tj) {
                    support = 1.0f;
                }
                if (support <= 0.0f) {
                    continue;
                }
                const size_t ij = static_cast<size_t>(i * qLen + j);
                const size_t ji = static_cast<size_t>(j * qLen + i);
                pairSupport[ij] += modelHit.weight * support;
                pairSupport[ji] += modelHit.weight * support;
            }
        }
    }

    if (totalWeight > 0.0f) {
        for (size_t i = 0; i < pairSupport.size(); ++i) {
            pairSupport[i] /= totalWeight;
        }
    }
    return pairSupport;
}

static std::vector<int> derivePartnersFromPairSupport(const std::vector<float> &pairSupport,
                                                      int len,
                                                      int minLoop,
                                                      float minPairSupport) {
    std::vector<int> partners(static_cast<size_t>(len), -1);
    if (len <= 1 || pairSupport.size() != static_cast<size_t>(len * len)) {
        return partners;
    }

    std::vector<float> dp(static_cast<size_t>(len * len), 0.0f);
    std::vector<uint8_t> trace(static_cast<size_t>(len * len), 0);
    std::vector<int> split(static_cast<size_t>(len * len), -1);
    const auto idx = [len](int i, int j) -> size_t {
        return static_cast<size_t>(i * len + j);
    };

    for (int span = 1; span < len; ++span) {
        for (int i = 0; i + span < len; ++i) {
            const int j = i + span;
            float best = dp[idx(i + 1, j)];
            uint8_t bestTrace = 1;
            int bestSplit = -1;
            if (dp[idx(i, j - 1)] > best) {
                best = dp[idx(i, j - 1)];
                bestTrace = 2;
            }
            const float pairScore = pairSupport[idx(i, j)];
            if (j - i > minLoop && pairScore >= minPairSupport) {
                const float inside = (i + 1 <= j - 1) ? dp[idx(i + 1, j - 1)] : 0.0f;
                const float cand = inside + pairScore;
                if (cand > best) {
                    best = cand;
                    bestTrace = 3;
                }
            }
            for (int k = i + 1; k < j; ++k) {
                const float cand = dp[idx(i, k)] + dp[idx(k + 1, j)];
                if (cand > best) {
                    best = cand;
                    bestTrace = 4;
                    bestSplit = k;
                }
            }
            dp[idx(i, j)] = best;
            trace[idx(i, j)] = bestTrace;
            split[idx(i, j)] = bestSplit;
        }
    }

    std::vector<std::pair<int, int>> stack;
    stack.push_back(std::make_pair(0, len - 1));
    while (!stack.empty()) {
        const int i = stack.back().first;
        const int j = stack.back().second;
        stack.pop_back();
        if (i >= j || i < 0 || j >= len) {
            continue;
        }
        const uint8_t action = trace[idx(i, j)];
        if (action == 1) {
            stack.push_back(std::make_pair(i + 1, j));
        } else if (action == 2) {
            stack.push_back(std::make_pair(i, j - 1));
        } else if (action == 3) {
            partners[static_cast<size_t>(i)] = j;
            partners[static_cast<size_t>(j)] = i;
            stack.push_back(std::make_pair(i + 1, j - 1));
        } else if (action == 4) {
            const int k = split[idx(i, j)];
            if (k >= i && k < j) {
                stack.push_back(std::make_pair(i, k));
                stack.push_back(std::make_pair(k + 1, j));
            }
        }
    }
    return partners;
}

static bool hasAnyPartners(const std::vector<int> &partners) {
    for (size_t i = 0; i < partners.size(); ++i) {
        if (partners[i] >= 0) {
            return true;
        }
    }
    return false;
}

static bool pairCrossesAccepted(int i,
                                int j,
                                const std::vector<std::pair<int, int>> &accepted) {
    for (size_t k = 0; k < accepted.size(); ++k) {
        const int a = accepted[k].first;
        const int b = accepted[k].second;
        if ((i < a && a < j && j < b) || (a < i && i < b && b < j)) {
            return true;
        }
    }
    return false;
}

static std::vector<int> derivePartnersFromDensePairMatrix(const std::vector<float> &pairMatrix,
                                                          int len,
                                                          int minLoop,
                                                          float minPairProb) {
    std::vector<int> partners(static_cast<size_t>(std::max(0, len)), -1);
    if (len <= 1 || pairMatrix.size() != static_cast<size_t>(len * len)) {
        return partners;
    }

    std::vector<float> rowMax(static_cast<size_t>(len), 0.0f);
    for (int i = 0; i < len; ++i) {
        for (int j = i + minLoop + 1; j < len; ++j) {
            const float p = pairMatrix[static_cast<size_t>(i * len + j)];
            if (p <= 0.0f) {
                continue;
            }
            rowMax[static_cast<size_t>(i)] = std::max(rowMax[static_cast<size_t>(i)], p);
            rowMax[static_cast<size_t>(j)] = std::max(rowMax[static_cast<size_t>(j)], p);
        }
    }

    struct PairCandidateLocal {
        int i;
        int j;
        float p;
    };
    std::vector<PairCandidateLocal> candidates;
    candidates.reserve(static_cast<size_t>(len));
    const float eps = 1e-6f;
    for (int i = 0; i < len; ++i) {
        for (int j = i + minLoop + 1; j < len; ++j) {
            const float p = pairMatrix[static_cast<size_t>(i * len + j)];
            if (p < minPairProb || p <= 0.0f) {
                continue;
            }
            if ((p + eps) < rowMax[static_cast<size_t>(i)] || (p + eps) < rowMax[static_cast<size_t>(j)]) {
                continue;
            }
            PairCandidateLocal cand;
            cand.i = i;
            cand.j = j;
            cand.p = p;
            candidates.push_back(cand);
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const PairCandidateLocal &a, const PairCandidateLocal &b) {
                  if (a.p != b.p) {
                      return a.p > b.p;
                  }
                  if (a.i != b.i) {
                      return a.i < b.i;
                  }
                  return a.j < b.j;
              });

    std::vector<std::pair<int, int>> accepted;
    accepted.reserve(static_cast<size_t>(len / 2));
    for (size_t idx = 0; idx < candidates.size(); ++idx) {
        const PairCandidateLocal &cand = candidates[idx];
        if (partners[static_cast<size_t>(cand.i)] != -1 || partners[static_cast<size_t>(cand.j)] != -1) {
            continue;
        }
        if (pairCrossesAccepted(cand.i, cand.j, accepted)) {
            continue;
        }
        partners[static_cast<size_t>(cand.i)] = cand.j;
        partners[static_cast<size_t>(cand.j)] = cand.i;
        accepted.push_back(std::make_pair(cand.i, cand.j));
    }
    return partners;
}

#ifdef HAVE_RNAFOLD_PREDICT
struct CmLitePlfoldAccumulator {
    std::vector<StructVec> *structVec;
    std::vector<float> *pairMatrix;
    int len;
};

static void accumulateCmLitePlfoldPairs(FLT_OR_DBL *pr,
                                        int prSize,
                                        int i,
                                        int /*max*/,
                                        unsigned int type,
                                        void *data) {
    if ((type & VRNA_PROBS_WINDOW_BPP) == 0U) {
        return;
    }
    CmLitePlfoldAccumulator *acc = static_cast<CmLitePlfoldAccumulator *>(data);
    if (acc == NULL || acc->structVec == NULL || acc->pairMatrix == NULL || i <= 0 || i > acc->len) {
        return;
    }
    std::vector<StructVec> &vec = *acc->structVec;
    std::vector<float> &pairMatrix = *acc->pairMatrix;
    const int len = acc->len;
    for (int j = i + 1; j <= prSize && j <= len; ++j) {
        const float p = static_cast<float>(pr[j]);
        if (p <= 0.0f || std::isnan(p)) {
            continue;
        }
        vec[static_cast<size_t>(i - 1)].r += p;
        vec[static_cast<size_t>(j - 1)].l += p;
        pairMatrix[static_cast<size_t>((i - 1) * len + (j - 1))] = p;
        pairMatrix[static_cast<size_t>((j - 1) * len + (i - 1))] = p;
    }
}

static bool rnaPlfoldPredictDense(const std::string &seq,
                                  int minLoop,
                                  int windowSize,
                                  int maxBpSpan,
                                  std::vector<StructVec> &structVec,
                                  std::vector<int> &partners,
                                  std::vector<float> &pairMatrix) {
    structVec.clear();
    partners.clear();
    pairMatrix.clear();
    if (seq.empty()) {
        return false;
    }
    const int len = static_cast<int>(seq.size());
    vrna_md_t md;
    vrna_md_set_default(&md);
    md.compute_bpp = 1;
    md.window_size = (windowSize > 0) ? std::min(windowSize, len) : len;
    md.max_bp_span = (maxBpSpan > 0) ? std::min(maxBpSpan, md.window_size) : md.window_size;

    vrna_fold_compound_t *fc = vrna_fold_compound(seq.c_str(), &md, VRNA_OPTION_PF | VRNA_OPTION_WINDOW);
    if (fc == NULL) {
        return false;
    }

    structVec.assign(static_cast<size_t>(len), StructVec(0.0f, 0.0f, 0.0f));
    pairMatrix.assign(static_cast<size_t>(len * len), 0.0f);
    CmLitePlfoldAccumulator acc = { &structVec, &pairMatrix, len };
    const int ok = vrna_probs_window(fc, 0, VRNA_PROBS_WINDOW_BPP, &accumulateCmLitePlfoldPairs, &acc);
    vrna_fold_compound_free(fc);
    if (!ok) {
        structVec.clear();
        pairMatrix.clear();
        return false;
    }

    for (size_t idx = 0; idx < structVec.size(); ++idx) {
        float l = std::max(0.0f, structVec[idx].l);
        float r = std::max(0.0f, structVec[idx].r);
        float paired = l + r;
        if (paired > 1.0f && paired > 0.0f) {
            const float scale = 1.0f / paired;
            l *= scale;
            r *= scale;
            paired = 1.0f;
        }
        structVec[idx].l = l;
        structVec[idx].r = r;
        structVec[idx].u = 1.0f - paired;
    }

    const float minPairProb = clampf(parseEnvFloat("MMSEQS_LOLALIGN_RNAPLFOLD_MIN_PAIR_PROB", 0.0f), 0.0f, 1.0f);
    partners = derivePartnersFromDensePairMatrix(pairMatrix, len, minLoop, minPairProb);
    if (hasAnyPartners(partners)) {
        structVec = pairMatrixToStructVec(pairMatrix, partners, len);
    }
    return true;
}
#endif

static std::string chooseStructureBackend(const StructureSettings &settings,
                                          int seqLen,
                                          bool hasUnknownResidues) {
    if (hasUnknownResidues) {
        return "canonical";
    }
    if (settings.backend != "auto") {
        return settings.backend;
    }
    if (seqLen >= settings.linearfoldMinLength) {
        return "linearfold";
    }
#ifdef HAVE_RNAFOLD_PREDICT
    return "rnafold";
#else
    return "canonical";
#endif
}

static FoldResult foldSequenceStructure(const std::string &seq,
                                        const StructureSettings &settings) {
    FoldResult out;
    out.backendUsed = "canonical";
    out.structVec.assign(seq.size(), StructVec());
    out.partners.assign(seq.size(), -1);
    if (seq.empty()) {
        return out;
    }

    std::string normalized(seq);
    for (size_t i = 0; i < normalized.size(); ++i) {
        normalized[i] = normalizeBase(normalized[i]);
    }
    const bool hasUnknown = std::find_if(normalized.begin(), normalized.end(),
                                         [](char c) {
                                             return c != 'A' && c != 'C' && c != 'G' && c != 'U';
                                         }) != normalized.end();
    const std::string backend = chooseStructureBackend(settings,
                                                       static_cast<int>(normalized.size()),
                                                       hasUnknown);
    if (backend == "rnaplfold" || backend == "plfold") {
#ifdef HAVE_RNAFOLD_PREDICT
        std::vector<StructVec> plfoldStruct;
        std::vector<int> plfoldPartners;
        std::vector<float> plfoldPairMatrix;
        if (rnaPlfoldPredictDense(normalized,
                                  settings.minLoop,
                                  settings.rnaPlfoldWindow,
                                  settings.rnaPlfoldSpan,
                                  plfoldStruct,
                                  plfoldPartners,
                                  plfoldPairMatrix)
            && plfoldStruct.size() == normalized.size()
            && plfoldPartners.size() == normalized.size()
            && plfoldPairMatrix.size() == normalized.size() * normalized.size()) {
            out.backendUsed = "rnaplfold";
            out.structVec.swap(plfoldStruct);
            out.partners.swap(plfoldPartners);
            if (settings.keepPairMatrix) {
                out.pairMatrix.swap(plfoldPairMatrix);
            }
            return out;
        }
#endif
    }
    if (backend == "linearfold" || backend == "partition") {
        std::string dot;
        std::vector<float> pairMatrix;
        if (rnaLinearPartitionPredict(normalized,
                                      settings.minLoop,
                                      settings.linearfoldBeamSize,
                                      dot,
                                      pairMatrix)
            && dot.size() == normalized.size()
            && pairMatrix.size() == normalized.size() * normalized.size()) {
            out.backendUsed = "linearfold";
            out.partners = dotBracketToPartners(dot);
            if (backend == "partition") {
                out.backendUsed = "partition";
                out.structVec.clear();
            } else if (settings.dotBracketOnly) {
                out.structVec = dotBracketToStructVec(dot);
            } else {
                out.structVec = pairMatrixToStructVec(pairMatrix,
                                                      out.partners,
                                                      static_cast<int>(normalized.size()));
            }
            if (settings.keepPairMatrix && !settings.dotBracketOnly) {
                if (settings.sparsePairMatrix) {
                    out.sparsePairMatrix = makeSparsePairMatrixFromDense(static_cast<int>(normalized.size()), pairMatrix);
                } else {
                    out.pairMatrix.swap(pairMatrix);
                }
            }
            return out;
        }
    }
    if (backend == "rnafold" || backend == "auto") {
        std::string dot;
        if (rnaFoldPredictDotBracket(normalized, dot, NULL) && dot.size() == normalized.size()) {
            out.backendUsed = "rnafold";
            out.partners = dotBracketToPartners(dot);
            out.structVec = dotBracketToStructVec(dot);
            return out;
        }
    }
    return out;
}


static unsigned long long lolalignAvailableMemoryBytes() {
    std::ifstream in("/proc/meminfo");
    if (in) {
        std::string key;
        unsigned long long value = 0;
        std::string unit;
        while (in >> key >> value >> unit) {
            if (key == "MemAvailable:") {
                return value * 1024ULL;
            }
        }
    }
#ifdef _SC_AVPHYS_PAGES
    const long pages = sysconf(_SC_AVPHYS_PAGES);
    const long pageSize = sysconf(_SC_PAGESIZE);
    if (pages > 0 && pageSize > 0) {
        return static_cast<unsigned long long>(pages) * static_cast<unsigned long long>(pageSize);
    }
#endif
    return 0ULL;
}

static unsigned long long saturatingAddBytes(unsigned long long a, unsigned long long b) {
    if (ULLONG_MAX - a < b) {
        return ULLONG_MAX;
    }
    return a + b;
}

static unsigned long long saturatingMulBytes(unsigned long long a, unsigned long long b) {
    if (a != 0ULL && b > ULLONG_MAX / a) {
        return ULLONG_MAX;
    }
    return a * b;
}

static unsigned long long scaleBytes(unsigned long long bytes, float scale) {
    if (bytes == 0ULL || scale <= 0.0f) {
        return 0ULL;
    }
    const long double scaled = static_cast<long double>(bytes) * static_cast<long double>(scale);
    if (scaled >= static_cast<long double>(ULLONG_MAX)) {
        return ULLONG_MAX;
    }
    return static_cast<unsigned long long>(scaled);
}

static bool lolalignDensePairRetained(const StructureSettings &settings) {
    return settings.keepPairMatrix && !settings.sparsePairMatrix && !settings.dotBracketOnly;
}

static bool lolalignDensePairScratchLikely(const StructureSettings &settings, size_t len) {
    if (len == 0) {
        return false;
    }
    const std::string &backend = settings.backend;
    if (backend == "canonical" || backend == "rnafold") {
        return false;
    }
    if (backend == "rnaplfold" || backend == "plfold") {
        return true;
    }
    if (backend == "linearfold" || backend == "partition") {
        return !(settings.keepPairMatrix && settings.sparsePairMatrix && !settings.dotBracketOnly);
    }
    if (backend == "auto") {
        return static_cast<int>(len) >= settings.linearfoldMinLength;
    }
    return true;
}

static unsigned long long lolalignLinearFoldBytesForLength(size_t len) {
    const unsigned long long linearBytes = static_cast<unsigned long long>(len)
        * static_cast<unsigned long long>(sizeof(StructVec) + sizeof(int) + 64U);
    return linearBytes;
}

static unsigned long long lolalignRetainedFoldBytes(size_t len, const StructureSettings &settings) {
    unsigned long long bytes = lolalignLinearFoldBytesForLength(len);
    if (lolalignDensePairRetained(settings)) {
        bytes = saturatingAddBytes(bytes, densePairMatrixBytesForLength(len));
    }
    return bytes;
}

static unsigned long long lolalignScratchFoldBytes(size_t len,
                                                   const StructureSettings &settings,
                                                   float scratchMultiplier) {
    unsigned long long bytes = lolalignLinearFoldBytesForLength(len);
    if (lolalignDensePairScratchLikely(settings, len)) {
        bytes = saturatingAddBytes(bytes, scaleBytes(densePairMatrixBytesForLength(len), scratchMultiplier));
    }
    return bytes;
}

static std::vector<StemPattern> buildStemPatterns(const std::string &querySeq,
                                                  const std::vector<int> &partners,
                                                  const std::vector<float> &pairSupport,
                                                  const std::vector<std::array<float, 4>> &profile) {
    std::vector<StemPattern> patterns;
    const int qLen = static_cast<int>(querySeq.size());
    if (qLen == 0 || partners.size() != static_cast<size_t>(qLen)) {
        return patterns;
    }

    const int cfgArmLen = std::max(2, parseEnvInt("MMSEQS_CMLITE_STEM_ARMLEN", 4));
    const int minStemLen = std::max(1, parseEnvInt("MMSEQS_CMLITE_STEM_MIN_LEN", 2));
    const int minSpan = std::max(cfgArmLen + 2, parseEnvInt("MMSEQS_CMLITE_STEM_MIN_SPAN", 8));
    const int maxPatterns = std::max(2, parseEnvInt("MMSEQS_CMLITE_MAX_PATTERNS", 16));

    for (int i = 0; i < qLen; ++i) {
        const int partner = partners[static_cast<size_t>(i)];
        if (partner <= i || partner >= qLen) {
            continue;
        }
        if (i > 0 && partners[static_cast<size_t>(i - 1)] == partner + 1) {
            continue;
        }

        int runLen = 1;
        while (i + runLen < qLen
               && partners[static_cast<size_t>(i + runLen)] == partner - runLen) {
            ++runLen;
        }
        if (runLen < minStemLen || partner - i < minSpan) {
            continue;
        }

        const int armLen = std::min(cfgArmLen, runLen);
        const int rightStart = partner - armLen + 1;
        if (rightStart <= i || i + armLen > qLen || rightStart < 0 || rightStart + armLen > qLen) {
            continue;
        }

        float conservation = 0.0f;
        float supportSum = 0.0f;
        for (int k = 0; k < armLen; ++k) {
            const int qLeft = i + k;
            const int qRight = partner - k;
            const int bLeft = encodeBase(querySeq[static_cast<size_t>(qLeft)]);
            const int bRight = encodeBase(querySeq[static_cast<size_t>(qRight)]);
            if (bLeft >= 0) {
                conservation += clampf((profile.empty() ? 0.75f : profile[static_cast<size_t>(qLeft)][static_cast<size_t>(bLeft)] / 8.0f), -0.5f, 1.5f);
            }
            if (bRight >= 0) {
                conservation += clampf((profile.empty() ? 0.75f : profile[static_cast<size_t>(qRight)][static_cast<size_t>(bRight)] / 8.0f), -0.5f, 1.5f);
            }
            if (pairSupport.size() == static_cast<size_t>(qLen * qLen)) {
                supportSum += pairSupport[static_cast<size_t>(qLeft * qLen + qRight)];
            }
        }
        conservation = std::max(0.25f, conservation / static_cast<float>(std::max(1, armLen * 2)));
        const float meanPairSupport = supportSum / static_cast<float>(std::max(1, armLen));

        StemPattern pat;
        pat.leftStart = i;
        pat.rightOuter = partner;
        pat.stemLen = runLen;
        pat.armLen = armLen;
        pat.span = partner - i;
        pat.pairSupport = meanPairSupport;
        pat.leftArm = querySeq.substr(static_cast<size_t>(i), static_cast<size_t>(armLen));
        pat.rightArm = querySeq.substr(static_cast<size_t>(rightStart), static_cast<size_t>(armLen));
        pat.priority = (0.45f * conservation + 0.55f * std::max(0.05f, meanPairSupport))
                     * static_cast<float>(runLen)
                     * (1.0f + 0.04f * static_cast<float>(pat.span));
        pat.leftProfile.reserve(static_cast<size_t>(armLen));
        pat.rightProfile.reserve(static_cast<size_t>(armLen));
        for (int k = 0; k < armLen; ++k) {
            pat.leftProfile.push_back(profile.empty()
                                          ? std::array<float, 4>{0.0f, 0.0f, 0.0f, 0.0f}
                                          : profile[static_cast<size_t>(i + k)]);
            pat.rightProfile.push_back(profile.empty()
                                           ? std::array<float, 4>{0.0f, 0.0f, 0.0f, 0.0f}
                                           : profile[static_cast<size_t>(rightStart + k)]);
        }
        patterns.push_back(pat);
    }

    std::sort(patterns.begin(), patterns.end(),
              [](const StemPattern &a, const StemPattern &b) {
                  if (a.priority != b.priority) {
                      return a.priority > b.priority;
                  }
                  return a.leftStart < b.leftStart;
              });
    if (static_cast<int>(patterns.size()) > maxPatterns) {
        patterns.resize(static_cast<size_t>(maxPatterns));
    }
    std::sort(patterns.begin(), patterns.end(),
              [](const StemPattern &a, const StemPattern &b) {
                  return a.leftStart < b.leftStart;
              });
    return patterns;
}

static int matchArmAt(const std::string &targetSeq,
                      int pos,
                      const std::string &pattern,
                      int maxMismatches) {
    if (pos < 0 || pos + static_cast<int>(pattern.size()) > static_cast<int>(targetSeq.size())) {
        return -1;
    }
    int mismatches = 0;
    for (size_t k = 0; k < pattern.size(); ++k) {
        const char q = normalizeBase(pattern[k]);
        const char t = normalizeBase(targetSeq[static_cast<size_t>(pos) + k]);
        if (q == 'N' || t == 'N') {
            return -1;
        }
        if (q != t) {
            ++mismatches;
            if (mismatches > maxMismatches) {
                return -1;
            }
        }
    }
    return mismatches;
}

static OrientedHit orientHit(const RnaMatcher::result_t &raw) {
    OrientedHit hit;
    hit.backtrace = raw.backtrace;

    const bool queryIsReversed = (raw.qStartPos > raw.qEndPos);
    const bool targetIsReversed = (raw.dbStartPos > raw.dbEndPos);
    int dbStart = raw.dbStartPos;
    int dbEnd = raw.dbEndPos;

    if (queryIsReversed && targetIsReversed) {
        std::swap(dbStart, dbEnd);
        std::reverse(hit.backtrace.begin(), hit.backtrace.end());
        hit.reverseStrand = false;
    } else if (queryIsReversed && !targetIsReversed) {
        std::swap(dbStart, dbEnd);
        std::reverse(hit.backtrace.begin(), hit.backtrace.end());
        hit.reverseStrand = true;
    } else if (!queryIsReversed && targetIsReversed) {
        hit.reverseStrand = true;
    } else {
        hit.reverseStrand = false;
    }

    hit.qStart = std::min(raw.qStartPos, raw.qEndPos);
    hit.qEnd = std::max(raw.qStartPos, raw.qEndPos);
    hit.dbStart = dbStart;
    hit.dbEnd = dbEnd;
    return hit;
}

static ProjectedAlignment projectAlignment(int qLen,
                                           const WindowInfo &window,
                                           const AlignmentResult &aln) {
    ProjectedAlignment proj;
    proj.row.assign(static_cast<size_t>(qLen), '-');
    proj.targetPosByQuery.assign(static_cast<size_t>(qLen), -1);
    proj.structByQuery.assign(static_cast<size_t>(qLen), StructVec());
    if (!aln.valid || aln.qStart < 0 || aln.dbStart < 0) {
        return proj;
    }

    int qPos = aln.qStart;
    int dbPos = aln.dbStart;
    for (size_t k = 0; k < aln.backtrace.size(); ++k) {
        const char op = aln.backtrace[k];
        if (op == 'M') {
            if (qPos >= 0 && qPos < qLen && dbPos >= 0 && dbPos < static_cast<int>(window.seq.size())) {
                proj.row[static_cast<size_t>(qPos)] = normalizeBase(window.seq[static_cast<size_t>(dbPos)]);
                proj.targetPosByQuery[static_cast<size_t>(qPos)] = dbPos;
                if (!window.structVec.empty()) {
                    proj.structByQuery[static_cast<size_t>(qPos)] = window.structVec[static_cast<size_t>(dbPos)];
                }
            }
            ++qPos;
            ++dbPos;
        } else if (op == 'I') {
            ++qPos;
        } else if (op == 'D') {
            ++dbPos;
        }
    }
    return proj;
}

static float structureSimilarity(const StructVec &a, const StructVec &b) {
    return a.u * b.u + a.l * b.l + a.r * b.r;
}

static AlignmentSummary summarizeAlignment(const std::string &querySeq,
                                           const std::vector<StructVec> &queryStruct,
                                           const StemGuide &stemGuide,
                                           const WindowInfo &window,
                                           const AlignmentResult &aln,
                                           unsigned int dbLen) {
    AlignmentSummary out;
    if (!aln.valid) {
        return out;
    }

    int qPos = aln.qStart;
    int dbPos = aln.dbStart;
    int identical = 0;
    float structSum = 0.0f;
    float stemSum = 0.0f;
    int stemCount = 0;
    for (size_t k = 0; k < aln.backtrace.size(); ++k) {
        const char op = aln.backtrace[k];
        if (op == 'M') {
            if (qPos >= 0 && qPos < static_cast<int>(querySeq.size())
                && dbPos >= 0 && dbPos < static_cast<int>(window.seq.size())) {
                if (normalizeBase(querySeq[static_cast<size_t>(qPos)])
                    == normalizeBase(window.seq[static_cast<size_t>(dbPos)])) {
                    ++identical;
                }
                if (qPos < static_cast<int>(queryStruct.size())
                    && dbPos < static_cast<int>(window.structVec.size())) {
                    structSum += structureSimilarity(queryStruct[static_cast<size_t>(qPos)],
                                                     window.structVec[static_cast<size_t>(dbPos)]);
                }
                if (qPos < static_cast<int>(stemGuide.anchorPosByQuery.size())
                    && stemGuide.anchorPosByQuery[static_cast<size_t>(qPos)] >= 0.0f
                    && stemGuide.supportByQuery[static_cast<size_t>(qPos)] > 0.0f) {
                    const float dist = std::fabs(static_cast<float>(dbPos) - stemGuide.anchorPosByQuery[static_cast<size_t>(qPos)]);
                    const float support = stemGuide.supportByQuery[static_cast<size_t>(qPos)];
                    const float radius = 8.0f;
                    stemSum += support * std::max(0.0f, 1.0f - dist / radius);
                    ++stemCount;
                }
            }
            ++out.queryResidues;
            ++out.targetResidues;
            ++out.matchColumns;
            ++qPos;
            ++dbPos;
        } else if (op == 'I') {
            ++out.queryResidues;
            ++qPos;
        } else if (op == 'D') {
            ++out.targetResidues;
            ++dbPos;
        }
    }

    const int qLen = static_cast<int>(querySeq.size());
    out.qcov = (qLen > 0) ? static_cast<float>(out.queryResidues) / static_cast<float>(qLen) : 0.0f;
    out.dbcov = (dbLen > 0) ? static_cast<float>(out.targetResidues) / static_cast<float>(dbLen) : 0.0f;
    out.seqId = (out.matchColumns > 0) ? static_cast<float>(identical) / static_cast<float>(out.matchColumns) : 0.0f;
    out.meanStructSim = (out.matchColumns > 0) ? structSum / static_cast<float>(out.matchColumns) : 0.0f;
    out.meanStemAnchor = (stemCount > 0) ? stemSum / static_cast<float>(stemCount) : 0.0f;
    out.consistency = clampf(0.52f * out.qcov
                           + 0.22f * out.seqId
                           + 0.16f * out.meanStructSim
                           + 0.10f * out.meanStemAnchor,
                             0.0f, 1.5f);
    return out;
}

static float modelWeightFromHit(const RnaMatcher::result_t &hit) {
    const float seqId = clampf(hit.seqId, 0.0f, 1.0f);
    const float qcov = clampf(hit.qcov, 0.0f, 1.0f);
    const float dbcov = clampf(hit.dbcov, 0.0f, 1.0f);
    const float spanBalance = std::min(qcov, dbcov) / std::max(0.10f, std::max(qcov, dbcov));
    const float quality = (0.20f + 0.80f * seqId)
                        * (0.15f + 0.85f * qcov)
                        * (0.35f + 0.65f * spanBalance);
    const float w = std::log1pf(std::max(1, hit.score)) * quality * quality;
    return std::max(0.05f, w);
}

static std::string buildQueryCentralSignature(const std::string &querySeq,
                                              const std::string &row) {
    const size_t len = std::min(querySeq.size(), row.size());
    std::string sig;
    sig.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        const char q = normalizeBase(querySeq[i]);
        const char r = normalizeBase(row[i]);
        if (encodeBase(r) < 0) {
            sig.push_back('-');
        } else if (r == q) {
            sig.push_back('=');
        } else {
            sig.push_back(r);
        }
    }
    return sig;
}

static void applyModelDiversityWeights(const std::string &querySeq,
                                       std::vector<ModelHit> &modelHits,
                                       const std::vector<CandidateHit> &candidates) {
    if (modelHits.size() < 2) {
        return;
    }

    std::vector<std::string> signatures(modelHits.size());
    std::unordered_map<std::string, size_t> signatureCounts;
    signatureCounts.reserve(modelHits.size() * 2);
    for (size_t i = 0; i < modelHits.size(); ++i) {
        const ModelHit &mh = modelHits[i];
        if (!mh.aln.valid) {
            continue;
        }
        ProjectedAlignment proj = projectAlignment(static_cast<int>(querySeq.size()),
                                                   candidates[mh.candidateIdx].window,
                                                   mh.aln);
        signatures[i] = buildQueryCentralSignature(querySeq, proj.row);
        if (!signatures[i].empty()) {
            ++signatureCounts[signatures[i]];
        }
    }

    for (size_t i = 0; i < modelHits.size(); ++i) {
        if (signatures[i].empty()) {
            continue;
        }
        const std::unordered_map<std::string, size_t>::const_iterator it = signatureCounts.find(signatures[i]);
        if (it == signatureCounts.end() || it->second == 0) {
            continue;
        }
        const float redundancy = static_cast<float>(it->second);
        modelHits[i].seedWeight /= redundancy;
        modelHits[i].weight = modelHits[i].seedWeight;
    }
}

static int effectiveDecodeSeqType(int dbtype, bool useDinucMapping) {
    if (!useDinucMapping) {
        return dbtype;
    }

    unsigned int ext = DBReader<unsigned int>::getExtendedDbtype(dbtype);
    ext |= Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE)) {
        return DBReader<unsigned int>::setExtendedDbtype(dbtype, ext);
    }
    return DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS, ext);
}

static std::string decodeMappedSequenceToRna(const Sequence &seq,
                                             const SubstitutionMatrix *subMat,
                                             const unsigned char *num2outputnum) {
    std::string out;
    out.reserve(static_cast<size_t>(seq.L));
    for (int i = 0; i < seq.L; ++i) {
        unsigned char code = seq.numSequence[i];
        if (num2outputnum != NULL) {
            code = num2outputnum[code];
        }
        char c = normalizeBase(subMat->num2aa[code]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'U') {
            c = 'N';
        }
        out.push_back(c);
    }
    return out;
}

static std::string decodeDbSequence(DBReader<unsigned int> &dbr,
                                    size_t id,
                                    unsigned int threadIdx,
                                    Sequence &mapper,
                                    bool useDinucMapping,
                                    bool isGpuDb,
                                    const SubstitutionMatrix *subMat) {
    const unsigned int seqLen = dbr.getSeqLen(id);
    if (isGpuDb) {
        const unsigned char *data = reinterpret_cast<const unsigned char *>(dbr.getDataUncompressed(id));
        mapper.mapSequence(id, dbr.getDbKey(id), std::make_pair(data, seqLen));
    } else {
        const char *data = dbr.getData(id, threadIdx);
        mapper.mapSequence(id, dbr.getDbKey(id), data, seqLen);
    }

    const Sequence::SeqAuxInfo *auxInfo = Sequence::getAuxInfo(mapper.getSeqType());
    const unsigned char *num2outputnum = (useDinucMapping && auxInfo != NULL)
        ? auxInfo->num2outputnum
        : NULL;
    return decodeMappedSequenceToRna(mapper, subMat, num2outputnum);
}

static void buildProfile(const std::string &querySeq,
                         const std::vector<ModelHit> &modelHits,
                         const std::vector<CandidateHit> &candidates,
                         const std::vector<StructVec> &baseQueryStruct,
                         float alpha,
                         std::vector<std::array<float, 4>> &profile,
                         std::vector<StructVec> &queryStruct) {
    const int qLen = static_cast<int>(querySeq.size());
    profile.assign(static_cast<size_t>(qLen), {0.0f, 0.0f, 0.0f, 0.0f});
    queryStruct = baseQueryStruct;

    std::vector<std::array<float, 4>> counts(static_cast<size_t>(qLen), {0.25f, 0.25f, 0.25f, 0.25f});
    std::vector<float> structWeight(static_cast<size_t>(qLen), 1.0f);

    for (int i = 0; i < qLen; ++i) {
        int b = encodeBase(querySeq[static_cast<size_t>(i)]);
        if (b >= 0) {
            counts[static_cast<size_t>(i)][static_cast<size_t>(b)] += 3.0f;
        }
    }

    for (size_t h = 0; h < modelHits.size(); ++h) {
        const ModelHit &modelHit = modelHits[h];
        if (!modelHit.aln.valid) {
            continue;
        }
        const CandidateHit &cand = candidates[modelHit.candidateIdx];
        ProjectedAlignment proj = projectAlignment(qLen, cand.window, modelHit.aln);
        for (int i = 0; i < qLen; ++i) {
            const char c = proj.row[static_cast<size_t>(i)];
            const int b = encodeBase(c);
            if (b < 0) {
                continue;
            }
            const float sim = structureSimilarity(queryStruct[static_cast<size_t>(i)],
                                                  proj.structByQuery[static_cast<size_t>(i)]);
            const float localWeight = modelHit.weight * (1.0f + alpha * sim);
            counts[static_cast<size_t>(i)][static_cast<size_t>(b)] += localWeight;
            queryStruct[static_cast<size_t>(i)].u += localWeight * proj.structByQuery[static_cast<size_t>(i)].u;
            queryStruct[static_cast<size_t>(i)].l += localWeight * proj.structByQuery[static_cast<size_t>(i)].l;
            queryStruct[static_cast<size_t>(i)].r += localWeight * proj.structByQuery[static_cast<size_t>(i)].r;
            structWeight[static_cast<size_t>(i)] += localWeight;
        }
    }

    for (int i = 0; i < qLen; ++i) {
        queryStruct[static_cast<size_t>(i)].u /= structWeight[static_cast<size_t>(i)];
        queryStruct[static_cast<size_t>(i)].l /= structWeight[static_cast<size_t>(i)];
        queryStruct[static_cast<size_t>(i)].r /= structWeight[static_cast<size_t>(i)];

        float sum = counts[static_cast<size_t>(i)][0] + counts[static_cast<size_t>(i)][1]
                  + counts[static_cast<size_t>(i)][2] + counts[static_cast<size_t>(i)][3];
        for (int b = 0; b < 4; ++b) {
            float prob = counts[static_cast<size_t>(i)][static_cast<size_t>(b)] / sum;
            float score = 5.0f * std::log2(prob / 0.25f);
            profile[static_cast<size_t>(i)][static_cast<size_t>(b)] = clampf(score, -10.0f, 10.0f);
        }
    }
}

static std::string consensusSequenceFromProfile(const std::string &querySeq,
                                                const std::vector<std::array<float, 4>> &profile) {
    static const char bases[4] = {'A', 'C', 'G', 'U'};
    std::string consensus;
    consensus.reserve(querySeq.size());
    for (size_t i = 0; i < querySeq.size(); ++i) {
        if (i >= profile.size()) {
            char c = normalizeBase(querySeq[i]);
            consensus.push_back((c == 'A' || c == 'C' || c == 'G' || c == 'U') ? c : 'N');
            continue;
        }
        int bestBase = 0;
        float bestScore = profile[i][0];
        for (int b = 1; b < 4; ++b) {
            if (profile[i][static_cast<size_t>(b)] > bestScore) {
                bestScore = profile[i][static_cast<size_t>(b)];
                bestBase = b;
            }
        }
        consensus.push_back(bases[bestBase]);
    }
    return consensus;
}

static RnaMatcher::result_t emitResultFromAlignment(const RnaMatcher::result_t &raw,
                                                    const WindowInfo &window,
                                                    const AlignmentResult &aln,
                                                    const AlignmentSummary &summary,
                                                    unsigned int qLen,
                                                    unsigned int dbLen) {
    if (!aln.valid) {
        return raw;
    }

    int dbStartFull = 0;
    int dbEndFull = 0;
    if (!window.reverseStrand) {
        dbStartFull = window.fullLo + aln.dbStart;
        dbEndFull = window.fullLo + aln.dbEnd;
    } else {
        dbStartFull = window.fullHi - aln.dbStart;
        dbEndFull = window.fullHi - aln.dbEnd;
    }

    const unsigned int alnLen = static_cast<unsigned int>(aln.backtrace.size());
    double eval = raw.eval;
    if (raw.score > 0) {
        double delta = static_cast<double>(aln.score) - static_cast<double>(raw.score);
        eval *= std::pow(2.0, -delta / 6.0);
    }
    eval = std::max(1e-300, std::min(DBL_MAX, eval));

    return RnaMatcher::result_t(raw.dbKey,
                                std::max(1, static_cast<int>(std::lrint(aln.score))),
                                summary.qcov,
                                summary.dbcov,
                                summary.seqId,
                                eval,
                                alnLen,
                                aln.qStart,
                                aln.qEnd,
                                qLen,
                                dbStartFull,
                                dbEndFull,
                                dbLen,
                                aln.backtrace);
}

static AlignmentResult coarseSeedAlignmentFromCandidate(const std::string &querySeq,
                                                        const CandidateHit &cand) {
    AlignmentResult out;
    if (cand.oriented.qStart < 0
        || cand.oriented.qEnd < cand.oriented.qStart
        || cand.window.localSeedStart < 0
        || cand.window.localSeedEnd < cand.window.localSeedStart
        || cand.window.localSeedEnd >= static_cast<int>(cand.window.seq.size())
        || cand.oriented.qEnd >= static_cast<int>(querySeq.size())) {
        return out;
    }

    const int qStart = cand.oriented.qStart;
    const int qEnd = cand.oriented.qEnd;
    const int tStart = cand.window.localSeedStart;
    const int tEnd = cand.window.localSeedEnd;
    const int qSpan = qEnd - qStart + 1;
    const int tSpan = tEnd - tStart + 1;
    if (qSpan <= 0 || tSpan <= 0) {
        return out;
    }

    const int gapPenalty = -2;
    const int matchScore = 2;
    const int mismatchScore = -1;
    std::vector<int> dp(static_cast<size_t>(qSpan + 1) * static_cast<size_t>(tSpan + 1), 0);
    std::vector<unsigned char> trace(static_cast<size_t>(qSpan + 1) * static_cast<size_t>(tSpan + 1), 0);
    const auto idx = [tSpan](int i, int j) {
        return static_cast<size_t>(i) * static_cast<size_t>(tSpan + 1) + static_cast<size_t>(j);
    };

    for (int i = 1; i <= qSpan; ++i) {
        dp[idx(i, 0)] = i * gapPenalty;
        trace[idx(i, 0)] = 'I';
    }
    for (int j = 1; j <= tSpan; ++j) {
        dp[idx(0, j)] = j * gapPenalty;
        trace[idx(0, j)] = 'D';
    }

    for (int i = 1; i <= qSpan; ++i) {
        const char qBase = normalizeBase(querySeq[static_cast<size_t>(qStart + i - 1)]);
        for (int j = 1; j <= tSpan; ++j) {
            const char tBase = normalizeBase(cand.window.seq[static_cast<size_t>(tStart + j - 1)]);
            const int diag = dp[idx(i - 1, j - 1)] + ((qBase == tBase && qBase != 'N') ? matchScore : mismatchScore);
            const int ins = dp[idx(i - 1, j)] + gapPenalty;
            const int del = dp[idx(i, j - 1)] + gapPenalty;
            if (diag >= ins && diag >= del) {
                dp[idx(i, j)] = diag;
                trace[idx(i, j)] = 'M';
            } else if (ins >= del) {
                dp[idx(i, j)] = ins;
                trace[idx(i, j)] = 'I';
            } else {
                dp[idx(i, j)] = del;
                trace[idx(i, j)] = 'D';
            }
        }
    }

    std::string bt;
    bt.reserve(static_cast<size_t>(qSpan + tSpan));
    int i = qSpan;
    int j = tSpan;
    int matches = 0;
    while (i > 0 || j > 0) {
        const unsigned char op = trace[idx(i, j)];
        if (op == 'M') {
            bt.push_back('M');
            const char qBase = normalizeBase(querySeq[static_cast<size_t>(qStart + i - 1)]);
            const char tBase = normalizeBase(cand.window.seq[static_cast<size_t>(tStart + j - 1)]);
            if (qBase == tBase && qBase != 'N') {
                ++matches;
            }
            --i;
            --j;
        } else if (op == 'I') {
            bt.push_back('I');
            --i;
        } else {
            bt.push_back('D');
            --j;
        }
    }
    std::reverse(bt.begin(), bt.end());
    if (bt.empty()) {
        return out;
    }

    out.valid = true;
    out.score = static_cast<float>(dp[idx(qSpan, tSpan)]);
    out.qStart = qStart;
    out.qEnd = qEnd;
    out.dbStart = tStart;
    out.dbEnd = tEnd;
    out.backtrace.swap(bt);
    out.matches = matches;
    return out;
}

static AlignmentResult rawAlignmentFromCandidate(const std::string &querySeq,
                                                 const CandidateHit &cand) {
    AlignmentResult out;
    out.valid = !cand.oriented.backtrace.empty()
             && cand.oriented.qStart >= 0
             && cand.oriented.qEnd >= cand.oriented.qStart
             && cand.window.localSeedStart >= 0
             && cand.window.localSeedEnd >= cand.window.localSeedStart;
    out.score = static_cast<float>(cand.raw.score);
    out.qStart = cand.oriented.qStart;
    out.qEnd = cand.oriented.qEnd;
    out.dbStart = cand.window.localSeedStart;
    out.dbEnd = cand.window.localSeedEnd;
    out.backtrace = cand.oriented.backtrace;
    if (out.valid) {
        for (size_t i = 0; i < out.backtrace.size(); ++i) {
            if (out.backtrace[i] == 'M') {
                ++out.matches;
            }
        }
        return out;
    }
    return coarseSeedAlignmentFromCandidate(querySeq, cand);
}

static WindowInfo buildWindowMetadataFromRaw(unsigned int tLen,
                                             const OrientedHit &hit,
                                             int qLen,
                                             float flankFrac) {
    WindowInfo w;
    w.reverseStrand = hit.reverseStrand;
    const int dbLo = std::min(hit.dbStart, hit.dbEnd);
    const int dbHi = std::max(hit.dbStart, hit.dbEnd);
    const int hitSpan = std::max(dbHi - dbLo + 1, hit.qEnd - hit.qStart + 1);
    int flank = (flankFrac > 0.0f)
        ? static_cast<int>(std::ceil(static_cast<float>(hitSpan) * flankFrac))
        : std::max(24, qLen / 4);
    flank = std::max(flank, 8);
    w.fullLo = std::max(0, dbLo - flank);
    w.fullHi = std::min(static_cast<int>(tLen) - 1, dbHi + flank);
    if (w.fullHi < w.fullLo) {
        w.fullLo = 0;
        w.fullHi = -1;
        return w;
    }

    if (!w.reverseStrand) {
        w.localSeedStart = hit.dbStart - w.fullLo;
        w.localSeedEnd = hit.dbEnd - w.fullLo;
    } else {
        w.localSeedStart = w.fullHi - hit.dbStart;
        w.localSeedEnd = w.fullHi - hit.dbEnd;
    }
    return w;
}

static void fillWindowSequenceFromTarget(WindowInfo &w, const char *targetData) {
    std::string().swap(w.seq);
    if (!windowBoundsValid(w)) {
        return;
    }

    const int winLen = w.fullHi - w.fullLo + 1;
    w.seq.reserve(static_cast<size_t>(winLen));
    if (!w.reverseStrand) {
        for (int pos = w.fullLo; pos <= w.fullHi; ++pos) {
            char c = normalizeBase(targetData[pos]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'U') {
                c = 'N';
            }
            w.seq.push_back(c);
        }
    } else {
        for (int pos = w.fullHi; pos >= w.fullLo; --pos) {
            w.seq.push_back(complementBase(targetData[pos]));
        }
    }
}

static WindowInfo buildWindowFromRaw(const char *targetData,
                                     unsigned int tLen,
                                     const OrientedHit &hit,
                                     int qLen,
                                     float flankFrac) {
    WindowInfo w = buildWindowMetadataFromRaw(tLen, hit, qLen, flankFrac);
    fillWindowSequenceFromTarget(w, targetData);
    return w;
}

static bool materializeCandidateWindow(CandidateHit &cand,
                                       DBReader<unsigned int> &tDbr,
                                       unsigned int threadIdx,
                                       Sequence &targetMapper,
                                       bool decodeTargetDinuc,
                                       bool targetGpuDb,
                                       const SubstitutionMatrix *subMat) {
    if (!cand.window.seq.empty()) {
        return true;
    }
    if (!windowBoundsValid(cand.window)) {
        return false;
    }
    const size_t tId = tDbr.getId(cand.raw.dbKey);
    if (tId == UINT_MAX) {
        return false;
    }
    const std::string targetSeq = decodeDbSequence(tDbr, tId, threadIdx, targetMapper, decodeTargetDinuc, targetGpuDb, subMat);
    fillWindowSequenceFromTarget(cand.window, targetSeq.c_str());
    return !cand.window.seq.empty();
}

static std::pair<int, int> chooseQuerySegment(const OrientedHit &hit,
                                              int qLen) {
    if (qLen <= 0) {
        return std::make_pair(0, -1);
    }
    const bool preferFullQuery = parseEnvInt("MMSEQS_CMLITE_QUERY_SEGMENT_FULL", 1) != 0;
    const int fullMaxLen = std::max(0, parseEnvInt("MMSEQS_CMLITE_QUERY_SEGMENT_FULL_MAX_LEN", 512));
    if (preferFullQuery && (fullMaxLen == 0 || qLen <= fullMaxLen)) {
        return std::make_pair(0, qLen - 1);
    }

    const int seedLo = std::max(0, std::min(hit.qStart, qLen - 1));
    const int seedHi = std::max(seedLo, std::min(hit.qEnd, qLen - 1));
    const int span = seedHi - seedLo + 1;
    const int minSlack = std::max(0, parseEnvInt("MMSEQS_CMLITE_QUERY_SEGMENT_MIN_SLACK", 0));
    const float slackFrac = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_QUERY_SEGMENT_SLACK_FRAC", 0.03f));
    const int slack = std::max(minSlack, static_cast<int>(std::ceil(static_cast<float>(span) * slackFrac)));

    int qLo = std::max(0, seedLo - slack);
    int qHi = std::min(qLen - 1, seedHi + slack);
    return std::make_pair(qLo, qHi);
}

static std::pair<int, int> chooseAnchoredQuerySegment(const std::vector<AnchorBlock> &blocks,
                                                      const OrientedHit &hit,
                                                      int qLen) {
    if (qLen <= 0) {
        return std::make_pair(0, -1);
    }
    const std::pair<int, int> seedSegment = chooseQuerySegment(hit, qLen);
    if (blocks.empty()) {
        return seedSegment;
    }

    int anchorLo = qLen - 1;
    int anchorHi = 0;
    for (size_t i = 0; i < blocks.size(); ++i) {
        anchorLo = std::min(anchorLo, blocks[i].qStart);
        anchorHi = std::max(anchorHi, blocks[i].qEnd);
    }
    anchorLo = std::max(0, std::min(anchorLo, qLen - 1));
    anchorHi = std::max(anchorLo, std::min(anchorHi, qLen - 1));

    const int span = anchorHi - anchorLo + 1;
    const int minSlack = std::max(0, parseEnvInt("MMSEQS_CMLITE_ANCHORED_QUERY_SEGMENT_MIN_SLACK", 2));
    const float slackFrac = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_ANCHORED_QUERY_SEGMENT_SLACK_FRAC", 0.03f));
    const int slack = std::max(minSlack, static_cast<int>(std::ceil(static_cast<float>(span) * slackFrac)));

    int qLo = std::max(0, anchorLo - slack);
    int qHi = std::min(qLen - 1, anchorHi + slack);

    qLo = std::max(seedSegment.first, qLo);
    qHi = std::min(seedSegment.second, qHi);
    if (qLo > qHi) {
        return seedSegment;
    }
    return std::make_pair(qLo, qHi);
}

static AlignmentResult alignProfileGlocal(const std::string &querySeq,
                                          const std::vector<std::array<float, 4>> &profile,
                                          const std::vector<StructVec> &queryStruct,
                                          const StemGuide &stemGuide,
                                          const std::string &targetSeq,
                                          const std::vector<StructVec> &targetStruct,
                                          int queryStart,
                                          int queryEnd,
                                          float structScoreWeight,
                                          float stemScoreWeight,
                                          int stemRadius,
                                          int diagCenter,
                                          int bandRadius) {
    AlignmentResult out;
    const int fullQueryLen = static_cast<int>(querySeq.size());
    const int qLo = std::max(0, std::min(queryStart, fullQueryLen - 1));
    const int qHi = std::max(qLo, std::min(queryEnd, fullQueryLen - 1));
    const int n = qHi - qLo + 1;
    const int m = static_cast<int>(targetSeq.size());
    if (n == 0 || m == 0) {
        return out;
    }

    const size_t cells = static_cast<size_t>(n + 1) * static_cast<size_t>(m + 1);
    std::vector<float> M(cells, 0.0f);
    std::vector<float> E(cells, NEG_INF_F);
    std::vector<float> F(cells, NEG_INF_F);
    std::vector<uint8_t> tM(cells, 0);
    std::vector<uint8_t> tE(cells, 0);
    std::vector<uint8_t> tF(cells, 0);

    const auto idx = [m](int i, int j) -> size_t {
        return static_cast<size_t>(i) * static_cast<size_t>(m + 1) + static_cast<size_t>(j);
    };

    const float gapOpen = 6.0f;
    const float gapExtend = 1.0f;
    M[idx(0, 0)] = 0.0f;
    F[idx(0, 0)] = 0.0f;
    for (int j = 1; j <= m; ++j) {
        M[idx(0, j)] = 0.0f;
        F[idx(0, j)] = 0.0f;
        tF[idx(0, j)] = 3;
    }
    for (int i = 1; i <= n; ++i) {
        const size_t cur = idx(i, 0);
        const size_t up = idx(i - 1, 0);
        float eFromM = M[up] - gapOpen;
        float eFromE = E[up] - gapExtend;
        if (eFromM >= eFromE) {
            E[cur] = eFromM;
            tE[cur] = 1;
        } else {
            E[cur] = eFromE;
            tE[cur] = 2;
        }
        M[cur] = NEG_INF_F;
        F[cur] = NEG_INF_F;
    }

    for (int i = 1; i <= n; ++i) {
        int jLo = 1;
        int jHi = m;
        if (bandRadius >= 0) {
            jLo = std::max(1, (i - 1) + diagCenter - bandRadius + 1);
            jHi = std::min(m, (i - 1) + diagCenter + bandRadius + 1);
            if (jLo > jHi) {
                continue;
            }
        }

        for (int j = jLo; j <= jHi; ++j) {
            const size_t cur = idx(i, j);
            const size_t up = idx(i - 1, j);
            const size_t left = idx(i, j - 1);
            const size_t diag = idx(i - 1, j - 1);

            float eFromM = M[up] - gapOpen;
            float eFromE = E[up] - gapExtend;
            if (eFromM >= eFromE) {
                E[cur] = eFromM;
                tE[cur] = 1;
            } else {
                E[cur] = eFromE;
                tE[cur] = 2;
            }

            float fFromM = M[left] - gapOpen;
            float fFromF = F[left] - gapExtend;
            if (fFromM >= fFromF) {
                F[cur] = fFromM;
                tF[cur] = 1;
            } else {
                F[cur] = fFromF;
                tF[cur] = 3;
            }

            const int qIdx = qLo + i - 1;
            int b = encodeBase(targetSeq[static_cast<size_t>(j - 1)]);
            float emit = (b >= 0 && qIdx < static_cast<int>(profile.size()))
                ? profile[static_cast<size_t>(qIdx)][static_cast<size_t>(b)]
                : -3.0f;
            if (qIdx < static_cast<int>(queryStruct.size()) && j - 1 < static_cast<int>(targetStruct.size())) {
                const float sim = structureSimilarity(queryStruct[static_cast<size_t>(qIdx)],
                                                      targetStruct[static_cast<size_t>(j - 1)]);
                emit += structScoreWeight * 4.0f * (sim - (1.0f / 3.0f));
            }
            if (qIdx < static_cast<int>(stemGuide.anchorPosByQuery.size())
                && stemGuide.anchorPosByQuery[static_cast<size_t>(qIdx)] >= 0.0f
                && stemGuide.supportByQuery[static_cast<size_t>(qIdx)] > 0.0f) {
                const float dist = std::fabs(static_cast<float>(j - 1) - stemGuide.anchorPosByQuery[static_cast<size_t>(qIdx)]);
                const float radius = static_cast<float>(std::max(1, stemRadius));
                const float anchorScore = std::max(0.0f, 1.0f - dist / radius);
                emit += stemScoreWeight * stemGuide.supportByQuery[static_cast<size_t>(qIdx)] * anchorScore;
            }

            float prev = M[diag];
            uint8_t prevState = 1;
            if (E[diag] > prev) {
                prev = E[diag];
                prevState = 2;
            }
            if (F[diag] > prev) {
                prev = F[diag];
                prevState = 3;
            }
            M[cur] = prev + emit;
            tM[cur] = prevState;
        }
    }

    float bestScore = NEG_INF_F;
    int bestJ = 0;
    uint8_t bestState = 0;
    for (int j = 0; j <= m; ++j) {
        const size_t cur = idx(n, j);
        if (M[cur] > bestScore) {
            bestScore = M[cur];
            bestJ = j;
            bestState = 1;
        }
        if (E[cur] > bestScore) {
            bestScore = E[cur];
            bestJ = j;
            bestState = 2;
        }
        if (F[cur] > bestScore) {
            bestScore = F[cur];
            bestJ = j;
            bestState = 3;
        }
    }

    if (bestScore <= NEG_INF_F / 2.0f || bestState == 0) {
        return out;
    }

    std::string bt;
    bt.reserve(static_cast<size_t>(n + m));
    int i = n;
    int j = bestJ;
    uint8_t state = bestState;
    int matches = 0;

    while (i > 0 && j >= 0 && state != 0) {
        const size_t cur = idx(i, j);
        if (state == 1) {
            if (j <= 0) {
                break;
            }
            bt.push_back('M');
            if (normalizeBase(querySeq[static_cast<size_t>(qLo + i - 1)]) == normalizeBase(targetSeq[static_cast<size_t>(j - 1)])) {
                ++matches;
            }
            state = tM[cur];
            --i;
            --j;
        } else if (state == 2) {
            bt.push_back('I');
            state = tE[cur];
            --i;
        } else {
            if (j <= 0) {
                break;
            }
            bt.push_back('D');
            state = tF[cur];
            --j;
        }
    }

    std::reverse(bt.begin(), bt.end());
    out.valid = !bt.empty();
    out.score = bestScore;
    out.qStart = qLo;
    out.qEnd = qHi;
    out.dbStart = j;
    out.dbEnd = bestJ - 1;
    out.matches = matches;
    out.backtrace.swap(bt);
    return out;
}

static AlignmentResult alignProfileLocalLocal(const std::string &querySeq,
                                              const std::vector<std::array<float, 4>> &profile,
                                              const std::vector<StructVec> &queryStruct,
                                              const StemGuide &stemGuide,
                                              const std::string &targetSeq,
                                              const std::vector<StructVec> &targetStruct,
                                              int queryStart,
                                              int queryEnd,
                                              float structScoreWeight,
                                              float stemScoreWeight,
                                              int stemRadius,
                                              int diagCenter,
                                              int bandRadius,
                                              const std::vector<float> *bonusMatrix) {
    AlignmentResult out;
    const int fullQueryLen = static_cast<int>(querySeq.size());
    const int qLo = std::max(0, std::min(queryStart, fullQueryLen - 1));
    const int qHi = std::max(qLo, std::min(queryEnd, fullQueryLen - 1));
    const int n = qHi - qLo + 1;
    const int m = static_cast<int>(targetSeq.size());
    if (n == 0 || m == 0) {
        return out;
    }

    const size_t cells = static_cast<size_t>(n + 1) * static_cast<size_t>(m + 1);
    std::vector<float> M(cells, 0.0f);
    std::vector<float> E(cells, 0.0f);
    std::vector<float> F(cells, 0.0f);
    std::vector<uint8_t> tM(cells, 0);
    std::vector<uint8_t> tE(cells, 0);
    std::vector<uint8_t> tF(cells, 0);

    const auto idx = [m](int i, int j) -> size_t {
        return static_cast<size_t>(i) * static_cast<size_t>(m + 1) + static_cast<size_t>(j);
    };

    const float gapOpen = 6.0f;
    const float gapExtend = 1.0f;
    float bestScore = 0.0f;
    int bestI = 0;
    int bestJ = 0;
    uint8_t bestState = 0;

    for (int i = 1; i <= n; ++i) {
        int jLo = 1;
        int jHi = m;
        if (bandRadius >= 0) {
            jLo = std::max(1, (i - 1) + diagCenter - bandRadius + 1);
            jHi = std::min(m, (i - 1) + diagCenter + bandRadius + 1);
            if (jLo > jHi) {
                continue;
            }
        }

        for (int j = jLo; j <= jHi; ++j) {
            const size_t cur = idx(i, j);
            const size_t up = idx(i - 1, j);
            const size_t left = idx(i, j - 1);
            const size_t diag = idx(i - 1, j - 1);

            float eFromM = M[up] - gapOpen;
            float eFromE = E[up] - gapExtend;
            float eBest = 0.0f;
            uint8_t eTrace = 0;
            if (eFromM >= eFromE && eFromM > 0.0f) {
                eBest = eFromM;
                eTrace = 1;
            } else if (eFromE > 0.0f) {
                eBest = eFromE;
                eTrace = 2;
            }
            E[cur] = eBest;
            tE[cur] = eTrace;

            float fFromM = M[left] - gapOpen;
            float fFromF = F[left] - gapExtend;
            float fBest = 0.0f;
            uint8_t fTrace = 0;
            if (fFromM >= fFromF && fFromM > 0.0f) {
                fBest = fFromM;
                fTrace = 1;
            } else if (fFromF > 0.0f) {
                fBest = fFromF;
                fTrace = 3;
            }
            F[cur] = fBest;
            tF[cur] = fTrace;

            const int qIdx = qLo + i - 1;
            const int b = encodeBase(targetSeq[static_cast<size_t>(j - 1)]);
            float emit = (b >= 0 && qIdx < static_cast<int>(profile.size()))
                ? profile[static_cast<size_t>(qIdx)][static_cast<size_t>(b)]
                : -3.0f;
            if (qIdx < static_cast<int>(queryStruct.size()) && j - 1 < static_cast<int>(targetStruct.size())) {
                const float sim = structureSimilarity(queryStruct[static_cast<size_t>(qIdx)],
                                                      targetStruct[static_cast<size_t>(j - 1)]);
                emit += structScoreWeight * 4.0f * (sim - (1.0f / 3.0f));
            }
            if (qIdx < static_cast<int>(stemGuide.anchorPosByQuery.size())
                && stemGuide.anchorPosByQuery[static_cast<size_t>(qIdx)] >= 0.0f
                && stemGuide.supportByQuery[static_cast<size_t>(qIdx)] > 0.0f) {
                const float dist = std::fabs(static_cast<float>(j - 1) - stemGuide.anchorPosByQuery[static_cast<size_t>(qIdx)]);
                const float radius = static_cast<float>(std::max(1, stemRadius));
                const float anchorScore = std::max(0.0f, 1.0f - dist / radius);
                emit += stemScoreWeight * stemGuide.supportByQuery[static_cast<size_t>(qIdx)] * anchorScore;
            }
            if (bonusMatrix != NULL && bonusMatrix->size() == static_cast<size_t>(n * m)) {
                emit += (*bonusMatrix)[static_cast<size_t>(i - 1) * static_cast<size_t>(m) + static_cast<size_t>(j - 1)];
            }

            float prev = 0.0f;
            uint8_t prevState = 0;
            if (M[diag] > prev) {
                prev = M[diag];
                prevState = 1;
            }
            if (E[diag] > prev) {
                prev = E[diag];
                prevState = 2;
            }
            if (F[diag] > prev) {
                prev = F[diag];
                prevState = 3;
            }

            const float mBest = prev + emit;
            if (mBest > 0.0f) {
                M[cur] = mBest;
                tM[cur] = prevState;
            } else {
                M[cur] = 0.0f;
                tM[cur] = 0;
            }

            if (M[cur] > bestScore) {
                bestScore = M[cur];
                bestI = i;
                bestJ = j;
                bestState = 1;
            }
            if (E[cur] > bestScore) {
                bestScore = E[cur];
                bestI = i;
                bestJ = j;
                bestState = 2;
            }
            if (F[cur] > bestScore) {
                bestScore = F[cur];
                bestI = i;
                bestJ = j;
                bestState = 3;
            }
        }
    }

    if (bestScore <= 0.0f || bestState == 0) {
        return out;
    }

    std::string bt;
    bt.reserve(static_cast<size_t>(n + m));
    int i = bestI;
    int j = bestJ;
    uint8_t state = bestState;
    int matches = 0;
    while (i > 0 && j >= 0 && state != 0) {
        const size_t cur = idx(i, j);
        if (state == 1) {
            if (j <= 0 || M[cur] <= 0.0f) {
                break;
            }
            bt.push_back('M');
            if (normalizeBase(querySeq[static_cast<size_t>(qLo + i - 1)]) == normalizeBase(targetSeq[static_cast<size_t>(j - 1)])) {
                ++matches;
            }
            state = tM[cur];
            --i;
            --j;
        } else if (state == 2) {
            if (E[cur] <= 0.0f) {
                break;
            }
            bt.push_back('I');
            state = tE[cur];
            --i;
        } else {
            if (j <= 0 || F[cur] <= 0.0f) {
                break;
            }
            bt.push_back('D');
            state = tF[cur];
            --j;
        }
    }

    if (bt.empty()) {
        return out;
    }

    std::reverse(bt.begin(), bt.end());
    out.valid = true;
    out.score = bestScore;
    out.qStart = qLo + i;
    out.qEnd = qLo + bestI - 1;
    out.dbStart = j;
    out.dbEnd = bestJ - 1;
    out.matches = matches;
    out.backtrace.swap(bt);
    return out;
}

static bool appendExactGlocalSegment(const std::string &querySeq,
                                     const std::vector<std::array<float, 4>> &profile,
                                     const std::vector<StructVec> &queryStruct,
                                     const std::vector<int> &partners,
                                     const std::string &targetSeq,
                                     const std::vector<StructVec> &targetStruct,
                                     int qLo,
                                     int qHi,
                                     int tLo,
                                     int tHi,
                                     float structScoreWeight,
                                     std::string &bt) {
    if (qLo > qHi) {
        for (int t = tLo; t <= tHi; ++t) {
            bt.push_back('D');
        }
        return true;
    }
    if (tLo > tHi) {
        for (int q = qLo; q <= qHi; ++q) {
            bt.push_back('I');
        }
        return true;
    }

    const std::string targetSub = targetSeq.substr(static_cast<size_t>(tLo),
                                                   static_cast<size_t>(tHi - tLo + 1));
    std::vector<StructVec> targetStructSub;
    if (!targetStruct.empty()) {
        targetStructSub.assign(targetStruct.begin() + tLo,
                               targetStruct.begin() + tHi + 1);
    }
    StemGuide neutralGuide;
    neutralGuide.partnerByQuery = partners;
    neutralGuide.anchorPosByQuery.assign(querySeq.size(), -1.0f);
    neutralGuide.supportByQuery.assign(querySeq.size(), 0.0f);
    AlignmentResult segAln = alignProfileGlocal(querySeq,
                                                profile,
                                                queryStruct,
                                                neutralGuide,
                                                targetSub,
                                                targetStructSub,
                                                qLo,
                                                qHi,
                                                structScoreWeight,
                                                0.0f,
                                                1,
                                                0,
                                                -1);
    if (!segAln.valid) {
        return false;
    }
    for (int d = 0; d < segAln.dbStart; ++d) {
        bt.push_back('D');
    }
    bt += segAln.backtrace;
    for (int d = segAln.dbEnd + 1; d < static_cast<int>(targetSub.size()); ++d) {
        bt.push_back('D');
    }
    return true;
}

struct PairAnchorCandidate {
    float score = NEG_INF_F;
    int ql = -1;
    int qr = -1;
    int tl = -1;
    int tr = -1;
};

struct SegmentStitchResult {
    bool valid = false;
    float score = NEG_INF_F;
    std::string bt;
};

static SegmentStitchResult buildExactGlocalSegment(const std::string &querySeq,
                                                   const std::vector<std::array<float, 4>> &profile,
                                                   const std::vector<StructVec> &queryStruct,
                                                   const std::vector<int> &partners,
                                                   const std::string &targetSeq,
                                                   const std::vector<StructVec> &targetStruct,
                                                   int qLo,
                                                   int qHi,
                                                   int tLo,
                                                   int tHi,
                                                   float structScoreWeight) {
    SegmentStitchResult out;
    if (qLo > qHi) {
        out.valid = true;
        out.score = -1.5f * static_cast<float>(std::max(0, tHi - tLo + 1));
        out.bt.assign(static_cast<size_t>(std::max(0, tHi - tLo + 1)), 'D');
        return out;
    }
    if (tLo > tHi) {
        out.valid = true;
        out.score = -1.5f * static_cast<float>(std::max(0, qHi - qLo + 1));
        out.bt.assign(static_cast<size_t>(std::max(0, qHi - qLo + 1)), 'I');
        return out;
    }

    const std::string targetSub = targetSeq.substr(static_cast<size_t>(tLo),
                                                   static_cast<size_t>(tHi - tLo + 1));
    std::vector<StructVec> targetStructSub;
    if (!targetStruct.empty()) {
        targetStructSub.assign(targetStruct.begin() + tLo,
                               targetStruct.begin() + tHi + 1);
    }
    StemGuide neutralGuide;
    neutralGuide.partnerByQuery = partners;
    neutralGuide.anchorPosByQuery.assign(querySeq.size(), -1.0f);
    neutralGuide.supportByQuery.assign(querySeq.size(), 0.0f);
    AlignmentResult segAln = alignProfileGlocal(querySeq,
                                                profile,
                                                queryStruct,
                                                neutralGuide,
                                                targetSub,
                                                targetStructSub,
                                                qLo,
                                                qHi,
                                                structScoreWeight,
                                                0.0f,
                                                1,
                                                0,
                                                -1);
    if (!segAln.valid) {
        return out;
    }
    out.valid = true;
    out.score = segAln.score
              - 0.75f * static_cast<float>(std::max(0, segAln.dbStart))
              - 0.75f * static_cast<float>(std::max(0, static_cast<int>(targetSub.size()) - 1 - segAln.dbEnd));
    out.bt.reserve(targetSub.size() + segAln.backtrace.size());
    for (int d = 0; d < segAln.dbStart; ++d) {
        out.bt.push_back('D');
    }
    out.bt += segAln.backtrace;
    for (int d = segAln.dbEnd + 1; d < static_cast<int>(targetSub.size()); ++d) {
        out.bt.push_back('D');
    }
    return out;
}

static float targetPairProbabilityAt(const WindowInfo *targetWindow,
                                     const std::vector<float> &targetPairMatrix,
                                     int targetLen,
                                     int i,
                                     int j) {
    if (targetWindow != NULL) {
        return windowPairMatrixAt(*targetWindow, i, j);
    }
    if (!targetPairMatrix.empty() && targetPairMatrix.size() == static_cast<size_t>(targetLen * targetLen)
        && i >= 0 && j >= 0 && i < targetLen && j < targetLen) {
        return targetPairMatrix[static_cast<size_t>(i * targetLen + j)];
    }
    return 0.0f;
}

static bool targetHasPairMatrix(const WindowInfo *targetWindow,
                                const std::vector<float> &targetPairMatrix) {
    return targetWindow != NULL ? windowHasPairMatrix(*targetWindow) : !targetPairMatrix.empty();
}

static std::vector<PairAnchorCandidate> collectPairAnchorCandidates(const std::string &querySeq,
                                                                    const std::vector<std::array<float, 4>> &profile,
                                                                    const std::vector<StructVec> &queryStruct,
                                                                    const std::vector<int> &queryPartners,
                                                                    const std::string &targetSeq,
                                                                    const std::vector<StructVec> &targetStruct,
                                                                    const std::vector<int> &targetPartners,
                                                                    const std::vector<float> &targetPairMatrix,
                                                                    const WindowInfo *targetWindow,
                                                                    int qLo,
                                                                    int qHi,
                                                                    int tLo,
                                                                    int tHi) {
    std::vector<PairAnchorCandidate> anchors;
    if (qLo >= qHi || tLo >= tHi || queryPartners.empty() || targetPartners.empty()) {
        return anchors;
    }

    const int minPairSpan = std::max(2, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_MIN_SPAN", 4));
    const int maxQSpan = std::max(minPairSpan, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_MAX_QUERY_SPAN", 192));
    const size_t maxQueryPairs = static_cast<size_t>(std::max(1, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_TOP_QPAIRS", 8)));
    const size_t maxTargetPairsPerQuery = static_cast<size_t>(std::max(1, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_TOP_TARGETS", 4)));
    const size_t maxAnchors = static_cast<size_t>(std::max(1, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_TOP_ANCHORS", 6)));
    const float pairStructWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_PAIR_STITCH_STRUCT_WEIGHT", 3.0f));
    const float pairSeqWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_PAIR_STITCH_SEQ_WEIGHT", 1.5f));
    const float pairCanonicalBonus = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_PAIR_STITCH_CANONICAL", 0.5f));
    const float minAnchorScore = parseEnvFloat("MMSEQS_CMLITE_PAIR_STITCH_MIN_SCORE", 1.5f);
    const int targetLen = static_cast<int>(targetSeq.size());

    struct QueryPairCand {
        float score;
        int ql;
        int qr;
    };
    std::vector<QueryPairCand> queryPairs;
    for (int ql = qLo; ql <= qHi; ++ql) {
        if (ql < 0 || ql >= static_cast<int>(queryPartners.size())) {
            continue;
        }
        const int qr = queryPartners[static_cast<size_t>(ql)];
        if (qr <= ql || qr < qLo || qr > qHi) {
            continue;
        }
        const int span = qr - ql;
        if (span < minPairSpan || span > maxQSpan) {
            continue;
        }
        float qScore = 0.0f;
        if (ql < static_cast<int>(queryStruct.size())) {
            qScore += 1.0f - queryStruct[static_cast<size_t>(ql)].u;
        }
        if (qr < static_cast<int>(queryStruct.size())) {
            qScore += 1.0f - queryStruct[static_cast<size_t>(qr)].u;
        }
        qScore += 0.02f * static_cast<float>(span);
        queryPairs.push_back(QueryPairCand{qScore, ql, qr});
    }
    std::sort(queryPairs.begin(), queryPairs.end(),
              [](const QueryPairCand &a, const QueryPairCand &b) {
                  if (a.score != b.score) {
                      return a.score > b.score;
                  }
                  return a.ql < b.ql;
              });
    if (queryPairs.size() > maxQueryPairs) {
        queryPairs.resize(maxQueryPairs);
    }

    for (size_t qi = 0; qi < queryPairs.size(); ++qi) {
        const int ql = queryPairs[qi].ql;
        const int qr = queryPairs[qi].qr;
        std::vector<PairAnchorCandidate> local;
        for (int tl = tLo; tl <= tHi; ++tl) {
            if (tl < 0 || tl >= static_cast<int>(targetPartners.size())) {
                continue;
            }
            const int tr = targetPartners[static_cast<size_t>(tl)];
            if (tr <= tl || tr < tLo || tr > tHi || tr >= targetLen) {
                continue;
            }
            const int tSpan = tr - tl;
            if (tSpan < minPairSpan) {
                continue;
            }
            const int qSpan = qr - ql;
            const float spanRatio = static_cast<float>(std::min(qSpan, tSpan))
                                  / static_cast<float>(std::max(qSpan, tSpan));
            if (spanRatio < 0.35f) {
                continue;
            }
            const int lb = encodeBase(targetSeq[static_cast<size_t>(tl)]);
            const int rb = encodeBase(targetSeq[static_cast<size_t>(tr)]);
            if (lb < 0 || rb < 0) {
                continue;
            }
            float score = queryPairs[qi].score;
            if (ql < static_cast<int>(profile.size())) {
                score += pairSeqWeight * basePreferenceFromProfile(profile[static_cast<size_t>(ql)], lb);
            }
            if (qr < static_cast<int>(profile.size())) {
                score += pairSeqWeight * basePreferenceFromProfile(profile[static_cast<size_t>(qr)], rb);
            }
            if (targetHasPairMatrix(targetWindow, targetPairMatrix)) {
                score += pairStructWeight * targetPairProbabilityAt(targetWindow, targetPairMatrix, targetLen, tl, tr);
            } else {
                score += pairStructWeight;
            }
            if (tl < static_cast<int>(targetStruct.size())) {
                score += 0.5f * (1.0f - targetStruct[static_cast<size_t>(tl)].u);
            }
            if (tr < static_cast<int>(targetStruct.size())) {
                score += 0.5f * (1.0f - targetStruct[static_cast<size_t>(tr)].u);
            }
            if (isCanonicalPair(targetSeq[static_cast<size_t>(tl)], targetSeq[static_cast<size_t>(tr)])) {
                score += pairCanonicalBonus;
            }
            if (score >= minAnchorScore) {
                PairAnchorCandidate cand;
                cand.score = score;
                cand.ql = ql;
                cand.qr = qr;
                cand.tl = tl;
                cand.tr = tr;
                local.push_back(cand);
            }
        }
        std::sort(local.begin(), local.end(),
                  [](const PairAnchorCandidate &a, const PairAnchorCandidate &b) {
                      if (a.score != b.score) {
                          return a.score > b.score;
                      }
                      if (a.ql != b.ql) {
                          return a.ql < b.ql;
                      }
                      return a.tl < b.tl;
                  });
        if (local.size() > maxTargetPairsPerQuery) {
            local.resize(maxTargetPairsPerQuery);
        }
        anchors.insert(anchors.end(), local.begin(), local.end());
    }

    std::sort(anchors.begin(), anchors.end(),
              [](const PairAnchorCandidate &a, const PairAnchorCandidate &b) {
                  if (a.score != b.score) {
                      return a.score > b.score;
                  }
                  if (a.ql != b.ql) {
                      return a.ql < b.ql;
                  }
                  return a.tl < b.tl;
              });
    std::vector<PairAnchorCandidate> filtered;
    filtered.reserve(std::min(maxAnchors, anchors.size()));
    for (size_t i = 0; i < anchors.size() && filtered.size() < maxAnchors; ++i) {
        bool nearDuplicate = false;
        for (size_t j = 0; j < filtered.size(); ++j) {
            if (std::abs(anchors[i].ql - filtered[j].ql) <= 1
                && std::abs(anchors[i].qr - filtered[j].qr) <= 1
                && std::abs(anchors[i].tl - filtered[j].tl) <= 2
                && std::abs(anchors[i].tr - filtered[j].tr) <= 2) {
                nearDuplicate = true;
                break;
            }
        }
        if (!nearDuplicate) {
            filtered.push_back(anchors[i]);
        }
    }
    return filtered;
}

static SegmentStitchResult buildPairAwareSegmentBest(const std::string &querySeq,
                                                     const std::vector<std::array<float, 4>> &profile,
                                                     const std::vector<StructVec> &queryStruct,
                                                     const std::vector<int> &queryPartners,
                                                     const std::string &targetSeq,
                                                     const std::vector<StructVec> &targetStruct,
                                                     const std::vector<int> &targetPartners,
                                                     const std::vector<float> &targetPairMatrix,
                                                     const WindowInfo *targetWindow,
                                                     int qLo,
                                                     int qHi,
                                                     int tLo,
                                                     int tHi,
                                                     float structScoreWeight,
                                                     int depth) {
    const int maxDepth = std::max(0, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_MAX_DEPTH", 8));
    const int maxSpan = std::max(8, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_MAX_SEGMENT", 192));
    const size_t maxAnchorBeam = static_cast<size_t>(std::max(1, parseEnvInt("MMSEQS_CMLITE_PAIR_STITCH_ANCHOR_BEAM", 3)));
    SegmentStitchResult fallback = buildExactGlocalSegment(querySeq,
                                                           profile,
                                                           queryStruct,
                                                           queryPartners,
                                                           targetSeq,
                                                           targetStruct,
                                                           qLo,
                                                           qHi,
                                                           tLo,
                                                           tHi,
                                                           structScoreWeight);
    if (depth >= maxDepth || (qHi - qLo + 1) > maxSpan || targetPartners.empty()) {
        return fallback;
    }

    const std::vector<PairAnchorCandidate> anchors = collectPairAnchorCandidates(querySeq,
                                                                                 profile,
                                                                                 queryStruct,
                                                                                 queryPartners,
                                                                                 targetSeq,
                                                                                 targetStruct,
                                                                                 targetPartners,
                                                                                 targetPairMatrix,
                                                                                 targetWindow,
                                                                                 qLo,
                                                                                 qHi,
                                                                                 tLo,
                                                                                 tHi);
    if (anchors.empty()) {
        return fallback;
    }

    SegmentStitchResult best = fallback;
    for (size_t ai = 0; ai < anchors.size() && ai < maxAnchorBeam; ++ai) {
        const PairAnchorCandidate &anchor = anchors[ai];
        SegmentStitchResult left = buildPairAwareSegmentBest(querySeq,
                                                             profile,
                                                             queryStruct,
                                                             queryPartners,
                                                             targetSeq,
                                                             targetStruct,
                                                             targetPartners,
                                                             targetPairMatrix,
                                                             targetWindow,
                                                             qLo,
                                                             anchor.ql - 1,
                                                             tLo,
                                                             anchor.tl - 1,
                                                             structScoreWeight,
                                                             depth + 1);
        if (!left.valid) {
            continue;
        }
        SegmentStitchResult mid = buildPairAwareSegmentBest(querySeq,
                                                            profile,
                                                            queryStruct,
                                                            queryPartners,
                                                            targetSeq,
                                                            targetStruct,
                                                            targetPartners,
                                                            targetPairMatrix,
                                                            targetWindow,
                                                            anchor.ql + 1,
                                                            anchor.qr - 1,
                                                            anchor.tl + 1,
                                                            anchor.tr - 1,
                                                            structScoreWeight,
                                                            depth + 1);
        if (!mid.valid) {
            continue;
        }
        SegmentStitchResult right = buildPairAwareSegmentBest(querySeq,
                                                              profile,
                                                              queryStruct,
                                                              queryPartners,
                                                              targetSeq,
                                                              targetStruct,
                                                              targetPartners,
                                                              targetPairMatrix,
                                                              targetWindow,
                                                              anchor.qr + 1,
                                                              qHi,
                                                              anchor.tr + 1,
                                                              tHi,
                                                              structScoreWeight,
                                                              depth + 1);
        if (!right.valid) {
            continue;
        }

        SegmentStitchResult combined;
        combined.valid = true;
        combined.score = left.score + anchor.score + mid.score + anchor.score + right.score;
        combined.bt.reserve(left.bt.size() + mid.bt.size() + right.bt.size() + 2);
        combined.bt += left.bt;
        combined.bt.push_back('M');
        combined.bt += mid.bt;
        combined.bt.push_back('M');
        combined.bt += right.bt;
        if (!best.valid || combined.score > best.score) {
            best = combined;
        }
    }

    return best;
}

static bool appendPairAwareSegment(const std::string &querySeq,
                                   const std::vector<std::array<float, 4>> &profile,
                                   const std::vector<StructVec> &queryStruct,
                                   const std::vector<int> &queryPartners,
                                   const std::string &targetSeq,
                                   const std::vector<StructVec> &targetStruct,
                                   const std::vector<int> &targetPartners,
                                   const std::vector<float> &targetPairMatrix,
                                   const WindowInfo *targetWindow,
                                   int qLo,
                                   int qHi,
                                   int tLo,
                                   int tHi,
                                   float structScoreWeight,
                                   int depth,
                                   std::string &bt) {
    SegmentStitchResult result = buildPairAwareSegmentBest(querySeq,
                                                           profile,
                                                           queryStruct,
                                                           queryPartners,
                                                           targetSeq,
                                                           targetStruct,
                                                           targetPartners,
                                                           targetPairMatrix,
                                                           targetWindow,
                                                           qLo,
                                                           qHi,
                                                           tLo,
                                                           tHi,
                                                           structScoreWeight,
                                                           depth);
    if (!result.valid) {
        return false;
    }
    bt += result.bt;
    return true;
}

static int estimateSeedBandRadius(const CandidateHit &cand,
                                  int bandExtra) {
    return bandExtra + std::abs((cand.window.localSeedEnd - cand.window.localSeedStart)
                              - (cand.oriented.qEnd - cand.oriented.qStart));
}

static void uniqueSorted(std::vector<int> &values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

static std::vector<StemMatch> findStemMatches(const std::vector<StemPattern> &patterns,
                                              const WindowInfo &window,
                                              int seedDiag,
                                              int seedBandRadius) {
    std::vector<StemMatch> matches;
    if (patterns.empty() || window.seq.empty()) {
        return matches;
    }

    const int spanSlack = std::max(2, parseEnvInt("MMSEQS_CMLITE_STEM_SPAN_SLACK", 12));
    const int diagSlack = std::max(8, parseEnvInt("MMSEQS_CMLITE_STEM_DIAG_SLACK", 16));
    const int candidateSlack = std::max(2, parseEnvInt("MMSEQS_CMLITE_STEM_PARTNER_SLACK", 6));
    const int maxPerPattern = std::max(4, parseEnvInt("MMSEQS_CMLITE_STEM_MAX_PER_PATTERN", 24));
    const float minSpanRatio = clampf(parseEnvFloat("MMSEQS_CMLITE_STEM_MIN_SPAN_RATIO", 0.55f), 0.0f, 1.0f);
    const float minArmSupport = clampf(parseEnvFloat("MMSEQS_CMLITE_STEM_MIN_ARM_SUPPORT", 0.40f), 0.0f, 1.0f);
    const int targetLen = static_cast<int>(window.seq.size());

    for (size_t p = 0; p < patterns.size(); ++p) {
        const StemPattern &pat = patterns[p];
        const int expectedRightStart = pat.rightOuter - pat.armLen + 1;
        std::vector<StemMatch> local;
        for (int leftPos = 0; leftPos + pat.armLen <= targetLen; ++leftPos) {
            if (seedBandRadius >= 0) {
                const int diag = leftPos - pat.leftStart;
                if (std::abs(diag - seedDiag) > seedBandRadius + diagSlack) {
                    continue;
                }
            }
            const float leftScore = pat.leftProfile.empty()
                ? ((matchArmAt(window.seq, leftPos, pat.leftArm, 1) >= 0) ? 1.0f : -1.0f)
                : scoreArmAt(window.seq, leftPos, pat.leftProfile, minArmSupport);
            if (leftScore < 0.0f) {
                continue;
            }

            std::vector<int> rightStarts;
            for (int k = 0; k < pat.armLen; ++k) {
                const int partner = (leftPos + k < static_cast<int>(window.partners.size()))
                    ? window.partners[static_cast<size_t>(leftPos + k)]
                    : -1;
                if (partner >= 0) {
                    const int rightStart = partner - (pat.armLen - 1 - k);
                    if (rightStart >= 0 && rightStart + pat.armLen <= targetLen) {
                        rightStarts.push_back(rightStart);
                    }
                }
            }
            const int expectedStart = leftPos + (expectedRightStart - pat.leftStart);
            if (rightStarts.empty()) {
                for (int rightStart = expectedStart - spanSlack;
                     rightStart <= expectedStart + spanSlack;
                     ++rightStart) {
                    if (rightStart >= 0 && rightStart + pat.armLen <= targetLen) {
                        rightStarts.push_back(rightStart);
                    }
                }
            } else {
                uniqueSorted(rightStarts);
                std::vector<int> expanded;
                for (size_t i = 0; i < rightStarts.size(); ++i) {
                    for (int delta = -candidateSlack; delta <= candidateSlack; ++delta) {
                        const int candidate = rightStarts[i] + delta;
                        if (candidate >= 0 && candidate + pat.armLen <= targetLen) {
                            expanded.push_back(candidate);
                        }
                    }
                }
                rightStarts.swap(expanded);
            }
            uniqueSorted(rightStarts);

            for (size_t r = 0; r < rightStarts.size(); ++r) {
                const int rightStart = rightStarts[r];
                const float rightScore = pat.rightProfile.empty()
                    ? ((matchArmAt(window.seq, rightStart, pat.rightArm, 1) >= 0) ? 1.0f : -1.0f)
                    : scoreArmAt(window.seq, rightStart, pat.rightProfile, minArmSupport);
                if (rightScore < 0.0f) {
                    continue;
                }
                const int rightOuter = rightStart + pat.armLen - 1;
                if (rightOuter <= leftPos) {
                    continue;
                }

                const int dbSpan = rightOuter - leftPos;
                const float spanRatio = static_cast<float>(std::min(pat.span, dbSpan))
                                      / static_cast<float>(std::max(pat.span, dbSpan));
                if (spanRatio < minSpanRatio) {
                    continue;
                }

                const int diagLeft = leftPos - pat.leftStart;
                const int diagRight = rightStart - expectedRightStart;
                const float drift = static_cast<float>(std::abs(diagLeft - diagRight));
                if (seedBandRadius >= 0 && drift > static_cast<float>(seedBandRadius + candidateSlack)) {
                    continue;
                }

                float pairSupport = 0.0f;
                for (int k = 0; k < pat.armLen; ++k) {
                    const int tl = leftPos + k;
                    const int tr = rightOuter - k;
                    const float pairProb = windowPairMatrixAt(window, tl, tr);
                    if (pairProb > 0.0f) {
                        pairSupport += pairProb;
                    } else if (window.partners.size() == static_cast<size_t>(targetLen)
                               && window.partners[static_cast<size_t>(tl)] == tr) {
                        pairSupport += 1.0f;
                    }
                }
                pairSupport /= static_cast<float>(std::max(1, pat.armLen));
                const int diagCenter = static_cast<int>(std::lrint(0.5f * static_cast<float>(diagLeft + diagRight)));
                const float seedCompat = (seedBandRadius >= 0)
                    ? std::max(0.0f, 1.0f - std::fabs(static_cast<float>(diagCenter - seedDiag))
                                        / static_cast<float>(seedBandRadius + diagSlack))
                    : 1.0f;
                const float moduleScore = scoreStemModuleAt(pat, window, leftPos, rightStart, leftScore, rightScore);
                if (moduleScore < 0.0f) {
                    continue;
                }
                const float support = clampf(0.30f * (leftScore + rightScore) * 0.5f
                                           + 0.45f * pairSupport
                                           + 0.15f * spanRatio
                                           + 0.10f * seedCompat, 0.0f, 1.5f);
                StemMatch match;
                match.patternIdx = static_cast<int>(p);
                match.leftQuery = pat.leftStart;
                match.rightQuery = pat.rightOuter;
                match.leftTarget = leftPos;
                match.rightTarget = rightOuter;
                match.diagCenter = diagCenter;
                match.support = support;
                match.score = pat.priority
                            * (0.35f + 0.65f * std::max(0.0f, pat.pairSupport))
                            * (0.25f + 0.75f * support)
                            * (0.30f + 0.70f * moduleScore)
                            / (1.0f + 0.08f * drift);
                local.push_back(match);
            }
        }

        std::sort(local.begin(), local.end(),
                  [](const StemMatch &a, const StemMatch &b) {
                      if (a.score != b.score) {
                          return a.score > b.score;
                      }
                      if (a.leftTarget != b.leftTarget) {
                          return a.leftTarget < b.leftTarget;
                      }
                      return a.rightTarget < b.rightTarget;
                  });
        if (static_cast<int>(local.size()) > maxPerPattern) {
            local.resize(static_cast<size_t>(maxPerPattern));
        }
        matches.insert(matches.end(), local.begin(), local.end());
    }

    std::sort(matches.begin(), matches.end(),
              [](const StemMatch &a, const StemMatch &b) {
                  if (a.leftQuery != b.leftQuery) {
                      return a.leftQuery < b.leftQuery;
                  }
                  if (a.leftTarget != b.leftTarget) {
                      return a.leftTarget < b.leftTarget;
                  }
                  return a.score > b.score;
              });
    return matches;
}

static float computeChainTransitionPenalty(const StemPattern &prevPattern,
                                           const StemPattern &currPattern,
                                           const StemMatch &prevMatch,
                                           const StemMatch &currMatch) {
    const float gapWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_CHAIN_GAP_WEIGHT", 0.12f));
    const float diagWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_CHAIN_DIAG_WEIGHT", 0.10f));
    const float overlapPenalty = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_CHAIN_OVERLAP_PENALTY", 1.0f));

    const int qGap = currPattern.leftStart - (prevPattern.leftStart + prevPattern.stemLen);
    const int tGap = currMatch.leftTarget - (prevMatch.leftTarget + prevPattern.stemLen);
    float penalty = 0.0f;
    if (qGap < -1 || tGap < -1) {
        penalty += overlapPenalty * static_cast<float>(std::max(-qGap, -tGap));
    } else {
        penalty += gapWeight * static_cast<float>(std::abs(tGap - qGap));
    }
    penalty += diagWeight * static_cast<float>(std::abs(currMatch.diagCenter - prevMatch.diagCenter));
    return penalty;
}

static bool areChainCompatible(const StemPattern &prevPattern,
                               const StemPattern &currPattern,
                               const StemMatch &prevMatch,
                               const StemMatch &currMatch) {
    if (currPattern.leftStart <= prevPattern.leftStart) {
        return false;
    }
    if (currMatch.leftTarget <= prevMatch.leftTarget) {
        return false;
    }
    if (currMatch.rightTarget <= prevMatch.rightTarget) {
        return false;
    }
    return true;
}

static void keepBestChainStates(std::vector<ChainState> &bucket,
                                const ChainState &candidate,
                                size_t keepLimit) {
    bucket.push_back(candidate);
    std::sort(bucket.begin(), bucket.end(),
              [](const ChainState &a, const ChainState &b) {
                  if (a.score != b.score) {
                      return a.score > b.score;
                  }
                  return a.hitIdx < b.hitIdx;
              });
    if (bucket.size() > keepLimit) {
        bucket.resize(keepLimit);
    }
}

static std::vector<std::vector<int>> buildTopChainsFromStemMatches(const std::vector<StemPattern> &patterns,
                                                                   const std::vector<StemMatch> &matches) {
    std::vector<std::vector<int>> chains;
    if (matches.empty()) {
        return chains;
    }

    const size_t keepPerEnd = static_cast<size_t>(std::max(2, parseEnvInt("MMSEQS_CMLITE_CHAIN_BEAM", 8)));
    const size_t maxChains = static_cast<size_t>(std::max(1, parseEnvInt("MMSEQS_CMLITE_TOP_CHAINS", 6)));
    const float chainMinScore = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_CHAIN_MIN_SCORE", 0.75f));

    std::vector<std::vector<ChainState>> bestByEnd(matches.size());
    for (size_t i = 0; i < matches.size(); ++i) {
        ChainState seed;
        seed.score = matches[i].score;
        seed.hitIdx = static_cast<int>(i);
        keepBestChainStates(bestByEnd[i], seed, keepPerEnd);
        for (size_t j = 0; j < i; ++j) {
            const StemPattern &prevPattern = patterns[static_cast<size_t>(matches[j].patternIdx)];
            const StemPattern &currPattern = patterns[static_cast<size_t>(matches[i].patternIdx)];
            if (!areChainCompatible(prevPattern, currPattern, matches[j], matches[i])) {
                continue;
            }
            const float penalty = computeChainTransitionPenalty(prevPattern, currPattern, matches[j], matches[i]);
            for (size_t rank = 0; rank < bestByEnd[j].size(); ++rank) {
                const ChainState &prev = bestByEnd[j][rank];
                ChainState next;
                next.score = prev.score + matches[i].score - penalty;
                next.hitIdx = static_cast<int>(i);
                next.prevHitIdx = static_cast<int>(j);
                next.prevRank = static_cast<int>(rank);
                keepBestChainStates(bestByEnd[i], next, keepPerEnd);
            }
        }
    }

    struct EndStateRef {
        float score;
        int endIdx;
        int rank;
    };
    std::vector<EndStateRef> finals;
    for (size_t i = 0; i < bestByEnd.size(); ++i) {
        for (size_t rank = 0; rank < bestByEnd[i].size(); ++rank) {
            if (bestByEnd[i][rank].score >= chainMinScore) {
                EndStateRef ref;
                ref.score = bestByEnd[i][rank].score;
                ref.endIdx = static_cast<int>(i);
                ref.rank = static_cast<int>(rank);
                finals.push_back(ref);
            }
        }
    }
    std::sort(finals.begin(), finals.end(),
              [](const EndStateRef &a, const EndStateRef &b) {
                  if (a.score != b.score) {
                      return a.score > b.score;
                  }
                  return a.endIdx < b.endIdx;
              });

    std::unordered_set<std::string> seen;
    for (size_t fi = 0; fi < finals.size() && chains.size() < maxChains; ++fi) {
        std::vector<int> chain;
        int endIdx = finals[fi].endIdx;
        int rank = finals[fi].rank;
        while (endIdx >= 0 && rank >= 0) {
            const ChainState &state = bestByEnd[static_cast<size_t>(endIdx)][static_cast<size_t>(rank)];
            chain.push_back(state.hitIdx);
            endIdx = state.prevHitIdx;
            rank = state.prevRank;
        }
        std::reverse(chain.begin(), chain.end());
        std::ostringstream oss;
        for (size_t ci = 0; ci < chain.size(); ++ci) {
            if (ci != 0) {
                oss << ',';
            }
            oss << chain[ci];
        }
        if (seen.insert(oss.str()).second) {
            chains.push_back(chain);
        }
    }

    return chains;
}

static bool buildGuideFromStemChain(int qLen,
                                    const std::vector<StemPattern> &patterns,
                                    const std::vector<StemMatch> &matches,
                                    int seedDiag,
                                    int seedBandRadius,
                                    const std::vector<int> &partners,
                                    const std::vector<int> &chain,
                                    StemGuide &guide,
                                    int &diagCenter,
                                    int &bandRadius,
                                    float &chainScore) {
    guide.partnerByQuery = partners;
    guide.anchorPosByQuery.assign(static_cast<size_t>(qLen), -1.0f);
    guide.supportByQuery.assign(static_cast<size_t>(qLen), 0.0f);
    diagCenter = seedDiag;
    bandRadius = seedBandRadius;
    chainScore = 0.0f;
    if (matches.empty() || chain.empty()) {
        return false;
    }

    std::vector<float> posSum(static_cast<size_t>(qLen), 0.0f);
    std::vector<float> weightSum(static_cast<size_t>(qLen), 0.0f);
    float diagSum = 0.0f;
    float diagWeightSum = 0.0f;
    int minDiag = INT_MAX;
    int maxDiag = INT_MIN;
    for (size_t ci = 0; ci < chain.size(); ++ci) {
        const StemMatch &match = matches[static_cast<size_t>(chain[ci])];
        const StemPattern &pat = patterns[static_cast<size_t>(match.patternIdx)];
        for (int k = 0; k < pat.stemLen; ++k) {
            const int qLeft = pat.leftStart + k;
            const int qRight = pat.rightOuter - k;
            const int tLeft = match.leftTarget + k;
            const int tRight = match.rightTarget - k;
            if (qLeft >= 0 && qLeft < qLen) {
                posSum[static_cast<size_t>(qLeft)] += match.support * static_cast<float>(tLeft);
                weightSum[static_cast<size_t>(qLeft)] += match.support;
            }
            if (qRight >= 0 && qRight < qLen) {
                posSum[static_cast<size_t>(qRight)] += match.support * static_cast<float>(tRight);
                weightSum[static_cast<size_t>(qRight)] += match.support;
            }
        }
        diagSum += match.support * static_cast<float>(match.diagCenter);
        diagWeightSum += match.support;
        minDiag = std::min(minDiag, match.diagCenter);
        maxDiag = std::max(maxDiag, match.diagCenter);
    }

    bool hasSupport = false;
    for (int i = 0; i < qLen; ++i) {
        if (weightSum[static_cast<size_t>(i)] > 0.0f) {
            guide.anchorPosByQuery[static_cast<size_t>(i)] = posSum[static_cast<size_t>(i)] / weightSum[static_cast<size_t>(i)];
            guide.supportByQuery[static_cast<size_t>(i)] = clampf(weightSum[static_cast<size_t>(i)] / static_cast<float>(std::max(1U, static_cast<unsigned int>(chain.size()))), 0.0f, 1.0f);
            hasSupport = true;
        }
    }
    if (!hasSupport) {
        return false;
    }

    const float bridgeSupportScale = clampf(parseEnvFloat("MMSEQS_CMLITE_GUIDE_BRIDGE_SCALE", 0.35f), 0.0f, 1.0f);
    const int maxBridgeGap = std::max(1, parseEnvInt("MMSEQS_CMLITE_GUIDE_BRIDGE_MAX_GAP", 24));
    int prevIdx = -1;
    for (int i = 0; i < qLen; ++i) {
        if (guide.anchorPosByQuery[static_cast<size_t>(i)] < 0.0f) {
            continue;
        }
        if (prevIdx >= 0) {
            const int qGap = i - prevIdx;
            if (qGap > 1 && qGap <= maxBridgeGap) {
                const float prevPos = guide.anchorPosByQuery[static_cast<size_t>(prevIdx)];
                const float currPos = guide.anchorPosByQuery[static_cast<size_t>(i)];
                const float prevSup = guide.supportByQuery[static_cast<size_t>(prevIdx)];
                const float currSup = guide.supportByQuery[static_cast<size_t>(i)];
                for (int mid = prevIdx + 1; mid < i; ++mid) {
                    const float frac = static_cast<float>(mid - prevIdx) / static_cast<float>(qGap);
                    guide.anchorPosByQuery[static_cast<size_t>(mid)] = prevPos + frac * (currPos - prevPos);
                    guide.supportByQuery[static_cast<size_t>(mid)] = bridgeSupportScale * std::min(prevSup, currSup);
                }
            }
        }
        prevIdx = i;
    }

    diagCenter = (diagWeightSum > 0.0f)
        ? static_cast<int>(std::lrint(diagSum / diagWeightSum))
        : seedDiag;
    const int inferredBand = std::max(8, (maxDiag >= minDiag) ? (maxDiag - minDiag + 10) : 16);
    bandRadius = (seedBandRadius >= 0) ? std::min(seedBandRadius, inferredBand) : inferredBand;
    float score = 0.0f;
    for (size_t ci = 0; ci < chain.size(); ++ci) {
        const StemMatch &curr = matches[static_cast<size_t>(chain[ci])];
        score += curr.score;
        if (ci > 0) {
            const StemMatch &prev = matches[static_cast<size_t>(chain[ci - 1])];
            const StemPattern &prevPattern = patterns[static_cast<size_t>(prev.patternIdx)];
            const StemPattern &currPattern = patterns[static_cast<size_t>(curr.patternIdx)];
            score -= computeChainTransitionPenalty(prevPattern, currPattern, prev, curr);
        }
    }
    chainScore = score;
    return true;
}

static StemGuide makeNeutralGuide(int qLen,
                                  const std::vector<int> &partners) {
    StemGuide guide;
    guide.partnerByQuery = partners;
    guide.anchorPosByQuery.assign(static_cast<size_t>(qLen), -1.0f);
    guide.supportByQuery.assign(static_cast<size_t>(qLen), 0.0f);
    return guide;
}

static float windowPairSupportAt(const WindowInfo &window,
                                 int a,
                                 int b) {
    if (a < 0 || b < 0 || a >= static_cast<int>(window.seq.size()) || b >= static_cast<int>(window.seq.size()) || a == b) {
        return 0.0f;
    }
    const int lo = std::min(a, b);
    const int hi = std::max(a, b);
    const int targetLen = static_cast<int>(window.seq.size());
    const float pairProb = windowPairMatrixAt(window, lo, hi);
    if (pairProb > 0.0f) {
        return pairProb;
    }
    if (window.partners.size() == static_cast<size_t>(targetLen) && window.partners[static_cast<size_t>(lo)] == hi) {
        return 1.0f;
    }
    return 0.0f;
}

static inline float logSumExp2(float a, float b) {
    if (a <= NEG_INF_F) return b;
    if (b <= NEG_INF_F) return a;
    const float mx = std::max(a, b);
    return mx + std::log(std::exp(a - mx) + std::exp(b - mx));
}

static float logSumExp4(float a, float b, float c, float d) {
    return logSumExp2(logSumExp2(a, b), logSumExp2(c, d));
}

static std::vector<float> buildLolScoreMatrix(const std::string &querySeq,
                                              const std::vector<std::array<float, 4>> &profile,
                                              const std::vector<StructVec> &queryStruct,
                                              const std::vector<int> &partners,
                                              const std::vector<float> &pairSupport,
                                              const StemGuide &stemGuide,
                                              const WindowInfo &window,
                                              const std::string &targetSeq,
                                              const std::vector<StructVec> &targetStruct,
                                              int queryStart,
                                              int queryEnd,
                                              float structScoreWeight,
                                              float stemScoreWeight,
                                              int stemRadius,
                                              int diagCenter,
                                              int bandRadius,
                                              const std::vector<float> *bonusMatrix) {
    const int fullQueryLen = static_cast<int>(querySeq.size());
    const int qLo = std::max(0, std::min(queryStart, fullQueryLen - 1));
    const int qHi = std::max(qLo, std::min(queryEnd, fullQueryLen - 1));
    const int n = qHi - qLo + 1;
    const int m = static_cast<int>(targetSeq.size());
    std::vector<float> score(static_cast<size_t>(std::max(0, n * m)), NEG_INF_F);
    if (n <= 0 || m <= 0) {
        return score;
    }
    const bool pairMatrixOnly = targetStruct.empty() && windowHasPairMatrix(window);
    const float cellPairWeight = pairMatrixOnly
        ? clampf(parseEnvFloat("MMSEQS_LOLALIGN_CELL_PAIR_WEIGHT", 0.60f), 0.0f, 4.0f)
        : 0.0f;
    const int pairAnchorRadius = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_CELL_PAIR_ANCHOR_RADIUS", 1));

    for (int i = 0; i < n; ++i) {
        const int qIdx = qLo + i;
        int jLo = 0;
        int jHi = m - 1;
        if (bandRadius >= 0) {
            jLo = std::max(0, i + diagCenter - bandRadius);
            jHi = std::min(m - 1, i + diagCenter + bandRadius);
            if (jLo > jHi) {
                continue;
            }
        }
        for (int j = jLo; j <= jHi; ++j) {
            const int b = encodeBase(targetSeq[static_cast<size_t>(j)]);
            float emit = (b >= 0 && qIdx < static_cast<int>(profile.size()))
                ? profile[static_cast<size_t>(qIdx)][static_cast<size_t>(b)]
                : -3.0f;
            if (!pairMatrixOnly
                && qIdx < static_cast<int>(queryStruct.size())
                && j < static_cast<int>(targetStruct.size())) {
                const float sim = structureSimilarity(queryStruct[static_cast<size_t>(qIdx)],
                                                      targetStruct[static_cast<size_t>(j)]);
                emit += structScoreWeight * 4.0f * (sim - (1.0f / 3.0f));
            }
            if (qIdx < static_cast<int>(stemGuide.anchorPosByQuery.size())
                && stemGuide.anchorPosByQuery[static_cast<size_t>(qIdx)] >= 0.0f
                && stemGuide.supportByQuery[static_cast<size_t>(qIdx)] > 0.0f) {
                const float dist = std::fabs(static_cast<float>(j) - stemGuide.anchorPosByQuery[static_cast<size_t>(qIdx)]);
                const float radius = static_cast<float>(std::max(1, stemRadius));
                const float anchorScore = std::max(0.0f, 1.0f - dist / radius);
                emit += stemScoreWeight * stemGuide.supportByQuery[static_cast<size_t>(qIdx)] * anchorScore;
            }
            if (cellPairWeight > 0.0f && qIdx < static_cast<int>(partners.size())) {
                const int partner = partners[static_cast<size_t>(qIdx)];
                if (partner >= 0
                    && partner < fullQueryLen
                    && partner < static_cast<int>(stemGuide.anchorPosByQuery.size())
                    && stemGuide.anchorPosByQuery[static_cast<size_t>(partner)] >= 0.0f) {
                    const int targetAnchor = static_cast<int>(std::lrint(stemGuide.anchorPosByQuery[static_cast<size_t>(partner)]));
                    float targetPair = 0.0f;
                    for (int da = -pairAnchorRadius; da <= pairAnchorRadius; ++da) {
                        targetPair = std::max(targetPair, windowPairSupportAt(window, j, targetAnchor + da));
                    }
                    float queryPair = 0.0f;
                    if (pairSupport.size() == static_cast<size_t>(fullQueryLen * fullQueryLen)) {
                        queryPair = std::max(0.0f, pairSupport[static_cast<size_t>(qIdx * fullQueryLen + partner)]);
                    } else {
                        queryPair = 1.0f;
                    }
                    const float supportScale = std::max(0.10f, queryPair);
                    emit += cellPairWeight * supportScale * (targetPair - 0.10f);
                }
            }
            if (bonusMatrix != NULL && bonusMatrix->size() == static_cast<size_t>(n * m)) {
                emit += (*bonusMatrix)[static_cast<size_t>(i) * static_cast<size_t>(m) + static_cast<size_t>(j)];
            }
            score[static_cast<size_t>(i) * static_cast<size_t>(m) + static_cast<size_t>(j)] = emit;
        }
    }
    return score;
}

static LolPosterior runLolForwardBackward(const std::vector<float> &scoreMatrix,
                                          int n,
                                          int m,
                                          float gapOpen,
                                          float gapExtend,
                                          float temperature) {
    LolPosterior out;
    out.n = n;
    out.m = m;
    if (n <= 0 || m <= 0 || scoreMatrix.size() != static_cast<size_t>(n * m)) {
        return out;
    }
    const float temp = std::max(0.25f, temperature);
    const float go = -std::fabs(gapOpen) / temp;
    const float ge = -std::fabs(gapExtend) / temp;
    const size_t cells = static_cast<size_t>(n + 1) * static_cast<size_t>(m + 1);
    const auto idx = [m](int i, int j) -> size_t {
        return static_cast<size_t>(i) * static_cast<size_t>(m + 1) + static_cast<size_t>(j);
    };

    std::vector<float> FM(cells, NEG_INF_F), FI(cells, NEG_INF_F), FD(cells, NEG_INF_F);
    std::vector<float> BM(cells, NEG_INF_F), BI(cells, NEG_INF_F), BD(cells, NEG_INF_F);
    float logZ = NEG_INF_F;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const float s = scoreMatrix[static_cast<size_t>(i - 1) * static_cast<size_t>(m) + static_cast<size_t>(j - 1)] / temp;
            if (s <= NEG_INF_F / 2.0f) {
                continue;
            }
            FM[idx(i, j)] = s + logSumExp4(0.0f,
                                           FM[idx(i - 1, j - 1)],
                                           FI[idx(i - 1, j - 1)],
                                           FD[idx(i - 1, j - 1)]);
            FI[idx(i, j)] = logSumExp2(FM[idx(i - 1, j)] + go,
                                       FI[idx(i - 1, j)] + ge);
            FD[idx(i, j)] = logSumExp2(FM[idx(i, j - 1)] + go,
                                       FD[idx(i, j - 1)] + ge);
            logZ = logSumExp4(logZ, FM[idx(i, j)], FI[idx(i, j)], FD[idx(i, j)]);
        }
    }
    if (logZ <= NEG_INF_F / 2.0f) {
        return out;
    }

    for (int i = n; i >= 1; --i) {
        for (int j = m; j >= 1; --j) {
            float bm = 0.0f;
            float bi = 0.0f;
            float bd = 0.0f;
            if (i < n && j < m) {
                const float sNext = scoreMatrix[static_cast<size_t>(i) * static_cast<size_t>(m) + static_cast<size_t>(j)] / temp;
                if (sNext > NEG_INF_F / 2.0f) {
                    bm = logSumExp2(bm, sNext + BM[idx(i + 1, j + 1)]);
                    bi = logSumExp2(bi, sNext + BM[idx(i + 1, j + 1)]);
                    bd = logSumExp2(bd, sNext + BM[idx(i + 1, j + 1)]);
                }
            }
            if (i < n) {
                bm = logSumExp2(bm, BI[idx(i + 1, j)] + go);
                bi = logSumExp2(bi, BI[idx(i + 1, j)] + ge);
            }
            if (j < m) {
                bm = logSumExp2(bm, BD[idx(i, j + 1)] + go);
                bd = logSumExp2(bd, BD[idx(i, j + 1)] + ge);
            }
            BM[idx(i, j)] = bm;
            BI[idx(i, j)] = bi;
            BD[idx(i, j)] = bd;
        }
    }

    out.posterior.assign(static_cast<size_t>(n * m), 0.0f);
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const float pLog = FM[idx(i, j)] + BM[idx(i, j)] - logZ;
            float p = 0.0f;
            if (pLog > NEG_INF_F / 2.0f) {
                p = std::exp(std::min(0.0f, pLog));
            }
            out.posterior[static_cast<size_t>(i - 1) * static_cast<size_t>(m) + static_cast<size_t>(j - 1)] = p;
            out.maxP = std::max(out.maxP, p);
        }
    }
    return out;
}

static std::vector<LolAnchor> extractLolAnchors(const LolPosterior &posterior,
                                                int queryOffset,
                                                float minSupport,
                                                int topK,
                                                int minSep) {
    std::vector<LolAnchor> candidates;
    if (posterior.n <= 0 || posterior.m <= 0 || posterior.posterior.empty()) {
        return candidates;
    }
    const float floor = std::max(minSupport, posterior.maxP * 0.35f);
    for (int i = 0; i < posterior.n; ++i) {
        for (int j = 0; j < posterior.m; ++j) {
            const float p = posterior.posterior[static_cast<size_t>(i) * static_cast<size_t>(posterior.m) + static_cast<size_t>(j)];
            if (p < floor) {
                continue;
            }
            LolAnchor a;
            a.q = queryOffset + i;
            a.t = j;
            a.support = p;
            candidates.push_back(a);
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const LolAnchor &a, const LolAnchor &b) {
                  return a.support > b.support;
              });
    std::vector<LolAnchor> anchors;
    for (size_t i = 0; i < candidates.size() && static_cast<int>(anchors.size()) < topK; ++i) {
        bool keep = true;
        for (size_t k = 0; k < anchors.size(); ++k) {
            if (std::abs(candidates[i].q - anchors[k].q) < minSep
                || std::abs(candidates[i].t - anchors[k].t) < minSep) {
                keep = false;
                break;
            }
        }
        if (keep) {
            anchors.push_back(candidates[i]);
        }
    }
    std::sort(anchors.begin(), anchors.end(),
              [](const LolAnchor &a, const LolAnchor &b) {
                  if (a.q != b.q) {
                      return a.q < b.q;
                  }
                  return a.t < b.t;
              });
    return anchors;
}

static StemGuide buildLolGuideFromAnchors(int qLen,
                                          const std::vector<int> &partners,
                                          const std::vector<LolAnchor> &anchors) {
    StemGuide guide = makeNeutralGuide(qLen, partners);
    for (size_t i = 0; i < anchors.size(); ++i) {
        if (anchors[i].q < 0 || anchors[i].q >= qLen) {
            continue;
        }
        guide.anchorPosByQuery[static_cast<size_t>(anchors[i].q)] = static_cast<float>(anchors[i].t);
        guide.supportByQuery[static_cast<size_t>(anchors[i].q)] =
            std::max(guide.supportByQuery[static_cast<size_t>(anchors[i].q)], anchors[i].support);
    }
    return guide;
}

static std::vector<float> buildLolBonusMatrixFromAnchors(const std::string &querySeq,
                                                         const std::vector<std::array<float, 4>> &profile,
                                                         const std::vector<int> &partners,
                                                         const std::vector<float> &pairSupport,
                                                         const WindowInfo &window,
                                                         const std::vector<LolAnchor> &anchors,
                                                         int queryStart,
                                                         int queryEnd) {
    const int qLen = static_cast<int>(querySeq.size());
    const int qLo = std::max(0, std::min(queryStart, qLen - 1));
    const int qHi = std::max(qLo, std::min(queryEnd, qLen - 1));
    const int n = qHi - qLo + 1;
    const int m = static_cast<int>(window.seq.size());
    std::vector<float> bonus(static_cast<size_t>(std::max(0, n * m)), 0.0f);
    if (n <= 0 || m <= 0 || anchors.empty()) {
        return bonus;
    }
    const float relWeight = std::max(0.0f, parseEnvFloat("MMSEQS_LOLALIGN_REL_WEIGHT", 2.0f));
    const float pairWeight = std::max(0.0f, parseEnvFloat("MMSEQS_LOLALIGN_PAIR_PROFILE_WEIGHT", 0.35f));
    const float seqWeight = std::max(0.0f, parseEnvFloat("MMSEQS_LOLALIGN_SEQDIST_WEIGHT", 0.40f));
    const int relRadius = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_REL_RADIUS", 16));

    for (size_t ai = 0; ai < anchors.size(); ++ai) {
        const int qa = anchors[ai].q;
        const int ta = anchors[ai].t;
        if (qa < 0 || qa >= qLen || ta < 0 || ta >= m) {
            continue;
        }
        for (int q = qLo; q <= qHi; ++q) {
            const float qs = (pairSupport.size() == static_cast<size_t>(qLen * qLen))
                ? std::max(0.0f, pairSupport[static_cast<size_t>(q * qLen + qa)])
                : 0.0f;
            const int qDelta = q - qa;
            const float qLogDist = std::copysign(1.0f, static_cast<float>(qDelta))
                                 * std::log(1.0f + std::fabs(static_cast<float>(qDelta)));
            for (int t = 0; t < m; ++t) {
                const float pairStruct = windowPairSupportAt(window, t, ta);
                const float canonical = isCanonicalPair(window.seq[static_cast<size_t>(t)],
                                                        window.seq[static_cast<size_t>(ta)]) ? 1.0f : 0.0f;
                float pairProfile = 0.0f;
                if (q < static_cast<int>(profile.size()) && qa < static_cast<int>(profile.size())) {
                    const int b1 = encodeBase(window.seq[static_cast<size_t>(t)]);
                    const int b2 = encodeBase(window.seq[static_cast<size_t>(ta)]);
                    if (b1 >= 0 && b2 >= 0) {
                        pairProfile = 0.5f * (basePreferenceFromProfile(profile[static_cast<size_t>(q)], b1)
                                            + basePreferenceFromProfile(profile[static_cast<size_t>(qa)], b2));
                    }
                }
                const int tDelta = t - ta;
                const float tLogDist = std::copysign(1.0f, static_cast<float>(tDelta))
                                     * std::log(1.0f + std::fabs(static_cast<float>(tDelta)));
                const float distScore = std::max(0.0f, 1.0f - std::fabs(qLogDist - tLogDist) / static_cast<float>(relRadius));
                const float relScore = 0.40f * pairStruct
                                     + 0.20f * canonical
                                     + pairWeight * pairProfile
                                     + seqWeight * distScore
                                     - 0.30f;
                const float anchorWeight = std::max(0.10f, std::max(qs, anchors[ai].support));
                bonus[static_cast<size_t>(q - qLo) * static_cast<size_t>(m) + static_cast<size_t>(t)]
                    += relWeight * anchorWeight * relScore;
            }
        }
    }
    return bonus;
}

static AlignmentResult tracebackLolPosteriorMAC(const LolPosterior &posterior,
                                                int queryStart,
                                                float mactPenalty) {
    AlignmentResult out;
    const int n = posterior.n;
    const int m = posterior.m;
    if (n <= 0 || m <= 0 || posterior.posterior.empty()) {
        return out;
    }
    const auto pidx = [m](int i, int j) -> size_t {
        return static_cast<size_t>(i - 1) * static_cast<size_t>(m) + static_cast<size_t>(j - 1);
    };
    const auto idx = [m](int i, int j) -> size_t {
        return static_cast<size_t>(i) * static_cast<size_t>(m + 1) + static_cast<size_t>(j);
    };
    std::vector<float> S(static_cast<size_t>(n + 1) * static_cast<size_t>(m + 1), 0.0f);
    std::vector<uint8_t> bt(static_cast<size_t>(n + 1) * static_cast<size_t>(m + 1), 0);
    float bestScore = NEG_INF_F;
    int bestI = 0;
    int bestJ = 0;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const float p = posterior.posterior[pidx(i, j)];
            const float stop = p - mactPenalty;
            const float mm = S[idx(i - 1, j - 1)] + p - mactPenalty;
            const float ins = S[idx(i, j - 1)] - 0.5f * mactPenalty;
            const float del = S[idx(i - 1, j)] - 0.5f * mactPenalty;
            float val = stop;
            uint8_t state = 0;
            if (mm > val) {
                val = mm;
                state = 1;
            }
            if (ins > val) {
                val = ins;
                state = 2;
            }
            if (del > val) {
                val = del;
                state = 3;
            }
            S[idx(i, j)] = val;
            bt[idx(i, j)] = state;
            if (val > bestScore) {
                bestScore = val;
                bestI = i;
                bestJ = j;
            }
        }
    }
    if (bestScore <= NEG_INF_F / 2.0f) {
        return out;
    }
    std::string backtrace;
    int i = bestI;
    int j = bestJ;
    int matches = 0;
    while (i > 0 && j > 0) {
        const uint8_t state = bt[idx(i, j)];
        if (state == 0) {
            break;
        }
        if (state == 1) {
            backtrace.push_back('M');
            ++matches;
            --i;
            --j;
        } else if (state == 2) {
            backtrace.push_back('D');
            --j;
        } else {
            backtrace.push_back('I');
            --i;
        }
    }
    while (!backtrace.empty() && backtrace.back() != 'M') {
        backtrace.pop_back();
    }
    if (backtrace.empty()) {
        return out;
    }
    std::reverse(backtrace.begin(), backtrace.end());
    out.valid = true;
    out.score = bestScore;
    out.backtrace = backtrace;
    out.qStart = queryStart + i;
    out.dbStart = j;
    out.qEnd = out.qStart - 1;
    out.dbEnd = out.dbStart - 1;
    out.matches = matches;
    for (size_t k = 0; k < backtrace.size(); ++k) {
        switch (backtrace[k]) {
            case 'M':
                ++out.qEnd;
                ++out.dbEnd;
                break;
            case 'I':
                ++out.qEnd;
                break;
            case 'D':
                ++out.dbEnd;
                break;
            default:
                break;
        }
    }
    return out;
}

static StemGuide buildLolGuideFromAlignment(int qLen,
                                            const std::vector<int> &partners,
                                            const std::vector<float> &pairSupport,
                                            const WindowInfo &window,
                                            const AlignmentResult &aln) {
    StemGuide guide = makeNeutralGuide(qLen, partners);
    if (!aln.valid) {
        return guide;
    }
    const ProjectedAlignment proj = projectAlignment(qLen, window, aln);
    for (int q = 0; q < qLen; ++q) {
        const int t = proj.targetPosByQuery[static_cast<size_t>(q)];
        if (t < 0) {
            continue;
        }
        float support = 0.15f;
        const int partner = (q < static_cast<int>(partners.size())) ? partners[static_cast<size_t>(q)] : -1;
        if (partner >= 0 && partner < qLen) {
            const int tp = proj.targetPosByQuery[static_cast<size_t>(partner)];
            if (tp >= 0) {
                const float qs = (pairSupport.size() == static_cast<size_t>(qLen * qLen))
                    ? pairSupport[static_cast<size_t>(q * qLen + partner)]
                    : 0.0f;
                const float ps = windowPairSupportAt(window, t, tp);
                support = clampf(0.15f + 0.45f * ps + 0.40f * std::max(0.0f, qs), 0.0f, 1.5f);
            }
        }
        guide.anchorPosByQuery[static_cast<size_t>(q)] = static_cast<float>(t);
        guide.supportByQuery[static_cast<size_t>(q)] = support;
    }
    return guide;
}

static std::vector<float> buildLolBonusMatrix(const std::string &querySeq,
                                              const std::vector<std::array<float, 4>> &profile,
                                              const std::vector<int> &partners,
                                              const std::vector<float> &pairSupport,
                                              const WindowInfo &window,
                                              const AlignmentResult &aln,
                                              int queryStart,
                                              int queryEnd) {
    const int qLen = static_cast<int>(querySeq.size());
    const int qLo = std::max(0, std::min(queryStart, qLen - 1));
    const int qHi = std::max(qLo, std::min(queryEnd, qLen - 1));
    const int n = qHi - qLo + 1;
    const int m = static_cast<int>(window.seq.size());
    std::vector<float> bonus(static_cast<size_t>(std::max(0, n * m)), 0.0f);
    if (!aln.valid || n <= 0 || m <= 0) {
        return bonus;
    }

    const ProjectedAlignment proj = projectAlignment(qLen, window, aln);
    const float relWeight = std::max(0.0f, parseEnvFloat("MMSEQS_LOLALIGN_REL_WEIGHT", 2.0f));
    const float pairWeight = std::max(0.0f, parseEnvFloat("MMSEQS_LOLALIGN_PAIR_PROFILE_WEIGHT", 0.35f));
    const int relRadius = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_REL_RADIUS", 16));

    for (int q = qLo; q <= qHi; ++q) {
        const int partner = (q < static_cast<int>(partners.size())) ? partners[static_cast<size_t>(q)] : -1;
        if (partner < 0 || partner >= qLen) {
            continue;
        }
        const int mappedPartner = proj.targetPosByQuery[static_cast<size_t>(partner)];
        const int currentTarget = proj.targetPosByQuery[static_cast<size_t>(q)];
        if (mappedPartner < 0) {
            continue;
        }
        const float qs = (pairSupport.size() == static_cast<size_t>(qLen * qLen))
            ? std::max(0.0f, pairSupport[static_cast<size_t>(q * qLen + partner)])
            : 0.0f;
        const float queryPairWeight = std::max(0.10f, qs);
        for (int t = 0; t < m; ++t) {
            if (t == mappedPartner) {
                continue;
            }
            const float pairStruct = windowPairSupportAt(window, t, mappedPartner);
            const float canonical = isCanonicalPair(window.seq[static_cast<size_t>(t)],
                                                    window.seq[static_cast<size_t>(mappedPartner)]) ? 1.0f : 0.0f;
            float pairProfile = 0.0f;
            if (q < static_cast<int>(profile.size()) && partner < static_cast<int>(profile.size())) {
                const int b1 = encodeBase(window.seq[static_cast<size_t>(t)]);
                const int b2 = encodeBase(window.seq[static_cast<size_t>(mappedPartner)]);
                if (b1 >= 0 && b2 >= 0) {
                    pairProfile = 0.5f * (basePreferenceFromProfile(profile[static_cast<size_t>(q)], b1)
                                        + basePreferenceFromProfile(profile[static_cast<size_t>(partner)], b2));
                }
            }
            const float distScore = (currentTarget >= 0)
                ? std::max(0.0f, 1.0f - std::fabs(static_cast<float>(t - currentTarget)) / static_cast<float>(relRadius))
                : 0.5f;
            const float relScore = 0.45f * pairStruct
                                 + 0.25f * canonical
                                 + pairWeight * pairProfile
                                 + 0.20f * distScore
                                 - 0.35f;
            bonus[static_cast<size_t>(q - qLo) * static_cast<size_t>(m) + static_cast<size_t>(t)]
                += relWeight * queryPairWeight * relScore;
        }
    }
    return bonus;
}

static AlignmentResult runLolAlignIteration(const std::string &querySeq,
                                            const std::vector<std::array<float, 4>> &profile,
                                            const std::vector<StructVec> &queryStruct,
                                            const std::vector<int> &partners,
                                            const std::vector<float> &pairSupport,
                                            const std::vector<StemPattern> &patterns,
                                            const WindowInfo &window,
                                            int queryStart,
                                            int queryEnd,
                                            int diagCenter,
                                            int bandRadius) {
    const int qLen = static_cast<int>(querySeq.size());
    const float structWeight = clampf(parseEnvFloat("MMSEQS_LOLALIGN_STRUCT_WEIGHT", 0.35f), 0.0f, 4.0f);
    const float guideWeight = clampf(parseEnvFloat("MMSEQS_LOLALIGN_GUIDE_WEIGHT", 0.75f), 0.0f, 6.0f);
    const int guideRadius = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_GUIDE_RADIUS", 8));
    const int iterations = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_ITERATIONS", 3));
    const float temperature = std::max(0.25f, parseEnvFloat("MMSEQS_LOLALIGN_TEMPERATURE", 1.5f));
    const float gapOpen = std::max(0.5f, parseEnvFloat("MMSEQS_LOLALIGN_GAP_OPEN", 6.0f));
    const float gapExtend = std::max(0.1f, parseEnvFloat("MMSEQS_LOLALIGN_GAP_EXTEND", 1.0f));
    const float mactPenalty = std::max(0.05f, parseEnvFloat("MMSEQS_LOLALIGN_MACT", 0.2f));
    const int topAnchors = std::max(2, parseEnvInt("MMSEQS_LOLALIGN_TOP_ANCHORS", 12));
    const int anchorSep = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_ANCHOR_MIN_SEP", 3));
    const float minPosterior = clampf(parseEnvFloat("MMSEQS_LOLALIGN_MIN_POSTERIOR", 0.05f), 0.0f, 1.0f);
    const bool useSectionCYK = parseEnvInt("MMSEQS_LOLALIGN_SECTION_CYK", 1) != 0;
    const int topChains = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_SECTION_TOP_CHAINS", 3));
    const int stemBandExtra = std::max(8, parseEnvInt("MMSEQS_LOLALIGN_SECTION_BAND_EXTRA",
                                                      std::max(12, bandRadius >= 0 ? bandRadius : 24)));

    const int qLo = std::max(0, std::min(queryStart, qLen - 1));
    const int qHi = std::max(qLo, std::min(queryEnd, qLen - 1));
    const int n = qHi - qLo + 1;
    const int m = static_cast<int>(window.seq.size());
    if (n <= 0 || m <= 0) {
        return AlignmentResult();
    }

    StemGuide guide = makeNeutralGuide(qLen, partners);
    std::vector<float> bonus;
    AlignmentResult best;
    float bestObjective = NEG_INF_F;

    for (int iter = 0; iter < iterations; ++iter) {
        const std::vector<float> scoreMatrix = buildLolScoreMatrix(querySeq,
                                                                   profile,
                                                                   queryStruct,
                                                                   partners,
                                                                   pairSupport,
                                                                   guide,
                                                                   window,
                                                                   window.seq,
                                                                   window.structVec,
                                                                   queryStart,
                                                                   queryEnd,
                                                                   structWeight,
                                                                   (iter == 0) ? 0.0f : guideWeight,
                                                                   guideRadius,
                                                                   diagCenter,
                                                                   bandRadius,
                                                                   bonus.empty() ? NULL : &bonus);
        const LolPosterior posterior = runLolForwardBackward(scoreMatrix, n, m, gapOpen, gapExtend, temperature);
        if (posterior.posterior.empty() || posterior.maxP <= 0.0f) {
            break;
        }
        const std::vector<LolAnchor> anchors = extractLolAnchors(posterior, qLo, minPosterior, topAnchors, anchorSep);
        const AlignmentResult current = tracebackLolPosteriorMAC(posterior, qLo, mactPenalty);
        if (current.valid) {
            const float objective = current.score + posterior.maxP * 4.0f + static_cast<float>(anchors.size()) * 0.25f;
            if (objective > bestObjective) {
                bestObjective = objective;
                best = current;
            }
        }
        if (anchors.empty()) {
            break;
        }
        guide = buildLolGuideFromAnchors(qLen, partners, anchors);
        bonus = buildLolBonusMatrixFromAnchors(querySeq,
                                               profile,
                                               partners,
                                               pairSupport,
                                               window,
                                               anchors,
                                               queryStart,
                                               queryEnd);
    }
    if (useSectionCYK && best.valid && !patterns.empty()) {
        const ProjectedAlignment proj = projectAlignment(qLen, window, best);
        StemGuide sectionGuide = buildLolGuideFromAlignment(qLen,
                                                            partners,
                                                            pairSupport,
                                                            window,
                                                            best);
        const int seedDiag = (best.dbStart >= 0 && best.qStart >= 0)
            ? (best.dbStart - best.qStart)
            : diagCenter;
        int sectionBand = stemBandExtra;
        if (best.dbEnd >= best.dbStart && best.qEnd >= best.qStart) {
            sectionBand = std::max(sectionBand,
                                   std::abs((best.dbEnd - best.dbStart) - (best.qEnd - best.qStart)) + stemBandExtra);
        }
        std::vector<StemMatch> sectionMatches = findStemMatches(patterns, window, seedDiag, sectionBand);
        if (!sectionMatches.empty()) {
            for (size_t mi = 0; mi < sectionMatches.size(); ++mi) {
                const StemMatch &m = sectionMatches[mi];
                const int leftQ = std::max(0, std::min(m.leftQuery, qLen - 1));
                const int rightQ = std::max(leftQ, std::min(m.rightQuery, qLen - 1));
                float mapped = 0.0f;
                float wsum = 0.0f;
                for (int q = leftQ; q <= rightQ; ++q) {
                    if (proj.targetPosByQuery[static_cast<size_t>(q)] >= 0) {
                        mapped += static_cast<float>(proj.targetPosByQuery[static_cast<size_t>(q)]);
                        wsum += 1.0f;
                    } else if (sectionGuide.anchorPosByQuery[static_cast<size_t>(q)] >= 0.0f) {
                        mapped += sectionGuide.anchorPosByQuery[static_cast<size_t>(q)];
                        wsum += 0.75f;
                    }
                }
                if (wsum > 0.0f) {
                    const float center = mapped / wsum;
                    const float dist = std::fabs(0.5f * static_cast<float>(m.leftTarget + m.rightTarget) - center);
                    const float compat = std::max(0.0f, 1.0f - dist / static_cast<float>(std::max(8, sectionBand)));
                    sectionMatches[mi].score *= (0.60f + 0.40f * compat);
                    sectionMatches[mi].support = clampf(0.70f * sectionMatches[mi].support + 0.30f * compat, 0.0f, 1.5f);
                }
            }
            const std::vector<std::vector<int>> chains = buildTopChainsFromStemMatches(patterns, sectionMatches);
            float bestSectionObjective = bestObjective;
            for (size_t ci = 0; ci < chains.size() && static_cast<int>(ci) < topChains; ++ci) {
                AlignmentResult sectionAln = stitchLolSectionConstrainedCYK(querySeq,
                                                                            profile,
                                                                            queryStruct,
                                                                            partners,
                                                                            window,
                                                                            patterns,
                                                                            sectionMatches,
                                                                            chains[ci],
                                                                            structWeight);
                if (!sectionAln.valid) {
                    continue;
                }
                const StemGuide cykGuide = buildLolGuideFromAlignment(qLen,
                                                                      partners,
                                                                      pairSupport,
                                                                      window,
                                                                      sectionAln);
                const AlignmentSummary cykSummary = summarizeAlignment(querySeq,
                                                                      queryStruct,
                                                                      cykGuide,
                                                                      window,
                                                                      sectionAln,
                                                                      static_cast<unsigned int>(window.seq.size()));
                const float objective = sectionAln.score
                                      + 2.5f * cykSummary.meanStructSim
                                      + 1.5f * cykSummary.meanStemAnchor
                                      + 1.5f * cykSummary.seqId;
                if (objective > bestSectionObjective) {
                    bestSectionObjective = objective;
                    best = sectionAln;
                }
            }
        }
    }
    return best;
}

static std::vector<AnchorBlock> buildAnchorBlocksFromChain(const std::vector<StemPattern> &patterns,
                                                           const std::vector<StemMatch> &matches,
                                                           const std::vector<int> &chain) {
    std::vector<AnchorBlock> blocks;
    blocks.reserve(chain.size() * 2);
    for (size_t i = 0; i < chain.size(); ++i) {
        const StemMatch &match = matches[static_cast<size_t>(chain[i])];
        const StemPattern &pat = patterns[static_cast<size_t>(match.patternIdx)];
        if (pat.armLen <= 0) {
            continue;
        }
        AnchorBlock left;
        left.qStart = pat.leftStart;
        left.qEnd = pat.leftStart + pat.armLen - 1;
        left.tStart = match.leftTarget;
        left.tEnd = match.leftTarget + pat.armLen - 1;
        left.support = match.support;
        blocks.push_back(left);

        AnchorBlock right;
        right.qStart = pat.rightOuter - pat.armLen + 1;
        right.qEnd = pat.rightOuter;
        right.tStart = match.rightTarget - pat.armLen + 1;
        right.tEnd = match.rightTarget;
        right.support = match.support;
        blocks.push_back(right);
    }

    std::sort(blocks.begin(), blocks.end(),
              [](const AnchorBlock &a, const AnchorBlock &b) {
                  if (a.qStart != b.qStart) {
                      return a.qStart < b.qStart;
                  }
                  return a.tStart < b.tStart;
              });

    std::vector<AnchorBlock> filtered;
    filtered.reserve(blocks.size());
    int lastQ = -1;
    int lastT = -1;
    for (size_t i = 0; i < blocks.size(); ++i) {
        const AnchorBlock &block = blocks[i];
        if (block.qStart <= lastQ || block.tStart <= lastT) {
            continue;
        }
        filtered.push_back(block);
        lastQ = block.qEnd;
        lastT = block.tEnd;
    }
    return filtered;
}

struct StockholmRow {
    std::string id;
    std::string aln;
};

static std::string buildStockholmTextLocal(const std::string &id,
                                           const std::vector<StockholmRow> &rows,
                                           const std::string &ssCons) {
    std::ostringstream out;
    out << "# STOCKHOLM 1.0\n";
    out << "#=GF ID " << (id.empty() ? "cmlite_module" : id) << "\n";
    for (size_t i = 0; i < rows.size(); ++i) {
        out << rows[i].id << " " << rows[i].aln << "\n";
    }
    if (!rows.empty() && !rows[0].aln.empty()) {
        out << "#=GC RF " << std::string(rows[0].aln.size(), 'x') << "\n";
    }
    if (!ssCons.empty()) {
        out << "#=GC SS_cons " << ssCons << "\n";
    }
    out << "//\n";
    return out.str();
}

static std::string buildModuleSsCons(const StemPattern &pat) {
    const int len = pat.span + 1;
    if (len <= 0) {
        return std::string();
    }
    std::string ss(static_cast<size_t>(len), '.');
    for (int k = 0; k < pat.stemLen; ++k) {
        const int left = k;
        const int right = len - 1 - k;
        if (left >= 0 && right >= 0 && left < len && right < len && left < right) {
            ss[static_cast<size_t>(left)] = '(';
            ss[static_cast<size_t>(right)] = ')';
        }
    }
    return ss;
}

static std::string normalizeAlignedRnaSlice(const std::string &row,
                                            int start,
                                            int len) {
    if (start < 0 || len <= 0 || start + len > static_cast<int>(row.size())) {
        return std::string();
    }
    std::string slice = row.substr(static_cast<size_t>(start), static_cast<size_t>(len));
    for (size_t i = 0; i < slice.size(); ++i) {
        char c = normalizeBase(slice[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'U' && c != '-') {
            c = 'N';
        }
        slice[i] = c;
    }
    return slice;
}

static std::string buildModuleStockholm(const std::string &querySeq,
                                        const StemPattern &pat,
                                        const std::vector<ModelHit> &modelHits,
                                        const std::vector<CandidateHit> &candidates) {
    const int start = pat.leftStart;
    const int len = pat.span + 1;
    if (start < 0 || len <= 0 || start + len > static_cast<int>(querySeq.size())) {
        return std::string();
    }

    const float minCov = clampf(parseEnvFloat("MMSEQS_CMLITECM_MODULE_MIN_ROW_COV", 0.55f), 0.0f, 1.0f);
    std::vector<StockholmRow> rows;
    rows.reserve(modelHits.size() + 1);
    StockholmRow queryRow;
    queryRow.id = "query";
    queryRow.aln = normalizeAlignedRnaSlice(querySeq, start, len);
    rows.push_back(queryRow);

    std::set<std::string> seenRows;
    seenRows.insert(queryRow.aln);
    for (size_t i = 0; i < modelHits.size(); ++i) {
        const ModelHit &mh = modelHits[i];
        if (!mh.aln.valid) {
            continue;
        }
        const CandidateHit &cand = candidates[mh.candidateIdx];
        ProjectedAlignment proj = projectAlignment(static_cast<int>(querySeq.size()), cand.window, mh.aln);
        std::string row = normalizeAlignedRnaSlice(proj.row, start, len);
        if (row.empty()) {
            continue;
        }
        int nongap = 0;
        for (size_t j = 0; j < row.size(); ++j) {
            if (row[j] != '-') {
                ++nongap;
            }
        }
        const float cov = row.empty() ? 0.0f : static_cast<float>(nongap) / static_cast<float>(row.size());
        if (cov < minCov) {
            continue;
        }
        if (!seenRows.insert(row).second) {
            continue;
        }
        StockholmRow seedRow;
        seedRow.id = "t" + SSTR(candidates[mh.candidateIdx].raw.dbKey);
        seedRow.aln.swap(row);
        rows.push_back(seedRow);
    }
    if (rows.size() < 2) {
        return std::string();
    }
    return buildStockholmTextLocal("stem_" + SSTR(pat.leftStart) + "_" + SSTR(pat.rightOuter),
                                   rows,
                                   buildModuleSsCons(pat));
}

static bool writeTextFile(const std::string &path,
                          const std::string &text) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out.good()) {
        return false;
    }
    out.write(text.data(), static_cast<std::streamsize>(text.size()));
    return out.good();
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

static void removePathRecursively(const std::string &path) {
    if (path.empty()) {
        return;
    }
    std::string cmd = "rm -rf '" + path + "'";
    std::system(cmd.c_str());
}

static int localToGlobalTargetPos(const CandidateHit &cand,
                                  int localPos) {
    if (!cand.window.reverseStrand) {
        return cand.window.fullLo + localPos;
    }
    return cand.window.fullHi - localPos;
}

static int pickIntegerQuantile(std::vector<int> values,
                               float q) {
    if (values.empty()) {
        return 0;
    }
    q = clampf(q, 0.0f, 1.0f);
    std::sort(values.begin(), values.end());
    const size_t idx = static_cast<size_t>(std::floor(q * static_cast<float>(values.size() - 1)));
    return values[idx];
}

static ModuleBandPrior computeModuleBandPrior(const StemPattern &pat,
                                              const std::vector<ModelHit> &modelHits,
                                              const std::vector<CandidateHit> &candidates,
                                              int qLen) {
    ModuleBandPrior prior;
    prior.expectedRightStart = pat.rightOuter - pat.armLen + 1;
    if (pat.armLen <= 0 || prior.expectedRightStart < 0 || pat.leftStart < 0 || pat.rightOuter >= qLen) {
        return prior;
    }

    std::vector<int> leftDiagOffsets;
    std::vector<int> rightDiagOffsets;
    std::vector<int> spans;
    leftDiagOffsets.reserve(modelHits.size());
    rightDiagOffsets.reserve(modelHits.size());
    spans.reserve(modelHits.size());

    for (size_t i = 0; i < modelHits.size(); ++i) {
        const ModelHit &mh = modelHits[i];
        if (!mh.aln.valid) {
            continue;
        }
        const CandidateHit &cand = candidates[mh.candidateIdx];
        if (cand.window.seq.empty()) {
            continue;
        }
        ProjectedAlignment proj = projectAlignment(qLen, cand.window, mh.aln);
        const int leftTarget = proj.targetPosByQuery[static_cast<size_t>(pat.leftStart)];
        const int rightStartTarget = proj.targetPosByQuery[static_cast<size_t>(prior.expectedRightStart)];
        const int rightOuterTarget = proj.targetPosByQuery[static_cast<size_t>(pat.rightOuter)];
        if (leftTarget < 0 || rightStartTarget < 0 || rightOuterTarget < 0 || rightOuterTarget <= leftTarget) {
            continue;
        }
        const int modelSeedDiag = cand.window.localSeedStart - cand.oriented.qStart;
        leftDiagOffsets.push_back((leftTarget - pat.leftStart) - modelSeedDiag);
        rightDiagOffsets.push_back((rightStartTarget - prior.expectedRightStart) - modelSeedDiag);
        spans.push_back(rightOuterTarget - leftTarget);
    }

    const int minSupport = std::max(3, parseEnvInt("MMSEQS_CMLITECM_PRIOR_MIN_SUPPORT", 8));
    if (static_cast<int>(spans.size()) < minSupport) {
        return prior;
    }

    const float centralMass = clampf(parseEnvFloat("MMSEQS_CMLITECM_PRIOR_CENTRAL_MASS", 0.90f), 0.50f, 0.999f);
    const float qLo = 0.5f * (1.0f - centralMass);
    const float qHi = 1.0f - qLo;
    const int pad = std::max(1, parseEnvInt("MMSEQS_CMLITECM_PRIOR_PAD", 3));

    prior.leftDiagCenter = static_cast<float>(pickIntegerQuantile(leftDiagOffsets, 0.5f));
    prior.rightDiagCenter = static_cast<float>(pickIntegerQuantile(rightDiagOffsets, 0.5f));
    prior.spanCenter = static_cast<float>(pickIntegerQuantile(spans, 0.5f));
    prior.leftDiagLo = pickIntegerQuantile(leftDiagOffsets, qLo) - pad;
    prior.leftDiagHi = pickIntegerQuantile(leftDiagOffsets, qHi) + pad;
    prior.rightDiagLo = pickIntegerQuantile(rightDiagOffsets, qLo) - pad;
    prior.rightDiagHi = pickIntegerQuantile(rightDiagOffsets, qHi) + pad;
    prior.spanLo = std::max(pat.armLen, pickIntegerQuantile(spans, qLo) - pad);
    prior.spanHi = pickIntegerQuantile(spans, qHi) + pad;
    prior.valid = true;
    return prior;
}

static float scorePriorBandAxis(int value,
                                int lo,
                                int hi,
                                float center,
                                int softPad) {
    const int paddedLo = lo - softPad;
    const int paddedHi = hi + softPad;
    if (value < paddedLo || value > paddedHi) {
        return -1.0f;
    }
    if (value < lo) {
        return std::max(0.0f, 1.0f - static_cast<float>(lo - value) / static_cast<float>(std::max(1, softPad)));
    }
    if (value > hi) {
        return std::max(0.0f, 1.0f - static_cast<float>(value - hi) / static_cast<float>(std::max(1, softPad)));
    }
    const float halfWidth = 0.5f * static_cast<float>(std::max(1, hi - lo + 1));
    return std::max(0.35f, 1.0f - std::fabs(static_cast<float>(value) - center) / std::max(1.0f, halfWidth));
}

static float scoreModulePriorCompatibility(const StemPattern &pat,
                                           const StemMatch &match,
                                           int seedDiag,
                                           const ModuleBandPrior &prior) {
    if (!prior.valid || pat.armLen <= 0) {
        return 1.0f;
    }
    const int rightStart = match.rightTarget - pat.armLen + 1;
    const int leftDiagOffset = (match.leftTarget - pat.leftStart) - seedDiag;
    const int rightDiagOffset = (rightStart - prior.expectedRightStart) - seedDiag;
    const int span = match.rightTarget - match.leftTarget;
    const int softPad = std::max(1, parseEnvInt("MMSEQS_CMLITECM_PRIOR_SOFT_PAD", 4));
    const float leftScore = scorePriorBandAxis(leftDiagOffset, prior.leftDiagLo, prior.leftDiagHi, prior.leftDiagCenter, softPad);
    const float rightScore = scorePriorBandAxis(rightDiagOffset, prior.rightDiagLo, prior.rightDiagHi, prior.rightDiagCenter, softPad);
    const float spanScore = scorePriorBandAxis(span, prior.spanLo, prior.spanHi, prior.spanCenter, softPad);
    if (leftScore < 0.0f || rightScore < 0.0f || spanScore < 0.0f) {
        return -1.0f;
    }
    return 0.35f * leftScore + 0.35f * rightScore + 0.30f * spanScore;
}

static bool buildModuleSeedResults(const StemPattern &pat,
                                   const ModuleBandPrior &prior,
                                   const std::vector<CandidateHit> &candidates,
                                   int numThreads,
                                   std::vector<RnaMatcher::result_t> &seedResults) {
    const int maxSeedsPerTarget = std::max(1, parseEnvInt("MMSEQS_CMLITECM_MAX_SEEDS_PER_TARGET", 2));
    const int moduleLen = pat.span + 1;
    if (moduleLen <= 0) {
        return false;
    }

    std::vector<StemPattern> onePattern(1, pat);
#ifdef OPENMP
    const int ompThreads = std::max(1, numThreads);
#pragma omp parallel num_threads(ompThreads)
    {
        std::vector<RnaMatcher::result_t> localRows;
        localRows.reserve(256);
#pragma omp for schedule(dynamic, 32)
        for (long ci = 0; ci < static_cast<long>(candidates.size()); ++ci) {
            const CandidateHit &cand = candidates[static_cast<size_t>(ci)];
            if (cand.window.seq.empty()) {
                continue;
            }
            const int seedDiag = cand.window.localSeedStart - cand.oriented.qStart;
            const int seedBandRadius = estimateSeedBandRadius(cand, std::max(8, parseEnvInt("MMSEQS_CMLITE_BAND_EXTRA", 32)));
            std::vector<StemMatch> matches = findStemMatches(onePattern, cand.window, seedDiag, seedBandRadius);
            if (matches.empty()) {
                continue;
            }
            for (size_t mi = 0; mi < matches.size(); ++mi) {
                const float priorCompat = scoreModulePriorCompatibility(pat, matches[mi], seedDiag, prior);
                if (prior.valid && priorCompat < 0.0f) {
                    matches[mi].score = NEG_INF_F;
                    continue;
                }
                matches[mi].support = clampf(0.70f * matches[mi].support
                                           + 0.30f * std::max(0.0f, priorCompat),
                                           0.0f, 1.5f);
                matches[mi].score *= (prior.valid ? (0.35f + 0.65f * std::max(0.0f, priorCompat)) : 1.0f);
            }
            std::sort(matches.begin(), matches.end(),
                      [](const StemMatch &a, const StemMatch &b) {
                          if (a.score != b.score) {
                              return a.score > b.score;
                          }
                          return a.leftTarget < b.leftTarget;
                      });
            const int clusterRadius = std::max(2, parseEnvInt("MMSEQS_CMLITECM_SEED_CLUSTER_RADIUS", std::max(3, pat.armLen)));
            std::vector<StemMatch> selected;
            selected.reserve(static_cast<size_t>(maxSeedsPerTarget));
            for (size_t mi = 0; mi < matches.size() && static_cast<int>(selected.size()) < maxSeedsPerTarget; ++mi) {
                const StemMatch &match = matches[mi];
                if (!(match.score > NEG_INF_F * 0.5f)) {
                    continue;
                }
                const int matchRightStart = match.rightTarget - pat.armLen + 1;
                bool nearExisting = false;
                for (size_t si = 0; si < selected.size(); ++si) {
                    const int selectedRightStart = selected[si].rightTarget - pat.armLen + 1;
                    if (std::abs(match.leftTarget - selected[si].leftTarget) <= clusterRadius
                        && std::abs(matchRightStart - selectedRightStart) <= clusterRadius) {
                        nearExisting = true;
                        break;
                    }
                }
                if (!nearExisting) {
                    selected.push_back(match);
                }
            }
            for (size_t mi = 0; mi < selected.size(); ++mi) {
                const StemMatch &match = selected[mi];
                const int globalLeft = localToGlobalTargetPos(cand, match.leftTarget);
                const int globalRight = localToGlobalTargetPos(cand, match.rightTarget);
                const int score = std::max(1, static_cast<int>(std::lrint(match.score * 100.0f)));
                RnaMatcher::result_t res(cand.raw.dbKey,
                                         score,
                                         1.0f,
                                         1.0f,
                                         clampf(match.support, 0.0f, 1.0f),
                                         std::max(1e-12, static_cast<double>(cand.raw.eval)),
                                         static_cast<unsigned int>(moduleLen),
                                         0,
                                         moduleLen - 1,
                                         static_cast<unsigned int>(moduleLen),
                                         globalLeft,
                                         globalRight,
                                         cand.raw.dbLen,
                                         "");
                localRows.push_back(res);
            }
        }
#pragma omp critical
        {
            seedResults.insert(seedResults.end(), localRows.begin(), localRows.end());
        }
    }
#else
    (void)numThreads;
    for (size_t ci = 0; ci < candidates.size(); ++ci) {
        const CandidateHit &cand = candidates[ci];
        if (cand.window.seq.empty()) {
            continue;
        }
        const int seedDiag = cand.window.localSeedStart - cand.oriented.qStart;
        const int seedBandRadius = estimateSeedBandRadius(cand, std::max(8, parseEnvInt("MMSEQS_CMLITE_BAND_EXTRA", 32)));
        std::vector<StemMatch> matches = findStemMatches(onePattern, cand.window, seedDiag, seedBandRadius);
        if (matches.empty()) {
            continue;
        }
        for (size_t mi = 0; mi < matches.size(); ++mi) {
            const float priorCompat = scoreModulePriorCompatibility(pat, matches[mi], seedDiag, prior);
            if (prior.valid && priorCompat < 0.0f) {
                matches[mi].score = NEG_INF_F;
                continue;
            }
            matches[mi].support = clampf(0.70f * matches[mi].support
                                       + 0.30f * std::max(0.0f, priorCompat),
                                       0.0f, 1.5f);
            matches[mi].score *= (prior.valid ? (0.35f + 0.65f * std::max(0.0f, priorCompat)) : 1.0f);
        }
        std::sort(matches.begin(), matches.end(),
                  [](const StemMatch &a, const StemMatch &b) {
                      if (a.score != b.score) {
                          return a.score > b.score;
                      }
                      return a.leftTarget < b.leftTarget;
                  });
        const int clusterRadius = std::max(2, parseEnvInt("MMSEQS_CMLITECM_SEED_CLUSTER_RADIUS", std::max(3, pat.armLen)));
        std::vector<StemMatch> selected;
        selected.reserve(static_cast<size_t>(maxSeedsPerTarget));
        for (size_t mi = 0; mi < matches.size() && static_cast<int>(selected.size()) < maxSeedsPerTarget; ++mi) {
            const StemMatch &match = matches[mi];
            if (!(match.score > NEG_INF_F * 0.5f)) {
                continue;
            }
            const int matchRightStart = match.rightTarget - pat.armLen + 1;
            bool nearExisting = false;
            for (size_t si = 0; si < selected.size(); ++si) {
                const int selectedRightStart = selected[si].rightTarget - pat.armLen + 1;
                if (std::abs(match.leftTarget - selected[si].leftTarget) <= clusterRadius
                    && std::abs(matchRightStart - selectedRightStart) <= clusterRadius) {
                    nearExisting = true;
                    break;
                }
            }
            if (!nearExisting) {
                selected.push_back(match);
            }
        }
        for (size_t mi = 0; mi < selected.size(); ++mi) {
            const StemMatch &match = selected[mi];
            const int globalLeft = localToGlobalTargetPos(cand, match.leftTarget);
            const int globalRight = localToGlobalTargetPos(cand, match.rightTarget);
            const int score = std::max(1, static_cast<int>(std::lrint(match.score * 100.0f)));
            RnaMatcher::result_t res(cand.raw.dbKey,
                                     score,
                                     1.0f,
                                     1.0f,
                                     clampf(match.support, 0.0f, 1.0f),
                                     std::max(1e-12, static_cast<double>(cand.raw.eval)),
                                     static_cast<unsigned int>(moduleLen),
                                     0,
                                     moduleLen - 1,
                                     static_cast<unsigned int>(moduleLen),
                                     globalLeft,
                                     globalRight,
                                     cand.raw.dbLen,
                                     "");
            seedResults.push_back(res);
        }
    }
#endif
    std::sort(seedResults.begin(), seedResults.end(), RnaMatcher::compareHitsByPosAndStrand);
    return !seedResults.empty();
}

static bool writeResultDb(const std::string &dbPath,
                          const std::vector<RnaMatcher::result_t> &rows) {
    DBWriter writer(dbPath.c_str(), (dbPath + ".index").c_str(), 1, false, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();
    writer.writeStart(0);
    char buffer[4096];
    for (size_t i = 0; i < rows.size(); ++i) {
        const size_t len = RnaMatcher::resultToBuffer(buffer, rows[i], false, true, false);
        writer.writeAdd(buffer, len, 0);
    }
    writer.writeEnd(0, 0);
    writer.close(true);
    return true;
}

static bool parseModuleCmsearchOutput(const std::string &resultDbPath,
                                      const StemPattern &pat,
                                      const std::vector<CandidateHit> &candidates,
                                      const std::unordered_map<unsigned int, size_t> &candidateIndexByTarget,
                                      std::unordered_map<unsigned int, std::vector<StemMatch>> &matchesByTarget) {
    DBReader<unsigned int> reader(resultDbPath.c_str(), (resultDbPath + ".index").c_str(), 1,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::NOSORT);
    const size_t rid = reader.getId(0);
    if (rid == UINT_MAX) {
        reader.close();
        return false;
    }

    std::vector<RnaMatcher::result_t> rows;
    RnaMatcher::readAlignmentResults(rows, reader.getData(rid, 0), false);
    reader.close();

    const int maxHitsPerTarget = std::max(1, parseEnvInt("MMSEQS_CMLITECM_MAX_REFINED_HITS_PER_TARGET", 4));
    for (size_t i = 0; i < rows.size(); ++i) {
        const RnaMatcher::result_t &raw = rows[i];
        const std::unordered_map<unsigned int, size_t>::const_iterator it = candidateIndexByTarget.find(raw.dbKey);
        if (it == candidateIndexByTarget.end()) {
            continue;
        }
        const CandidateHit &cand = candidates[it->second];
        int localLeft = -1;
        int localRight = -1;
        if (!cand.window.reverseStrand) {
            localLeft = raw.dbStartPos - cand.window.fullLo;
            localRight = raw.dbEndPos - cand.window.fullLo;
        } else {
            localLeft = cand.window.fullHi - raw.dbStartPos;
            localRight = cand.window.fullHi - raw.dbEndPos;
        }
        if (localLeft > localRight) {
            std::swap(localLeft, localRight);
        }
        if (localLeft < 0 || localRight < localLeft || localRight >= static_cast<int>(cand.window.seq.size())) {
            continue;
        }

        StemMatch match;
        match.patternIdx = -1;
        match.leftQuery = pat.leftStart;
        match.rightQuery = pat.rightOuter;
        match.leftTarget = localLeft;
        match.rightTarget = localRight;
        const int expectedRightStart = pat.rightOuter - pat.armLen + 1;
        const int rightStart = match.rightTarget - pat.armLen + 1;
        const int diagLeft = match.leftTarget - pat.leftStart;
        const int diagRight = rightStart - expectedRightStart;
        match.diagCenter = static_cast<int>(std::lrint(0.5f * static_cast<float>(diagLeft + diagRight)));
        const float scoreNorm = clampf(static_cast<float>(std::max(0, raw.score))
                                       / static_cast<float>(std::max(4, pat.stemLen * 8)),
                                       0.0f,
                                       2.0f);
        match.support = clampf(0.45f * raw.seqId + 0.35f * raw.qcov + 0.20f * scoreNorm, 0.0f, 2.0f);
        match.score = pat.priority * (0.35f + 0.65f * match.support) * (0.50f + 0.50f * scoreNorm);
        std::vector<StemMatch> &bucket = matchesByTarget[raw.dbKey];
        if (static_cast<int>(bucket.size()) < maxHitsPerTarget) {
            bucket.push_back(match);
        } else {
            size_t worstIdx = 0;
            for (size_t j = 1; j < bucket.size(); ++j) {
                if (bucket[j].score < bucket[worstIdx].score) {
                    worstIdx = j;
                }
            }
            if (match.score > bucket[worstIdx].score) {
                bucket[worstIdx] = match;
            }
        }
    }
    return !matchesByTarget.empty();
}

static bool refineStemMatchesWithRealCm(const std::string &execPath,
#ifdef HAVE_INFERNAL_BRIDGE
                                        InfernalBridge::WorkerPool *workerPool,
#endif
                                        int moduleThreads,
                                        const std::string &querySeq,
                                        const StemPattern &pat,
                                        const std::vector<ModelHit> &modelHits,
                                        const std::vector<CandidateHit> &candidates,
                                        const std::unordered_map<unsigned int, size_t> &candidateIndexByTarget,
                                        const LocalParameters &par,
                                        std::unordered_map<unsigned int, std::vector<StemMatch>> &matchesByTarget) {
#ifndef HAVE_INFERNAL_BRIDGE
    (void)execPath;
    (void)querySeq;
    (void)pat;
    (void)modelHits;
    (void)candidates;
    (void)candidateIndexByTarget;
    (void)par;
    (void)matchesByTarget;
    return false;
#else
    if (workerPool == NULL || execPath.empty()) {
        return false;
    }

    const std::string stoText = buildModuleStockholm(querySeq, pat, modelHits, candidates);
    if (stoText.empty()) {
        return false;
    }

    std::vector<RnaMatcher::result_t> seedRows;
    const ModuleBandPrior prior = computeModuleBandPrior(pat, modelHits, candidates, static_cast<int>(querySeq.size()));
    if (!buildModuleSeedResults(pat, prior, candidates, std::max(1, moduleThreads), seedRows)) {
        return false;
    }

    std::string cmText;
    std::string infernalErr;
    if (!InfernalBridge::buildCmFromStockholmText(workerPool, stoText, cmText, infernalErr, false)) {
        Debug(Debug::WARNING) << "cmlitecm: module cmbuild failed for stem "
                              << pat.leftStart << "-" << pat.rightOuter
                              << ": " << infernalErr << "\n";
        return false;
    }

    const std::string tempDir = createTempDir("riboseek_cmlitecm");
    if (tempDir.empty()) {
        return false;
    }
    const std::string cmPath = tempDir + "/module.cm";
    const std::string resultIn = tempDir + "/module_result";
    const std::string resultOut = tempDir + "/module_aln";
    bool ok = writeTextFile(cmPath, cmText) && writeResultDb(resultIn, seedRows);
    if (!ok) {
        removePathRecursively(tempDir);
        return false;
    }

    std::vector<std::string> args;
    args.push_back(execPath);
    args.push_back("cmsearch");
    args.push_back(cmPath);
    args.push_back(par.db2);
    args.push_back(resultIn);
    args.push_back(resultOut);
    args.push_back("-e");
    args.push_back("inf");
    args.push_back("--cm-region");
    args.push_back(SSTR((par.cmRegionFlanking > 0.0f) ? par.cmRegionFlanking : 1.5f));
    args.push_back("--threads");
    args.push_back(SSTR(std::max(1, moduleThreads)));
    args.push_back("-v");
    args.push_back("0");

    std::string runErr;
    std::vector<std::pair<std::string, std::string>> env;
    env.push_back(std::make_pair("MMSEQS_CMSCAN_ALLOW_MULTIHIT", "1"));
    if (!runExternalCommand(args, env, runErr)) {
        Debug(Debug::WARNING) << "cmlitecm: module cmsearch failed for stem "
                              << pat.leftStart << "-" << pat.rightOuter
                              << ": " << runErr << "\n";
        removePathRecursively(tempDir);
        return false;
    }

    ok = parseModuleCmsearchOutput(resultOut, pat, candidates, candidateIndexByTarget, matchesByTarget);
    removePathRecursively(tempDir);
    return ok;
#endif
}

static void appendTargetAdvance(std::string &bt,
                                int &targetPos,
                                int targetGoal) {
    while (targetPos < targetGoal) {
        bt.push_back('D');
        ++targetPos;
    }
}

static AlignmentResult stitchCmLiteAlignment(const std::string &querySeq,
                                             const std::vector<std::array<float, 4>> &profile,
                                             const std::vector<StructVec> &queryStruct,
                                             const std::vector<int> &partners,
                                             const std::vector<AnchorBlock> &blocks,
                                             const std::string &targetSeq,
                                             const std::vector<StructVec> &targetStruct,
                                             const std::vector<int> &targetPartners,
                                             const std::vector<float> &targetPairMatrix,
                                             const WindowInfo *targetWindow,
                                             int queryStart,
                                             int queryEnd,
                                             float structScoreWeight) {
    AlignmentResult out;
    const int qLen = static_cast<int>(querySeq.size());
    const int tLen = static_cast<int>(targetSeq.size());
    if (qLen <= 0 || tLen <= 0) {
        return out;
    }

    const int qLo = std::max(0, std::min(queryStart, qLen - 1));
    const int qHi = std::max(qLo, std::min(queryEnd, qLen - 1));
    std::string bt;
    int cursorQ = qLo;
    int cursorT = 0;

    for (size_t bi = 0; bi <= blocks.size(); ++bi) {
        const bool hasBlock = (bi < blocks.size());
        const int nextQ = hasBlock ? blocks[bi].qStart : (qHi + 1);
        const int nextT = hasBlock ? blocks[bi].tStart : tLen;

        if (nextQ > cursorQ) {
            const int segQStart = cursorQ;
            const int segQEnd = nextQ - 1;
            const int segTStart = cursorT;
            const int segTEnd = std::max(segTStart - 1, nextT - 1);
            if (!appendPairAwareSegment(querySeq,
                                        profile,
                                        queryStruct,
                                        partners,
                                        targetSeq,
                                        targetStruct,
                                        targetPartners,
                                        targetPairMatrix,
                                        targetWindow,
                                        segQStart,
                                        segQEnd,
                                        segTStart,
                                        segTEnd,
                                        structScoreWeight,
                                        0,
                                        bt)) {
                return AlignmentResult();
            }
            cursorQ = segQEnd + 1;
            cursorT = nextT;
        }

        if (!hasBlock) {
            break;
        }

        appendTargetAdvance(bt, cursorT, blocks[bi].tStart);
        for (int q = blocks[bi].qStart, t = blocks[bi].tStart;
             q <= blocks[bi].qEnd && t <= blocks[bi].tEnd;
             ++q, ++t) {
            bt.push_back('M');
        }
        cursorQ = blocks[bi].qEnd + 1;
        cursorT = blocks[bi].tEnd + 1;
    }

    while (cursorQ <= qHi) {
        bt.push_back('I');
        ++cursorQ;
    }

    if (bt.empty()) {
        return out;
    }

    int firstConsumedT = -1;
    int firstMatchedT = -1;
    int lastMatchedT = -1;
    int scanQ = qLo;
    int scanT = 0;
    for (size_t i = 0; i < bt.size(); ++i) {
        const char op = bt[i];
        if (op == 'M') {
            if (firstConsumedT < 0) {
                firstConsumedT = scanT;
            }
            if (firstMatchedT < 0) {
                firstMatchedT = scanT;
            }
            lastMatchedT = scanT;
            ++scanQ;
            ++scanT;
        } else if (op == 'I') {
            ++scanQ;
        } else if (op == 'D') {
            if (firstConsumedT < 0) {
                firstConsumedT = scanT;
            }
            ++scanT;
        }
    }
    if (firstConsumedT < 0 || firstMatchedT < 0 || lastMatchedT < firstMatchedT) {
        return out;
    }

    out.valid = true;
    out.qStart = qLo;
    out.qEnd = qHi;
    out.dbStart = firstConsumedT;
    out.dbEnd = lastMatchedT;
    out.backtrace.swap(bt);
    out.score = 0.0f;
    out.matches = 0;
    int qPos = out.qStart;
    int tPos = out.dbStart;
    bool started = false;
    for (size_t i = 0; i < out.backtrace.size(); ++i) {
        const char op = out.backtrace[i];
        if (op == 'M') {
            if (!started) {
                started = true;
            }
            const int b = (tPos >= 0 && tPos < tLen) ? encodeBase(targetSeq[static_cast<size_t>(tPos)]) : -1;
            out.score += (b >= 0 && qPos >= 0 && qPos < qLen) ? profile[static_cast<size_t>(qPos)][static_cast<size_t>(b)] : -3.0f;
            if (qPos >= 0 && qPos < qLen && tPos >= 0 && tPos < tLen
                && normalizeBase(querySeq[static_cast<size_t>(qPos)]) == normalizeBase(targetSeq[static_cast<size_t>(tPos)])) {
                ++out.matches;
            }
            ++qPos;
            ++tPos;
        } else if (op == 'I') {
            out.score -= 6.0f;
            ++qPos;
        } else if (op == 'D') {
            out.score -= 6.0f;
            ++tPos;
        }
    }
    return out;
}

static AlignmentResult stitchPairAwareLocalWindow(const std::string &querySeq,
                                                  const std::vector<std::array<float, 4>> &profile,
                                                  const std::vector<StructVec> &queryStruct,
                                                  const std::vector<int> &partners,
                                                  const std::string &targetSeq,
                                                  const std::vector<StructVec> &targetStruct,
                                                  const std::vector<int> &targetPartners,
                                                  const std::vector<float> &targetPairMatrix,
                                                  const WindowInfo *targetWindow,
                                                  int qLo,
                                                  int qHi,
                                                  int tLo,
                                                  int tHi,
                                                  float structScoreWeight) {
    AlignmentResult out;
    const int qLen = static_cast<int>(querySeq.size());
    const int tLen = static_cast<int>(targetSeq.size());
    if (qLen <= 0 || tLen <= 0 || qLo < 0 || qHi < qLo || tLo < 0 || tHi < tLo || tHi >= tLen) {
        return out;
    }

    SegmentStitchResult seg = buildPairAwareSegmentBest(querySeq,
                                                        profile,
                                                        queryStruct,
                                                        partners,
                                                        targetSeq,
                                                        targetStruct,
                                                        targetPartners,
                                                        targetPairMatrix,
                                                        targetWindow,
                                                        qLo,
                                                        qHi,
                                                        tLo,
                                                        tHi,
                                                        structScoreWeight,
                                                        0);
    if (!seg.valid || seg.bt.empty()) {
        return out;
    }

    out.valid = true;
    out.qStart = qLo;
    out.qEnd = qHi;
    out.dbStart = tLo;
    out.dbEnd = tHi;
    out.backtrace.swap(seg.bt);
    out.score = seg.score;
    out.matches = 0;
    int qPos = qLo;
    int tPos = tLo;
    int firstMatchedT = -1;
    int lastMatchedT = -1;
    for (size_t i = 0; i < out.backtrace.size(); ++i) {
        const char op = out.backtrace[i];
        if (op == 'M') {
            if (firstMatchedT < 0) {
                firstMatchedT = tPos;
            }
            lastMatchedT = tPos;
            if (qPos >= 0 && qPos < qLen && tPos >= 0 && tPos < tLen
                && normalizeBase(querySeq[static_cast<size_t>(qPos)]) == normalizeBase(targetSeq[static_cast<size_t>(tPos)])) {
                ++out.matches;
            }
            ++qPos;
            ++tPos;
        } else if (op == 'I') {
            ++qPos;
        } else if (op == 'D') {
            ++tPos;
        }
    }
    if (firstMatchedT >= 0 && lastMatchedT >= firstMatchedT) {
        out.dbStart = firstMatchedT;
        out.dbEnd = lastMatchedT;
    }
    return out;
}

static SegmentStitchResult buildLocalLocalSectionSegment(const std::string &querySeq,
                                                         const std::vector<std::array<float, 4>> &profile,
                                                         const std::vector<StructVec> &queryStruct,
                                                         const std::vector<int> &partners,
                                                         const std::string &targetSeq,
                                                         const std::vector<StructVec> &targetStruct,
                                                         int qLo,
                                                         int qHi,
                                                         int tLo,
                                                         int tHi,
                                                         float structScoreWeight) {
    SegmentStitchResult out;
    if (qLo > qHi) {
        out.valid = true;
        out.score = -0.50f * static_cast<float>(std::max(0, tHi - tLo + 1));
        out.bt.assign(static_cast<size_t>(std::max(0, tHi - tLo + 1)), 'D');
        return out;
    }
    if (tLo > tHi) {
        out.valid = true;
        out.score = -0.50f * static_cast<float>(std::max(0, qHi - qLo + 1));
        out.bt.assign(static_cast<size_t>(std::max(0, qHi - qLo + 1)), 'I');
        return out;
    }

    const std::string targetSub = targetSeq.substr(static_cast<size_t>(tLo),
                                                   static_cast<size_t>(tHi - tLo + 1));
    std::vector<StructVec> targetStructSub;
    if (!targetStruct.empty()) {
        targetStructSub.assign(targetStruct.begin() + tLo,
                               targetStruct.begin() + tHi + 1);
    }
    StemGuide neutralGuide = makeNeutralGuide(static_cast<int>(querySeq.size()), partners);
    AlignmentResult aln = alignProfileLocalLocal(querySeq,
                                                 profile,
                                                 queryStruct,
                                                 neutralGuide,
                                                 targetSub,
                                                 targetStructSub,
                                                 qLo,
                                                 qHi,
                                                 structScoreWeight,
                                                 0.0f,
                                                 1,
                                                 0,
                                                 -1,
                                                 NULL);
    if (!aln.valid) {
        return out;
    }

    out.valid = true;
    out.score = aln.score
              - 0.40f * static_cast<float>(std::max(0, aln.qStart - qLo))
              - 0.40f * static_cast<float>(std::max(0, qHi - aln.qEnd))
              - 0.25f * static_cast<float>(std::max(0, aln.dbStart))
              - 0.25f * static_cast<float>(std::max(0, static_cast<int>(targetSub.size()) - 1 - aln.dbEnd));
    out.bt.reserve(static_cast<size_t>(qHi - qLo + 1 + tHi - tLo + 1));
    for (int q = qLo; q < aln.qStart; ++q) {
        out.bt.push_back('I');
    }
    for (int d = 0; d < aln.dbStart; ++d) {
        out.bt.push_back('D');
    }
    out.bt += aln.backtrace;
    for (int q = aln.qEnd + 1; q <= qHi; ++q) {
        out.bt.push_back('I');
    }
    for (int d = aln.dbEnd + 1; d < static_cast<int>(targetSub.size()); ++d) {
        out.bt.push_back('D');
    }
    return out;
}

static std::vector<LolSectionBlock> buildLolSectionBlocks(const std::vector<StemPattern> &patterns,
                                                          const std::vector<StemMatch> &matches,
                                                          const std::vector<int> &chain,
                                                          int qLen) {
    std::vector<LolSectionBlock> blocks;
    if (qLen <= 0) {
        return blocks;
    }
    int qCursor = 0;
    int tCursor = 0;
    for (size_t ci = 0; ci < chain.size(); ++ci) {
        const StemMatch &match = matches[static_cast<size_t>(chain[ci])];
        if (match.patternIdx < 0 || match.patternIdx >= static_cast<int>(patterns.size())) {
            continue;
        }
        const StemPattern &pat = patterns[static_cast<size_t>(match.patternIdx)];
        if (pat.leftStart > qCursor) {
            LolSectionBlock loop;
            loop.stem = false;
            loop.qLo = qCursor;
            loop.qHi = pat.leftStart - 1;
            loop.tLo = tCursor;
            loop.tHi = std::max(tCursor - 1, match.leftTarget - 1);
            loop.support = match.support * 0.5f;
            blocks.push_back(loop);
        }
        LolSectionBlock stem;
        stem.stem = true;
        stem.qLo = pat.leftStart;
        stem.qHi = pat.rightOuter;
        stem.tLo = match.leftTarget;
        stem.tHi = match.rightTarget;
        stem.support = match.support;
        blocks.push_back(stem);
        qCursor = pat.rightOuter + 1;
        tCursor = match.rightTarget + 1;
    }
    if (qCursor < qLen) {
        LolSectionBlock tail;
        tail.stem = false;
        tail.qLo = qCursor;
        tail.qHi = qLen - 1;
        tail.tLo = tCursor;
        tail.tHi = std::max(tCursor - 1, tCursor + (tail.qHi - tail.qLo));
        tail.support = 0.25f;
        blocks.push_back(tail);
    }
    return blocks;
}

static AlignmentResult stitchLolSectionConstrainedCYK(const std::string &querySeq,
                                                      const std::vector<std::array<float, 4>> &profile,
                                                      const std::vector<StructVec> &queryStruct,
                                                      const std::vector<int> &partners,
                                                      const WindowInfo &window,
                                                      const std::vector<StemPattern> &patterns,
                                                      const std::vector<StemMatch> &matches,
                                                      const std::vector<int> &chain,
                                                      float structScoreWeight) {
    AlignmentResult out;
    if (chain.empty()) {
        return out;
    }
    std::vector<LolSectionBlock> blocks = buildLolSectionBlocks(patterns,
                                                                matches,
                                                                chain,
                                                                static_cast<int>(querySeq.size()));
    if (blocks.empty()) {
        return out;
    }

    const int targetLen = static_cast<int>(window.seq.size());
    const int sectionSlack = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_SECTION_SLACK", 6));
    SegmentStitchResult total;
    total.valid = true;
    total.score = 0.0f;
    int globalQStart = -1;
    int globalQEnd = -1;
    int globalTStart = -1;
    int globalTEnd = -1;

    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        LolSectionBlock block = blocks[bi];
        block.qLo = std::max(0, block.qLo);
        block.qHi = std::min(static_cast<int>(querySeq.size()) - 1, block.qHi);
        block.tLo = std::max(0, block.tLo - (block.stem ? sectionSlack : 0));
        if (block.stem) {
            block.tHi = std::min(targetLen - 1, block.tHi + sectionSlack);
        } else {
            const int qSpan = std::max(0, block.qHi - block.qLo + 1);
            block.tHi = std::min(targetLen - 1, std::max(block.tHi, block.tLo + qSpan + sectionSlack));
        }
        if (block.qLo > block.qHi || block.tLo > block.tHi || block.tLo >= targetLen) {
            return AlignmentResult();
        }

        SegmentStitchResult seg;
        if (block.stem) {
            seg = buildPairAwareSegmentBest(querySeq,
                                           profile,
                                           queryStruct,
                                           partners,
                                           window.seq,
                                           window.structVec,
                                           window.partners,
                                           window.pairMatrix,
                                           &window,
                                           block.qLo,
                                           block.qHi,
                                           block.tLo,
                                           block.tHi,
                                           structScoreWeight,
                                           0);
        } else {
            seg = buildLocalLocalSectionSegment(querySeq,
                                               profile,
                                               queryStruct,
                                               partners,
                                               window.seq,
                                               window.structVec,
                                               block.qLo,
                                               block.qHi,
                                               block.tLo,
                                               block.tHi,
                                               structScoreWeight);
        }
        if (!seg.valid || seg.bt.empty()) {
            return AlignmentResult();
        }
        total.score += seg.score + (block.stem ? (1.5f * block.support) : (0.25f * block.support));
        total.bt += seg.bt;
        if (globalQStart < 0) {
            globalQStart = block.qLo;
            globalTStart = block.tLo;
        }
        globalQEnd = block.qHi;
        globalTEnd = block.tHi;
    }

    if (total.bt.empty()) {
        return out;
    }
    out.valid = true;
    out.score = total.score;
    out.qStart = std::max(0, globalQStart);
    out.qEnd = std::max(out.qStart, globalQEnd);
    out.dbStart = std::max(0, globalTStart);
    out.dbEnd = std::max(out.dbStart, globalTEnd);
    out.backtrace.swap(total.bt);
    out.matches = 0;
    int qPos = out.qStart;
    int tPos = out.dbStart;
    int firstMatchedT = -1;
    int lastMatchedT = -1;
    for (size_t i = 0; i < out.backtrace.size(); ++i) {
        const char op = out.backtrace[i];
        if (op == 'M') {
            if (firstMatchedT < 0) {
                firstMatchedT = tPos;
            }
            lastMatchedT = tPos;
            if (qPos >= 0 && qPos < static_cast<int>(querySeq.size())
                && tPos >= 0 && tPos < targetLen
                && normalizeBase(querySeq[static_cast<size_t>(qPos)]) == normalizeBase(window.seq[static_cast<size_t>(tPos)])) {
                ++out.matches;
            }
            ++qPos;
            ++tPos;
        } else if (op == 'I') {
            ++qPos;
        } else if (op == 'D') {
            ++tPos;
        }
    }
    if (firstMatchedT >= 0) {
        out.dbStart = firstMatchedT;
        out.dbEnd = std::max(firstMatchedT, lastMatchedT);
    }
    return out;
}

} // namespace

static int runCmliteImpl(int argc,
                         const char **argv,
                         const Command &command,
                         bool useRealModuleCm) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    par.evalThr = DBL_MAX;
    if (MMseqsMPI::isMaster() == false) {
        return EXIT_SUCCESS;
    }

    const float alpha = clampf(parseEnvFloat("MMSEQS_CMLITE_STRUCT_ALPHA", 0.35f), 0.0f, 2.0f);
    const float stemWeight = clampf(parseEnvFloat("MMSEQS_CMLITE_STEM_WEIGHT", 2.0f), 0.0f, 6.0f);
    const int stemRadius = std::max(1, parseEnvInt("MMSEQS_CMLITE_STEM_RADIUS", 6));
    const int bandExtra = std::max(8, parseEnvInt("MMSEQS_CMLITE_BAND_EXTRA", 32));
    const bool useFullQueryRefine = parseEnvInt("MMSEQS_CMLITE_FULLREFINE", 1) != 0;
    const float fullRefineStemWeight = clampf(parseEnvFloat("MMSEQS_CMLITE_FULLREFINE_STEM_WEIGHT", 1.25f), 0.0f, 6.0f);
    const int fullRefineBandSlack = std::max(4, parseEnvInt("MMSEQS_CMLITE_FULLREFINE_BAND_SLACK", 24));
    const float fullRefineCoverageWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_FULLREFINE_COVERAGE_WEIGHT", 4.0f));
    const int fullRefineWindowSlack = std::max(0, parseEnvInt("MMSEQS_CMLITE_FULLREFINE_WINDOW_SLACK", 6));
    const bool acceptAllRefinements = parseEnvInt("MMSEQS_CMLITE_ACCEPT_ALL", 0) != 0;
    const bool useRowDedup = parseEnvInt("MMSEQS_CMLITE_DEDUP_ROWS", 1) != 0;
    const bool useDinuc = par.emsearchDinucleotide || (parseEnvInt("MMSEQS_CMLITE_DINUC", 0) != 0);
    const double msaEvalThr = par.cmliteMsaEvalThr;
    const bool uncappedModelHits = std::isfinite(msaEvalThr);
    const int maxModelHits = std::max(4, parseEnvInt("MMSEQS_CMLITE_MAX_MODEL_HITS", 64));
    const int querySupportMaxLen = std::max(0, parseEnvInt("MMSEQS_CMLITE_QUERY_SUPPORT_MAX_LEN", 0));
    const float flankFrac = (par.cmRegionFlanking > 0.0f)
        ? par.cmRegionFlanking
        : parseEnvFloat("MMSEQS_CMLITE_FLANK_FRAC", 0.5f);
#ifdef HAVE_INFERNAL_BRIDGE
    InfernalBridge::WorkerPool *workerPool = NULL;
    std::string execPath;
    if (useRealModuleCm) {
        workerPool = InfernalBridge::startWorkerPool(std::max(1, par.threads));
        execPath = getSelfExecutablePath();
        if (workerPool == NULL || execPath.empty()) {
            Debug(Debug::ERROR) << "cmlitecm: failed to initialize Infernal bridge or executable path\n";
            if (workerPool != NULL) {
                InfernalBridge::stopWorkerPool(workerPool);
            }
            return EXIT_FAILURE;
        }
    }
#else
    if (useRealModuleCm) {
        Debug(Debug::ERROR) << "cmlitecm requires Infernal bridge support\n";
        return EXIT_FAILURE;
    }
#endif
    StructureSettings structureSettings;
    structureSettings.backend = std::getenv("MMSEQS_CMLITE_STRUCTURE_BACKEND")
        ? std::getenv("MMSEQS_CMLITE_STRUCTURE_BACKEND")
        : "auto";
    structureSettings.minLoop = std::max(0, parseEnvInt("MMSEQS_CMLITE_MIN_LOOP", 3));
    structureSettings.linearfoldMinLength = std::max(0, parseEnvInt("MMSEQS_CMLITE_LINEARFOLD_MIN_LENGTH", 128));
    structureSettings.linearfoldBeamSize = std::max(1, parseEnvInt("MMSEQS_CMLITE_LINEARFOLD_BEAM", 100));
    structureSettings.rnaPlfoldWindow = std::max(0, parseEnvInt("MMSEQS_CMLITE_RNAPLFOLD_WINDOW", 0));
    structureSettings.rnaPlfoldSpan = std::max(0, parseEnvInt("MMSEQS_CMLITE_RNAPLFOLD_SPAN", 0));
    structureSettings.sparsePairMatrix = parseEnvInt("MMSEQS_CMLITE_SPARSE_PAIR_MATRIX", 0) != 0;

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0f);

    Debug(Debug::INFO) << (useRealModuleCm ? "Running CM-lite realignment with real module CMs\n"
                                           : "Running CM-lite realignment\n");
    Debug(Debug::INFO) << "  refinement mode: ignoring alignment E-value threshold, candidate set comes from input result DB\n";
    std::ostringstream msaEvalStream;
    msaEvalStream << std::scientific << msaEvalThr;
    Debug(Debug::INFO) << "  model hits: "
                       << (uncappedModelHits
                               ? "all deduplicated hits with E-value <= threshold"
                               : ("top " + SSTR(maxModelHits) + " deduplicated hits (threshold is inf)"))
                       << ", stem weight: " << stemWeight
                       << ", stem radius: " << stemRadius
                       << ", query support max length: " << querySupportMaxLen
                       << ", flank frac: " << flankFrac
                       << ", msa e-value threshold: " << msaEvalStream.str()
                       << ", structure backend: " << structureSettings.backend
                       << ", linearfold min length: " << structureSettings.linearfoldMinLength
                       << ", linearfold beam: " << structureSettings.linearfoldBeamSize
                       << ", rnaplfold window: " << structureSettings.rnaPlfoldWindow
                       << ", rnaplfold span: " << structureSettings.rnaPlfoldSpan
                       << ", keep pair matrix: " << (structureSettings.keepPairMatrix ? "yes" : "no")
                       << ", sparse pair matrix: " << (structureSettings.sparsePairMatrix ? "yes" : "no")
                       << ", dot-bracket only: " << (structureSettings.dotBracketOnly ? "yes" : "no")
                       << "\n";

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qDbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> tDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tDbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> resReader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                     DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const unsigned int queryExt = DBReader<unsigned int>::getExtendedDbtype(qDbr.getDbtype());
    const unsigned int targetExt = DBReader<unsigned int>::getExtendedDbtype(tDbr.getDbtype());
    const bool queryGpuDb = (DBReader<unsigned int>::getExtendedDbtype(qDbr.getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    const bool targetGpuDb = (DBReader<unsigned int>::getExtendedDbtype(tDbr.getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    const bool decodeQueryDinuc = useDinuc || ((queryExt & Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0) || queryGpuDb;
    const bool decodeTargetDinuc = useDinuc || ((targetExt & Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0) || targetGpuDb;
    const int querySeqType = effectiveDecodeSeqType(qDbr.getDbtype(), decodeQueryDinuc);
    const int targetSeqType = effectiveDecodeSeqType(tDbr.getDbtype(), decodeTargetDinuc);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), par.threads,
                 par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    Debug::Progress progress(resReader.getSize());
    const bool parallelQueries = (!useRealModuleCm && resReader.getSize() > 1 && par.threads > 1);
    const bool parallelCandidates = (!parallelQueries && par.threads > 1);
    Debug(Debug::INFO) << "  query parallelism: " << (parallelQueries ? "outer" : "serial") << "\n";
    Debug(Debug::INFO) << "  target parallelism: " << (parallelCandidates ? "inner" : "shared across queries") << "\n";

    const auto processQuery = [&](size_t rid, unsigned int threadIdx) {
        progress.updateProgress();

        Sequence queryMapper(qDbr.getMaxSeqLen(), querySeqType, &subMat, 0, false, false);
        Sequence targetMapper(tDbr.getMaxSeqLen(), targetSeqType, &subMat, 0, false, false);

        const unsigned int queryKey = resReader.getDbKey(rid);
        const size_t qId = qDbr.getId(queryKey);
        if (qId == UINT_MAX) {
            dbw.writeData("", 0, queryKey, threadIdx);
            return;
        }

        const unsigned int qLen = qDbr.getSeqLen(qId);
        std::string querySeq = decodeDbSequence(qDbr, qId, threadIdx, queryMapper, decodeQueryDinuc, queryGpuDb, &subMat);
        FoldResult baseQueryFold = foldSequenceStructure(querySeq, structureSettings);
        if (baseQueryFold.structVec.empty()) {
            baseQueryFold.structVec.assign(querySeq.size(), StructVec());
        }

        std::vector<RnaMatcher::result_t> rawHits;
        RnaMatcher::readAlignmentResults(rawHits, resReader.getData(rid, threadIdx), false);
        if (rawHits.empty()) {
            dbw.writeData("", 0, queryKey, threadIdx);
            return;
        }

        std::unordered_map<unsigned int, size_t> bestByTarget;
        std::vector<CandidateHit> candidates;
        candidates.reserve(rawHits.size());
        for (size_t i = 0; i < rawHits.size(); ++i) {
            const RnaMatcher::result_t &raw = rawHits[i];
            std::unordered_map<unsigned int, size_t>::iterator it = bestByTarget.find(raw.dbKey);
            if (it != bestByTarget.end()) {
                if (!RnaMatcher::compareHits(raw, candidates[it->second].raw)) {
                    continue;
                }
                candidates[it->second].raw = raw;
                continue;
            }
            CandidateHit cand;
            cand.raw = raw;
            candidates.push_back(cand);
            bestByTarget[raw.dbKey] = candidates.size() - 1;
        }

        for (size_t i = 0; i < candidates.size(); ++i) {
            CandidateHit &cand = candidates[i];
            cand.oriented = orientHit(cand.raw);
            const size_t tId = tDbr.getId(cand.raw.dbKey);
            if (tId == UINT_MAX) {
                continue;
            }
            std::string targetSeq = decodeDbSequence(tDbr, tId, threadIdx, targetMapper, decodeTargetDinuc, targetGpuDb, &subMat);
            cand.window = buildWindowFromRaw(targetSeq.c_str(),
                                             static_cast<unsigned int>(targetSeq.size()),
                                             cand.oriented,
                                             static_cast<int>(qLen),
                                             flankFrac);
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const CandidateHit &a, const CandidateHit &b) {
                      return RnaMatcher::compareHits(a.raw, b.raw);
                  });
        std::unordered_map<unsigned int, size_t> candidateIndexByTarget;
        for (size_t i = 0; i < candidates.size(); ++i) {
            candidateIndexByTarget[candidates[i].raw.dbKey] = i;
        }

        const int maxFoldHits = std::max(64, parseEnvInt("MMSEQS_CMLITE_MAX_FOLD_HITS", 512));
        const size_t foldLimit = std::min(candidates.size(), static_cast<size_t>(maxFoldHits));
#pragma omp parallel for schedule(dynamic, 64) if(parallelCandidates)
        for (size_t i = 0; i < foldLimit; ++i) {
            if (!candidates[i].window.seq.empty()) {
                FoldResult fold = foldSequenceStructure(candidates[i].window.seq, structureSettings);
                candidates[i].window.structVec.swap(fold.structVec);
                candidates[i].window.partners.swap(fold.partners);
                candidates[i].window.pairMatrix.swap(fold.pairMatrix);
                candidates[i].window.sparsePairMatrix = fold.sparsePairMatrix;
            }
        }

        std::vector<ModelHit> modelHits;
        modelHits.reserve(uncappedModelHits ? candidates.size() : static_cast<size_t>(maxModelHits));
        for (size_t i = 0; i < candidates.size(); ++i) {
            if (!uncappedModelHits && static_cast<int>(modelHits.size()) >= maxModelHits) {
                break;
            }
            if (candidates[i].window.seq.empty()) {
                continue;
            }
            const AlignmentResult seedAln = rawAlignmentFromCandidate(querySeq, candidates[i]);
            if (!seedAln.valid) {
                continue;
            }
            if (candidates[i].raw.eval > msaEvalThr) {
                continue;
            }
            ModelHit mh;
            mh.candidateIdx = i;
            mh.seedWeight = modelWeightFromHit(candidates[i].raw);
            mh.weight = mh.seedWeight;
            mh.aln = seedAln;
            modelHits.push_back(mh);
        }
        applyModelDiversityWeights(querySeq, modelHits, candidates);

        std::vector<std::array<float, 4>> profile;
        std::vector<StructVec> queryStruct = baseQueryFold.structVec;
        buildProfile(querySeq, modelHits, candidates, baseQueryFold.structVec, alpha, profile, queryStruct);
        std::vector<float> pairSupport;
        std::vector<int> msaPartners = baseQueryFold.partners;
        const bool buildQuerySupport = querySupportMaxLen <= 0 || static_cast<int>(qLen) <= querySupportMaxLen;
        if (buildQuerySupport) {
            pairSupport = buildPairSupportMatrix(static_cast<int>(qLen), modelHits, candidates);
            const float minPairSupport = clampf(parseEnvFloat("MMSEQS_CMLITE_MIN_PAIR_SUPPORT", 0.18f), 0.0f, 1.0f);
            std::vector<int> supportPartners = derivePartnersFromPairSupport(pairSupport,
                                                                              static_cast<int>(qLen),
                                                                              structureSettings.minLoop,
                                                                              minPairSupport);
            if (hasAnyPartners(supportPartners)) {
                msaPartners.swap(supportPartners);
                queryStruct = pairMatrixToStructVec(pairSupport, msaPartners, static_cast<int>(qLen));
            }
        }
        const std::string consensusSeq = consensusSequenceFromProfile(querySeq, profile);
        std::vector<StemPattern> patterns = buildStemPatterns(consensusSeq, msaPartners, pairSupport, profile);
        if (patterns.empty()) {
            const std::vector<float> basePairMatrix = baseQueryFold.pairMatrix.empty() && baseQueryFold.sparsePairMatrix
                ? sparsePairMatrixToDense(*baseQueryFold.sparsePairMatrix)
                : baseQueryFold.pairMatrix;
            patterns = buildStemPatterns(querySeq, baseQueryFold.partners, basePairMatrix, profile);
        }

        std::unordered_map<unsigned int, std::vector<StemMatch>> refinedModuleMatchesByTarget;
        if (useRealModuleCm && !patterns.empty()) {
            const int requestedModuleJobs = parseEnvInt("MMSEQS_CMLITECM_MODULE_JOBS", 1);
            const int defaultThreadsPerModule = std::max(1, parseEnvInt("MMSEQS_CMLITECM_THREADS_PER_MODULE", std::max(1, par.threads)));
            int moduleJobs = requestedModuleJobs;
            if (moduleJobs <= 0) {
                moduleJobs = std::max(1, par.threads / defaultThreadsPerModule);
            }
            moduleJobs = std::min(moduleJobs, static_cast<int>(patterns.size()));
            moduleJobs = std::max(1, moduleJobs);
            const int moduleThreads = std::max(1, par.threads / moduleJobs);
            Debug(Debug::INFO) << "  query " << queryKey
                               << ": cmlitecm module jobs=" << moduleJobs
                               << ", threads/module=" << moduleThreads
                               << ", stems=" << patterns.size()
                               << "\n";
#pragma omp parallel for schedule(dynamic, 1) num_threads(moduleJobs) if(moduleJobs > 1)
            for (long pli = 0; pli < static_cast<long>(patterns.size()); ++pli) {
                const size_t pi = static_cast<size_t>(pli);
                std::unordered_map<unsigned int, std::vector<StemMatch>> moduleMatches;
                if (!refineStemMatchesWithRealCm(execPath,
#ifdef HAVE_INFERNAL_BRIDGE
                                                workerPool,
#endif
                                                moduleThreads,
                                                querySeq,
                                                patterns[pi],
                                                modelHits,
                                                candidates,
                                                candidateIndexByTarget,
                                                par,
                                                moduleMatches)) {
                    continue;
                }
                for (std::unordered_map<unsigned int, std::vector<StemMatch>>::iterator it = moduleMatches.begin();
                     it != moduleMatches.end();
                     ++it) {
                    for (size_t mi = 0; mi < it->second.size(); ++mi) {
                        it->second[mi].patternIdx = static_cast<int>(pi);
                    }
                }
                #pragma omp critical
                {
                    for (std::unordered_map<unsigned int, std::vector<StemMatch>>::iterator it = moduleMatches.begin();
                         it != moduleMatches.end();
                         ++it) {
                        std::vector<StemMatch> &bucket = refinedModuleMatchesByTarget[it->first];
                        bucket.insert(bucket.end(), it->second.begin(), it->second.end());
                    }
                }
            }
            for (std::unordered_map<unsigned int, std::vector<StemMatch>>::iterator it = refinedModuleMatchesByTarget.begin();
                 it != refinedModuleMatchesByTarget.end();
                 ++it) {
                std::sort(it->second.begin(), it->second.end(),
                          [](const StemMatch &a, const StemMatch &b) {
                              if (a.leftQuery != b.leftQuery) {
                                  return a.leftQuery < b.leftQuery;
                              }
                              if (a.leftTarget != b.leftTarget) {
                                  return a.leftTarget < b.leftTarget;
                              }
                              return a.score > b.score;
                          });
            }
            Debug(Debug::INFO) << "  query " << queryKey
                               << ": real module CM stems=" << patterns.size()
                               << ", refined targets=" << refinedModuleMatchesByTarget.size()
                               << "\n";
        }

        std::vector<RnaMatcher::result_t> refined(candidates.size());
        std::vector<std::string> rowSignatures(candidates.size());
        size_t validRawAlignmentCount = 0;
        size_t rawFallbackCount = 0;
        size_t chainHitCount = 0;
        size_t acceptedRefinementCount = 0;
#pragma omp parallel for schedule(dynamic, 32) if(parallelCandidates)
        for (size_t i = 0; i < candidates.size(); ++i) {
            const CandidateHit &cand = candidates[i];
            if (cand.window.seq.empty()) {
                refined[i] = cand.raw;
                continue;
            }
            const AlignmentResult rawAln = rawAlignmentFromCandidate(querySeq, cand);
            if (rawAln.valid) {
#pragma omp atomic
                ++validRawAlignmentCount;
            }
            const StemGuide rawGuide = makeNeutralGuide(static_cast<int>(qLen), msaPartners);
            const AlignmentSummary rawSummary = summarizeAlignment(querySeq, queryStruct, rawGuide, cand.window, rawAln, cand.raw.dbLen);
            const auto emitRawFallback = [&]() {
#pragma omp atomic
                ++rawFallbackCount;
                if (rawAln.valid) {
                    refined[i] = emitResultFromAlignment(cand.raw, cand.window, rawAln, rawSummary, qLen, cand.raw.dbLen);
                } else {
                    refined[i] = cand.raw;
                }
                if (useRowDedup && rawAln.valid) {
                    ProjectedAlignment proj = projectAlignment(static_cast<int>(querySeq.size()), cand.window, rawAln);
                    rowSignatures[i] = proj.row;
                }
            };

            const int seedDiag = cand.window.localSeedStart - cand.oriented.qStart;
            const int seedBandRadius = estimateSeedBandRadius(cand, bandExtra);
            std::vector<StemMatch> matches;
            if (useRealModuleCm) {
                const std::unordered_map<unsigned int, std::vector<StemMatch>>::const_iterator hitIt =
                    refinedModuleMatchesByTarget.find(cand.raw.dbKey);
                if (hitIt != refinedModuleMatchesByTarget.end()) {
                    matches = hitIt->second;
                }
            }
            if (matches.empty()) {
                matches = findStemMatches(patterns, cand.window, seedDiag, seedBandRadius);
            }
            const std::vector<std::vector<int>> chains = buildTopChainsFromStemMatches(patterns, matches);
            if (!chains.empty()) {
#pragma omp atomic
                ++chainHitCount;
            }
            if (chains.empty()) {
                emitRawFallback();
                continue;
            }

            const size_t maxEvaluatedChains = static_cast<size_t>(std::max(1, parseEnvInt("MMSEQS_CMLITE_EVAL_CHAINS", 4)));
            const float chainScoreWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_CHAIN_SCORE_WEIGHT", 0.10f));
            const float chainSupportWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_CHAIN_SUPPORT_WEIGHT", 0.12f));
            const float chainCoverageWeight = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_CHAIN_COVERAGE_WEIGHT", 0.0f));
            bool foundCandidate = false;
            float bestObjective = NEG_INF_F;
            float bestChainScore = 0.0f;
            float bestAnchoredFrac = 0.0f;
            StemGuide bestGuide;
            AlignmentResult bestAln;
            AlignmentSummary bestSummary;
            for (size_t ci = 0; ci < chains.size() && ci < maxEvaluatedChains; ++ci) {
                StemGuide chainGuide;
                int diagCenter = seedDiag;
                int bandRadius = seedBandRadius;
                float chainScore = 0.0f;
                if (!buildGuideFromStemChain(static_cast<int>(qLen),
                                             patterns,
                                             matches,
                                             seedDiag,
                                             seedBandRadius,
                                             msaPartners,
                                             chains[ci],
                                             chainGuide,
                                             diagCenter,
                                             bandRadius,
                                             chainScore)) {
                    continue;
                }
                const std::vector<AnchorBlock> blocks = buildAnchorBlocksFromChain(patterns, matches, chains[ci]);
                const std::pair<int, int> qSegment = chooseAnchoredQuerySegment(blocks, cand.oriented, static_cast<int>(qLen));
                AlignmentResult aln = stitchCmLiteAlignment(querySeq,
                                                            profile,
                                                            queryStruct,
                                                            msaPartners,
                                                            blocks,
                                                            cand.window.seq,
                                                            cand.window.structVec,
                                                            cand.window.partners,
                                                            cand.window.pairMatrix,
                                                            &cand.window,
                                                            qSegment.first,
                                                            qSegment.second,
                                                            alpha);
                if (!aln.valid) {
                    continue;
                }
                AlignmentSummary summary = summarizeAlignment(querySeq, queryStruct, chainGuide, cand.window, aln, cand.raw.dbLen);
                const float anchored = std::count_if(chainGuide.supportByQuery.begin(),
                                                     chainGuide.supportByQuery.end(),
                                                     [](float x) { return x > 0.0f; });
                const float anchoredFrac = (qLen > 0)
                    ? anchored / static_cast<float>(qLen)
                    : 0.0f;
                const auto computeObjective = [&](const AlignmentResult &candAln,
                                                  const AlignmentSummary &candSummary) {
                    return candAln.score
                         + chainScoreWeight * chainScore
                         + chainSupportWeight * candSummary.meanStemAnchor * static_cast<float>(std::max<size_t>(1, chains[ci].size()))
                         + chainCoverageWeight * anchoredFrac * static_cast<float>(qLen);
                };
                float objective = computeObjective(aln, summary);

                if (useFullQueryRefine) {
                    int fullDiagCenter = diagCenter;
                    if (aln.qStart >= 0 && aln.qEnd >= aln.qStart && aln.dbStart >= 0 && aln.dbEnd >= aln.dbStart) {
                        const int qMid = (aln.qStart + aln.qEnd) / 2;
                        const int tMid = (aln.dbStart + aln.dbEnd) / 2;
                        fullDiagCenter = tMid - qMid;
                    }
                    const int fullBandRadius = std::max(bandRadius,
                                                        (aln.dbEnd >= aln.dbStart)
                                                            ? ((aln.dbEnd - aln.dbStart + 1) + fullRefineBandSlack)
                                                            : (bandRadius + fullRefineBandSlack));
                    AlignmentResult fullAln = alignProfileLocalLocal(querySeq,
                                                                     profile,
                                                                     queryStruct,
                                                                     chainGuide,
                                                                     cand.window.seq,
                                                                     cand.window.structVec,
                                                                     0,
                                                                     static_cast<int>(qLen) - 1,
                                                                     alpha,
                                                                     fullRefineStemWeight,
                                                                     stemRadius,
                                                                     fullDiagCenter,
                                                                     fullBandRadius,
                                                                     NULL);
                    if (fullAln.valid) {
                        const AlignmentSummary fullSummary = summarizeAlignment(querySeq, queryStruct, chainGuide, cand.window, fullAln, cand.raw.dbLen);
                        const float fullObjective = computeObjective(fullAln, fullSummary)
                                                 + fullRefineCoverageWeight * fullSummary.qcov * static_cast<float>(qLen);
                        if (fullObjective > objective) {
                            aln = fullAln;
                            summary = fullSummary;
                            objective = fullObjective;
                        }

                        const int pairQLo = std::max(0, fullAln.qStart - fullRefineWindowSlack);
                        const int pairQHi = std::min(static_cast<int>(qLen) - 1, fullAln.qEnd + fullRefineWindowSlack);
                        const int pairTLo = std::max(0, fullAln.dbStart - fullRefineWindowSlack);
                        const int pairTHi = std::min(static_cast<int>(cand.window.seq.size()) - 1, fullAln.dbEnd + fullRefineWindowSlack);
                        AlignmentResult pairWindowAln = stitchPairAwareLocalWindow(querySeq,
                                                                                   profile,
                                                                                   queryStruct,
                                                                                   msaPartners,
                                                                                   cand.window.seq,
                                                                                   cand.window.structVec,
                                                                                   cand.window.partners,
                                                                                   cand.window.pairMatrix,
                                                                                   &cand.window,
                                                                                   pairQLo,
                                                                                   pairQHi,
                                                                                   pairTLo,
                                                                                   pairTHi,
                                                                                   alpha);
                        if (pairWindowAln.valid) {
                            const AlignmentSummary pairSummary = summarizeAlignment(querySeq, queryStruct, chainGuide, cand.window, pairWindowAln, cand.raw.dbLen);
                            const float pairObjective = computeObjective(pairWindowAln, pairSummary)
                                                     + fullRefineCoverageWeight * pairSummary.qcov * static_cast<float>(qLen);
                            if (pairObjective > objective) {
                                aln = pairWindowAln;
                                summary = pairSummary;
                                objective = pairObjective;
                            }
                        }
                    }
                }

                if (!foundCandidate || objective > bestObjective) {
                    foundCandidate = true;
                    bestObjective = objective;
                    bestGuide = chainGuide;
                    bestAln = aln;
                    bestSummary = summary;
                    bestChainScore = chainScore;
                    bestAnchoredFrac = anchoredFrac;
                }
            }
            if (!foundCandidate) {
                emitRawFallback();
                continue;
            }
            const float minConsistencyGain = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_MIN_CONSISTENCY_GAIN", 0.01f));
            const float minSeqIdGain = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_MIN_SEQID_GAIN", 0.005f));
            const float minChainScore = std::max(0.0f, parseEnvFloat("MMSEQS_CMLITE_ACCEPT_MIN_CHAIN_SCORE", 0.75f));
            float minAcceptedQcov = std::max(0.20f, rawSummary.qcov * 0.7f);
            if (useFullQueryRefine) {
                const float anchoredQcovFloor = std::max(0.10f, bestAnchoredFrac * 0.80f);
                minAcceptedQcov = std::min(minAcceptedQcov, anchoredQcovFloor);
            }
            if (!acceptAllRefinements
                && (bestSummary.qcov < minAcceptedQcov
                    || bestSummary.consistency + minConsistencyGain < rawSummary.consistency
                    || bestSummary.seqId + minSeqIdGain < rawSummary.seqId
                    || bestChainScore < minChainScore)) {
                emitRawFallback();
                continue;
            }

            refined[i] = emitResultFromAlignment(cand.raw, cand.window, bestAln, bestSummary, qLen, cand.raw.dbLen);
            #pragma omp atomic
            ++acceptedRefinementCount;
            if (useRowDedup) {
                ProjectedAlignment proj = projectAlignment(static_cast<int>(querySeq.size()), cand.window, bestAln);
                rowSignatures[i] = proj.row;
            }
        }

        refined.erase(std::remove_if(refined.begin(), refined.end(),
                                     [](const RnaMatcher::result_t &r) {
                                         return r.qStartPos < 0 || r.dbStartPos < 0 || r.backtrace.empty();
                                     }),
                      refined.end());
        std::sort(refined.begin(), refined.end(), RnaMatcher::compareHits);

        std::string out;
        out.reserve(refined.size() * 64);
        char buffer[4096];
        std::unordered_set<std::string> seenRows;
        for (size_t i = 0; i < refined.size(); ++i) {
            if (useRowDedup) {
                const std::unordered_map<unsigned int, size_t>::const_iterator it = candidateIndexByTarget.find(refined[i].dbKey);
                if (it != candidateIndexByTarget.end()) {
                    const std::string &sig = rowSignatures[it->second];
                    if (!sig.empty()) {
                        std::pair<std::unordered_set<std::string>::iterator, bool> inserted = seenRows.insert(sig);
                        if (!inserted.second) {
                            continue;
                        }
                    }
                }
            }
            const size_t written = RnaMatcher::resultToBuffer(buffer, refined[i], par.addBacktrace, true, false);
            out.append(buffer, written);
        }
        Debug(Debug::INFO) << "  query " << queryKey
                           << ": raw hits=" << rawHits.size()
                           << ", candidates=" << candidates.size()
                           << ", model hits=" << modelHits.size()
                           << ", valid raw alignments=" << validRawAlignmentCount
                           << ", candidates with chains=" << chainHitCount
                           << ", accepted refinements=" << acceptedRefinementCount
                           << ", raw fallbacks=" << rawFallbackCount
                           << ", final rows=" << refined.size()
                           << ", output bytes=" << out.size()
                           << "\n";
        dbw.writeData(out.c_str(), out.size(), queryKey, threadIdx);
    };

    if (parallelQueries) {
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t rid = 0; rid < resReader.getSize(); ++rid) {
            unsigned int threadIdx = 0;
#ifdef OPENMP
            threadIdx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            processQuery(rid, threadIdx);
        }
    } else {
        for (size_t rid = 0; rid < resReader.getSize(); ++rid) {
            processQuery(rid, 0);
        }
    }

    dbw.close(true);
    resReader.close();
    tDbr.close();
    qDbr.close();
#ifdef HAVE_INFERNAL_BRIDGE
    if (workerPool != NULL) {
        InfernalBridge::stopWorkerPool(workerPool);
    }
#endif
    return EXIT_SUCCESS;
}

static int runLolalignImpl(int argc,
                           const char **argv,
                           const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    par.evalThr = DBL_MAX;
    if (MMseqsMPI::isMaster() == false) {
        return EXIT_SUCCESS;
    }

    const int bandExtra = std::max(8, parseEnvInt("MMSEQS_LOLALIGN_BAND_EXTRA", 32));
    const bool useDinuc = par.emsearchDinucleotide || (parseEnvInt("MMSEQS_LOLALIGN_DINUC", 0) != 0);
    const bool acceptAllRefinements = parseEnvInt("MMSEQS_LOLALIGN_ACCEPT_ALL", 0) != 0;
    const bool useRowDedup = parseEnvInt("MMSEQS_LOLALIGN_DEDUP_ROWS", 0) != 0;
    const double msaEvalThr = par.lolalignMsaEvalThr;
    const bool uncappedModelHits = std::isfinite(msaEvalThr);
    const int maxModelHits = std::max(4, parseEnvInt("MMSEQS_LOLALIGN_MAX_MODEL_HITS", 64));
    const int querySupportMaxLen = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_QUERY_SUPPORT_MAX_LEN", 0));
    const float flankFrac = (par.cmRegionFlanking > 0.0f)
        ? par.cmRegionFlanking
        : parseEnvFloat("MMSEQS_LOLALIGN_FLANK_FRAC", 1.5f);
    const int candidateBatchSize = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_BATCH_SIZE", 2048));
    const float minConsistencyGain = parseEnvFloat("MMSEQS_LOLALIGN_MIN_CONSISTENCY_GAIN", 0.0f);
    const float minSeqIdGain = parseEnvFloat("MMSEQS_LOLALIGN_MIN_SEQID_GAIN", 0.0f);

    StructureSettings structureSettings;
    const char *lolalignBackendEnv = std::getenv("MMSEQS_LOLALIGN_STRUCTURE_BACKEND");
    structureSettings.backend = lolalignBackendEnv != NULL
        ? lolalignBackendEnv
        : "linearfold";
    structureSettings.minLoop = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_MIN_LOOP", 3));
    structureSettings.linearfoldMinLength = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_LINEARFOLD_MIN_LENGTH", 0));
    structureSettings.linearfoldBeamSize = std::max(1, parseEnvInt("MMSEQS_LOLALIGN_LINEARFOLD_BEAM", 100));
    structureSettings.rnaPlfoldWindow = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_RNAPLFOLD_WINDOW", 0));
    structureSettings.rnaPlfoldSpan = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_RNAPLFOLD_SPAN", 0));
    const char *keepPairMatrixEnv = std::getenv("MMSEQS_LOLALIGN_KEEP_PAIR_MATRIX");
    structureSettings.keepPairMatrix = keepPairMatrixEnv != NULL
        ? std::atoi(keepPairMatrixEnv) != 0
        : true;
    structureSettings.sparsePairMatrix = parseEnvInt("MMSEQS_LOLALIGN_SPARSE_PAIR_MATRIX", 0) != 0;
    structureSettings.dotBracketOnly = parseEnvInt("MMSEQS_LOLALIGN_DOTBRACKET_ONLY", 0) != 0;
    if (structureSettings.dotBracketOnly) {
        structureSettings.keepPairMatrix = false;
        structureSettings.sparsePairMatrix = false;
    }

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0f);

    std::ostringstream msaEvalStream;
    msaEvalStream << std::scientific << msaEvalThr;
    Debug(Debug::INFO) << "Running LoL-style RNA realignment\n";
    Debug(Debug::INFO) << "  refinement mode: ignoring alignment E-value threshold, candidate set comes from input result DB\n";
    Debug(Debug::INFO) << "  model hits: "
                       << (uncappedModelHits
                               ? "all deduplicated hits with E-value <= threshold"
                               : ("top " + SSTR(maxModelHits) + " deduplicated hits (threshold is inf)"))
                       << ", query support max length: " << querySupportMaxLen
                       << ", flank frac: " << flankFrac
                       << ", msa e-value threshold: " << msaEvalStream.str()
                       << ", structure backend: " << structureSettings.backend
                       << ", candidate batch size: " << candidateBatchSize
                       << ", linearfold min length: " << structureSettings.linearfoldMinLength
                       << ", linearfold beam: " << structureSettings.linearfoldBeamSize
                       << ", rnaplfold window: " << structureSettings.rnaPlfoldWindow
                       << ", rnaplfold span: " << structureSettings.rnaPlfoldSpan
                       << ", keep pair matrix: " << (structureSettings.keepPairMatrix ? "yes" : "no")
                       << ", sparse pair matrix: " << (structureSettings.sparsePairMatrix ? "yes" : "no")
                       << "\n";

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qDbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> tDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tDbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> resReader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                     DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const unsigned int queryExt = DBReader<unsigned int>::getExtendedDbtype(qDbr.getDbtype());
    const unsigned int targetExt = DBReader<unsigned int>::getExtendedDbtype(tDbr.getDbtype());
    const bool queryGpuDb = (DBReader<unsigned int>::getExtendedDbtype(qDbr.getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    const bool targetGpuDb = (DBReader<unsigned int>::getExtendedDbtype(tDbr.getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    const bool decodeQueryDinuc = useDinuc || ((queryExt & Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0) || queryGpuDb;
    const bool decodeTargetDinuc = useDinuc || ((targetExt & Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0) || targetGpuDb;
    const int querySeqType = effectiveDecodeSeqType(qDbr.getDbtype(), decodeQueryDinuc);
    const int targetSeqType = effectiveDecodeSeqType(tDbr.getDbtype(), decodeTargetDinuc);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), par.threads,
                 par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    Debug::Progress progress(resReader.getSize());
    const bool parallelCandidates = (par.threads > 1);
    Debug(Debug::INFO) << "  query parallelism: serial\n";
    Debug(Debug::INFO) << "  target parallelism: " << (parallelCandidates ? "inner" : "serial") << "\n";

    const auto processQuery = [&](size_t rid, unsigned int threadIdx) {
        progress.updateProgress();

        Sequence queryMapper(qDbr.getMaxSeqLen(), querySeqType, &subMat, 0, false, false);
        Sequence targetMapper(tDbr.getMaxSeqLen(), targetSeqType, &subMat, 0, false, false);

        const unsigned int queryKey = resReader.getDbKey(rid);
        const size_t qId = qDbr.getId(queryKey);
        if (qId == UINT_MAX) {
            dbw.writeData("", 0, queryKey, threadIdx);
            return;
        }

        const unsigned int qLen = qDbr.getSeqLen(qId);
        std::string querySeq = decodeDbSequence(qDbr, qId, threadIdx, queryMapper, decodeQueryDinuc, queryGpuDb, &subMat);
        FoldResult baseQueryFold = foldSequenceStructure(querySeq, structureSettings);
        if (baseQueryFold.structVec.empty()) {
            baseQueryFold.structVec.assign(querySeq.size(), StructVec());
        }

        std::vector<RnaMatcher::result_t> rawHits;
        RnaMatcher::readAlignmentResults(rawHits, resReader.getData(rid, threadIdx), false);
        if (rawHits.empty()) {
            dbw.writeData("", 0, queryKey, threadIdx);
            return;
        }

        std::unordered_map<unsigned int, size_t> bestByTarget;
        std::vector<CandidateHit> candidates;
        candidates.reserve(rawHits.size());
        for (size_t i = 0; i < rawHits.size(); ++i) {
            const RnaMatcher::result_t &raw = rawHits[i];
            std::unordered_map<unsigned int, size_t>::iterator it = bestByTarget.find(raw.dbKey);
            if (it != bestByTarget.end()) {
                if (!RnaMatcher::compareHits(raw, candidates[it->second].raw)) {
                    continue;
                }
                candidates[it->second].raw = raw;
                continue;
            }
            CandidateHit cand;
            cand.raw = raw;
            bestByTarget[raw.dbKey] = candidates.size();
            candidates.push_back(cand);
        }

        for (size_t i = 0; i < candidates.size(); ++i) {
            CandidateHit &cand = candidates[i];
            cand.oriented = orientHit(cand.raw);
            const size_t tId = tDbr.getId(cand.raw.dbKey);
            if (tId == UINT_MAX) {
                continue;
            }
            cand.window = buildWindowMetadataFromRaw(static_cast<unsigned int>(tDbr.getSeqLen(tId)),
                                                     cand.oriented,
                                                     static_cast<int>(qLen),
                                                     flankFrac);
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const CandidateHit &a, const CandidateHit &b) {
                      return RnaMatcher::compareHits(a.raw, b.raw);
                  });

        size_t maxWindowLen = 0;
        unsigned long long foldCacheDenseBytes = 0;
        for (size_t i = 0; i < candidates.size(); ++i) {
            const size_t windowLen = windowBoundsValid(candidates[i].window)
                ? static_cast<size_t>(candidates[i].window.fullHi - candidates[i].window.fullLo + 1)
                : 0u;
            maxWindowLen = std::max(maxWindowLen, windowLen);
            foldCacheDenseBytes = saturatingAddBytes(foldCacheDenseBytes, densePairMatrixBytesForLength(windowLen));
        }

        const unsigned long long MiB = 1024ULL * 1024ULL;
        const int denseCapMb = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_DENSE_PAIR_MATRIX_MAX_MB", 1024));
        const int denseCapMinLen = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_DENSE_PAIR_MATRIX_MIN_LEN", 2048));
        unsigned long long denseCapBytes = static_cast<unsigned long long>(denseCapMb) * MiB;
        const bool memoryGuardEnabled = parseEnvInt("MMSEQS_LOLALIGN_MEMORY_GUARD", 1) != 0;
        const unsigned long long explicitBudgetMb = parseEnvUnsignedLongLong("MMSEQS_LOLALIGN_MEMORY_BUDGET_MB", 0ULL);
        const float memoryBudgetFrac = clampf(parseEnvFloat("MMSEQS_LOLALIGN_MEMORY_FRACTION", 0.8f), 0.05f, 1.0f);
        const unsigned long long detectedAvailableBytes = lolalignAvailableMemoryBytes();
        const unsigned long long memoryBudgetBytes = explicitBudgetMb > 0ULL
            ? explicitBudgetMb * MiB
            : scaleBytes(detectedAvailableBytes, memoryBudgetFrac);
        const float cacheBudgetFrac = clampf(parseEnvFloat("MMSEQS_LOLALIGN_MEMORY_CACHE_FRACTION", 0.50f), 0.0f, 0.90f);
        const bool denseRetained = lolalignDensePairRetained(structureSettings);
        if (memoryGuardEnabled && memoryBudgetBytes > 0ULL && denseRetained) {
            const unsigned long long guardCacheBytes = std::max(1ULL, scaleBytes(memoryBudgetBytes, cacheBudgetFrac));
            denseCapBytes = (denseCapBytes > 0ULL) ? std::min(denseCapBytes, guardCacheBytes) : guardCacheBytes;
        }
        const bool denseCapActive = denseCapBytes > 0ULL
            && (maxWindowLen >= static_cast<size_t>(denseCapMinLen)
                || (memoryGuardEnabled && memoryBudgetBytes > 0ULL))
            && denseRetained;
        int effectiveBatchSize = candidateBatchSize;
        if (denseCapActive) {
            const unsigned long long perWindowBytes = std::max(1ULL, densePairMatrixBytesForLength(maxWindowLen));
            const unsigned long long capped = std::max(1ULL, denseCapBytes / perWindowBytes);
            effectiveBatchSize = std::max(1, std::min(candidateBatchSize, static_cast<int>(std::min<unsigned long long>(capped, static_cast<unsigned long long>(INT_MAX)))));
        }

        const bool foldCacheRequestedForMemory = parseEnvInt("MMSEQS_LOLALIGN_FOLD_CACHE", 1) != 0;
        unsigned long long estimatedPeakBytes = 0ULL;
        unsigned long long estimatedCacheBytes = 0ULL;
        unsigned long long perWindowRetainedBytes = 0ULL;
        unsigned long long perConcurrentFoldBytes = 0ULL;
        if (memoryGuardEnabled && memoryBudgetBytes > 0ULL) {
            const unsigned long long baseOverheadBytes = parseEnvUnsignedLongLong("MMSEQS_LOLALIGN_MEMORY_BASE_MB", 512ULL) * MiB;
            const float scratchMultiplier = std::max(1.0f, parseEnvFloat("MMSEQS_LOLALIGN_MEMORY_SCRATCH_MULTIPLIER", 6.0f));
            const int workerThreads = std::max(1, par.threads);
            const bool buildQuerySupportForMemory = querySupportMaxLen <= 0 || static_cast<int>(qLen) <= querySupportMaxLen;
            unsigned long long queryBytes = lolalignRetainedFoldBytes(querySeq.size(), structureSettings);
            if (buildQuerySupportForMemory) {
                queryBytes = saturatingAddBytes(queryBytes, densePairMatrixBytesForLength(querySeq.size()));
            }
            perWindowRetainedBytes = lolalignRetainedFoldBytes(maxWindowLen, structureSettings);
            perConcurrentFoldBytes = lolalignScratchFoldBytes(maxWindowLen, structureSettings, scratchMultiplier);
            estimatedCacheBytes = (foldCacheRequestedForMemory && denseRetained)
                ? std::min(foldCacheDenseBytes, denseCapActive ? denseCapBytes : foldCacheDenseBytes)
                : 0ULL;

            const unsigned long long fixedBytes = saturatingAddBytes(baseOverheadBytes,
                saturatingAddBytes(queryBytes, estimatedCacheBytes));
            const auto estimateBatchPeak = [&](int batchSize) -> unsigned long long {
                const unsigned long long batchBytes = saturatingMulBytes(perWindowRetainedBytes,
                    static_cast<unsigned long long>(std::max(1, batchSize)));
                const unsigned long long concurrentBytes = saturatingMulBytes(perConcurrentFoldBytes,
                    static_cast<unsigned long long>(std::min(workerThreads, std::max(1, batchSize))));
                return saturatingAddBytes(fixedBytes, saturatingAddBytes(batchBytes, concurrentBytes));
            };

            int lo = 1;
            int hi = effectiveBatchSize;
            int best = 1;
            while (lo <= hi) {
                const int mid = lo + (hi - lo) / 2;
                const unsigned long long peak = estimateBatchPeak(mid);
                if (peak <= memoryBudgetBytes) {
                    best = mid;
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }
            effectiveBatchSize = std::max(1, std::min(effectiveBatchSize, best));
            estimatedPeakBytes = estimateBatchPeak(effectiveBatchSize);
            if (effectiveBatchSize < candidateBatchSize || parseEnvInt("MMSEQS_LOLALIGN_MEMORY_VERBOSE", 0) != 0) {
                Debug(Debug::INFO) << "  lolalign memory estimate: query=" << qLen
                                   << ", max target window=" << maxWindowLen
                                   << ", threads=" << workerThreads
                                   << ", budget MB=" << (memoryBudgetBytes / MiB)
                                   << ", estimated peak MB=" << (estimatedPeakBytes / MiB)
                                   << ", per-window retained MB=" << (perWindowRetainedBytes / MiB)
                                   << ", per-concurrent-fold MB=" << (perConcurrentFoldBytes / MiB)
                                   << ", cache MB=" << (estimatedCacheBytes / MiB)
                                   << ", batch=" << effectiveBatchSize
                                   << "\n";
            }
        }

        std::vector<ModelHit> modelHits;
        modelHits.reserve(uncappedModelHits ? candidates.size() : static_cast<size_t>(maxModelHits));
        for (size_t i = 0; i < candidates.size(); ++i) {
            if (!uncappedModelHits && static_cast<int>(modelHits.size()) >= maxModelHits) {
                break;
            }
            if (candidates[i].raw.eval > msaEvalThr || !windowBoundsValid(candidates[i].window)) {
                continue;
            }
            if (!materializeCandidateWindow(candidates[i], tDbr, threadIdx, targetMapper,
                                            decodeTargetDinuc, targetGpuDb, &subMat)) {
                continue;
            }
            const AlignmentResult seedAln = rawAlignmentFromCandidate(querySeq, candidates[i]);
            if (!seedAln.valid) {
                clearWindowMaterialized(candidates[i].window);
                continue;
            }
            ModelHit mh;
            mh.candidateIdx = i;
            mh.seedWeight = modelWeightFromHit(candidates[i].raw);
            mh.weight = mh.seedWeight;
            mh.aln = seedAln;
            modelHits.push_back(mh);
        }
        applyModelDiversityWeights(querySeq, modelHits, candidates);

        std::vector<char> modelCandidate(candidates.size(), 0);
        for (size_t i = 0; i < modelHits.size(); ++i) {
            modelCandidate[modelHits[i].candidateIdx] = 1;
        }

        bool useFoldCache = parseEnvInt("MMSEQS_LOLALIGN_FOLD_CACHE", 1) != 0;
        const int foldCacheMaxUnique = std::max(0, parseEnvInt("MMSEQS_LOLALIGN_FOLD_CACHE_MAX_UNIQUE", 0));
        if (effectiveBatchSize != candidateBatchSize) {
            Debug(Debug::INFO) << "  candidate batch size capped for lolalign memory: requested="
                               << candidateBatchSize
                               << ", effective=" << effectiveBatchSize
                               << ", dense cap MB=" << (denseCapBytes / MiB)
                               << ", max window length=" << maxWindowLen
                               << ", min length=" << denseCapMinLen
                               << "\n";
        }
        struct CachedFold {
            std::vector<StructVec> structVec;
            std::vector<int> partners;
            std::shared_ptr<const std::vector<float>> pairMatrix;
            std::shared_ptr<const SparsePairMatrix> sparsePairMatrix;
        };
        std::vector<int> foldCacheIdx(candidates.size(), -1);
        std::vector<std::string> foldCacheSeqs;
        std::vector<CachedFold> foldCache;
        if (useFoldCache) {
            std::unordered_map<std::string, int> foldIndexBySeq;
            foldIndexBySeq.reserve(candidates.size());
            foldCacheSeqs.reserve(candidates.size());
            unsigned long long cachedDenseBytes = 0;
            size_t uniqueWindowCount = 0;
            size_t uncachedUniqueWindowCount = 0;
            const size_t cacheEntryLimit = foldCacheMaxUnique > 0
                ? static_cast<size_t>(foldCacheMaxUnique)
                : static_cast<size_t>(std::numeric_limits<size_t>::max());
            for (size_t i = 0; i < candidates.size(); ++i) {
                if (!windowBoundsValid(candidates[i].window)) {
                    continue;
                }
                if (!materializeCandidateWindow(candidates[i], tDbr, threadIdx, targetMapper,
                                                decodeTargetDinuc, targetGpuDb, &subMat)) {
                    continue;
                }
                const std::string seqKey = candidates[i].window.seq;
                std::unordered_map<std::string, int>::iterator it = foldIndexBySeq.find(seqKey);
                if (it != foldIndexBySeq.end()) {
                    foldCacheIdx[i] = it->second;
                    clearWindowMaterialized(candidates[i].window);
                    continue;
                }

                ++uniqueWindowCount;
                int cacheIdx = -1;
                const unsigned long long windowBytes = denseRetained
                    ? densePairMatrixBytesForLength(seqKey.size())
                    : 0ULL;
                const bool fitsByteBudget = !denseCapActive || saturatingAddBytes(cachedDenseBytes, windowBytes) <= denseCapBytes;
                const bool fitsEntryBudget = foldCacheSeqs.size() < cacheEntryLimit;
                if (fitsByteBudget && fitsEntryBudget) {
                    cacheIdx = static_cast<int>(foldCacheSeqs.size());
                    foldCacheSeqs.push_back(seqKey);
                    cachedDenseBytes = saturatingAddBytes(cachedDenseBytes, windowBytes);
                } else {
                    ++uncachedUniqueWindowCount;
                }
                foldIndexBySeq.insert(std::make_pair(seqKey, cacheIdx));
                foldCacheIdx[i] = cacheIdx;
                clearWindowMaterialized(candidates[i].window);
            }

            foldCache.resize(foldCacheSeqs.size());
#pragma omp parallel for schedule(dynamic, 32) if(parallelCandidates)
            for (size_t ci = 0; ci < foldCacheSeqs.size(); ++ci) {
                FoldResult fold = foldSequenceStructure(foldCacheSeqs[ci], structureSettings);
                foldCache[ci].structVec.swap(fold.structVec);
                foldCache[ci].partners.swap(fold.partners);
                if (!fold.pairMatrix.empty()) {
                    foldCache[ci].pairMatrix = std::make_shared<std::vector<float>>(std::move(fold.pairMatrix));
                }
                foldCache[ci].sparsePairMatrix = fold.sparsePairMatrix;
            }
            Debug(Debug::INFO) << "  fold cache: cached windows=" << foldCacheSeqs.size()
                               << ", uncached unique windows=" << uncachedUniqueWindowCount
                               << ", unique windows=" << uniqueWindowCount
                               << ", candidates=" << candidates.size()
                               << ", dense cache MB=" << (cachedDenseBytes / (1024ULL * 1024ULL))
                               << ", dense cap MB=" << (denseCapActive ? (denseCapBytes / MiB) : 0)
                               << ", entry cap=" << foldCacheMaxUnique
                               << "\n";
        }

        const auto assignCachedFold = [&](size_t idx) -> bool {
            const int cacheIdx = (idx < foldCacheIdx.size()) ? foldCacheIdx[idx] : -1;
            if (cacheIdx < 0 || static_cast<size_t>(cacheIdx) >= foldCache.size()) {
                return false;
            }
            const CachedFold &fold = foldCache[static_cast<size_t>(cacheIdx)];
            candidates[idx].window.structVec = fold.structVec;
            candidates[idx].window.partners = fold.partners;
            candidates[idx].window.pairMatrix.clear();
            candidates[idx].window.pairMatrixRef = fold.pairMatrix;
            candidates[idx].window.sparsePairMatrix = fold.sparsePairMatrix;
            return true;
        };

        for (size_t h = 0; h < modelHits.size(); ++h) {
            const size_t idx = modelHits[h].candidateIdx;
            if (idx < candidates.size()) {
                materializeCandidateWindow(candidates[idx], tDbr, threadIdx, targetMapper,
                                           decodeTargetDinuc, targetGpuDb, &subMat);
            }
        }

#pragma omp parallel for schedule(dynamic, 32) if(parallelCandidates)
        for (size_t i = 0; i < candidates.size(); ++i) {
            if (!modelCandidate[i] || candidates[i].window.seq.empty()) {
                continue;
            }
            if (!useFoldCache || !assignCachedFold(i)) {
                FoldResult fold = foldSequenceStructure(candidates[i].window.seq, structureSettings);
                candidates[i].window.structVec.swap(fold.structVec);
                candidates[i].window.partners.swap(fold.partners);
                candidates[i].window.pairMatrix.swap(fold.pairMatrix);
                candidates[i].window.sparsePairMatrix = fold.sparsePairMatrix;
            }
        }

        std::vector<std::array<float, 4>> profile;
        std::vector<StructVec> queryStruct = baseQueryFold.structVec;
        buildProfile(querySeq, modelHits, candidates, baseQueryFold.structVec, 0.35f, profile, queryStruct);
        std::vector<float> pairSupport;
        std::vector<int> msaPartners = baseQueryFold.partners;
        const bool buildQuerySupport = querySupportMaxLen <= 0 || static_cast<int>(qLen) <= querySupportMaxLen;
        if (buildQuerySupport) {
            pairSupport = buildPairSupportMatrix(static_cast<int>(qLen), modelHits, candidates);
            const float minPairSupport = clampf(parseEnvFloat("MMSEQS_LOLALIGN_MIN_PAIR_SUPPORT", 0.18f), 0.0f, 1.0f);
            std::vector<int> supportPartners = derivePartnersFromPairSupport(pairSupport,
                                                                              static_cast<int>(qLen),
                                                                              structureSettings.minLoop,
                                                                              minPairSupport);
            if (hasAnyPartners(supportPartners)) {
                msaPartners.swap(supportPartners);
                queryStruct = pairMatrixToStructVec(pairSupport, msaPartners, static_cast<int>(qLen));
            }
        }
        std::string lolConsensusSeq = consensusSequenceFromProfile(querySeq, profile);
        std::vector<StemPattern> lolPatterns = buildStemPatterns(lolConsensusSeq, msaPartners, pairSupport, profile);
        if (lolPatterns.empty()) {
            lolPatterns = buildStemPatterns(querySeq, msaPartners, pairSupport, profile);
        }

        for (size_t i = 0; i < candidates.size(); ++i) {
            clearWindowMaterialized(candidates[i].window);
        }
        Debug(Debug::INFO) << "  model/profile windows released before batch refinement: query length="
                           << qLen
                           << ", candidates=" << candidates.size()
                           << "\n";

        size_t validRawAlignmentCount = 0;
        size_t acceptedRefinementCount = 0;
        size_t rawFallbackCount = 0;

        std::string out;
        out.reserve(candidates.size() * 128);
        std::unordered_set<std::string> seenRows;
        size_t emittedRowCount = 0;
        if (useRowDedup) {
            seenRows.reserve(candidates.size() * 2);
        }
        char buffer[4096];
        for (size_t batchStart = 0; batchStart < candidates.size(); batchStart += static_cast<size_t>(effectiveBatchSize)) {
            const size_t batchEnd = std::min(candidates.size(), batchStart + static_cast<size_t>(effectiveBatchSize));
            const size_t batchLen = batchEnd - batchStart;

#pragma omp parallel if(parallelCandidates)
            {
                unsigned int localThreadIdx = threadIdx;
#ifdef OPENMP
                localThreadIdx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                Sequence batchTargetMapper(tDbr.getMaxSeqLen(), targetSeqType, &subMat, 0, false, false);
#pragma omp for schedule(dynamic, 32)
                for (size_t bi = 0; bi < batchLen; ++bi) {
                    const size_t idx = batchStart + bi;
                    if (!windowBoundsValid(candidates[idx].window)) {
                        continue;
                    }
                    if (!materializeCandidateWindow(candidates[idx], tDbr, localThreadIdx, batchTargetMapper,
                                                    decodeTargetDinuc, targetGpuDb, &subMat)) {
                        continue;
                    }
                    if (!candidates[idx].window.structVec.empty()
                        && candidates[idx].window.structVec.size() == candidates[idx].window.seq.size()) {
                        continue;
                    }
                    if (!useFoldCache || !assignCachedFold(idx)) {
                        FoldResult fold = foldSequenceStructure(candidates[idx].window.seq, structureSettings);
                        candidates[idx].window.structVec.swap(fold.structVec);
                        candidates[idx].window.partners.swap(fold.partners);
                        candidates[idx].window.pairMatrix.swap(fold.pairMatrix);
                        candidates[idx].window.sparsePairMatrix = fold.sparsePairMatrix;
                    }
                }
            }

            std::vector<RnaMatcher::result_t> refined(batchLen);
            std::vector<std::string> rowSignatures(batchLen);

#pragma omp parallel for schedule(dynamic, 32) if(parallelCandidates)
            for (size_t bi = 0; bi < batchLen; ++bi) {
                const size_t idx = batchStart + bi;
                const CandidateHit &cand = candidates[idx];
                if (cand.window.seq.empty()) {
                    refined[bi] = cand.raw;
                    continue;
                }
                const AlignmentResult rawAln = rawAlignmentFromCandidate(querySeq, cand);
                if (rawAln.valid) {
#pragma omp atomic
                    ++validRawAlignmentCount;
                }
                const StemGuide rawGuide = makeNeutralGuide(static_cast<int>(qLen), msaPartners);
                const AlignmentSummary rawSummary = summarizeAlignment(querySeq, queryStruct, rawGuide, cand.window, rawAln, cand.raw.dbLen);
                const auto emitRawFallback = [&]() {
#pragma omp atomic
                    ++rawFallbackCount;
                    if (rawAln.valid) {
                        refined[bi] = emitResultFromAlignment(cand.raw, cand.window, rawAln, rawSummary, qLen, cand.raw.dbLen);
                    } else {
                        refined[bi] = cand.raw;
                    }
                    if (useRowDedup && rawAln.valid) {
                        ProjectedAlignment proj = projectAlignment(static_cast<int>(querySeq.size()), cand.window, rawAln);
                        rowSignatures[bi] = proj.row;
                    }
                };

                const int seedDiag = cand.window.localSeedStart - cand.oriented.qStart;
                const int seedBandRadius = estimateSeedBandRadius(cand, bandExtra);
                const AlignmentResult lolAln = runLolAlignIteration(querySeq,
                                                                    profile,
                                                                    queryStruct,
                                                                    msaPartners,
                                                                    pairSupport,
                                                                    lolPatterns,
                                                                    cand.window,
                                                                    0,
                                                                    static_cast<int>(qLen) - 1,
                                                                    seedDiag,
                                                                    seedBandRadius);
                if (!lolAln.valid) {
                    emitRawFallback();
                    continue;
                }

                const StemGuide lolGuide = buildLolGuideFromAlignment(static_cast<int>(qLen),
                                                                      msaPartners,
                                                                      pairSupport,
                                                                      cand.window,
                                                                      lolAln);
                const AlignmentSummary lolSummary = summarizeAlignment(querySeq, queryStruct, lolGuide, cand.window, lolAln, cand.raw.dbLen);
                const bool accept = acceptAllRefinements
                    || (!rawAln.valid)
                    || (lolSummary.qcov >= std::max(0.10f, rawSummary.qcov * 0.60f)
                        && lolSummary.consistency + minConsistencyGain >= rawSummary.consistency
                        && lolSummary.seqId + minSeqIdGain >= rawSummary.seqId);
                if (!accept) {
                    emitRawFallback();
                    continue;
                }

#pragma omp atomic
                ++acceptedRefinementCount;
                refined[bi] = emitResultFromAlignment(cand.raw, cand.window, lolAln, lolSummary, qLen, cand.raw.dbLen);
                if (useRowDedup) {
                    ProjectedAlignment proj = projectAlignment(static_cast<int>(querySeq.size()), cand.window, lolAln);
                    rowSignatures[bi] = proj.row;
                }
            }

            for (size_t bi = 0; bi < batchLen; ++bi) {
                if (useRowDedup && !rowSignatures[bi].empty()) {
                    if (!seenRows.insert(rowSignatures[bi]).second) {
                        continue;
                    }
                }
                const size_t written = RnaMatcher::resultToBuffer(buffer, refined[bi], par.addBacktrace, true, false);
                out.append(buffer, written);
                ++emittedRowCount;
            }

            for (size_t bi = 0; bi < batchLen; ++bi) {
                const size_t idx = batchStart + bi;
                clearWindowMaterialized(candidates[idx].window);
            }
        }
        for (size_t i = 0; i < candidates.size(); ++i) {
            clearWindowMaterialized(candidates[i].window);
        }
        Debug(Debug::INFO) << "  query " << queryKey
                           << ": raw hits=" << rawHits.size()
                           << ", candidates=" << candidates.size()
                           << ", model hits=" << modelHits.size()
                           << ", valid raw alignments=" << validRawAlignmentCount
                           << ", accepted refinements=" << acceptedRefinementCount
                           << ", raw fallbacks=" << rawFallbackCount
                           << ", final rows=" << emittedRowCount
                           << ", output bytes=" << out.size()
                           << "\n";
        dbw.writeData(out.c_str(), out.size(), queryKey, threadIdx);
    };

    for (size_t rid = 0; rid < resReader.getSize(); ++rid) {
        processQuery(rid, 0);
    }

    dbw.close(true);
    resReader.close();
    tDbr.close();
    qDbr.close();
    return EXIT_SUCCESS;
}

int cmlite(int argc, const char **argv, const Command &command) {
    return runCmliteImpl(argc, argv, command, false);
}

int cmlitecm(int argc, const char **argv, const Command &command) {
    return runCmliteImpl(argc, argv, command, true);
}

int lolalign(int argc, const char **argv, const Command &command) {
    return runLolalignImpl(argc, argv, command);
}
