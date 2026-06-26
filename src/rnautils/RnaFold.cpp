#include "LocalCommandDeclarations.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "RnaTurnerParams.h"
#include "RnaFold.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cctype>
#include <limits>
#include <string>
#include <vector>

namespace {

static const int MIN_HAIRPIN_LOOP = 3;
static const int INF = 1000000000;
static const int INF_CLAMP = INF / 2;

static std::string normalizeRnaSequence(const char *data, size_t len) {
    std::string seq;
    seq.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        char c = static_cast<char>(std::toupper(static_cast<unsigned char>(data[i])));
        if (c == 'T') {
            c = 'U';
        }
        switch (c) {
            case 'A':
            case 'C':
            case 'G':
            case 'U':
            case 'N':
                seq.push_back(c);
                break;
            default:
                seq.push_back('N');
                break;
        }
    }
    return seq;
}

static int pairType(char left, char right) {
    if (left == 'C' && right == 'G') {
        return 1;
    }
    if (left == 'G' && right == 'C') {
        return 2;
    }
    if (left == 'G' && right == 'U') {
        return 3;
    }
    if (left == 'U' && right == 'G') {
        return 4;
    }
    if (left == 'A' && right == 'U') {
        return 5;
    }
    if (left == 'U' && right == 'A') {
        return 6;
    }
    return 0;
}

static bool canPair(char left, char right) {
    return pairType(left, right) > 0;
}

static int terminalPairPenalty(int type) {
    return (type > 2) ? 50 : 0;
}


static int ntIndex(char c) {
    switch (c) {
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'U': return 4;
        default: return 0;
    }
}

static int hairpinMismatchEnergy(int type, char leftMismatch, char rightMismatch) {
    static const int mismatchH37[8][5][5] = {
        {{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}},
        {{-80, -100, -110, -100, -80}, {-140, -150, -150, -140, -150}, {-80, -100, -110, -100, -80}, {-150, -230, -150, -240, -150}, {-100, -100, -140, -100, -210}},
        {{-50, -110, -70, -110, -50}, {-110, -110, -150, -130, -150}, {-50, -110, -70, -110, -50}, {-150, -250, -150, -220, -150}, {-100, -110, -100, -110, -160}},
        {{20, 20, -20, -10, -20}, {20, 20, -50, -30, -50}, {-10, -10, -20, -10, -20}, {-50, -100, -50, -110, -50}, {-10, -10, -30, -10, -100}},
        {{0, -20, -10, -20, 0}, {-30, -50, -30, -60, -30}, {0, -20, -10, -20, 0}, {-30, -90, -30, -110, -30}, {-10, -20, -10, -20, -90}},
        {{-10, -10, -20, -10, -20}, {-30, -30, -50, -30, -50}, {-10, -10, -20, -10, -20}, {-50, -120, -50, -110, -50}, {-10, -10, -30, -10, -120}},
        {{0, -20, -10, -20, 0}, {-30, -50, -30, -50, -30}, {0, -20, -10, -20, 0}, {-30, -150, -30, -150, -30}, {-10, -20, -10, -20, -90}},
        {{20, 20, -10, -10, 0}, {20, 20, -30, -30, -30}, {0, -10, -10, -10, 0}, {-30, -90, -30, -110, -30}, {-10, -10, -10, -10, -90}}
    };
    if (type <= 0) {
        return 0;
    }
    return mismatchH37[type][ntIndex(leftMismatch)][ntIndex(rightMismatch)];
}

static int dangle5Energy(int type, char base) {
    static const int dangle5_37[8][5] = {
        {INF, INF, INF, INF, INF},
        {-10, -50, -30, -20, -10},
        {0, -20, -30, 0, 0},
        {-20, -30, -30, -40, -20},
        {-10, -30, -10, -20, -20},
        {-20, -30, -30, -40, -20},
        {-10, -30, -10, -20, -20},
        {0, -20, -10, 0, 0}
    };
    if (type <= 0) {
        return 0;
    }
    return dangle5_37[type][ntIndex(base)];
}

static int mismatchMEnergy(int type, char leftMismatch, char rightMismatch) {
    static const int mismatchM37[8][5][5] = {
        {{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}},
        {{-50, -110, -50, -140, -70}, {-110, -110, -110, -160, -110}, {-70, -150, -70, -150, -100}, {-110, -130, -110, -140, -110}, {-50, -150, -50, -150, -70}},
        {{-80, -140, -80, -140, -100}, {-100, -150, -100, -140, -100}, {-110, -150, -110, -150, -140}, {-100, -140, -100, -160, -100}, {-80, -150, -80, -150, -120}},
        {{-50, -80, -50, -50, -50}, {-50, -100, -70, -50, -70}, {-60, -80, -60, -80, -60}, {-70, -110, -70, -80, -70}, {-50, -80, -50, -80, -50}},
        {{-30, -30, -60, -60, -60}, {-30, -30, -60, -60, -60}, {-70, -100, -70, -100, -80}, {-60, -80, -60, -80, -60}, {-60, -100, -70, -100, -60}},
        {{-50, -80, -50, -80, -50}, {-70, -100, -70, -110, -70}, {-60, -80, -60, -80, -60}, {-70, -110, -70, -120, -70}, {-50, -80, -50, -80, -50}},
        {{-60, -80, -60, -80, -60}, {-60, -80, -60, -80, -60}, {-70, -100, -70, -100, -80}, {-60, -80, -60, -80, -60}, {-70, -100, -70, -100, -80}},
        {{-30, -30, -50, -50, -50}, {-30, -30, -60, -50, -60}, {-60, -80, -60, -80, -60}, {-60, -80, -60, -80, -60}, {-50, -80, -50, -80, -50}}
    };
    if (type <= 0) {
        return 0;
    }
    return mismatchM37[type][ntIndex(leftMismatch)][ntIndex(rightMismatch)];
}

static int dangle3Energy(int type, char base) {
    static const int dangle3_37[8][5] = {
        {INF, INF, INF, INF, INF},
        {-40, -110, -40, -130, -60},
        {-80, -170, -80, -170, -120},
        {-10, -70, -10, -70, -10},
        {-50, -80, -50, -80, -60},
        {-10, -70, -10, -70, -10},
        {-50, -80, -50, -80, -60},
        {-10, -70, -10, -70, -10}
    };
    if (type <= 0) {
        return 0;
    }
    return dangle3_37[type][ntIndex(base)];
}

static int exteriorStemEnergy(int type, char leftDangle, char rightDangle, bool hasLeftDangle, bool hasRightDangle) {
    int energy = terminalPairPenalty(type);
    if (hasLeftDangle && hasRightDangle) {
        energy += mismatchMEnergy(type, leftDangle, rightDangle);
    } else {
        if (hasLeftDangle) {
            energy += dangle5Energy(type, leftDangle);
        }
        if (hasRightDangle) {
            energy += dangle3Energy(type, rightDangle);
        }
    }
    return energy;
}

static int multibranchStemEnergy(int type, char leftDangle, char rightDangle, bool hasLeftDangle, bool hasRightDangle) {
    int energy = -90;
    if (hasLeftDangle && hasRightDangle) {
        energy += mismatchMEnergy(type, leftDangle, rightDangle);
    } else if (hasLeftDangle) {
        energy += dangle5Energy(type, leftDangle);
    } else if (hasRightDangle) {
        energy += dangle3Energy(type, rightDangle);
    }
    return energy + terminalPairPenalty(type);
}

static int specialHairpinPenalty(const std::string &seq, size_t i, size_t loopSize) {
    static const char *triloops[] = {
        "CAACG", "GUUAC"
    };
    static const int triloop37[] = {680, 690};

    static const char *tetraloops[] = {
        "CAACGG", "CCAAGG", "CCACGG", "CCCAGG", "CCGAGG", "CCGCGG", "CCUAGG", "CCUCGG",
        "CUAAGG", "CUACGG", "CUCAGG", "CUCCGG", "CUGCGG", "CUUAGG", "CUUCGG", "CUUUGG"
    };
    static const int tetraloop37[] = {550, 330, 370, 340, 350, 360, 370, 250, 360, 280, 370, 270, 280, 350, 370, 370};

    static const char *hexaloops[] = {
        "ACAGUACU", "ACAGUGAU", "ACAGUGCU", "ACAGUGUU"
    };
    static const int hexaloop37[] = {280, 360, 290, 180};

    if (loopSize == 3) {
        for (size_t k = 0; k < sizeof(triloop37) / sizeof(triloop37[0]); ++k) {
            if (seq.compare(i, 5, triloops[k]) == 0) {
                return triloop37[k];
            }
        }
    } else if (loopSize == 4) {
        for (size_t k = 0; k < sizeof(tetraloop37) / sizeof(tetraloop37[0]); ++k) {
            if (seq.compare(i, 6, tetraloops[k]) == 0) {
                return tetraloop37[k];
            }
        }
    } else if (loopSize == 6) {
        for (size_t k = 0; k < sizeof(hexaloop37) / sizeof(hexaloop37[0]); ++k) {
            if (seq.compare(i, 8, hexaloops[k]) == 0) {
                return hexaloop37[k];
            }
        }
    }
    return INF;
}

static int hairpinPenalty(const std::string &seq, size_t i, size_t j, int closingType) {
    const size_t loopSize = j - i - 1;
    static const int hairpin37[31] = {
        INF, INF, INF, 540, 560, 570, 540, 600, 550, 640, 650,
        660, 670, 680, 690, 690, 700, 710, 710, 720, 720, 730,
        730, 740, 740, 750, 750, 750, 760, 760, 770
    };
    int energy = 0;
    if (loopSize <= 30) {
        energy = hairpin37[loopSize];
    } else {
        energy = hairpin37[30] + static_cast<int>(107.856 * std::log(static_cast<double>(loopSize) / 30.0));
    }

    const int specialEnergy = specialHairpinPenalty(seq, i, loopSize);
    if (specialEnergy < INF) {
        return specialEnergy;
    }

    if (loopSize == 3) {
        return energy + terminalPairPenalty(closingType);
    }
    if (loopSize > 3) {
        energy += hairpinMismatchEnergy(closingType, seq[i + 1], seq[j - 1]);
    }
    return energy;
}

static int stackEnergy(int outerType, int innerType) {
    static const int stack37[8][8] = {
        {INF, INF, INF, INF, INF, INF, INF, INF},
        {INF, -240, -330, -210, -140, -210, -210, -140},
        {INF, -330, -340, -250, -150, -220, -240, -150},
        {INF, -210, -250,  130,  -50, -140, -130,  130},
        {INF, -140, -150,  -50,   30,  -60, -100,   30},
        {INF, -210, -220, -140,  -60, -110,  -90,  -60},
        {INF, -210, -240, -130, -100,  -90, -130,  -90},
        {INF, -140, -150,  130,   30,  -60,  -90,  130}
    };
    if (outerType <= 0 || innerType <= 0) {
        return INF;
    }
    return stack37[outerType][innerType];
}

static int internalLoopPenalty(const std::string &seq, size_t i, size_t j, size_t p, size_t q, size_t leftLoop, size_t rightLoop, int outerType, int innerType) {
    static const int bulge37[31] = {
        INF, 380, 280, 320, 360, 400, 440, 460, 470, 480, 490,
        500, 510, 520, 530, 540, 540, 550, 550, 560, 570, 570,
        580, 580, 580, 590, 590, 600, 600, 600, 610
    };
    static const int internalLoop37[31] = {
        INF, INF, 100, 100, 110, 200, 200, 210, 230, 240, 250,
        260, 270, 280, 290, 290, 300, 310, 310, 320, 330, 330,
        340, 340, 350, 350, 350, 360, 360, 370, 370
    };

    const size_t loopSize = leftLoop + rightLoop;
    if (loopSize == 0) {
        return stackEnergy(outerType, innerType);
    }

    if (leftLoop == 1 && rightLoop == 1) {
        return RnaTurnerParams::int11_37[outerType][innerType][ntIndex(seq[i + 1])][ntIndex(seq[j - 1])];
    }
    if (leftLoop == 1 && rightLoop == 2) {
        return RnaTurnerParams::int21_37[outerType][innerType][ntIndex(seq[i + 1])][ntIndex(seq[q + 1])][ntIndex(seq[j - 1])];
    }
    if (leftLoop == 2 && rightLoop == 1) {
        return RnaTurnerParams::int21_37[innerType][outerType][ntIndex(seq[q + 1])][ntIndex(seq[i + 1])][ntIndex(seq[p - 1])];
    }
    if (leftLoop == 2 && rightLoop == 2) {
        return RnaTurnerParams::int22_37[outerType][innerType][ntIndex(seq[i + 1])][ntIndex(seq[p - 1])][ntIndex(seq[q + 1])][ntIndex(seq[j - 1])];
    }

    if (leftLoop == 0 || rightLoop == 0) {
        int energy = 0;
        if (loopSize <= 30) {
            energy = bulge37[loopSize];
        } else {
            energy = bulge37[30] + static_cast<int>(107.856 * std::log(static_cast<double>(loopSize) / 30.0));
        }
        if (loopSize == 1) {
            energy += stackEnergy(outerType, innerType);
        } else {
            energy += terminalPairPenalty(outerType) + terminalPairPenalty(innerType);
        }
        return energy;
    }

    const size_t nl = std::max(leftLoop, rightLoop);
    const size_t ns = std::min(leftLoop, rightLoop);
    int energy = 0;
    if (ns == 1) {
        if (loopSize <= 30) {
            energy = internalLoop37[loopSize];
        } else {
            energy = internalLoop37[30] + static_cast<int>(107.856 * std::log(static_cast<double>(loopSize) / 30.0));
        }
        energy += static_cast<int>(std::min<size_t>((nl - ns) * 60, 300));
        energy += RnaTurnerParams::mismatch1nI37[outerType][ntIndex(seq[i + 1])][ntIndex(seq[j - 1])];
        energy += RnaTurnerParams::mismatch1nI37[innerType][ntIndex(seq[q + 1])][ntIndex(seq[p - 1])];
        return energy;
    }
    if (ns == 2 && nl == 3) {
        energy = internalLoop37[5] + 60;
        energy += RnaTurnerParams::mismatch23I37[outerType][ntIndex(seq[i + 1])][ntIndex(seq[j - 1])];
        energy += RnaTurnerParams::mismatch23I37[innerType][ntIndex(seq[q + 1])][ntIndex(seq[p - 1])];
        return energy;
    }
    if (loopSize <= 30) {
        energy = internalLoop37[loopSize];
    } else {
        energy = internalLoop37[30] + static_cast<int>(107.856 * std::log(static_cast<double>(loopSize) / 30.0));
    }
    energy += static_cast<int>(std::min<size_t>((nl - ns) * 60, 300));
    energy += RnaTurnerParams::mismatchI37[outerType][ntIndex(seq[i + 1])][ntIndex(seq[j - 1])];
    energy += RnaTurnerParams::mismatchI37[innerType][ntIndex(seq[q + 1])][ntIndex(seq[p - 1])];
    return energy;
}

struct FoldMatrices {
    size_t n;
    std::vector<size_t> rowOffset;
    std::vector<size_t> colOffset;
    std::vector<int> c;
    std::vector<int> mRow;
    std::vector<int> mCol;
    std::vector<int> f5;

    explicit FoldMatrices(size_t n_)
        : n(n_),
          rowOffset(n_),
          colOffset(n_),
          c(packedSize(n_), INF),
          mRow(packedSize(n_), INF),
          mCol(packedSize(n_), INF),
          f5(n_, 0) {
        size_t row = 0;
        size_t col = 0;
        for (size_t i = 0; i < n; ++i) {
            rowOffset[i] = row;
            row += n - i;
            colOffset[i] = col;
            col += i + 1;
        }
    }

    static size_t packedSize(size_t n) {
        return n * (n + 1) / 2;
    }

    size_t idx(size_t i, size_t j) const {
        return rowOffset[i] + (j - i);
    }

    size_t colIdx(size_t i, size_t j) const {
        return colOffset[j] + i;
    }
};

static int getC(const FoldMatrices &mx, size_t i, size_t j) {
    if (i >= mx.n || j >= mx.n || i > j) {
        return INF;
    }
    return mx.c[mx.idx(i, j)];
}

static int getM(const FoldMatrices &mx, size_t i, size_t j) {
    if (i >= mx.n || j >= mx.n || i > j) {
        return INF;
    }
    return mx.mRow[mx.idx(i, j)];
}

static int getM2(const FoldMatrices &mx, size_t i, size_t j) {
    if (i >= mx.n || j >= mx.n || i > j) {
        return INF;
    }
    int best = INF;
    size_t leftIndex = mx.idx(i, i + 1);
    size_t rightIndex = mx.colIdx(i + 2, j);
    for (size_t k = i + 1; k + 1 < j; ++k, ++leftIndex, ++rightIndex) {
        const int energy = mx.mRow[leftIndex] + mx.mCol[rightIndex];
        if (energy < best) {
            best = energy;
        }
    }
    return (best >= INF_CLAMP) ? INF : best;
}

static void tracebackC(const FoldMatrices &mx, const std::string &seq, size_t i, size_t j, std::string &structure);

static void tracebackM(const FoldMatrices &mx, const std::string &seq, size_t i, size_t j, std::string &structure) {
    if (i >= mx.n || j >= mx.n || i > j) {
        return;
    }

    while (i <= j) {
        const int current = getM(mx, i, j);
        if (j > 0 && current == getM(mx, i, j - 1)) {
            --j;
            continue;
        }
        if (i + 1 < mx.n && current == getM(mx, i + 1, j)) {
            ++i;
            continue;
        }
        break;
    }

    if (i >= mx.n || j >= mx.n || i > j) {
        return;
    }

    const int current = getM(mx, i, j);
    const int branchC = getC(mx, i, j);
    if (branchC < INF) {
        const bool hasLeftMismatch = i > 0;
        const bool hasRightMismatch = j + 1 < mx.n;
        const int type = pairType(seq[i], seq[j]);
        const int branch = branchC + multibranchStemEnergy(type,
                                                           hasLeftMismatch ? seq[i - 1] : 'N',
                                                           hasRightMismatch ? seq[j + 1] : 'N',
                                                           hasLeftMismatch,
                                                           hasRightMismatch);
        if (current == branch) {
            tracebackC(mx, seq, i, j, structure);
            return;
        }
    }

    for (size_t k = i + 1; k + 1 < j; ++k) {
        const int left = getM(mx, i, k);
        const int right = getM(mx, k + 1, j);
        if (left < INF && right < INF && current == left + right) {
            tracebackM(mx, seq, i, k, structure);
            tracebackM(mx, seq, k + 1, j, structure);
            return;
        }
    }
}

static void tracebackM2(const FoldMatrices &mx, const std::string &seq, size_t i, size_t j, std::string &structure) {
    if (i >= mx.n || j >= mx.n || i > j) {
        return;
    }
    const int current = getM2(mx, i, j);
    for (size_t k = i + 1; k + 1 < j; ++k) {
        const int left = getM(mx, i, k);
        const int right = getM(mx, k + 1, j);
        if (left < INF && right < INF && current == left + right) {
            tracebackM(mx, seq, i, k, structure);
            tracebackM(mx, seq, k + 1, j, structure);
            return;
        }
    }
}

static void tracebackC(const FoldMatrices &mx, const std::string &seq, size_t i, size_t j, std::string &structure) {
    if (i >= mx.n || j >= mx.n || i >= j) {
        return;
    }

    const int type = pairType(seq[i], seq[j]);
    const int current = getC(mx, i, j);
    if (type <= 0 || current >= INF) {
        return;
    }

    structure[i] = '(';
    structure[j] = ')';

    if (current == hairpinPenalty(seq, i, j, type)) {
        return;
    }

    const size_t maxLoop = 30;
    for (size_t p = i + 1; p <= std::min(j - 2, i + maxLoop + 1); ++p) {
        const size_t leftLoop = p - i - 1;
        const size_t minQForLoopBudget = (j > maxLoop - leftLoop + 1) ? (j - (maxLoop - leftLoop) - 1) : 0;
        const size_t minQ = std::max(p + 1, minQForLoopBudget);
        for (size_t q = j - 1; q >= minQ; --q) {
            const int inner = getC(mx, p, q);
            if (inner < INF) {
                const size_t rightLoop = j - q - 1;
                const int innerType = pairType(seq[q], seq[p]);
                const int energy = inner + internalLoopPenalty(seq, i, j, p, q, leftLoop, rightLoop, type, innerType);
                if (current == energy) {
                    tracebackC(mx, seq, p, q, structure);
                    return;
                }
            }
            if (q == 0) {
                break;
            }
        }
    }

    if (i + 1 < j - 1) {
        const int enclosed = getM2(mx, i + 1, j - 1);
        if (enclosed < INF) {
            const int multibranch = 930 + enclosed + multibranchStemEnergy(pairType(seq[j], seq[i]), seq[j - 1], seq[i + 1], true, true);
            if (current == multibranch) {
                tracebackM2(mx, seq, i + 1, j - 1, structure);
                return;
            }
        }
    }
}


static void computeExteriorPrefix(const std::string &seq, FoldMatrices &mx) {
    for (size_t j = 0; j < mx.n; ++j) {
        int best = (j == 0) ? 0 : mx.f5[j - 1];
        for (size_t i = 0; i <= j; ++i) {
            const int direct = getC(mx, i, j);
            if (direct >= INF) {
                continue;
            }
            const int prefix = (i == 0) ? 0 : mx.f5[i - 1];
            const bool hasLeftDangle = i > 0;
            const bool hasRightDangle = j + 1 < mx.n;
            const int type = pairType(seq[i], seq[j]);
            const int energy = prefix + direct + exteriorStemEnergy(type,
                                                                    hasLeftDangle ? seq[i - 1] : 'N',
                                                                    hasRightDangle ? seq[j + 1] : 'N',
                                                                    hasLeftDangle,
                                                                    hasRightDangle);
            if (energy < best) {
                best = energy;
            }
        }
        mx.f5[j] = best;
    }
}

static void tracebackExteriorPrefix(const FoldMatrices &mx, const std::string &seq, size_t j, std::string &structure) {
    if (j >= mx.n) {
        return;
    }

    while (j > 0 && mx.f5[j] == mx.f5[j - 1]) {
        --j;
    }
    if (j == 0) {
        return;
    }

    const int current = mx.f5[j];
    for (size_t i = j; i > 0; --i) {
        const int direct = getC(mx, i, j);
        if (direct >= INF) {
            continue;
        }
        const int prefix = mx.f5[i - 1];
        const bool hasRightDangle = j + 1 < mx.n;
        const int type = pairType(seq[i], seq[j]);
        const int energy = prefix + direct + exteriorStemEnergy(type,
                                                                seq[i - 1],
                                                                hasRightDangle ? seq[j + 1] : 'N',
                                                                true,
                                                                hasRightDangle);
        if (current == energy) {
            tracebackExteriorPrefix(mx, seq, i - 1, structure);
            tracebackC(mx, seq, i, j, structure);
            return;
        }
    }

    const int wholeC = getC(mx, 0, j);
    if (wholeC < INF) {
        const bool hasRightDangle = j + 1 < mx.n;
        const int type = pairType(seq[0], seq[j]);
        const int energy = wholeC + exteriorStemEnergy(type,
                                                       'N',
                                                       hasRightDangle ? seq[j + 1] : 'N',
                                                       false,
                                                       hasRightDangle);
        if (current == energy) {
            tracebackC(mx, seq, 0, j, structure);
            return;
        }
    }
}


static std::string predictDotBracketThermo(const std::string &seq) {
    const size_t n = seq.size();
    std::string structure(n, '.');
    if (n <= MIN_HAIRPIN_LOOP + 1) {
        return structure;
    }

    FoldMatrices mx(n);
    int *const cData = mx.c.data();
    int *const mRowData = mx.mRow.data();
    int *const mColData = mx.mCol.data();
    const size_t *const rowOffset = mx.rowOffset.data();
    const size_t *const colOffset = mx.colOffset.data();
    std::vector<std::vector<int> > m2BySpan(3, std::vector<int>(n, INF));

    for (size_t span = 1; span < n; ++span) {
        std::vector<int> &currentM2 = m2BySpan[span % 3];
        for (size_t i = 0; i + span < n; ++i) {
            const size_t j = i + span;
            const size_t rowBase = rowOffset[i];
            const size_t rowMinusI = rowBase - i;
            const int type = pairType(seq[i], seq[j]);

            int bestC = INF;
            if (type > 0 && j > i + MIN_HAIRPIN_LOOP) {
                bestC = hairpinPenalty(seq, i, j, type);

                if (i + 1 < j) {
                    const int innerType = pairType(seq[j - 1], seq[i + 1]);
                    const int inner = cData[rowOffset[i + 1] + (j - i - 2)];
                    if (innerType > 0 && inner < INF) {
                        const int energy = inner + stackEnergy(type, innerType);
                        if (energy < bestC) {
                            bestC = energy;
                        }
                    }
                }

                const size_t maxLoop = 30;
                for (size_t p = i + 1; p < j; ++p) {
                    const size_t leftLoop = p - i - 1;
                    if (leftLoop > maxLoop) {
                        break;
                    }
                    const size_t maxRightLoop = maxLoop - leftLoop;
                    const size_t minQForHairpin = p + MIN_HAIRPIN_LOOP + 1;
                    const size_t minQForLoopBudget = (j > maxRightLoop + 1) ? (j - maxRightLoop - 1) : 0;
                    const size_t qStart = std::max(minQForHairpin, minQForLoopBudget);
                    for (size_t q = qStart; q < j; ++q) {
                        const size_t rightLoop = j - q - 1;
                        if (leftLoop + rightLoop == 0) {
                            continue;
                        }
                        const int inner = cData[rowOffset[p] + (q - p)];
                        if (inner >= INF) {
                            continue;
                        }
                        const int innerType = pairType(seq[q], seq[p]);
                        const int energy = inner + internalLoopPenalty(seq, i, j, p, q, leftLoop, rightLoop, type, innerType);
                        if (energy < bestC) {
                            bestC = energy;
                        }
                    }
                }

                if (i + 1 < j - 1) {
                    const int enclosed = m2BySpan[(span - 2) % 3][i + 1];
                    if (enclosed < INF) {
                        const int multibranch = 930 + enclosed + multibranchStemEnergy(pairType(seq[j], seq[i]), seq[j - 1], seq[i + 1], true, true);
                        if (multibranch < bestC) {
                            bestC = multibranch;
                        }
                    }
                }
            }
            cData[rowMinusI + j] = bestC;

            int bestM2 = INF;
            size_t leftIndex = rowBase + 1;
            size_t rightIndex = colOffset[j] + i + 2;
            for (size_t k = i + 1; k + 1 < j; ++k, ++leftIndex, ++rightIndex) {
                const int energy = mRowData[leftIndex] + mColData[rightIndex];
                if (energy < bestM2) {
                    bestM2 = energy;
                }
            }
            if (bestM2 >= INF_CLAMP) {
                bestM2 = INF;
            }
            currentM2[i] = bestM2;

            int bestM = bestM2;
            if (i + 1 <= j) {
                const int unpairedM = mRowData[rowOffset[i + 1] + (j - i - 1)];
                if (unpairedM < bestM) {
                    bestM = unpairedM;
                }
            }
            if (i <= j - 1) {
                const int rightUnpairedM = mRowData[rowMinusI + j - 1];
                if (rightUnpairedM < bestM) {
                    bestM = rightUnpairedM;
                }
            }
            if (bestC < INF) {
                const bool hasLeftMismatch = i > 0;
                const bool hasRightMismatch = j + 1 < n;
                const int branch = bestC + multibranchStemEnergy(type,
                                                                 hasLeftMismatch ? seq[i - 1] : 'N',
                                                                 hasRightMismatch ? seq[j + 1] : 'N',
                                                                 hasLeftMismatch,
                                                                 hasRightMismatch);
                if (branch < bestM) {
                    bestM = branch;
                }
            }
            mRowData[rowMinusI + j] = bestM;
            mColData[colOffset[j] + i] = bestM;

        }
    }

    computeExteriorPrefix(seq, mx);


    tracebackExteriorPrefix(mx, seq, n - 1, structure);
    return structure;
}

}

namespace RiboseekRnaFold {

bool predictDotBracket(const std::string &rna, std::string &dotBracket, double *scoreOut) {
    dotBracket.clear();
    if (scoreOut != 0) {
        *scoreOut = 0.0;
    }
    if (rna.empty()) {
        return false;
    }
    const std::string seq = normalizeRnaSequence(rna.data(), rna.size());
    dotBracket = predictDotBracketThermo(seq);
    return dotBracket.size() == seq.size();
}

}

int rnafold(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    queryReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter structureWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed,
                             Parameters::DBTYPE_GENERIC_DB);
    structureWriter.open();

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min(static_cast<size_t>(par.threads), queryReader.getSize()), static_cast<size_t>(1));
#endif

    Debug::Progress progress(queryReader.getSize());
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int threadIdx = 0;
#ifdef OPENMP
        threadIdx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < queryReader.getSize(); ++id) {
            progress.updateProgress();
            const char *data = queryReader.getData(id, threadIdx);
            const size_t seqLen = queryReader.getSeqLen(id);
            std::string dotBracket;
            RiboseekRnaFold::predictDotBracket(std::string(data, seqLen), dotBracket, 0);
            structureWriter.writeData(dotBracket.c_str(), dotBracket.size(), queryReader.getDbKey(id), threadIdx);
        }
    }
    Debug(Debug::INFO) << "\n";

    structureWriter.close();
    queryReader.close();
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, (DBFiles::Files)(DBFiles::LOOKUP | DBFiles::SOURCE));
    return EXIT_SUCCESS;
}
