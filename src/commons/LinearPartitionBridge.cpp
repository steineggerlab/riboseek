#include "RNAFoldBridge.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sys/time.h>

#if __has_include("../../lib/LinearTurboFold/src/LinearPartition/src/LinearPartition.h")

class LinearTurboFold {
public:
    double get_folding_extrinsic_information(int, int, int) {
        return 0.0;
    }
};

#define LINEAR_TURBOFOLD_H
#define RNA_FOLD_H

#define base_pair ltf_lp_base_pair
#define internal_1x1_nucleotides ltf_lp_internal_1x1_nucleotides
#define helix_stacking ltf_lp_helix_stacking
#define terminal_mismatch ltf_lp_terminal_mismatch
#define bulge_0x1_nucleotides ltf_lp_bulge_0x1_nucleotides
#define helix_closing ltf_lp_helix_closing
#define dangle_left ltf_lp_dangle_left
#define dangle_right ltf_lp_dangle_right
#define internal_explicit ltf_lp_internal_explicit
#define hairpin_length ltf_lp_hairpin_length
#define bulge_length ltf_lp_bulge_length
#define internal_length ltf_lp_internal_length
#define internal_symmetric_length ltf_lp_internal_symmetric_length
#define internal_asymmetry ltf_lp_internal_asymmetry
#define hairpin_length_at_least ltf_lp_hairpin_length_at_least
#define bulge_length_at_least ltf_lp_bulge_length_at_least
#define internal_length_at_least ltf_lp_internal_length_at_least
#define internal_symmetric_length_at_least ltf_lp_internal_symmetric_length_at_least
#define internal_asymmetry_at_least ltf_lp_internal_asymmetry_at_least
#define multi_base ltf_lp_multi_base
#define multi_unpaired ltf_lp_multi_unpaired
#define multi_paired ltf_lp_multi_paired
#define external_unpaired ltf_lp_external_unpaired
#define external_paired ltf_lp_external_paired

#define lxc37 ltf_lp_lxc37
#define ML_intern37 ltf_lp_ML_intern37
#define ML_closing37 ltf_lp_ML_closing37
#define ML_BASE37 ltf_lp_ML_BASE37
#define MAX_NINIO ltf_lp_MAX_NINIO
#define ninio37 ltf_lp_ninio37
#define TerminalAU37 ltf_lp_TerminalAU37
#define Triloops ltf_lp_Triloops
#define Triloop37 ltf_lp_Triloop37
#define Tetraloops ltf_lp_Tetraloops
#define Tetraloop37 ltf_lp_Tetraloop37
#define Hexaloops ltf_lp_Hexaloops
#define Hexaloop37 ltf_lp_Hexaloop37
#define stack37 ltf_lp_stack37
#define hairpin37 ltf_lp_hairpin37
#define bulge37 ltf_lp_bulge37
#define internal_loop37 ltf_lp_internal_loop37
#define mismatchI37 ltf_lp_mismatchI37
#define mismatchH37 ltf_lp_mismatchH37
#define mismatchM37 ltf_lp_mismatchM37
#define mismatch1nI37 ltf_lp_mismatch1nI37
#define mismatch23I37 ltf_lp_mismatch23I37
#define mismatchExt37 ltf_lp_mismatchExt37
#define dangle5_37 ltf_lp_dangle5_37
#define dangle3_37 ltf_lp_dangle3_37
#define int11_37 ltf_lp_int11_37
#define int21_37 ltf_lp_int21_37
#define int22_37 ltf_lp_int22_37

#include "../../lib/LinearTurboFold/src/LinearPartition/src/LinearPartition.h"
#include "../../lib/LinearTurboFold/src/LinearPartition/src/LinearPartition.cpp"

#undef base_pair
#undef internal_1x1_nucleotides
#undef helix_stacking
#undef terminal_mismatch
#undef bulge_0x1_nucleotides
#undef helix_closing
#undef dangle_left
#undef dangle_right
#undef internal_explicit
#undef hairpin_length
#undef bulge_length
#undef internal_length
#undef internal_symmetric_length
#undef internal_asymmetry
#undef hairpin_length_at_least
#undef bulge_length_at_least
#undef internal_length_at_least
#undef internal_symmetric_length_at_least
#undef internal_asymmetry_at_least
#undef multi_base
#undef multi_unpaired
#undef multi_paired
#undef external_unpaired
#undef external_paired
#undef lxc37
#undef ML_intern37
#undef ML_closing37
#undef ML_BASE37
#undef MAX_NINIO
#undef ninio37
#undef TerminalAU37
#undef Triloops
#undef Triloop37
#undef Tetraloops
#undef Tetraloop37
#undef Hexaloops
#undef Hexaloop37
#undef stack37
#undef hairpin37
#undef bulge37
#undef internal_loop37
#undef mismatchI37
#undef mismatchH37
#undef mismatchM37
#undef mismatch1nI37
#undef mismatch23I37
#undef mismatchExt37
#undef dangle5_37
#undef dangle3_37
#undef int11_37
#undef int21_37
#undef int22_37

namespace {

struct PairCandidate {
    int i = -1;
    int j = -1;
    float p = 0.0f;
};

static inline double linearGetBasePairProb(int i,
                                           int j,
                                           std::unordered_map<int, State> *pfscore,
                                           std::unordered_map<int, ExtValue> *extInfo,
                                           double viterbiScore) {
    if (i <= 0 || j <= 0 || i >= j) {
        return 0.0;
    }
    if (pfscore[j - 1].find(i - 1) == pfscore[j - 1].end()) {
        return 0.0;
    }
    if (extInfo[j][i].score < TO_XLOG(EPSILON)) {
        return 0.0;
    }

    double prob = DIV(PROD(pfscore[j - 1][i - 1].alpha, pfscore[j - 1][i - 1].beta),
                      PROD(viterbiScore POWSCALING2(2), extInfo[j][i].score));
    if (prob > float(-9.91152)) {
        prob = TO_LINEAR(prob);
        if (prob > 1.0) {
            prob = 1.0;
        }
        return prob;
    }
    return 0.0;
}

std::string normalizeRnaSequence(const std::string &rna) {
    std::string normalized(rna);
    for (size_t i = 0; i < normalized.size(); ++i) {
        char c = static_cast<char>(std::toupper(static_cast<unsigned char>(normalized[i])));
        normalized[i] = (c == 'T') ? 'U' : c;
    }
    return normalized;
}

bool crossesAnyAccepted(int i, int j, const std::vector<std::pair<int, int> > &accepted) {
    for (size_t idx = 0; idx < accepted.size(); ++idx) {
        const int a = accepted[idx].first;
        const int b = accepted[idx].second;
        if ((i < a && a < j && j < b) || (a < i && i < b && b < j)) {
            return true;
        }
    }
    return false;
}

void buildDotBracketFromProbabilities(const std::vector<float> &pairMatrix,
                                      int len,
                                      int minLoop,
                                      std::string &dotBracket) {
    dotBracket.assign(static_cast<size_t>(len), '.');
    if (len <= 0) {
        return;
    }

    std::vector<float> rowMax(static_cast<size_t>(len), 0.0f);
    for (int i = 0; i < len; ++i) {
        for (int j = i + minLoop + 1; j < len; ++j) {
            const float p = pairMatrix[static_cast<size_t>(i * len + j)];
            rowMax[static_cast<size_t>(i)] = std::max(rowMax[static_cast<size_t>(i)], p);
            rowMax[static_cast<size_t>(j)] = std::max(rowMax[static_cast<size_t>(j)], p);
        }
    }

    std::vector<PairCandidate> candidates;
    candidates.reserve(static_cast<size_t>(len * 4));
    const float eps = 1e-6f;
    for (int i = 0; i < len; ++i) {
        for (int j = i + minLoop + 1; j < len; ++j) {
            const float p = pairMatrix[static_cast<size_t>(i * len + j)];
            if (p <= 0.0f) {
                continue;
            }
            if ((p + eps) < rowMax[static_cast<size_t>(i)] || (p + eps) < rowMax[static_cast<size_t>(j)]) {
                continue;
            }
            PairCandidate cand;
            cand.i = i;
            cand.j = j;
            cand.p = p;
            candidates.push_back(cand);
        }
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const PairCandidate &lhs, const PairCandidate &rhs) {
                  if (lhs.p != rhs.p) {
                      return lhs.p > rhs.p;
                  }
                  if (lhs.i != rhs.i) {
                      return lhs.i < rhs.i;
                  }
                  return lhs.j < rhs.j;
              });

    std::vector<int> partner(static_cast<size_t>(len), -1);
    std::vector<std::pair<int, int> > accepted;
    accepted.reserve(static_cast<size_t>(len / 2));
    for (size_t idx = 0; idx < candidates.size(); ++idx) {
        const PairCandidate &cand = candidates[idx];
        if (partner[static_cast<size_t>(cand.i)] != -1 || partner[static_cast<size_t>(cand.j)] != -1) {
            continue;
        }
        if (crossesAnyAccepted(cand.i, cand.j, accepted)) {
            continue;
        }
        partner[static_cast<size_t>(cand.i)] = cand.j;
        partner[static_cast<size_t>(cand.j)] = cand.i;
        accepted.push_back(std::make_pair(cand.i, cand.j));
    }

    for (int i = 0; i < len; ++i) {
        const int j = partner[static_cast<size_t>(i)];
        if (j > i) {
            dotBracket[static_cast<size_t>(i)] = '(';
            dotBracket[static_cast<size_t>(j)] = ')';
        }
    }
}

} // namespace

bool rnaLinearPartitionPredict(const std::string &rna,
                               int minLoop,
                               int beamSize,
                               std::string &dotBracket,
                               std::vector<float> &pairMatrix) {
    dotBracket.clear();
    pairMatrix.clear();
    if (rna.empty()) {
        return false;
    }

    const std::string normalized = normalizeRnaSequence(rna);
    const int len = static_cast<int>(normalized.size());
    pairMatrix.assign(static_cast<size_t>(len * len), 0.0f);

    std::unordered_map<int, State> *pfscore = nullptr;
    std::unordered_map<int, ExtValue> *extInfo = nullptr;
    try {
        pfscore = new std::unordered_map<int, State>[static_cast<size_t>(len)];
        extInfo = new std::unordered_map<int, ExtValue>[static_cast<size_t>(len + 1)];

        BeamCKYParser parser(std::max(1, beamSize), true, false);
        std::string parserSeq(normalized);
        const double viterbiScore = parser.parse(0, 0, nullptr, parserSeq, pfscore, extInfo, "", "");

        for (int i = 1; i < len; ++i) {
            for (int j = i + minLoop + 1; j <= len; ++j) {
                const float p = static_cast<float>(linearGetBasePairProb(i, j, pfscore, extInfo, viterbiScore));
                pairMatrix[static_cast<size_t>((i - 1) * len + (j - 1))] = p;
                pairMatrix[static_cast<size_t>((j - 1) * len + (i - 1))] = p;
            }
        }

        buildDotBracketFromProbabilities(pairMatrix, len, minLoop, dotBracket);
        delete[] extInfo;
        delete[] pfscore;
        return dotBracket.size() == normalized.size();
    } catch (...) {
        delete[] extInfo;
        delete[] pfscore;
        dotBracket.assign(normalized.size(), '.');
        pairMatrix.assign(static_cast<size_t>(len * len), 0.0f);
        return false;
    }
}

#else

bool rnaLinearPartitionPredict(const std::string &,
                               int,
                               int,
                               std::string &dotBracket,
                               std::vector<float> &pairMatrix) {
    dotBracket.clear();
    pairMatrix.clear();
    return false;
}

#endif
