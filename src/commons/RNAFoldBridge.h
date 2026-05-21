#ifndef MMSEQS_RNAFOLD_BRIDGE_H
#define MMSEQS_RNAFOLD_BRIDGE_H

#include <string>
#include <vector>

#ifdef HAVE_RNAFOLD_PREDICT

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <algorithm>

extern "C" {
#include <ViennaRNA/mfe/global.h>
#include <ViennaRNA/fold.h>
}

static inline bool rnaFoldPredictDotBracket(const std::string &rna,
                                            std::string &dotBracket,
                                            double *scoreOut) {
    dotBracket.clear();
    if (scoreOut != nullptr) {
        *scoreOut = 0.0;
    }
    if (rna.empty()) {
        return false;
    }

    // Prepare sequence: uppercase, T->U
    std::string seq(rna);
    for (size_t i = 0; i < seq.size(); ++i) {
        seq[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[i])));
        if (seq[i] == 'T') {
            seq[i] = 'U';
        }
    }

    // Allocate structure string
    char *structure = (char *)calloc(seq.size() + 1, sizeof(char));
    if (structure == nullptr) {
        return false;
    }

    // Run MFE prediction
    float mfe = vrna_fold(seq.c_str(), structure);

    dotBracket.assign(structure, seq.size());
    if (scoreOut != nullptr) {
        *scoreOut = static_cast<double>(mfe);
    }

    free(structure);
    return true;
}

// Consensus dot-bracket via single-seq vrna_fold on the QUERY sequence
// (rows[0], ungapped), with the predicted structure mapped back into
// alignment-column space (gap columns receive '.'). Mirrors the reference
// workflow: reformat MSA → RNAfold(query) → splice as #=GC SS_cons →
// cmbuild --hand. vrna_alifold is intentionally not used.
static inline bool rnaFoldAlifoldDotBracket(const std::vector<std::string> &rows,
                                            std::string &dotBracket,
                                            double *scoreOut) {
    dotBracket.clear();
    if (scoreOut != nullptr) *scoreOut = 0.0;
    if (rows.empty()) return false;

    const size_t alnLen = rows[0].size();
    if (alnLen == 0) return false;

    std::string ungapped;
    std::vector<bool> keep(alnLen, false);
    ungapped.reserve(alnLen);
    for (size_t col = 0; col < alnLen; ++col) {
        char c = static_cast<char>(std::toupper(static_cast<unsigned char>(rows[0][col])));
        if (c == '-' || c == '.') continue;
        if (c == 'T') c = 'U';
        keep[col] = true;
        ungapped.push_back(c);
    }
    if (ungapped.empty()) return false;

    std::string compactSs;
    if (!rnaFoldPredictDotBracket(ungapped, compactSs, scoreOut)) return false;
    if (compactSs.size() != ungapped.size()) return false;

    dotBracket.assign(alnLen, '.');
    size_t j = 0;
    for (size_t col = 0; col < alnLen; ++col) {
        if (keep[col]) {
            dotBracket[col] = compactSs[j++];
        }
    }
    return true;
}

#else

static inline bool rnaFoldPredictDotBracket(const std::string & /*rna*/,
                                            std::string &dotBracket,
                                            double *scoreOut) {
    dotBracket.clear();
    if (scoreOut != nullptr) *scoreOut = 0.0;
    return false;
}

static inline bool rnaFoldAlifoldDotBracket(const std::vector<std::string> & /*rows*/,
                                            std::string &dotBracket,
                                            double *scoreOut) {
    dotBracket.clear();
    if (scoreOut != nullptr) *scoreOut = 0.0;
    return false;
}

#endif

bool rnaLinearPartitionPredict(const std::string &rna,
                               int minLoop,
                               int beamSize,
                               std::string &dotBracket,
                               std::vector<float> &pairMatrix);

#endif
