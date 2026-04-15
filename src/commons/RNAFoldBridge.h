#ifndef MMSEQS_RNAFOLD_BRIDGE_H
#define MMSEQS_RNAFOLD_BRIDGE_H

#include <string>

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

#else

static inline bool rnaFoldPredictDotBracket(const std::string & /*rna*/,
                                            std::string &dotBracket,
                                            double *scoreOut) {
    dotBracket.clear();
    if (scoreOut != nullptr) {
        *scoreOut = 0.0;
    }
    return false;
}

#endif

#endif
