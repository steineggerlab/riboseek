#ifndef MMSEQS_LINEARFOLD_BRIDGE_H
#define MMSEQS_LINEARFOLD_BRIDGE_H

#include <string>

#ifdef HAVE_LINEARFOLD_PREDICT
bool linearFoldPredictDotBracket(const std::string &rna,
                                 std::string &dotBracket,
                                 double *scoreOut);
#else
static inline bool linearFoldPredictDotBracket(const std::string & /*rna*/,
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
