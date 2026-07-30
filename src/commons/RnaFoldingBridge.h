#ifndef MMSEQS_RNA_FOLDING_BRIDGE_H
#define MMSEQS_RNA_FOLDING_BRIDGE_H

#include <string>
#include <vector>

#include <cctype>
#include <tornadofold.h>

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

    tornadofold::TornadoFold folder;
    int mfe = folder.fold(rna);
    dotBracket = folder.traceback(mfe);
    if (dotBracket.size() != rna.size()) {
        dotBracket.clear();
        return false;
    }
    if (scoreOut != nullptr) {
        *scoreOut = static_cast<double>(mfe) / 100.0;
    }
    return true;
}

// Consensus dot-bracket via single-seq folding of the QUERY sequence (rows[0],
// ungapped), with the predicted structure mapped back into alignment-column
// space (gap columns receive '.')
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

#endif
