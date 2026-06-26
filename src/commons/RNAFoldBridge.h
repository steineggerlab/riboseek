#ifndef MMSEQS_RNAFOLD_BRIDGE_H
#define MMSEQS_RNAFOLD_BRIDGE_H

#include "rnautils/RnaFold.h"

#include <cctype>
#include <string>
#include <vector>

static inline bool rnaFoldPredictDotBracket(const std::string &rna,
                                            std::string &dotBracket,
                                            double *scoreOut) {
    return RiboseekRnaFold::predictDotBracket(rna, dotBracket, scoreOut);
}

// Consensus dot-bracket via single-seq internal fold on the QUERY sequence
// (rows[0], ungapped), with the predicted structure mapped back into
// alignment-column space (gap columns receive '.'). Mirrors the reference
// workflow: reformat MSA -> fold(query) -> splice as #=GC SS_cons ->
// cmbuild --hand. Comparative/alifold folding is intentionally not used.
static inline bool rnaFoldAlifoldDotBracket(const std::vector<std::string> &rows,
                                            std::string &dotBracket,
                                            double *scoreOut) {
    dotBracket.clear();
    if (scoreOut != 0) {
        *scoreOut = 0.0;
    }
    if (rows.empty()) {
        return false;
    }

    const size_t alnLen = rows[0].size();
    if (alnLen == 0) {
        return false;
    }

    std::string ungapped;
    std::vector<bool> keep(alnLen, false);
    ungapped.reserve(alnLen);
    for (size_t col = 0; col < alnLen; ++col) {
        char c = static_cast<char>(std::toupper(static_cast<unsigned char>(rows[0][col])));
        if (c == '-' || c == '.') {
            continue;
        }
        if (c == 'T') {
            c = 'U';
        }
        keep[col] = true;
        ungapped.push_back(c);
    }
    if (ungapped.empty()) {
        return false;
    }

    std::string compactSs;
    if (!rnaFoldPredictDotBracket(ungapped, compactSs, scoreOut)) {
        return false;
    }
    if (compactSs.size() != ungapped.size()) {
        return false;
    }

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
