#ifndef RIBOSEEK_RNAREALIGN_H
#define RIBOSEEK_RNAREALIGN_H

#include <string>

// Thin POD boundary around riboseek's dinucleotide Smith-Waterman.
// RnaSmithWaterman.h and mmseqs' SmithWaterman.h both define s_align/aln_t, so a
// translation unit cannot include Matcher.h and RnaMatcher.h at the same time. Callers on
// the mmseqs side (FragmentMerge.cpp) talk to the aligner through this interface instead.

struct RnaRealignHit {
    int qStartPos;
    int qEndPos;
    int dbStartPos;
    int dbEndPos;
    unsigned int alnLength;
    int score;
    float qcov;
    float dbcov;
    float seqId;
    double eval;
};

// Align querySeq (raw nucleotide characters, reverse-complemented internally when
// reverse == true) against targetWindow. Coordinates come back relative to the window and
// to the forward query. Returns false when no alignment was produced.
bool rnaRealignWindow(const char *querySeq, unsigned int queryLen, bool reverse,
                      const char *targetWindow, unsigned int windowLen,
                      unsigned int queryKey, unsigned int targetKey,
                      size_t targetTotalResidues, int threadIdx,
                      RnaRealignHit &out, std::string &backtrace);

// Largest window this aligner will accept.
unsigned int rnaRealignMaxWindow();

#endif
