#ifndef RIBOSEEK_RNAREALIGN_H
#define RIBOSEEK_RNAREALIGN_H

#include <string>

template <typename T> class DBReader;

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

// Align one query against targetWindow. Coordinates come back relative to the window and to the
// forward query. Returns false when no alignment was produced. targetDb is the ORIGINAL (unsplit)
// target database: its GC statistic gives the realigned hit the same per-query E-value model as
// every other hit.
//
// queryDbType >= 0: queryData is entry queryEntryKey of the iteration's strand-split
// query DB (forward copy at 2k, reverse at 2k+1) of that dbtype -- the sequence in the first
// iteration, a profile in later ones -- and is mapped exactly as rnaalign maps it, so alignment
// and E-value are rnaalign's. queryDbType < 0: legacy, queryData is the raw forward sequence and
// the reverse strand is built here.
bool rnaRealignWindow(const char *queryData, unsigned int queryLen, int queryDbType,
                      unsigned int queryEntryKey, bool reverse,
                      const char *targetWindow, unsigned int windowLen,
                      unsigned int queryKey, unsigned int targetKey,
                      DBReader<unsigned int> *targetDb, int threadIdx,
                      RnaRealignHit &out, std::string &backtrace);

// Largest window this aligner will accept.
unsigned int rnaRealignMaxWindow();

#endif
