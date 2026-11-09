#ifndef PAIRSW_H
#define PAIRSW_H

// Batched pairwise local alignment scores (Smith-Waterman, affine gaps) for sparse (query, target) pairs, e.g. the
// k-mer pairs of kmermatcher: one group of threads per pair, many pairs per launch.  Scores only (int32, no
// saturation), with the query profile of MMseqs2's SmithWaterman (SUBSTITUTIONMATRIX mode):
//   s(i, j) = subMat[target_j][query_i] + bias_i,
// gaps cost gapOpen for the first position and gapExtend for each further one (H - gapOpen, E - gapExtend), as
// StripedSmithWaterman.  Any alphabet up to 32 letters (riboseek dinucleotides: 25).

#include <cstddef>
#include <cstdint>

class PairSW {
public:
    struct Pair {
        uint64_t qOff, tOff;  // offsets of the query / target letters in the batch buffer
        uint32_t qLen, tLen;
    };

    // subMat: alphabetSize x alphabetSize, row = target letter, column = query letter.  letterMap (256 entries, or
    // NULL = letters are already encoded): byte -> letter, applied on the device to the raw sequence bytes.
    // maxLetters / maxPairs bound one batch; device buffers are allocated once here.
    PairSW(const int8_t *subMat, int alphabetSize, int gapOpen, int gapExtend, size_t maxLetters, size_t maxPairs,
           const uint8_t *letterMap = nullptr);
    ~PairSW();

    // letters: sequences (raw bytes when a letterMap is set); bias: per-letter composition bias (read at query
    // offsets only).  Encoded letters must be < alphabetSize.  Not thread-safe: one call at a time.  scores[k] = best local score of pairs[k].
    void score(const uint8_t *letters, const int8_t *bias, size_t nLetters, const Pair *pairs, size_t nPairs,
               int32_t *scores);

    static size_t deviceBytes(size_t maxLetters, size_t maxPairs);
    static int deviceCount();

private:
    int alphabetSize, gapOpen, gapExtend;
    size_t maxLetters, maxPairs;
    void *dLetters, *dBias, *dPairs, *dScores, *dOrder, *dCounter, *dTemp, *dMat, *dMap;
    size_t tempBytes;
    void *stream;
};

#endif
