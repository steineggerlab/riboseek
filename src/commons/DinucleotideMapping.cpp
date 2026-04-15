#include "DinucleotideMapping.h"
#include "Sequence.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "MathUtil.h"
#include "Util.h"

#include "SubstitutionMatrix.h"

#include <algorithm>
#include <climits>
#include <cctype>
#include <cstring>
#include <vector>
#include <cstdio>
#include <mutex>

static void dinucDbgDump(const char *tag, unsigned int dbKey, const unsigned char *num, unsigned int L) {
    if (dbKey > 2) return;
    static std::mutex m;
    std::lock_guard<std::mutex> lk(m);
    FILE *f = fopen("/tmp/riboseek_dbg.txt", "a");
    if (!f) return;
    fprintf(f, "%s key=%u L=%u |", tag, dbKey, L);
    unsigned int n = L < 64 ? L : 64;
    for (unsigned int i = 0; i < n; i++) fprintf(f, " %02x", num[i]);
    fprintf(f, "\n");
    fclose(f);
}

// Static lookup tables populated by dinucSetupMatrix
static unsigned char dinucHead[256];      // dinuc code -> XN sentinel (first nuc as head)
static unsigned char dinucTail[256];      // dinuc code -> NX sentinel (second nuc as tail)
static unsigned char dinucToNucTbl[256];  // dinuc code -> second nucleotide (numeric)
static bool dinucTablesReady = false;
static BaseMatrix *dinucAlignSubMat = NULL;  // bitFactor=2.0 for profile_for_alignment

// Local calcLocalAaBiasCorrection matching reference (windowSize=50, reverse support)
static void dinucCalcLocalAaBiasCorrection(const BaseMatrix *m,
                                            unsigned char *int_sequence,
                                            const int N,
                                            float *compositionBias,
                                            float scale,
                                            bool reverse) {
    if (reverse) {
        std::reverse(int_sequence, int_sequence + N);
    }
    const int windowSize = 50;
    for (int i = 0; i < N; i++) {
        const int minPos = std::max(0, (i - windowSize / 2));
        const int maxPos = std::min(N, (i + (windowSize+1) / 2));
        const int windowLength = maxPos - minPos;

        int sumSubScores = 0;
        short *subMat = m->subMatrix[int_sequence[i]];
        for (int j = minPos; j < maxPos; j++) {
            sumSubScores += subMat[int_sequence[j]];
        }

        sumSubScores -= subMat[int_sequence[i]];
        float deltaS_i = (float) sumSubScores;
        deltaS_i /= -1.0 * static_cast<float>(windowLength);
        for (int a = 0; a < m->alphabetSize; a++) {
            deltaS_i += m->pBack[a] * static_cast<float>(subMat[a]);
        }
        compositionBias[i] = scale * deltaS_i;
    }
    if (reverse) {
        std::reverse(compositionBias, compositionBias + N);
        std::reverse(int_sequence, int_sequence + N);
    }
}

static void dinucSetupMatrix(BaseMatrix *mat) {
    unsigned char *aa2num = mat->aa2num;

    // Restore single-char aa2num entries from num2aa (setupLetterMapping corrupts B/J/O/U)
    for (int i = 0; i < mat->alphabetSize; i++) {
        aa2num[(unsigned char)mat->num2aa[i]] = i;
        aa2num[(unsigned char)tolower(mat->num2aa[i])] = i;
    }

    // Populate the 16-bit dinucleotide pair entries in aa2num (which is USHRT_MAX entries).
    // For each pair (char1, char2), aa2num[(char1 << 8) | char2] = dinuc amino acid index.
    // Fill all 16-bit entries with X first, then overwrite known pairs.
    for (int i = 256; i < USHRT_MAX; i++) {
        aa2num[i] = aa2num[(int)'X'];
    }
    for (int i = 1; i < UCHAR_MAX; i++) {
        unsigned short upperLetter1 = (static_cast<unsigned short>(toupper(static_cast<unsigned char>(i))) << 8);
        unsigned short first = (static_cast<unsigned short>(i) << 8);
        for (int j = 0; j < UCHAR_MAX; j++) {
            unsigned short dinucleotide = first | static_cast<unsigned char>(j);
            unsigned char upperLetter2 = toupper(static_cast<unsigned char>(j));
            unsigned short upperLetters = upperLetter1 | upperLetter2;
            switch(upperLetters) {
                case 0b0100000101000001: aa2num[dinucleotide] = aa2num[(int)'C']; break; // AA
                case 0b0100000101000011: aa2num[dinucleotide] = aa2num[(int)'G']; break; // AC
                case 0b0100000101000111: aa2num[dinucleotide] = aa2num[(int)'L']; break; // AG
                case 0b0100000101010100: aa2num[dinucleotide] = aa2num[(int)'Q']; break; // AT
                case 0b0100000101010101: aa2num[dinucleotide] = aa2num[(int)'Q']; break; // AU
                case 0b0100001101000001: aa2num[dinucleotide] = aa2num[(int)'D']; break; // CA
                case 0b0100001101000011: aa2num[dinucleotide] = aa2num[(int)'F']; break; // CC
                case 0b0100001101000111: aa2num[dinucleotide] = aa2num[(int)'R']; break; // CG
                case 0b0100001101010100: aa2num[dinucleotide] = aa2num[(int)'K']; break; // CT
                case 0b0100001101010101: aa2num[dinucleotide] = aa2num[(int)'K']; break; // CU
                case 0b0100011101000001: aa2num[dinucleotide] = aa2num[(int)'M']; break; // GA
                case 0b0100011101000011: aa2num[dinucleotide] = aa2num[(int)'A']; break; // GC
                case 0b0100011101000111: aa2num[dinucleotide] = aa2num[(int)'P']; break; // GG
                case 0b0100011101010100: aa2num[dinucleotide] = aa2num[(int)'I']; break; // GT
                case 0b0100011101010101: aa2num[dinucleotide] = aa2num[(int)'I']; break; // GU
                case 0b0101010001000001: aa2num[dinucleotide] = aa2num[(int)'E']; break; // TA
                case 0b0101010001000011: aa2num[dinucleotide] = aa2num[(int)'N']; break; // TC
                case 0b0101010001000111: aa2num[dinucleotide] = aa2num[(int)'H']; break; // TG
                case 0b0101010001010100: aa2num[dinucleotide] = aa2num[(int)'S']; break; // TT
                case 0b0101010001010101: aa2num[dinucleotide] = aa2num[(int)'S']; break; // TU
                case 0b0101010101000001: aa2num[dinucleotide] = aa2num[(int)'E']; break; // UA
                case 0b0101010101000011: aa2num[dinucleotide] = aa2num[(int)'N']; break; // UC
                case 0b0101010101000111: aa2num[dinucleotide] = aa2num[(int)'H']; break; // UG
                case 0b0101010101010100: aa2num[dinucleotide] = aa2num[(int)'S']; break; // UT
                case 0b0101010101010101: aa2num[dinucleotide] = aa2num[(int)'S']; break; // UU
                default: {
                    unsigned char up1 = static_cast<char>(upperLetter1 >> 8);
                    unsigned char up2 = upperLetter2;
                    bool firstCanonical  = (up1 == 'A' || up1 == 'C' || up1 == 'G' || up1 == 'T' || up1 == 'U');
                    bool secondCanonical = (up2 == 'A' || up2 == 'C' || up2 == 'G' || up2 == 'T' || up2 == 'U');
                    if (firstCanonical && !secondCanonical) {
                        switch (up1) {
                            case 'A': aa2num[dinucleotide] = aa2num[(int)'T']; break;
                            case 'C': aa2num[dinucleotide] = aa2num[(int)'V']; break;
                            case 'G': aa2num[dinucleotide] = aa2num[(int)'W']; break;
                            case 'T': case 'U': aa2num[dinucleotide] = aa2num[(int)'Y']; break;
                            default: aa2num[dinucleotide] = aa2num[(int)'X']; break;
                        }
                    } else if (!firstCanonical && secondCanonical) {
                        switch (up2) {
                            case 'A': aa2num[dinucleotide] = aa2num[(int)'B']; break;
                            case 'C': aa2num[dinucleotide] = aa2num[(int)'J']; break;
                            case 'G': aa2num[dinucleotide] = aa2num[(int)'O']; break;
                            case 'T': case 'U': aa2num[dinucleotide] = aa2num[(int)'U']; break;
                            default: aa2num[dinucleotide] = aa2num[(int)'X']; break;
                        }
                    } else {
                        aa2num[dinucleotide] = aa2num[(int)'X'];
                    }
                    break;
                }
            }
        }
    }

    // num2revcompnum — 25-entry dinucleotide reverse complement
    unsigned char *rc = mat->num2revcompnum;
    rc[aa2num[(int)'C']] = aa2num[(int)'S'];
    rc[aa2num[(int)'G']] = aa2num[(int)'I'];
    rc[aa2num[(int)'L']] = aa2num[(int)'K'];
    rc[aa2num[(int)'Q']] = aa2num[(int)'Q'];
    rc[aa2num[(int)'D']] = aa2num[(int)'H'];
    rc[aa2num[(int)'F']] = aa2num[(int)'P'];
    rc[aa2num[(int)'R']] = aa2num[(int)'R'];
    rc[aa2num[(int)'K']] = aa2num[(int)'L'];
    rc[aa2num[(int)'M']] = aa2num[(int)'N'];
    rc[aa2num[(int)'A']] = aa2num[(int)'A'];
    rc[aa2num[(int)'P']] = aa2num[(int)'F'];
    rc[aa2num[(int)'I']] = aa2num[(int)'G'];
    rc[aa2num[(int)'E']] = aa2num[(int)'E'];
    rc[aa2num[(int)'N']] = aa2num[(int)'M'];
    rc[aa2num[(int)'H']] = aa2num[(int)'D'];
    rc[aa2num[(int)'S']] = aa2num[(int)'C'];
    rc[aa2num[(int)'T']] = aa2num[(int)'U'];
    rc[aa2num[(int)'V']] = aa2num[(int)'O'];
    rc[aa2num[(int)'W']] = aa2num[(int)'J'];
    rc[aa2num[(int)'Y']] = aa2num[(int)'B'];
    rc[aa2num[(int)'B']] = aa2num[(int)'Y'];
    rc[aa2num[(int)'J']] = aa2num[(int)'W'];
    rc[aa2num[(int)'O']] = aa2num[(int)'V'];
    rc[aa2num[(int)'U']] = aa2num[(int)'T'];
    rc[aa2num[(int)'X']] = aa2num[(int)'X'];

    // head — dinuc code -> XN sentinel (first nucleotide as head dinucleotide)
    dinucHead[aa2num[(int)'C']] = aa2num[(int)'B'];
    dinucHead[aa2num[(int)'G']] = aa2num[(int)'B'];
    dinucHead[aa2num[(int)'L']] = aa2num[(int)'B'];
    dinucHead[aa2num[(int)'Q']] = aa2num[(int)'B'];
    dinucHead[aa2num[(int)'D']] = aa2num[(int)'J'];
    dinucHead[aa2num[(int)'F']] = aa2num[(int)'J'];
    dinucHead[aa2num[(int)'R']] = aa2num[(int)'J'];
    dinucHead[aa2num[(int)'K']] = aa2num[(int)'J'];
    dinucHead[aa2num[(int)'M']] = aa2num[(int)'O'];
    dinucHead[aa2num[(int)'A']] = aa2num[(int)'O'];
    dinucHead[aa2num[(int)'P']] = aa2num[(int)'O'];
    dinucHead[aa2num[(int)'I']] = aa2num[(int)'O'];
    dinucHead[aa2num[(int)'E']] = aa2num[(int)'U'];
    dinucHead[aa2num[(int)'N']] = aa2num[(int)'U'];
    dinucHead[aa2num[(int)'H']] = aa2num[(int)'U'];
    dinucHead[aa2num[(int)'S']] = aa2num[(int)'U'];
    dinucHead[aa2num[(int)'T']] = aa2num[(int)'B'];
    dinucHead[aa2num[(int)'V']] = aa2num[(int)'J'];
    dinucHead[aa2num[(int)'W']] = aa2num[(int)'O'];
    dinucHead[aa2num[(int)'Y']] = aa2num[(int)'U'];
    dinucHead[aa2num[(int)'B']] = aa2num[(int)'X'];
    dinucHead[aa2num[(int)'J']] = aa2num[(int)'X'];
    dinucHead[aa2num[(int)'O']] = aa2num[(int)'X'];
    dinucHead[aa2num[(int)'U']] = aa2num[(int)'X'];
    dinucHead[aa2num[(int)'X']] = aa2num[(int)'X'];

    // tail — dinuc code -> NX sentinel (second nucleotide as tail dinucleotide)
    dinucTail[aa2num[(int)'C']] = aa2num[(int)'T'];
    dinucTail[aa2num[(int)'G']] = aa2num[(int)'V'];
    dinucTail[aa2num[(int)'L']] = aa2num[(int)'W'];
    dinucTail[aa2num[(int)'Q']] = aa2num[(int)'Y'];
    dinucTail[aa2num[(int)'D']] = aa2num[(int)'T'];
    dinucTail[aa2num[(int)'F']] = aa2num[(int)'V'];
    dinucTail[aa2num[(int)'R']] = aa2num[(int)'W'];
    dinucTail[aa2num[(int)'K']] = aa2num[(int)'Y'];
    dinucTail[aa2num[(int)'M']] = aa2num[(int)'T'];
    dinucTail[aa2num[(int)'A']] = aa2num[(int)'V'];
    dinucTail[aa2num[(int)'P']] = aa2num[(int)'W'];
    dinucTail[aa2num[(int)'I']] = aa2num[(int)'Y'];
    dinucTail[aa2num[(int)'E']] = aa2num[(int)'T'];
    dinucTail[aa2num[(int)'N']] = aa2num[(int)'V'];
    dinucTail[aa2num[(int)'H']] = aa2num[(int)'W'];
    dinucTail[aa2num[(int)'S']] = aa2num[(int)'Y'];
    dinucTail[aa2num[(int)'B']] = aa2num[(int)'T'];
    dinucTail[aa2num[(int)'J']] = aa2num[(int)'V'];
    dinucTail[aa2num[(int)'O']] = aa2num[(int)'W'];
    dinucTail[aa2num[(int)'U']] = aa2num[(int)'Y'];
    dinucTail[aa2num[(int)'T']] = aa2num[(int)'X'];
    dinucTail[aa2num[(int)'V']] = aa2num[(int)'X'];
    dinucTail[aa2num[(int)'W']] = aa2num[(int)'X'];
    dinucTail[aa2num[(int)'Y']] = aa2num[(int)'X'];
    dinucTail[aa2num[(int)'X']] = aa2num[(int)'X'];

    // dinucToNuc — maps dinucleotide encoding back to second nucleotide
    memset(dinucToNucTbl, 0, sizeof(dinucToNucTbl));
    // AA=C, AC=G, AG=L, AU=Q -> A
    dinucToNucTbl[aa2num[(int)'C']] = aa2num[(int)'A'];
    dinucToNucTbl[aa2num[(int)'G']] = aa2num[(int)'A'];
    dinucToNucTbl[aa2num[(int)'L']] = aa2num[(int)'A'];
    dinucToNucTbl[aa2num[(int)'Q']] = aa2num[(int)'A'];
    // CA=D, CC=F, CG=R, CU=K -> C
    dinucToNucTbl[aa2num[(int)'D']] = aa2num[(int)'C'];
    dinucToNucTbl[aa2num[(int)'F']] = aa2num[(int)'C'];
    dinucToNucTbl[aa2num[(int)'R']] = aa2num[(int)'C'];
    dinucToNucTbl[aa2num[(int)'K']] = aa2num[(int)'C'];
    // GA=M, GC=A, GG=P, GU=I -> G
    dinucToNucTbl[aa2num[(int)'M']] = aa2num[(int)'G'];
    dinucToNucTbl[aa2num[(int)'A']] = aa2num[(int)'G'];
    dinucToNucTbl[aa2num[(int)'P']] = aa2num[(int)'G'];
    dinucToNucTbl[aa2num[(int)'I']] = aa2num[(int)'G'];
    // UA=E, UC=N, UG=H, UU=S -> U
    dinucToNucTbl[aa2num[(int)'E']] = aa2num[(int)'U'];
    dinucToNucTbl[aa2num[(int)'N']] = aa2num[(int)'U'];
    dinucToNucTbl[aa2num[(int)'H']] = aa2num[(int)'U'];
    dinucToNucTbl[aa2num[(int)'S']] = aa2num[(int)'U'];
    // XN sentinels (first nuc known, second unknown)
    dinucToNucTbl[aa2num[(int)'T']] = aa2num[(int)'A'];
    dinucToNucTbl[aa2num[(int)'V']] = aa2num[(int)'C'];
    dinucToNucTbl[aa2num[(int)'W']] = aa2num[(int)'G'];
    dinucToNucTbl[aa2num[(int)'Y']] = aa2num[(int)'U'];
    // NX sentinels (first nuc unknown)
    dinucToNucTbl[aa2num[(int)'B']] = aa2num[(int)'X'];
    dinucToNucTbl[aa2num[(int)'J']] = aa2num[(int)'X'];
    dinucToNucTbl[aa2num[(int)'O']] = aa2num[(int)'X'];
    dinucToNucTbl[aa2num[(int)'U']] = aa2num[(int)'X'];
    dinucToNucTbl[aa2num[(int)'X']] = aa2num[(int)'X'];

    dinucTablesReady = true;

    // Create alignSubMat (bitFactor=2.0) for profile_for_alignment if not yet created
    // Guard against recursion: SubstitutionMatrix constructor calls all setupMatrixFn callbacks
    static bool creatingAlignSubMat = false;
    if (dinucAlignSubMat == NULL && !creatingAlignSubMat) {
        creatingAlignSubMat = true;
        char *serialized = BaseMatrix::serialize(mat->matrixName, mat->matrixData);
        dinucAlignSubMat = new SubstitutionMatrix(serialized, 2.0, 0.0f);
        free(serialized);
        creatingAlignSubMat = false;
    }
}

// Build profile arrays from an already-encoded numSequence.
// Generates the profile on-the-fly without writing to disk (no resolution loss).
// Matches reference extractqueryprofiles + prefilter remap behavior.
static void dinucBuildProfileFromSequence(Sequence *seq) {
    if (seq->profile_score == NULL) {
        return;  // profile arrays not allocated
    }
    BaseMatrix *subMat = seq->subMat;
    BaseMatrix *alignMat = dinucAlignSubMat ? dinucAlignSubMat : subMat;
    const bool aaBiasCorrection = (Parameters::getInstance().compBiasCorrection != 0);
    const bool reverse = (seq->getDbKey() % 2 == 1);
    const int seqLen = seq->L;
    const float scale = Parameters::getInstance().compBiasCorrectionScale;

    // For reverse-strand queries, reverse the dinucleotide position order.
    // This matches the old system's extractqueryprofiles which reverses numSequence
    // before writing to disk. The profile columns are then swapped via revcomp.
    if (reverse) {
        std::reverse(seq->numSequence, seq->numSequence + seqLen);
    }

    // Step 1: Build profile_score from kmerSubMat + composition bias.
    // Composition bias IS applied to profile_score because QueryMatcher sets
    // compositionBias=0 for HMM_PROFILE queries, so bias must be baked in here.
    // This matches the old system's forceCompBias behavior.
    float *kmerBias = new float[seqLen + 1];
    memset(kmerBias, 0, (seqLen + 1) * sizeof(float));
    if (aaBiasCorrection) {
        dinucCalcLocalAaBiasCorrection(subMat, seq->numSequence, seqLen,
                                       kmerBias, scale, reverse);
    }
    for (int pos = 0; pos < seqLen; pos++) {
        unsigned char aaIdx = seq->numSequence[pos];
        short bias = static_cast<short>((kmerBias[pos] < 0.0f) ? (kmerBias[pos] - 0.5f) : (kmerBias[pos] + 0.5f));
        if (reverse) {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                seq->profile_score[pos * seq->profile_row_size + aa] =
                    subMat->subMatrix[aaIdx][subMat->num2revcompnum[aa]] + bias;
            }
        } else {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                seq->profile_score[pos * seq->profile_row_size + aa] =
                    subMat->subMatrix[aaIdx][aa] + bias;
            }
        }
    }
    delete[] kmerBias;

    // consensus = sequence for single sequences
    memcpy(seq->numConsensusSequence, seq->numSequence, sizeof(unsigned char) * seqLen);

    // neffM = 1.0 for single sequences
    for (int i = 0; i < seqLen; i++) {
        seq->neffM[i] = 1.0f;
    }

    // Step 2: Build profile_for_alignment from alignSubMat + composition bias.
    // Composition bias is applied here for ungapped diagonal scoring in prefilter.
    float *compositionBias = new float[seqLen + 1];
    memset(compositionBias, 0, (seqLen + 1) * sizeof(float));
    if (aaBiasCorrection) {
        dinucCalcLocalAaBiasCorrection(alignMat, seq->numSequence, seqLen,
                                       compositionBias, scale, reverse);
    }
    for (int i = 0; i < seqLen; i++) {
        unsigned char queryLetter = seq->numSequence[i];
        char bias = static_cast<char>((compositionBias[i] < 0.0f) ? (compositionBias[i] - 0.5f) : (compositionBias[i] + 0.5f));
        if (reverse) {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                seq->profile_for_alignment[aa * seqLen + i] =
                    alignMat->subMatrix[queryLetter][alignMat->num2revcompnum[aa]] + bias;
            }
        } else {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                seq->profile_for_alignment[aa * seqLen + i] =
                    alignMat->subMatrix[queryLetter][aa] + bias;
            }
        }
    }
    // Set X value to 0
    if (subMat->alphabetSize - Sequence::PROFILE_AA_SIZE != 0) {
        memset(&seq->profile_for_alignment[(subMat->alphabetSize - 1) * seqLen], 0, seqLen);
    }
    delete[] compositionBias;

    // Step 3: Suppress non-canonical dinucleotides (16+) for k-mer generation
    if (seq->getKmerSize() != 0) {
        for (int i = 0; i < seqLen; i++) {
            unsigned int indexArray[Sequence::PROFILE_AA_SIZE] = {
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
            for (size_t nc = 16; nc < Sequence::PROFILE_AA_SIZE; nc++) {
                seq->profile_score[i * seq->profile_row_size + nc] = -SHRT_MAX;
            }
            Util::rankedDescSort20(&seq->profile_score[i * seq->profile_row_size],
                                   (unsigned int *)&indexArray);
            memcpy(&seq->profile_index[i * seq->profile_row_size], &indexArray,
                   Sequence::PROFILE_AA_SIZE * sizeof(int));
        }
    }
}

static void dinucMapSequencePadded(Sequence *seq, const char *sequence, unsigned int dataLen) {
    const unsigned char *lookup = seq->subMat->aa2num;
    unsigned char *numSequence = seq->numSequence;
    unsigned int i = 0;
    for (; i + 4 <= dataLen; i += 4) {
        numSequence[i]   = lookup[(unsigned char)sequence[i]];
        numSequence[i+1] = lookup[(unsigned char)sequence[i+1]];
        numSequence[i+2] = lookup[(unsigned char)sequence[i+2]];
        numSequence[i+3] = lookup[(unsigned char)sequence[i+3]];
    }
    for (; i < dataLen; i++) {
        numSequence[i] = lookup[(unsigned char)sequence[i]];
    }
    seq->L = i;
}

static void dinucMapSequenceReverse(Sequence *seq, const char *sequence, unsigned int dataLen);

static void dinucMapSequence(Sequence *seq, const char *sequence, unsigned int dataLen) {
    // Strand-split DBs use even keys for forward, odd keys for the reverse-strand
    // companion entry. Match MMseqs2's mapSequenceReverse layout (X-as-head sentinel)
    // for odd keys so prefilter/align see identical numSequence to the upstream pipeline.
    // [TEMP] X-head dispatch off for diagnostic
    unsigned char *numSequence = seq->numSequence;
    const unsigned char *aa2num = seq->subMat->aa2num;
    const unsigned char *__restrict p = (const unsigned char *)sequence;
    unsigned int i = 0;
    unsigned int s_1 = p[0];
    for (; i + 4 < dataLen; i += 4) {
        unsigned int s_2 = p[i+1], s_3 = p[i+2], s_4 = p[i+3], s_5 = p[i+4];
        numSequence[i]   = aa2num[(s_1 << 8) | s_2];
        numSequence[i+1] = aa2num[(s_2 << 8) | s_3];
        numSequence[i+2] = aa2num[(s_3 << 8) | s_4];
        numSequence[i+3] = aa2num[(s_4 << 8) | s_5];
        s_1 = s_5;
    }
    for (; i < dataLen - 1; i++) {
        unsigned int s_2 = p[i+1];
        numSequence[i] = aa2num[(s_1 << 8) | s_2];
        s_1 = s_2;
    }
    numSequence[i] = aa2num[(s_1 << 8) | (unsigned int)'X'];
    seq->L = i + 1;
    dinucDbgDump("RBS mapSeq", seq->getDbKey(), numSequence, seq->L);
}

static void dinucMapSequenceReverse(Sequence *seq, const char *sequence, unsigned int dataLen) {
    unsigned char *numSequence = seq->numSequence;
    const unsigned char *aa2num = seq->subMat->aa2num;
    const unsigned char *__restrict p = (const unsigned char *)sequence;
    unsigned int i = 0;
    unsigned int s_1 = (unsigned int)'X';
    for (; i + 3 < dataLen; i += 4) {
        unsigned int s_2 = p[i], s_3 = p[i+1], s_4 = p[i+2], s_5 = p[i+3];
        numSequence[i]   = aa2num[(s_1 << 8) | s_2];
        numSequence[i+1] = aa2num[(s_2 << 8) | s_3];
        numSequence[i+2] = aa2num[(s_3 << 8) | s_4];
        numSequence[i+3] = aa2num[(s_4 << 8) | s_5];
        s_1 = s_5;
    }
    for (; i < dataLen; i++) {
        unsigned int s_2 = p[i];
        numSequence[i] = aa2num[(s_1 << 8) | s_2];
        s_1 = s_2;
    }
    seq->L = i;
    dinucDbgDump("RBS mapSeqRev", seq->getDbKey(), numSequence, seq->L);
}

static void dinucMapProfile(Sequence *seq, const char *profileData, unsigned int seqLen) {
    unsigned int extDbtype = DBReader<unsigned int>::getExtendedDbtype(seq->getSeqType());
    bool srcIsSequence = (extDbtype & Parameters::DBTYPE_EXTENDED_SRC_SEQUENCE) != 0;

    if (srcIsSequence) {
        // Data on disk is raw nucleotide characters (from splitstrand on a sequence DB).
        // Dinuc-encode the sequence, then build profile from substitution matrix.
        dinucMapSequence(seq, profileData, seqLen);
        dinucBuildProfileFromSequence(seq);
    } else {
        // Data on disk is a PSSM profile already in dinucleotide space
        // (result2profile built MSA with dinuc-encoded sequences).
        // splitstrand reversed entry order for reverse strand but didn't complement.
        seq->mapProfile(profileData, seqLen);

        // Only apply reverse complement if the profile went through splitstrand,
        // which renumbers keys so odd = reverse strand.
        // result2profile output keeps original keys — no reverse complement needed.
        bool isStrandSplit = (extDbtype & LocalParameters::DBTYPE_EXTENDED_STRAND_SPLIT) != 0;
        const bool reverse = isStrandSplit && (seq->getDbKey() % 2 == 1);
        BaseMatrix *subMat = seq->subMat;
        const unsigned char *revcomp = subMat->num2revcompnum;
        const int L = seq->L;

        // For reverse strand: complement dinucleotide codes in profile.
        // splitstrand on profile DBs already reversed entry order on disk;
        // here we permute scores per position to apply revcomp.
        if (reverse) {
            for (int i = 0; i < L; i++) {
                // Permute profile_score: score at dinuc j moves to revcomp(j)
                short tmpScores[Sequence::PROFILE_AA_SIZE];
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    tmpScores[aa] = seq->profile_score[i * seq->profile_row_size + aa];
                }
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    seq->profile_score[i * seq->profile_row_size + revcomp[aa]] = tmpScores[aa];
                }
                // Complement the query letter
                seq->numSequence[i] = revcomp[seq->numSequence[i]];
            }
        }

        // consensus and neffM (neffM already loaded by mapProfile)
        memcpy(seq->numConsensusSequence, seq->numSequence, sizeof(unsigned char) * L);

        // profile_for_alignment: profile_score / 4 to match reference behavior
        // Must happen BEFORE non-canonical clamping
        for (int i = 0; i < L; i++) {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                seq->profile_for_alignment[aa * L + i] = seq->profile_score[i * seq->profile_row_size + aa] / 4;
            }
        }
        if (subMat->alphabetSize - Sequence::PROFILE_AA_SIZE != 0) {
            memset(&seq->profile_for_alignment[(subMat->alphabetSize - 1) * L], 0, L);
        }

        // Sort profile_score for k-mer generation
        if (seq->getKmerSize() != 0) {
            // Clamp non-canonical entries (16+) to -SHRT_MAX for k-mer generation
            for (int i = 0; i < L; i++) {
                for (size_t nc = 16; nc < Sequence::PROFILE_AA_SIZE; nc++) {
                    seq->profile_score[i * seq->profile_row_size + nc] = -SHRT_MAX;
                }
            }
            for (int i = 0; i < L; i++) {
                unsigned int indexArray[Sequence::PROFILE_AA_SIZE] = {
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
                Util::rankedDescSort20(&seq->profile_score[i * seq->profile_row_size],
                                       (unsigned int *)&indexArray);
                memcpy(&seq->profile_index[i * seq->profile_row_size], &indexArray,
                       Sequence::PROFILE_AA_SIZE * sizeof(int));
            }
        }
    }
}

static void dinucReverseComplement(unsigned char *numSequence, unsigned int L,
                                    const BaseMatrix *subMat) {
    std::reverse(numSequence, numSequence + L - 1); // reverse all except last
    const unsigned char *__restrict complement = subMat->num2revcompnum;
    unsigned int i = 0;
    for (; i + 4 < L - 1; i += 4) {
        numSequence[i]     = complement[numSequence[i]];
        numSequence[i + 1] = complement[numSequence[i + 1]];
        numSequence[i + 2] = complement[numSequence[i + 2]];
        numSequence[i + 3] = complement[numSequence[i + 3]];
    }
    for (; i < L - 1; i++) {
        numSequence[i] = complement[numSequence[i]];
    }
    // Last character: derive from the preceding dinucleotide
    numSequence[i] = dinucTail[numSequence[i - 1]];
}

const unsigned char* getDinucToNucTable() {
    return dinucToNucTbl;
}

void dinucEncodeReverse(Sequence *seq) {
    unsigned char *numSequence = seq->numSequence;
    const char *sequence = seq->getSeqData();
    const unsigned int dataLen = seq->L;
    const unsigned char *aa2num = seq->subMat->aa2num;
    const unsigned char *__restrict p = (const unsigned char *)sequence;
    unsigned int i = 0;
    unsigned int s_1 = (unsigned int)'X';
    for (; i + 3 < dataLen; i += 4) {
        unsigned int s_2 = p[i], s_3 = p[i+1], s_4 = p[i+2], s_5 = p[i+3];
        numSequence[i]   = aa2num[(s_1 << 8) | s_2];
        numSequence[i+1] = aa2num[(s_2 << 8) | s_3];
        numSequence[i+2] = aa2num[(s_3 << 8) | s_4];
        numSequence[i+3] = aa2num[(s_4 << 8) | s_5];
        s_1 = s_5;
    }
    for (; i < dataLen; i++) {
        unsigned int s_2 = p[i];
        numSequence[i] = aa2num[(s_1 << 8) | s_2];
        s_1 = s_2;
    }
    seq->L = i;
}

void registerDinucleotideMapping() {
    Sequence::registerAuxSplit(
        LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE,  // extFlag — dinucleotide encoding
        NULL, NULL, NULL, 0, 16,          // canonicalAlphabetSize=16 (dinuc canonical pairs); non-canonical clamped to ANY in MSA
        dinucMapSequence,
        dinucMapSequenceReverse,
        dinucMapProfile,
        dinucReverseComplement,
        dinucSetupMatrix,
        dinucToNucTbl
    );
}
