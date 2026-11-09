// redkmermatcher: linclust k-mer matching on a reduced alphabet built from the nucleotides of a dinucleotide db, both
// strands at once.
//
// Input, either
//   a nucleotide db (riboseek createdb; one entry per sequence): the pairs are written for rnaalign, query db =
//     splitstrand --strand 2 of the input (key 2i forward, 2i + 1 reverse complement of the i-th entry by key),
//     target db = the input; same-strand members under the forward query, opposite-strand ones under the reverse one
//   a dinucleotide-encoded db (amino acid typed, e.g. igr_cluster.sh's subdb "+-") with both strands as entries
//     (lookup names ending in '+' / '-'): only the '+' entries are read (first nucleotide of each dinucleotide letter;
//     sentinels / X = unknown) and the pairs are written over the input keys for align / gpupairfilter on the input
// The nucleotides are turned into k-mers over a reduced alphabet (--kmer-alph):
//   n4   one letter per nucleotide (A C G U), k nucleotides per k-mer
//   d16  one letter per dinucleotide (the 16 letters of the dinucleotide alphabet), k + 1 nucleotides per k-mer
//   d12  d16 with CA=CG, GC=GU, AC=AU, CC=UC merged (MI-greedy reduction of dinuc.out), k + 1 nucleotides
//   d8   further merged to 8 classes, k + 1 nucleotides
// Every k-mer is compared with the k-mer of the reverse-complement strand over the same nucleotides and the smaller one
// is kept (canonical k-mer, MMseqs2's nucleotide linclust convention: forward k-mers carry bit 63, reverse ones do
// not). Words wider than 63 bits (d12 / d16 k > 17, d8 k > 21) are replaced by a 63-bit hash of the exact word, so k is
// not limited by the 64-bit k-mer index. Selection (kmer-per-seq + scale x length lowest-hash k-mers), grouping, center
// choice (longest, adjacency, count table) and the merge are those of kmermatcher (linclust 2): only the k-mer
// extraction is replaced, the rest calls the vendored MMseqs2 code.
//
// A nucleotide db whose dbtype carries DBTYPE_EXTENDED_DINUCLEOTIDE (riboseek linclust on nucleotides, converted to
// dinucleotides on the fly) gets the pairs for align2clust instead: keyed by the center's input key, members by their
// input key, an opposite-strand member with a negative score and the diagonal of the forward center against the
// member's reverse complement (DBTYPE_PREFILTER_REV_RES, MMseqs2's convention for reverse-strand prefilter hits).
//
// Each pair is written once (canonical strands). Dinucleotide input: keyed by the '+' entry of the center, a same-strand
// member as its '+' entry, an opposite-strand member as its '-' entry (diagonal moved to the '+' query / '-' target
// frame), a plain prefilter db over the input keys. Nucleotide input: keyed by the center's forward (same-strand
// members) or reverse (opposite-strand members) splitstrand query, targets are input keys.
#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "QueryMatcher.h"
#include "SequenceWeights.h"
#include "Timer.h"
#include "Util.h"
#include "kmermatcher.h"
#include "CommandDeclarations.h"
#include "LocalCommandDeclarations.h"
#include "xxhash.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

// vendored kmermatcher.cpp (explicitly instantiated there)
template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
size_t assignGroup(KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair, KmerPosition<T, false, IncludeSeqLen> *writeSeqPair,
                   bool includeOnlyExtendable, int covMode, float covThr,
                   SequenceWeights *sequenceWeights, float weightThr, int threads,
                   std::vector<size_t> &threadOffsets, BaseMatrix *subMat,
                   AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template <typename T, bool includeAdjacency, bool IncludeSeqLen>
KmerPosition<T, includeAdjacency, IncludeSeqLen> *initKmerPositionMemory(size_t size);
uint64_t hashUInt64(uint64_t in, uint64_t seed);
void setLinearFilterDefault(Parameters *p);

namespace {

__extension__ typedef unsigned __int128 u128;
const unsigned char NT_UNK = 4;
const uint64_t MASK63 = (uint64_t(1) << 63) - 1;

struct RedAlphabet {
    std::string name;
    int letters = 0;          // classes; the X letter is index 'letters'
    bool dinuc = false;       // letter of position p from nucleotides p, p + 1
    unsigned char cls[16];    // dinucleotide (4 x first + second) or nucleotide -> class
    unsigned char first[256]; // dinucleotide-db letter -> first nucleotide (0..3) or NT_UNK
};

// dinucleotide letters by (first, second) nucleotide, A C G U order (igrtool encode / dinuc.out)
const char *DI = "CGLQ" "DFRK" "MAPI" "ENHS";
const char *HEAD = "TVWY";   // first nucleotide known, second missing

bool setupAlphabet(const std::string &name, bool nucleotideInput, RedAlphabet &a) {
    a.name = name;
    for (int i = 0; i < 256; i++) a.first[i] = NT_UNK;
    if (nucleotideInput) {
        const char *NT = "ACGU";
        for (int f = 0; f < 4; f++) {
            a.first[(unsigned char) NT[f]] = f;
            a.first[(unsigned char) tolower(NT[f])] = f;
        }
        a.first[(unsigned char) 'T'] = 3;
        a.first[(unsigned char) 't'] = 3;
    }
    for (int f = 0; f < 4 && nucleotideInput == false; f++) {
        for (int s = 0; s < 4; s++) {
            a.first[(unsigned char) DI[4 * f + s]] = f;
            a.first[(unsigned char) tolower(DI[4 * f + s])] = f;
        }
        a.first[(unsigned char) HEAD[f]] = f;
        a.first[(unsigned char) tolower(HEAD[f])] = f;
    }
    // classes as groups of dinucleotide letters (mi_reduce.py)
    std::vector<std::string> groups;
    if (name == "n4") {
        a.dinuc = false; a.letters = 4;
        for (int i = 0; i < 4; i++) a.cls[i] = i;
        return true;
    } else if (name == "d16") {
        for (int i = 0; i < 16; i++) groups.push_back(std::string(1, DI[i]));
    } else if (name == "d12") {
        groups = {"C", "L", "K", "M", "P", "E", "H", "S", "DR", "AI", "GQ", "FN"};
    } else if (name == "d8") {
        groups = {"H", "AI", "GQ", "FN", "KS", "EDR", "MP", "CL"};
    } else {
        return false;
    }
    a.dinuc = true; a.letters = groups.size();
    for (size_t g = 0; g < groups.size(); g++) {
        for (char c : groups[g]) {
            a.cls[strchr(DI, c) - DI] = g;
        }
    }
    return true;
}

// a minimal matrix over the class letters for assignGroup's adjacency scores (+1 / -2, X -1)
class ClassMatrix : public BaseMatrix {
public:
    explicit ClassMatrix(int letters) {
        alphabetSize = letters + 1;
        initMatrixMemory(alphabetSize);
        for (int i = 0; i < alphabetSize; i++) {
            num2aa[i] = (i == letters) ? 'X' : static_cast<char>('A' + i);
            for (int j = 0; j < alphabetSize; j++) {
                subMatrix[i][j] = (i == letters || j == letters) ? -1 : (i == j ? 1 : -2);
            }
        }
    }
};

struct StrandInfo {
    std::vector<size_t> plusIds;           // reader ids of the '+' entries
    std::vector<unsigned int> partner;     // key -> key of the other strand (UINT_MAX: none)
    bool stranded = false;
    bool nucleotide = false;               // nucleotide input: output for splitstrand queries / rnaalign
    bool clusterPairs = false;             // nucleotide input flagged dinucleotide: center-keyed pairs for align2clust
    std::vector<unsigned int> rank;        // nucleotide input: key -> index by key (splitstrand query 2 x rank)
};

StrandInfo nucleotideStrands(DBReader<unsigned int> &seqDbr) {
    StrandInfo s;
    s.nucleotide = true;
    std::vector<unsigned int> keys(seqDbr.getSize());
    for (size_t id = 0; id < seqDbr.getSize(); id++) {
        keys[id] = seqDbr.getDbKey(id);
        s.plusIds.push_back(id);
    }
    std::sort(keys.begin(), keys.end());
    s.rank.assign(size_t(seqDbr.getLastKey()) + 1, UINT_MAX);
    for (size_t i = 0; i < keys.size(); i++) s.rank[keys[i]] = i;
    return s;
}

// '+' / '-' entries from the lookup (names ending in '+' / '-'); without strand suffixes every entry is '+'
StrandInfo readStrands(DBReader<unsigned int> &seqDbr, const std::string &lookupFile) {
    StrandInfo s;
    unsigned int lastKey = seqDbr.getLastKey();
    s.partner.assign(size_t(lastKey) + 1, UINT_MAX);
    std::vector<char> isPlus(size_t(lastKey) + 1, 1);
    struct NameKey { uint64_t h; unsigned int key; char strand; };
    std::vector<NameKey> names;
    if (FileUtil::fileExists(lookupFile.c_str())) {
        std::ifstream in(lookupFile);
        std::string line;
        while (std::getline(in, line)) {
            size_t t1 = line.find('\t');
            if (t1 == std::string::npos) continue;
            size_t t2 = line.find('\t', t1 + 1);
            std::string name = line.substr(t1 + 1, (t2 == std::string::npos ? line.size() : t2) - t1 - 1);
            unsigned int key = strtoul(line.c_str(), NULL, 10);
            if (name.empty() || key > lastKey) continue;
            char st = name.back();
            if (st != '+' && st != '-') continue;
            names.push_back({XXH64(name.data(), name.size() - 1, 0), key, st});
        }
    }
    if (names.empty() == false) {
        s.stranded = true;
        SORT_PARALLEL(names.begin(), names.end(), [](const NameKey &a, const NameKey &b) {
            return a.h < b.h || (a.h == b.h && a.strand < b.strand);
        });
        for (size_t i = 0; i < names.size(); i++) {
            isPlus[names[i].key] = (names[i].strand == '+');
            if (i + 1 < names.size() && names[i].h == names[i + 1].h && names[i].strand != names[i + 1].strand) {
                s.partner[names[i].key] = names[i + 1].key;
                s.partner[names[i + 1].key] = names[i].key;
            }
        }
    }
    for (size_t id = 0; id < seqDbr.getSize(); id++) {
        if (isPlus[seqDbr.getDbKey(id)]) s.plusIds.push_back(id);
    }
    return s;
}

// per-thread k-mer extraction state of one sequence
struct SeqKmers {
    std::vector<unsigned char> nt, sym, rsym;
    std::vector<u128> fwd, rev;
    std::vector<char> valid;
    size_t L = 0, S = 0;      // sequence length, letters
};

void computeKmers(const char *data, size_t L, const RedAlphabet &a, int k, SeqKmers &s) {
    s.L = L;
    s.nt.resize(L + 1);
    for (size_t p = 0; p < L; p++) s.nt[p] = a.first[(unsigned char) data[p]];
    const unsigned char X = a.letters;
    s.S = a.dinuc ? (L > 0 ? L - 1 : 0) : L;
    const size_t S = s.S;
    s.sym.resize(S + 1); s.rsym.resize(S + 1);
    for (size_t p = 0; p < S; p++) {
        if (a.dinuc) {
            unsigned char f = s.nt[p], n = s.nt[p + 1];
            s.sym[p] = (f < 4 && n < 4) ? a.cls[4 * f + n] : X;
            // reverse-complement strand letter q covers nucleotides L-1-q, L-2-q
            unsigned char rf = s.nt[L - 1 - p], rn = s.nt[L - 2 - p];
            s.rsym[p] = (rf < 4 && rn < 4) ? a.cls[4 * (3 - rf) + (3 - rn)] : X;
        } else {
            unsigned char f = s.nt[p], r = s.nt[L - 1 - p];
            s.sym[p] = (f < 4) ? a.cls[f] : X;
            s.rsym[p] = (r < 4) ? a.cls[3 - r] : X;
        }
    }
    s.fwd.resize(S + 1); s.rev.resize(S + 1); s.valid.assign(S + 1, 0);
    if (S < static_cast<size_t>(k)) return;
    u128 top = 1;
    for (int i = 0; i < k - 1; i++) top *= a.letters;
    for (int strand = 0; strand < 2; strand++) {
        const std::vector<unsigned char> &sy = strand ? s.rsym : s.sym;
        std::vector<u128> &out = strand ? s.rev : s.fwd;
        u128 code = 0; size_t run = 0;
        for (size_t p = 0; p < S; p++) {
            unsigned char c = sy[p];
            if (c == X) { run = 0; code = 0; continue; }
            if (run >= static_cast<size_t>(k)) code -= sy[p - k] * top;
            code = code * a.letters + c;
            run++;
            if (run >= static_cast<size_t>(k)) {
                out[p + 1 - k] = code;
                if (strand == 0) s.valid[p + 1 - k] = 1;
            }
        }
    }
}

struct Options {
    RedAlphabet alph;
    bool fits63 = true;
};

// kmermatcher's fillKmerPositionArray (nucleotide branch) with the reduced canonical k-mers of the '+' entries
template <typename T, bool includeAdjacency>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition<T, includeAdjacency> *kmerArray, size_t kmerArraySize,
                                                DBReader<unsigned int> &seqDbr, const StrandInfo &strands,
                                                Parameters &par, const Options &opt,
                                                size_t hashStartRange, size_t hashEndRange, size_t *hashDistribution) {
    size_t offset = 0;
    const int k = par.kmerSize;
    const unsigned char xIndex = opt.alph.letters;
    const size_t nSeq = strands.plusIds.size();
    Debug::Progress progress(nSeq);
#pragma omp parallel num_threads(par.threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        unsigned short *scoreDist = new unsigned short[65536];
        unsigned int *hierarchicalScoreDist = new unsigned int[128];
        const unsigned int BUFFER_SIZE = 1048576;
        size_t bufferPos = 0;
        KmerPosition<T, includeAdjacency> *threadKmerBuffer = new KmerPosition<T, includeAdjacency>[BUFFER_SIZE];
        std::vector<SequencePosition> kmers;
        SeqKmers sk;

        auto flush = [&]() {
            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
            if (writeOffset + bufferPos < kmerArraySize) {
                if (kmerArray != NULL) {
                    memcpy(kmerArray + writeOffset, threadKmerBuffer, sizeof(KmerPosition<T, includeAdjacency>) * bufferPos);
                }
            } else {
                Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset=" << writeOffset
                                    << ", kmerBufferPos=" << bufferPos << ", kmerArraySize=" << kmerArraySize << ".\n";
                EXIT(EXIT_FAILURE);
            }
            bufferPos = 0;
        };

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < nSeq; i++) {
            progress.updateProgress();
            const size_t id = strands.plusIds[i];
            const unsigned int seqId = seqDbr.getDbKey(id);
            const size_t L = seqDbr.getSeqLen(id);
            computeKmers(seqDbr.getData(id, thread_idx), L, opt.alph, k, sk);
            const size_t S = sk.S;
            memset(scoreDist, 0, sizeof(unsigned short) * 65536);
            memset(hierarchicalScoreDist, 0, sizeof(unsigned int) * 128);

            // identity k-mer (hash of the whole nucleotide sequence), forward
            size_t seqHash = hashUInt64(Util::hash(sk.nt.data(), L), par.hashShift);
            seqHash = BIT_SET(seqHash, 63);
            if (seqHash == SIZE_T_MAX) seqHash--;

            kmers.clear();
            if (S >= static_cast<size_t>(k)) {
                for (size_t w = 0; w + k <= S; w++) {
                    if (sk.valid[w] == 0) continue;
                    const size_t q = S - k - w;
                    const u128 f = sk.fwd[w], r = sk.rev[q];
                    // forward and reverse identical: the strand of a match is undefined
                    if (f == r) continue;
                    const bool pickReverse = r < f;
                    const u128 c = pickReverse ? r : f;
                    uint64_t kmerIdx;
                    if (opt.fits63) {
                        kmerIdx = static_cast<uint64_t>(c);
                    } else {
                        kmerIdx = XXH64(&c, sizeof(c), 0x5eed) & MASK63;
                        if (kmerIdx == MASK63) kmerIdx--;
                    }
                    const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);
                    SequencePosition sp;
                    sp.kmer = pickReverse ? BIT_CLEAR(kmerIdx, 63) : BIT_SET(kmerIdx, 63);
                    sp.pos = pickReverse ? q : w;
                    sp.score = hash;
                    kmers.push_back(sp);
                    scoreDist[hash]++;
                    hierarchicalScoreDist[hash >> 9]++;
                }
            }
            const size_t seqKmerCount = kmers.size();
            float kmersPerSequenceScale = par.kmersPerSequenceScale.values.nucleotide();
            size_t kmerConsidered = std::min(static_cast<size_t>(par.kmersPerSequence - 1 + (kmersPerSequenceScale * L)), seqKmerCount);

            unsigned int threshold = 0;
            size_t kmerInBins = 0;
            if (seqKmerCount > 0) {
                size_t hierarchicaThreshold = 0;
                for (hierarchicaThreshold = 0; hierarchicaThreshold < 128 && kmerInBins < kmerConsidered; hierarchicaThreshold++) {
                    kmerInBins += hierarchicalScoreDist[hierarchicaThreshold];
                }
                hierarchicaThreshold -= (hierarchicaThreshold > 0) ? 1 : 0;
                kmerInBins -= hierarchicalScoreDist[hierarchicaThreshold];
                for (threshold = hierarchicaThreshold * 512; threshold <= USHRT_MAX && kmerInBins < kmerConsidered; threshold++) {
                    kmerInBins += scoreDist[threshold];
                }
            }
            int tooMuchElemInLastBin = (kmerInBins - kmerConsidered);

            if (static_cast<unsigned short>(seqHash) >= hashStartRange && static_cast<unsigned short>(seqHash) <= hashEndRange) {
                if (hashDistribution != NULL) {
                    __sync_fetch_and_add(&hashDistribution[static_cast<unsigned short>(seqHash)], 1);
                } else {
                    threadKmerBuffer[bufferPos].kmer = seqHash;
                    threadKmerBuffer[bufferPos].id = seqId;
                    threadKmerBuffer[bufferPos].pos = 0;
                    threadKmerBuffer[bufferPos].sl.setSeqLen(static_cast<T>(L));
                    if (includeAdjacency) {
                        for (size_t j = 0; j < 6; j++) threadKmerBuffer[bufferPos].setAdjacentSeq(j, xIndex);
                    }
                    bufferPos++;
                    if (bufferPos >= BUFFER_SIZE) flush();
                }
            }

            size_t selectedKmer = 0;
            for (size_t kmerIdx = 0; kmerIdx < seqKmerCount && selectedKmer < kmerConsidered; kmerIdx++) {
                SequencePosition &sp = kmers[kmerIdx];
                if (sp.score < threshold) {
                    if (sp.score == (threshold - 1) && tooMuchElemInLastBin) {
                        tooMuchElemInLastBin--;
                        threshold -= (tooMuchElemInLastBin == 0) ? 1 : 0;
                    }
                    selectedKmer++;
                    if (sp.score >= hashStartRange && sp.score <= hashEndRange) {
                        if (hashDistribution != NULL) {
                            __sync_fetch_and_add(&hashDistribution[sp.score], 1);
                            continue;
                        }
                        threadKmerBuffer[bufferPos].kmer = sp.kmer;
                        threadKmerBuffer[bufferPos].id = seqId;
                        threadKmerBuffer[bufferPos].pos = sp.pos;
                        threadKmerBuffer[bufferPos].sl.setSeqLen(static_cast<T>(L));
                        if (includeAdjacency) {
                            // three letters on each side of the k-mer, on the strand it was taken from
                            const bool isRev = BIT_CHECK(sp.kmer, 63) == false;
                            const std::vector<unsigned char> &sy = isRev ? sk.rsym : sk.sym;
                            const long start = sp.pos, end = static_cast<long>(sp.pos) + k - 1;
                            for (int j = 0; j < 3; j++) {
                                long before = start - 3 + j, after = end + 1 + j;
                                threadKmerBuffer[bufferPos].setAdjacentSeq(j, before >= 0 ? sy[before] : xIndex);
                                threadKmerBuffer[bufferPos].setAdjacentSeq(3 + j, after < static_cast<long>(S) ? sy[after] : xIndex);
                            }
                        }
                        bufferPos++;
                        if (bufferPos >= BUFFER_SIZE) flush();
                    }
                }
            }
        }
        if (bufferPos > 0) {
            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
            if (kmerArray != NULL) {
                memcpy(kmerArray + writeOffset, threadKmerBuffer, sizeof(KmerPosition<T, includeAdjacency>) * bufferPos);
            }
        }
        delete[] threadKmerBuffer;
        delete[] hierarchicalScoreDist;
        delete[] scoreDist;
    }
    return std::make_pair(offset, static_cast<size_t>(par.kmerSize));
}

// kmermatcher.cpp's runIteration, nucleotide branch
template <typename T, bool includeAdjacency>
void runIteration(AssignGroupMask mask, int &iteration, size_t &writePos, size_t hashEndRange, const std::string &splitFile,
                  KmerPosition<T, includeAdjacency> *hashSeqPair, KmerPosition<T, false> *writeSeqPair,
                  size_t nSeqForThreads, Parameters &par, BaseMatrix *subMat, short *countTable,
                  std::vector<size_t> &threadOffsets, Timer &timer) {
    timer.reset();
    writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, false>(
        hashSeqPair, writeSeqPair, par.includeOnlyExtendable, par.covMode, par.covThr, NULL, par.weightThr,
        par.threads, threadOffsets, subMat, mask, ComputationPhase::Main, countTable);
    Debug(Debug::INFO) << "Time for assign: " << timer.lap() << "\n";
    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    if (par.needWriteBuffer) {
        SORT_PARALLEL(writeSeqPair, writeSeqPair + writePos, KmerPosition<T, false>::compareRepSequenceAndIdAndDiagReverse);
    } else {
        SORT_PARALLEL(hashSeqPair, hashSeqPair + writePos, KmerPosition<T, includeAdjacency>::compareRepSequenceAndIdAndDiagReverse);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";
    if (hashEndRange != SIZE_T_MAX || par.needWriteBuffer) {
        std::vector<size_t> threadQueryOffsets(par.threads + 1);
        size_t qSplitSize = nSeqForThreads / par.threads;
        threadQueryOffsets[0] = 0;
        // split the written pairs on rep key = thread x (key range / threads); the rep keys carry the strand in bit
        // 63, so compare without it (kmermatcher's plain comparison can cut a rep's pairs across two threads)
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (int thread = 1; thread < par.threads; thread++) {
            size_t startqid = qSplitSize * thread;
            if (par.needWriteBuffer) {
                KmerPosition<T, false> *it = std::lower_bound(writeSeqPair, writeSeqPair + writePos, startqid,
                    [](const KmerPosition<T, false> &elem, size_t key) { return (elem.kmer & MASK63) < key; });
                threadQueryOffsets[thread] = it - writeSeqPair;
            } else {
                KmerPosition<T, includeAdjacency> *it = std::lower_bound(hashSeqPair, hashSeqPair + writePos, startqid,
                    [](const KmerPosition<T, includeAdjacency> &elem, size_t key) { return (elem.kmer & MASK63) < key; });
                threadQueryOffsets[thread] = it - hashSeqPair;
            }
        }
        threadQueryOffsets[par.threads] = writePos;
        timer.reset();
        if (par.needWriteBuffer) {
            writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, false, false>(
                splitFile, writeSeqPair, writePos + 1, par.threads, &threadQueryOffsets, iteration);
        } else {
            writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, includeAdjacency, false>(
                splitFile, hashSeqPair, writePos + 1, par.threads, &threadQueryOffsets, iteration);
        }
        Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
    }
    iteration++;
}

// kmermatcher.cpp's doComputation, nucleotide branch
template <typename T, bool includeAdjacency>
void doComputation(size_t totalKmers, size_t hashStartRange, size_t hashEndRange, const std::string &splitFile,
                   AssignGroupMask assignGroupMask, ComputationPhase phase, DBReader<unsigned int> &seqDbr,
                   const StrandInfo &strands, Parameters &par, const Options &opt, BaseMatrix *subMat, short *countTable) {
    KmerPosition<T, includeAdjacency> *hashSeqPair = initKmerPositionMemory<T, includeAdjacency, false>(totalKmers);
    KmerPosition<T, false> *writeSeqPair = NULL;
    if (phase == ComputationPhase::Main && par.needWriteBuffer) {
        writeSeqPair = initKmerPositionMemory<T, false, false>(totalKmers);
    }
    std::pair<size_t, size_t> ret = fillKmerPositionArray<T, includeAdjacency>(
        hashSeqPair, totalKmers, seqDbr, strands, par, opt, hashStartRange, hashEndRange, NULL);
    size_t elementsToSort = ret.first;
    Debug(Debug::INFO) << "\nSort kmer ";
    Timer timer;
    SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition<T, includeAdjacency>::compareRepSequenceAndIdAndPosReverse);
    Debug(Debug::INFO) << timer.lap() << "\n";

    std::vector<size_t> threadOffsets;
    size_t splitSize = elementsToSort / par.threads;
    threadOffsets.push_back(0);
    for (int thread = 1; thread < par.threads; thread++) {
        if (!par.useParallelism) {
            threadOffsets.push_back(elementsToSort);
            continue;
        }
        size_t prevHash = hashSeqPair[thread * splitSize].kmer;
        if (prevHash == SIZE_T_MAX) {
            for (int i = thread; i < par.threads; i++) threadOffsets.push_back(elementsToSort);
            break;
        }
        prevHash = BIT_SET(prevHash, 63);
        bool wasSet = false;
        for (size_t pos = thread * splitSize; pos < elementsToSort; pos++) {
            size_t currKmer = BIT_SET(hashSeqPair[pos].kmer, 63);
            if (prevHash != currKmer) {
                wasSet = true;
                threadOffsets.push_back(pos);
                break;
            }
        }
        if (wasSet == false) {
            for (int i = thread; i < par.threads; i++) threadOffsets.push_back(elementsToSort);
            break;
        }
    }
    threadOffsets.push_back(elementsToSort);

    if (phase == ComputationPhase::SetupCountTable) {
        assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, false>(
            hashSeqPair, NULL, par.includeOnlyExtendable, par.covMode, par.covThr, NULL, par.weightThr,
            par.threads, threadOffsets, subMat, AssignGroupFeature::Default, ComputationPhase::SetupCountTable, countTable);
        delete[] hashSeqPair;
        delete[] writeSeqPair;
        return;
    }
    const size_t nSeqForThreads = seqDbr.getLastKey() + 1;
    int iteration = 0;
    size_t writePos = 0;
    runIteration<T, includeAdjacency>(AssignGroupFeature::Default, iteration, writePos, hashEndRange, splitFile,
                                      hashSeqPair, writeSeqPair, nSeqForThreads, par, subMat, countTable, threadOffsets, timer);
    if (hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq)) {
        for (int iter = 0; iter < par.adjIteration; iter++) {
            runIteration<T, includeAdjacency>(AssignGroupFeature::AdjacentSeq, iteration, writePos, hashEndRange, splitFile,
                                              hashSeqPair, writeSeqPair, nSeqForThreads, par, subMat, countTable, threadOffsets, timer);
        }
    }
    if (hasFeature(assignGroupMask, AssignGroupFeature::CountTable)) {
        for (int iter = 0; iter < par.countTableIteration; iter++) {
            runIteration<T, includeAdjacency>(AssignGroupFeature::CountTable, iteration, writePos, hashEndRange, splitFile,
                                              hashSeqPair, writeSeqPair, nSeqForThreads, par, subMat, countTable, threadOffsets, timer);
        }
    }
    std::string splitDone = splitFile + ".done";
    FILE *done = FileUtil::openFileOrDie(splitDone.c_str(), "w", false);
    if (fclose(done) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << splitDone << "\n";
        EXIT(EXIT_FAILURE);
    }
    delete[] hashSeqPair;
    delete[] writeSeqPair;
}

// hash-range splits (kmermatcher.cpp's setupKmerSplits / setupCountTable over the '+' entries)
template <typename T, bool includeAdjacency>
std::vector<std::pair<size_t, size_t>> hashRangesFor(DBReader<unsigned int> &seqDbr, const StrandInfo &strands,
                                                     Parameters &par, const Options &opt, size_t perSplit,
                                                     size_t countTableTotal, bool countTable) {
    std::vector<std::pair<size_t, size_t>> hashRanges;
    size_t *hashDist = new size_t[USHRT_MAX + 1];
    memset(hashDist, 0, sizeof(size_t) * (USHRT_MAX + 1));
    fillKmerPositionArray<T, includeAdjacency>(NULL, SIZE_T_MAX, seqDbr, strands, par, opt, 0, SIZE_T_MAX, hashDist);
    size_t maxBucketSize = 0;
    for (size_t i = 0; i < (USHRT_MAX + 1); i++) maxBucketSize = std::max(maxBucketSize, hashDist[i]);
    if (maxBucketSize > perSplit) {
        Debug(Debug::ERROR) << "Not enough memory to run redkmermatcher. Minimum is at least "
                            << maxBucketSize * sizeof(KmerPosition<T, includeAdjacency>) << " bytes\n";
        EXIT(EXIT_FAILURE);
    }
    size_t currBucketSize = 0, currBucketStart = 0, totalBucketSize = 0;
    for (size_t i = 0; i < (USHRT_MAX + 1); i++) {
        if (countTable && totalBucketSize + hashDist[i] >= countTableTotal) {
            hashRanges.emplace_back(currBucketStart, i - 1);
            delete[] hashDist;
            return hashRanges;
        }
        if (currBucketSize + hashDist[i] >= perSplit) {
            hashRanges.emplace_back(currBucketStart, i - 1);
            currBucketSize = 0;
            currBucketStart = i;
        }
        currBucketSize += hashDist[i];
        totalBucketSize += hashDist[i];
    }
    if (countTable == false) hashRanges.emplace_back(currBucketStart, (USHRT_MAX + 1));
    delete[] hashDist;
    return hashRanges;
}

// center-keyed pairs with a reverse flag (negative score) -> plain prefilter pairs (see the header for the two layouts)
void strandPairs(const std::string &tmpDb, DBReader<unsigned int> &seqDbr, const StrandInfo &strands, Parameters &par) {
    std::string tmpIndex = tmpDb + ".index";
    DBReader<unsigned int> pref(tmpDb.c_str(), tmpIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    pref.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed,
                 strands.clusterPairs ? Parameters::DBTYPE_PREFILTER_REV_RES : Parameters::DBTYPE_PREFILTER_RES);
    dbw.open();
    std::vector<int> lenByKey(size_t(seqDbr.getLastKey()) + 1, 0);
    for (size_t id = 0; id < seqDbr.getSize(); id++) lenByKey[seqDbr.getDbKey(id)] = seqDbr.getSeqLen(id);
    size_t dropped = 0;
    Debug::Progress progress(pref.getSize());
#pragma omp parallel num_threads(par.threads) reduction(+:dropped)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string out;
        char buffer[100];
#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < pref.getSize(); i++) {
            progress.updateProgress();
            const unsigned int qKey = pref.getDbKey(i);
            const int qLen = lenByKey[qKey];
            char *data = pref.getData(i, thread_idx);
            out.clear();
            if (strands.clusterPairs) {
                // center key, member keys; reverse: the stored diagonal is reverse center vs forward member, written
                // as forward center vs reverse-complement member, the score stays negative
                while (*data != '\0') {
                    hit_t h = QueryMatcher::parsePrefilterHit(data);
                    data = Util::skipLine(data);
                    if (h.prefScore < 0) {
                        if (h.seqId == qKey) continue;
                        const int diag = (qLen - lenByKey[h.seqId]) - static_cast<short>(h.diagonal);
                        h.diagonal = static_cast<unsigned short>(static_cast<short>(diag));
                    }
                    size_t len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                    out.append(buffer, len);
                }
                dbw.writeData(out.c_str(), out.size(), qKey, thread_idx);
                continue;
            }
            if (strands.nucleotide) {
                // forward query 2r: same-strand members; reverse query 2r + 1: opposite-strand members (the stored
                // diagonal already is reverse query vs forward target; rnaalign takes the strand from the query key)
                std::string outRev;
                while (*data != '\0') {
                    hit_t h = QueryMatcher::parsePrefilterHit(data);
                    data = Util::skipLine(data);
                    bool rev = h.prefScore < 0;
                    if (rev && h.seqId == qKey) continue;
                    h.prefScore = std::abs(h.prefScore);
                    size_t len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                    (rev ? outRev : out).append(buffer, len);
                }
                const unsigned int r = strands.rank[qKey];
                if (out.empty() == false) dbw.writeData(out.c_str(), out.size(), 2 * r, thread_idx);
                if (outRev.empty() == false) dbw.writeData(outRev.c_str(), outRev.size(), 2 * r + 1, thread_idx);
                continue;
            }
            while (*data != '\0') {
                hit_t h = QueryMatcher::parsePrefilterHit(data);
                data = Util::skipLine(data);
                if (h.prefScore < 0) {
                    // reverse: the stored diagonal is query '-' vs target '+'; write query '+' vs target '-'
                    const unsigned int tMinus = strands.partner[h.seqId];
                    if (tMinus == UINT_MAX) { dropped++; continue; }
                    // an inverted repeat inside one entry pairs it with its own reverse complement
                    if (h.seqId == qKey) continue;
                    const int tLen = lenByKey[h.seqId];
                    const int diag = (qLen - tLen) - static_cast<short>(h.diagonal);
                    h.seqId = tMinus;
                    h.prefScore = -h.prefScore;
                    h.diagonal = static_cast<unsigned short>(static_cast<short>(diag));
                }
                size_t len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                out.append(buffer, len);
            }
            dbw.writeData(out.c_str(), out.size(), qKey, thread_idx);
        }
    }
    dbw.close();
    pref.close();
    if (dropped > 0) {
        Debug(Debug::WARNING) << dropped << " opposite-strand pairs dropped: no '-' entry for the member\n";
    }
}

template <typename T, bool includeAdjacency>
int redKmerMatcherInner(Parameters &par, DBReader<unsigned int> &seqDbr, const StrandInfo &strands, const Options &opt) {
    const size_t dbKeySize = seqDbr.getLastKey() + 1;
    ClassMatrix subMat(opt.alph.letters);
    size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);

    // k-mers kept per '+' entry (computeKmerCount)
    const float scale = par.kmersPerSequenceScale.values.nucleotide();
    size_t totalKmers = 0;
    for (size_t id : strands.plusIds) {
        int seqLen = static_cast<int>(seqDbr.getSeqLen(id));
        int kmerAdjustedSeqLen = std::max(1, seqLen - static_cast<int>(par.kmerSize) + 2);
        totalKmers += std::min(kmerAdjustedSeqLen, static_cast<int>(par.kmersPerSequence + (scale * seqLen)));
    }
    const size_t perKmer = sizeof(KmerPosition<T, includeAdjacency>) + (par.needWriteBuffer ? sizeof(KmerPosition<T, false>) : 0);
    size_t totalSizeNeeded = totalKmers * perKmer;
    size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
    size_t totalKmersPerSplit = std::max(static_cast<size_t>(1024 + 1), static_cast<size_t>(std::min(totalSizeNeeded, memoryLimit) / perKmer) + 1);
    Debug(Debug::INFO) << strands.plusIds.size() << " '+' entries, " << totalKmers << " k-mers\n";

    T *seqkey_to_len = new (std::nothrow) T[dbKeySize + 1];
    Util::checkAllocation(seqkey_to_len, "Can not allocate seqkey_to_len memory");
    memset(seqkey_to_len, 0, sizeof(T) * (dbKeySize + 1));
    for (size_t id = 0; id < seqDbr.getSize(); id++) {
        seqkey_to_len[seqDbr.getDbKey(id)] = static_cast<T>(seqDbr.getSeqLen(id));
    }
    SeqLenData<T, false>::seqkey_to_len = seqkey_to_len;

    std::vector<short> countTable;
    if (par.includeCountTable) {
        countTable.assign(dbKeySize, 0);
        size_t countTableTotalKmers = static_cast<size_t>(totalKmers * par.countTableScale);
        size_t countTableTotalSizeNeeded = sizeof(KmerPosition<T, false>) * countTableTotalKmers * 2;
        size_t availableMemory = memoryLimit - countTable.size() * sizeof(short);
        size_t countTableKmersPerSplit = std::max(static_cast<size_t>(1024 + 1),
            static_cast<size_t>(std::min(countTableTotalSizeNeeded, availableMemory) / (sizeof(KmerPosition<T, false>) * 2)) + 1);
        Debug(Debug::INFO) << "Initiating count table\n";
        std::vector<std::pair<size_t, size_t>> ranges = hashRangesFor<T, false>(
            seqDbr, strands, par, opt, countTableKmersPerSplit, countTableTotalKmers, true);
        for (size_t split = 0; split < ranges.size(); split++) {
            Debug(Debug::INFO) << "Fill count table for " << (split + 1) << " split\n";
            doComputation<T, false>(countTableKmersPerSplit, ranges[split].first, ranges[split].second, "COUNT_TABLE",
                                    AssignGroupFeature::Default, ComputationPhase::SetupCountTable,
                                    seqDbr, strands, par, opt, &subMat, countTable.data());
        }
    }
    AssignGroupMask assignGroupMask = AssignGroupFeature::Default;
    if (par.includeAdjacency) assignGroupMask |= AssignGroupFeature::AdjacentSeq;
    if (par.includeCountTable) assignGroupMask |= AssignGroupFeature::CountTable;

    std::vector<std::pair<size_t, size_t>> hashRanges;
    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split\n";
        hashRanges = hashRangesFor<T, includeAdjacency>(seqDbr, strands, par, opt, totalKmersPerSplit, 0, false);
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    } else {
        hashRanges.emplace_back(0, SIZE_T_MAX);
    }
    // the pairs go to a '+'-keyed temporary db first (strandPairs writes the result)
    const std::string tmpDb = par.db2 + "_strand";
    std::vector<std::string> splitFiles;
    for (size_t split = 0; split < hashRanges.size(); split++) {
        std::string splitFileName = par.db2 + "_split_" + SSTR(split);
        Debug(Debug::INFO) << "Generate k-mers list for " << (split + 1) << " split\n";
        std::string splitFileNameDone = splitFileName + ".done";
        if (FileUtil::fileExists(splitFileNameDone.c_str()) == false) {
            doComputation<T, includeAdjacency>(totalKmersPerSplit, hashRanges[split].first, hashRanges[split].second,
                                               splitFileName, assignGroupMask, ComputationPhase::Main,
                                               seqDbr, strands, par, opt, &subMat, countTable.empty() ? NULL : countTable.data());
        }
        splitFiles.push_back(splitFileName);
    }
    if (splits <= 1 && par.needWriteBuffer == false) {
        Debug(Debug::ERROR) << "redkmermatcher needs the merge path (linclust 2 with adjacency or count table)\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<char> repSequence(dbKeySize, false);
    {
        std::string tmpIndex = tmpDb + ".index";
        DBWriter dbw(tmpDb.c_str(), tmpIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_PREFILTER_REV_RES);
        dbw.open();
        Timer timer;
        int maxIter = 1;
        if (par.includeAdjacency) maxIter += par.adjIteration;
        if (par.includeCountTable) maxIter += par.countTableIteration;
        mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, includeAdjacency>(dbw, splitFiles, repSequence, par.threads, maxIter);
        for (int iter = 0; iter < maxIter; ++iter) {
            for (size_t i = 0; i < splitFiles.size(); ++i) {
                std::string doneFilePath = splitFiles[i] + "_iter_" + std::to_string(iter) + ".done";
                if (FileUtil::fileExists(doneFilePath.c_str())) FileUtil::remove(doneFilePath.c_str());
            }
        }
        for (size_t i = 0; i < splitFiles.size(); ++i) {
            std::string splitDone = splitFiles[i] + ".done";
            if (FileUtil::fileExists(splitDone.c_str())) FileUtil::remove(splitDone.c_str());
        }
        Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";
        // '+' entries that are in no group: self hit only
        char buffer[100];
        for (size_t id : strands.plusIds) {
            unsigned int dbKey = seqDbr.getDbKey(id);
            if (repSequence[dbKey] == false) {
                hit_t h;
                h.prefScore = 0;
                h.diagonal = 0;
                h.seqId = dbKey;
                int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                dbw.writeData(buffer, len, dbKey, 0);
            }
        }
        dbw.close(false, false);
    }
    delete[] seqkey_to_len;
    SeqLenData<T, false>::seqkey_to_len = NULL;

    Debug(Debug::INFO) << "Write both-strand pairs\n";
    strandPairs(tmpDb, seqDbr, strands, par);
    DBReader<unsigned int>::removeDb(tmpDb);
    return EXIT_SUCCESS;
}

}  // namespace

int redkmermatcher(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setLinearFilterDefault(&par);
    par.kmerSize = 0;
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    const bool nucleotideInput = Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    if (nucleotideInput == false && Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) == false) {
        Debug(Debug::ERROR) << "redkmermatcher needs a nucleotide or a dinucleotide-encoded (amino acid typed) database\n";
        return EXIT_FAILURE;
    }
    Options opt;
    if (setupAlphabet(par.kmerAlph, nucleotideInput, opt.alph) == false) {
        Debug(Debug::ERROR) << "--kmer-alph must be one of n4, d16, d12, d8\n";
        return EXIT_FAILURE;
    }
    if (par.maskLowerCaseMode) {
        // --mask-lower-case: soft-masked (lowercase) letters start or extend no k-mer, as in kmermatcher
        for (int c = 'a'; c <= 'z'; c++) opt.alph.first[c] = NT_UNK;
    }
    StrandInfo strands = nucleotideInput ? nucleotideStrands(seqDbr) : readStrands(seqDbr, par.db1 + ".lookup");
    strands.clusterPairs = nucleotideInput && (DBReader<unsigned int>::getExtendedDbtype(seqDbr.getDbtype())
                                               & LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0;
    if (nucleotideInput == false && strands.stranded == false) {
        Debug(Debug::WARNING) << "no '+' / '-' names in the lookup: every entry is read as '+', opposite-strand pairs are dropped\n";
    }
    if (par.kmerSize == 0) {
        // MMseqs2's nucleotide rule (one random k-mer per database) in nucleotides; a dinucleotide k-mer spans k + 1
        size_t residues = 0;
        for (size_t id : strands.plusIds) residues += seqDbr.getSeqLen(id);
        int nt = std::max(17, static_cast<int>(log(static_cast<double>(std::max<size_t>(residues, 2))) / log(4.0)));
        par.kmerSize = opt.alph.dinuc ? nt - 1 : nt;
    }
    const int maxK = opt.alph.dinuc ? static_cast<int>(127.0 / std::log2(static_cast<double>(opt.alph.letters))) : 63;
    if (par.kmerSize < 1 || par.kmerSize > maxK) {
        Debug(Debug::ERROR) << "-k must be in 1.." << maxK << " for --kmer-alph " << opt.alph.name << "\n";
        return EXIT_FAILURE;
    }
    opt.fits63 = (par.kmerSize * std::log2(static_cast<double>(opt.alph.letters))) < 63.0;

    std::vector<MMseqsParameter *> *params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    Debug(Debug::INFO) << "Database size: " << seqDbr.getSize() << ", alphabet " << opt.alph.name << " (" << opt.alph.letters
                       << " letters), k " << par.kmerSize << (opt.fits63 ? "" : " (hashed words)") << "\n";

    if (par.includeCountTable || par.includeAdjacency) {
        par.needWriteBuffer = true;
        par.useParallelism = true;
    } else {
        par.needWriteBuffer = par.useParallelism;
    }
    if (par.includeAdjacency && par.adjIteration == 0) {
        Debug(Debug::ERROR) << "Adjacent Iteration must be greater than 0 when include Adjacent is true\n";
        return EXIT_FAILURE;
    }
    if (par.includeCountTable && par.countTableIteration == 0) {
        Debug(Debug::ERROR) << "CountTable Iteration must be greater than 0 when include CountTable is true\n";
        return EXIT_FAILURE;
    }
    int ret;
    if (seqDbr.getMaxSeqLen() < SHRT_MAX) {
        ret = par.includeAdjacency ? redKmerMatcherInner<short, true>(par, seqDbr, strands, opt)
                                   : redKmerMatcherInner<short, false>(par, seqDbr, strands, opt);
    } else {
        ret = par.includeAdjacency ? redKmerMatcherInner<int, true>(par, seqDbr, strands, opt)
                                   : redKmerMatcherInner<int, false>(par, seqDbr, strands, opt);
    }
    seqDbr.close();
    return ret;
}

// kmermatcher as riboseek linclust calls it (both levels of linclust 2) runs redkmermatcher when
//   the input is a nucleotide db flagged DBTYPE_EXTENDED_DINUCLEOTIDE (n4 k-mers, -k in nucleotides, pairs for
//     align2clust over both strands), or
//   RIBOSEEK_KMER_ALPH is set (riboseek linclust --kmer-alph: that alphabet, also on a dinucleotide-encoded db);
// otherwise the MMseqs2 kmermatcher, unchanged ("base:kmermatcher" always reaches it)
int riboseekKmermatcher(int argc, const char **argv, const Command &command) {
    const char *alph = getenv("RIBOSEEK_KMER_ALPH");
    bool flagged = false;
    if (argc > 0) {
        const std::string dbtypeFile = std::string(argv[0]) + ".dbtype";
        if (FileUtil::fileExists(dbtypeFile.c_str())) {
            const int dbtype = FileUtil::parseDbType(argv[0]);
            flagged = Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_NUCLEOTIDES)
                      && (DBReader<unsigned int>::getExtendedDbtype(dbtype) & LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0;
        }
    }
    const bool alphSet = (alph != NULL && *alph != '\0');
    if (flagged == false && alphSet == false) {
        return kmermatcher(argc, argv, command);
    }
    std::vector<const char *> args(argv, argv + argc);
    if (alphSet) {
        args.push_back("--kmer-alph");
        args.push_back(alph);
    }
    return redkmermatcher(static_cast<int>(args.size()), args.data(), command);
}
