// gpupairfilter: score every (query, target) pair of a prefilter-format result db (e.g. kmermatcher pairs) with
// Smith-Waterman on the GPU (libmarv PairSW: int32, no saturation) and keep only the pairs whose E-value could pass
// `align -e`.  `align` on the kept pairs then gives the same result as on all pairs: the GPU score uses align's
// query profile (substitution matrix + rounded local composition bias, StripedSmithWaterman SUBSTITUTIONMATRIX
// mode) and gap costs, and is never below align's score (align's 16-bit pass saturates; the GPU does not).
// Kept: E-value <= -e x 10 (slack against rounding), and the query itself (align scores identities separately).
// Takes the align parameters; only the amino-acid (riboseek dinucleotide) path is implemented.
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "EvalueComputation.h"
#include "Parameters.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Timer.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef HAVE_CUDA
#include "pairsw.h"
#endif

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <cmath>
#include <string>
#include <vector>

int gpupairfilter(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
#ifndef HAVE_CUDA
    Debug(Debug::ERROR) << "gpupairfilter needs a build with ENABLE_CUDA=1\n";
    EXIT(EXIT_FAILURE);
#else
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    const bool same = (par.db1 == par.db2);
    DBReader<unsigned int> *tdbr = &qdbr;
    if (!same) {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        tdbr->open(DBReader<unsigned int>::NOSORT);
    }
    DBReader<unsigned int> pref(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    pref.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    const int qtype = qdbr.getDbtype(), ttype = tdbr->getDbtype();
    if (Parameters::isEqualDbtype(qtype, Parameters::DBTYPE_NUCLEOTIDES) || Parameters::isEqualDbtype(qtype, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "gpupairfilter: only amino-acid typed (dinucleotide) sequence dbs are supported\n";
        EXIT(EXIT_FAILURE);
    }
    SubstitutionMatrix m(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    const int gapOpen = par.gapOpen.values.aminoacid(), gapExtend = par.gapExtend.values.aminoacid();
    const int alph = m.alphabetSize;
    std::vector<int8_t> mat(size_t(alph) * alph);
    for (int i = 0; i < alph; i++)
        for (int j = 0; j < alph; j++) mat[size_t(i) * alph + j] = int8_t(m.subMatrix[i][j]);
    EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), &m, gapOpen, gapExtend);
    const double evalKeep = par.evalThr * 10.0;
    // GPUPAIRFILTER_BIAS=upper: per-position composition bias max(0, rounded bias), so the score bounds every scoring
    // align2clust uses (ungapped without bias, banded with truncated bias, Smith-Waterman with rounded bias): the
    // filter then keeps every pair align2clust could accept (linclust.sh sets it)
    const bool biasUpper = getenv("GPUPAIRFILTER_BIAS") && std::string(getenv("GPUPAIRFILTER_BIAS")) == "upper";
    if (biasUpper) Debug(Debug::INFO) << "gpupairfilter: upper-bound composition bias (align2clust mode)\n";

    const size_t maxLetters = size_t(1) << 30, maxPairs = size_t(8) << 20;
    if (PairSW::deviceCount() < 1) {
        Debug(Debug::ERROR) << "gpupairfilter: no CUDA device\n";
        EXIT(EXIT_FAILURE);
    }
    // byte -> letter, as Sequence::mapSequence does it (one letter at a time), applied on the GPU to raw bytes
    uint8_t letterMap[256] = {0};
    {
        Sequence s(par.maxSeqLen, ttype, &m, 0, false, false);
        for (int c = 1; c < 256; c++) {
            const char ch = char(c);
            s.mapSequence(0, 0, &ch, 1);
            letterMap[c] = s.numSequence[0];
        }
    }
    PairSW sw(mat.data(), alph, gapOpen, gapExtend, maxLetters, maxPairs, letterMap);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, pref.getDbtype());
    writer.open();

    Timer tPhase;
    // key -> id of the target db (direct array instead of a binary search per pair)
    unsigned int maxKey = 0;
    for (size_t i = 0; i < tdbr->getSize(); i++) maxKey = std::max(maxKey, tdbr->getDbKey(i));
    std::vector<uint32_t> tIdOf(size_t(maxKey) + 1, UINT32_MAX);
    for (size_t i = 0; i < tdbr->getSize(); i++) tIdOf[tdbr->getDbKey(i)] = uint32_t(i);
    auto tid = [&](const char *key) -> size_t {
        const unsigned long k = strtoul(key, NULL, 10);
        return (k <= maxKey && tIdOf[k] != UINT32_MAX) ? size_t(tIdOf[k]) : size_t(UINT_MAX);
    };
    // sizes per query entry: letters (query + all its targets) and pairs
    const size_t nQ = pref.getSize();
    std::vector<uint64_t> qLetters(nQ, 0);
    std::vector<uint32_t> qPairs(nQ, 0);
#pragma omp parallel for schedule(dynamic, 1024) num_threads(par.threads)
    for (size_t k = 0; k < nQ; k++) {
        unsigned int thread = 0;
#ifdef OPENMP
        thread = (unsigned int) omp_get_thread_num();
#endif
        char *data = pref.getData(k, thread);
        if (*data == '\0') continue;
        const unsigned int qKey = pref.getDbKey(k);
        const size_t qId = qdbr.getId(qKey);
        if (qId == UINT_MAX) { Debug(Debug::ERROR) << "query " << qKey << " not in the query db\n"; EXIT(EXIT_FAILURE); }
        uint64_t nl = qdbr.getSeqLen(qId);
        uint32_t np = 0;
        char key[256];
        while (*data != '\0') {
            Util::parseKey(data, key);
            const size_t tId = tid(key);
            if (tId == UINT_MAX) { Debug(Debug::ERROR) << "target " << key << " not in the target db\n"; EXIT(EXIT_FAILURE); }
            nl += tdbr->getSeqLen(tId); np++;
            data = Util::skipLine(data);
        }
        qLetters[k] = nl; qPairs[k] = np;
    }

    const double secSizes = tPhase.getTimediff();
    double secFill = 0, secWrite = 0, secWait = 0;
    struct Batch {
        size_t first = 0, n = 0;
        std::vector<size_t> letterOff, pairOff;
        std::vector<uint8_t> letters;
        std::vector<int8_t> bias;
        std::vector<PairSW::Pair> pairs;
        std::vector<int32_t> scores;
    };
    Batch buf[2];
    for (Batch &B : buf) { B.letters.resize(maxLetters); B.bias.resize(maxLetters); B.pairs.resize(maxPairs); B.scores.resize(maxPairs); }

    // fill: raw query / target bytes (mapped on the GPU), query composition bias from the mapped query
    auto fill = [&](Batch &B) {
        B.letterOff.assign(B.n + 1, 0); B.pairOff.assign(B.n + 1, 0);
        for (size_t k = 0; k < B.n; k++) {
            B.letterOff[k + 1] = B.letterOff[k] + qLetters[B.first + k];
            B.pairOff[k + 1] = B.pairOff[k] + qPairs[B.first + k];
        }
#pragma omp parallel num_threads(par.threads)
        {
            unsigned int thread = 0;
#ifdef OPENMP
            thread = (unsigned int) omp_get_thread_num();
#endif
            Sequence qSeq(par.maxSeqLen, qtype, &m, 0, false, par.compBiasCorrection);
            std::vector<float> cb(par.maxSeqLen + 1);
            char key[256];
#pragma omp for schedule(dynamic, 64)
            for (size_t k = 0; k < B.n; k++) {
                if (qPairs[B.first + k] == 0) continue;
                char *data = pref.getData(B.first + k, thread);
                const unsigned int qKey = pref.getDbKey(B.first + k);
                const size_t qId = qdbr.getId(qKey);
                const char *qData = qdbr.getData(qId, thread);
                const unsigned int qLen = qdbr.getSeqLen(qId);
                size_t o = B.letterOff[k];
                const size_t qOff = o;
                memcpy(&B.letters[o], qData, qLen);
                if (par.compBiasCorrection) {
                    qSeq.mapSequence(qId, qKey, qData, qLen);
                    SubstitutionMatrix::calcLocalAaBiasCorrection(&m, qSeq.numSequence, qSeq.L, cb.data(), par.compBiasCorrectionScale);
                    for (unsigned int i = 0; i < qLen; i++) {
                        const int8_t r = (int8_t) ((cb[i] < 0.0) ? cb[i] - 0.5 : cb[i] + 0.5);
                        B.bias[o + i] = biasUpper ? std::max<int8_t>(0, r) : r;
                    }
                } else {
                    memset(&B.bias[o], 0, qLen);
                }
                o += qLen;
                size_t p = B.pairOff[k];
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    const size_t tId = tid(key);
                    const unsigned int tLen = tdbr->getSeqLen(tId);
                    memcpy(&B.letters[o], tdbr->getData(tId, thread), tLen);
                    B.pairs[p++] = {qOff, o, qLen, tLen};
                    o += tLen;
                    data = Util::skipLine(data);
                }
            }
        }
    };
    size_t totalPairs = 0, keptPairs = 0;
    double cells = 0, gpuSec = 0;
    Debug::Progress progress(nQ);
    // write the kept lines of a scored batch
    auto write = [&](Batch &B) {
        size_t kept = 0;
        double c = 0;
#pragma omp parallel num_threads(par.threads) reduction(+ : kept, c)
        {
            unsigned int thread = 0;
#ifdef OPENMP
            thread = (unsigned int) omp_get_thread_num();
#endif
            std::string out;
            char key[256];
#pragma omp for schedule(dynamic, 64)
            for (size_t k = 0; k < B.n; k++) {
                const unsigned int qKey = pref.getDbKey(B.first + k);
                char *data = pref.getData(B.first + k, thread);
                out.clear();
                size_t p = B.pairOff[k];
                while (*data != '\0') {
                    const char *line = data;
                    Util::parseKey(data, key);
                    const unsigned int tKey = (unsigned int) strtoul(key, NULL, 10);
                    data = Util::skipLine(data);
                    const PairSW::Pair &P = B.pairs[p];
                    c += double(P.qLen) * double(P.tLen);
                    const bool keep = (tKey == qKey) || evaluer.computeEvalue(B.scores[p], P.qLen) <= evalKeep;
                    p++;
                    if (keep) { out.append(line, data - line); kept++; }
                }
                writer.writeData(out.c_str(), out.size(), qKey, thread);
                progress.updateProgress();
            }
        }
        totalPairs += B.pairOff[B.n]; keptPairs += kept; cells += c;
    };
    // pipeline: fill batch i while the GPU scores batch i - 1, then write batch i - 1 while it scores batch i
    std::thread gpu;
    bool pending = false;
    int cur = 0;
    size_t b0 = 0;
    while (b0 < nQ) {
        size_t b1 = b0, nl = 0, np = 0;
        while (b1 < nQ && (b1 == b0 || (nl + qLetters[b1] <= maxLetters && np + qPairs[b1] <= maxPairs))) {
            nl += qLetters[b1]; np += qPairs[b1]; b1++;
        }
        if (nl > maxLetters || np > maxPairs) {
            Debug(Debug::ERROR) << "gpupairfilter: query " << pref.getDbKey(b0) << " alone exceeds a batch\n";
            EXIT(EXIT_FAILURE);
        }
        Batch &B = buf[cur];
        B.first = b0; B.n = b1 - b0;
        { Timer t; fill(B); secFill += t.getTimediff(); }
        { Timer t; if (pending) gpu.join(); secWait += t.getTimediff(); }
        gpu = std::thread([&sw, &B, &gpuSec]() {
            // CUDA errors end the process inside PairSW (this build has no exceptions)
            Timer gt;
            sw.score(B.letters.data(), B.bias.data(), B.letterOff[B.n], B.pairs.data(), B.pairOff[B.n], B.scores.data());
            gpuSec += gt.getTimediff();
        });
        if (pending) { Timer t; write(buf[1 - cur]); secWrite += t.getTimediff(); }
        pending = true; cur = 1 - cur; b0 = b1;
    }
    if (pending) {
        gpu.join();
        write(buf[1 - cur]);
    }
    writer.close();
    Debug(Debug::INFO) << totalPairs << " pairs, " << keptPairs << " kept (E <= " << evalKeep << " or identity), "
                       << "GPU " << gpuSec << " s, sizes " << secSizes << " s, fill " << secFill << " s, wait " << secWait << " s, write " << secWrite << " s, " << (gpuSec > 0 ? cells / gpuSec / 1e9 : 0) << " GCUPS\n";
    pref.close();
    if (!same) { tdbr->close(); delete tdbr; }
    qdbr.close();
    return EXIT_SUCCESS;
#endif
}
