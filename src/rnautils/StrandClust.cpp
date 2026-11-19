#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Matcher.h"
#include "Util.h"
#include "itoa.h"

#include <algorithm>
#include <climits>
#include <unordered_map>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

// riboseek cluster on both strands (cascaded_clustering.sh with STRAND_DB): every step aligns both strands of its
// sequences (dinucdb --strand 2: key 2k forward, 2k+1 reverse complement) against their forward strands (keys k).
// mergestrandaln turns that into one alignment list per sequence for clust, addclustrand writes each member's strand
// relative to its representative into the cluster db ("key\t1": reverse complement), which mergeclusters carries
// through the cascade.

// the better of two hits of one target: higher score, then lower E-value (the forward strand wins ties, it is read first)
static bool betterStrandHit(const Matcher::result_t &a, const Matcher::result_t &b) {
    if (a.score != b.score) {
        return a.score > b.score;
    }
    return a.eval < b.eval;
}

// mergestrandaln <i:sequenceDB> <i:strandAlnDB> <o:alnDB>: per sequence k, the hits of queries 2k and 2k+1 merged,
// one per target (the better strand), E-values doubled (both strands searched) and filtered by -e. Reverse strand hits
// get their query coordinates on the forward strand (start > end) and a reversed backtrace. The self hit comes first
// (added when missing, so clust sees every sequence), the rest are sorted as align sorts them.
int mergestrandaln(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> seqDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    seqDb.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> alnDb(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                 DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    alnDb.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    Debug::Progress progress(seqDb.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::vector<Matcher::result_t> hits;
        std::unordered_map<unsigned int, size_t> byTarget;
        std::string result;
        std::vector<char> buffer(4096);

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < seqDb.getSize(); ++i) {
            progress.updateProgress();
            const unsigned int key = seqDb.getDbKey(i);
            hits.clear();
            byTarget.clear();
            for (unsigned int reverse = 0; reverse < 2; ++reverse) {
                const size_t id = alnDb.getId(2 * key + reverse);
                if (id >= alnDb.getSize()) { // missing key: UINT_MAX or DB_ENTRY_NOT_FOUND, depending on the MMseqs2 version
                    continue;
                }
                char *data = alnDb.getData(id, thread_idx);
                while (*data != '\0') {
                    Matcher::result_t hit = Matcher::parseAlignmentRecord(data, false);
                    data = Util::skipLine(data);
                    hit.eval *= 2.0;
                    if (hit.eval > par.evalThr || (reverse == 1 && hit.dbKey == key)) {
                        continue;
                    }
                    if (reverse == 1) {
                        hit.qStartPos = static_cast<int>(hit.qLen) - 1 - hit.qStartPos;
                        hit.qEndPos = static_cast<int>(hit.qLen) - 1 - hit.qEndPos;
                        std::reverse(hit.backtrace.begin(), hit.backtrace.end());
                    }
                    std::unordered_map<unsigned int, size_t>::iterator it = byTarget.find(hit.dbKey);
                    if (it == byTarget.end()) {
                        byTarget.emplace(hit.dbKey, hits.size());
                        hits.push_back(hit);
                    } else if (betterStrandHit(hit, hits[it->second])) {
                        hits[it->second] = hit;
                    }
                }
            }

            int maxScore = 1;
            for (size_t j = 0; j < hits.size(); ++j) {
                maxScore = std::max(maxScore, hits[j].score);
            }
            if (byTarget.find(key) == byTarget.end()) {
                const unsigned int len = static_cast<unsigned int>(seqDb.getSeqLen(i));
                const int last = static_cast<int>(len) - 1;
                hits.emplace_back(key, maxScore, 1.0f, 1.0f, 1.0f, 0.0, len, 0, last, len, 0, last, len,
                                  -1, -1, -1, -1, "");
            }
            std::sort(hits.begin(), hits.end(), [key](const Matcher::result_t &a, const Matcher::result_t &b) {
                if ((a.dbKey == key) != (b.dbKey == key)) {
                    return a.dbKey == key;
                }
                return Matcher::compareHits(a, b);
            });

            result.clear();
            for (size_t j = 0; j < hits.size(); ++j) {
                const bool backtrace = hits[j].backtrace.empty() == false;
                if (buffer.size() < 1024 + 2 * hits[j].backtrace.size()) {
                    buffer.resize(1024 + 2 * hits[j].backtrace.size());
                }
                const size_t len = Matcher::resultToBuffer(buffer.data(), hits[j], backtrace, true);
                result.append(buffer.data(), len);
            }
            writer.writeData(result.c_str(), result.length(), key, thread_idx);
        }
    }
    writer.close();
    alnDb.close();
    seqDb.close();
    return EXIT_SUCCESS;
}

// strand of the best alignment of query (both strand entries) to target in the strand alignment db: 0 forward,
// 1 reverse complement, -1 none
static int strandOf(DBReader<unsigned int> &alnDb, unsigned int query, unsigned int target, int thread_idx) {
    int strand = -1;
    int bestScore = INT_MIN;
    double bestEval = 0.0;
    for (unsigned int reverse = 0; reverse < 2; ++reverse) {
        const size_t id = alnDb.getId(2 * query + reverse);
        if (id >= alnDb.getSize()) { // missing key: UINT_MAX or DB_ENTRY_NOT_FOUND, depending on the MMseqs2 version
            continue;
        }
        char *data = alnDb.getData(id, thread_idx);
        char keyBuffer[255];
        while (*data != '\0') {
            Util::parseKey(data, keyBuffer);
            if (Util::fast_atoi<unsigned int>(keyBuffer) == target) {
                Matcher::result_t hit = Matcher::parseAlignmentRecord(data, true);
                if (strand == -1 || hit.score > bestScore || (hit.score == bestScore && hit.eval < bestEval)) {
                    strand = static_cast<int>(reverse);
                    bestScore = hit.score;
                    bestEval = hit.eval;
                }
            }
            data = Util::skipLine(data);
        }
    }
    return strand;
}

// addclustrand <i:clusterDB> <i:strandAlnDB> <o:clusterDB>: each member line gets its strand relative to the
// representative, from the representative's alignments to it (both strands), else from the member's to the representative
int addclustrand(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> cluDb(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                 DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    cluDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBReader<unsigned int> alnDb(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                 DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    alnDb.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_CLUSTER_RES);
    writer.open();

    size_t missing = 0;
    Debug::Progress progress(cluDb.getSize());
#pragma omp parallel reduction(+:missing)
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string result;
        char keyBuffer[255];
        char numBuffer[32];

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < cluDb.getSize(); ++i) {
            progress.updateProgress();
            const unsigned int rep = cluDb.getDbKey(i);
            char *data = cluDb.getData(i, thread_idx);
            result.clear();
            while (*data != '\0') {
                Util::parseKey(data, keyBuffer);
                const unsigned int member = Util::fast_atoi<unsigned int>(keyBuffer);
                data = Util::skipLine(data);
                int strand = 0;
                if (member != rep) {
                    strand = strandOf(alnDb, rep, member, thread_idx);
                    if (strand == -1) {
                        strand = strandOf(alnDb, member, rep, thread_idx);
                    }
                    if (strand == -1) {
                        missing++;
                        strand = 0;
                    }
                }
                char *end = Itoa::u32toa_sse2(member, numBuffer);
                result.append(numBuffer, end - numBuffer - 1);
                result.append(strand == 1 ? "\t1\n" : "\t0\n");
            }
            writer.writeData(result.c_str(), result.length(), rep, thread_idx);
        }
    }
    writer.close();
    alnDb.close();
    cluDb.close();
    if (missing > 0) {
        Debug(Debug::WARNING) << missing << " members without an alignment to their representative, kept as forward\n";
    }
    return EXIT_SUCCESS;
}

// splitclustrand <i:clusterDB> <o:clusterDB>: a cluster db as align input on the both-strand db, every member of
// representative r listed under both of r's strand entries 2r and 2r+1 (reassignment aligns each member to its
// representative in both orientations and keeps the better one, mergestrandaln, instead of trusting the strand the
// cascade composed through its intermediate representatives)
int splitclustrand(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> cluDb(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                 DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    cluDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_CLUSTER_RES);
    writer.open();

    Debug::Progress progress(cluDb.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string members;
        char keyBuffer[255];

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < cluDb.getSize(); ++i) {
            progress.updateProgress();
            const unsigned int rep = cluDb.getDbKey(i);
            char *data = cluDb.getData(i, thread_idx);
            members.clear();
            while (*data != '\0') {
                Util::parseKey(data, keyBuffer);
                members.append(keyBuffer);
                members.push_back('\n');
                data = Util::skipLine(data);
            }
            writer.writeData(members.c_str(), members.length(), 2 * rep, thread_idx);
            writer.writeData(members.c_str(), members.length(), 2 * rep + 1, thread_idx);
        }
    }
    writer.close();
    cluDb.close();
    return EXIT_SUCCESS;
}
