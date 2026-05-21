#include "LocalCommandDeclarations.h"
#include "CommandDeclarations.h"
#include "Debug.h"
#include "DBReader.h"
#include "MMseqsMPI.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "RNAFoldBridge.h"
#include "Matcher.h"
#include "Orf.h"
#include "LocalParameters.h"
#include "MsaFilter.h"
#include "MultipleAlignment.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "DinucleotideMapping.h"
#ifdef HAVE_INFERNAL_BRIDGE
#include "infernal/InfernalBridge.h"
#endif

#ifdef OPENMP
#include <omp.h>
#endif

#include <cctype>
#include <cstdlib>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Native generatecm command:
// maps DB -> DB covariance-model generation to MMseqs profile construction.
int generatecm(int argc, const char **argv, const Command &command) {
    return result2profile(argc, argv, command);
}

namespace {

struct AlnSeq {
    std::string id;
    std::string aln;
};

static inline char normalizeBase(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (c == 'T') {
        return 'U';
    }
    return c;
}

static int effectiveDecodeSeqType(int dbtype, bool useDinucMapping) {
    if (!useDinucMapping) {
        return dbtype;
    }
    unsigned int ext = DBReader<unsigned int>::getExtendedDbtype(dbtype);
    ext |= Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE)) {
        return DBReader<unsigned int>::setExtendedDbtype(dbtype, ext);
    }
    return DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS, ext);
}

static std::string decodeMappedSequenceToRna(const Sequence &seq,
                                             const SubstitutionMatrix &subMat,
                                             const unsigned char *num2outputnum) {
    std::string out;
    out.reserve(static_cast<size_t>(seq.L));
    for (int i = 0; i < seq.L; ++i) {
        unsigned char code = seq.numSequence[i];
        if (num2outputnum != NULL) {
            code = num2outputnum[code];
        }
        char c = normalizeBase(subMat.num2aa[code]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'U') {
            c = 'N';
        }
        out.push_back(c);
    }
    return out;
}

static std::string decodeDbSequence(DBReader<unsigned int> &dbr,
                                    size_t id,
                                    unsigned int threadIdx,
                                    Sequence &mapper,
                                    bool useDinucMapping,
                                    bool isGpuDb,
                                    const SubstitutionMatrix &subMat) {
    const unsigned int seqLen = dbr.getSeqLen(id);
    if (isGpuDb) {
        const unsigned char *data =
            reinterpret_cast<const unsigned char *>(dbr.getDataUncompressed(id));
        mapper.mapSequence(id, dbr.getDbKey(id), std::make_pair(data, seqLen));
    } else {
        const char *data = dbr.getData(id, threadIdx);
        mapper.mapSequence(id, dbr.getDbKey(id), data, seqLen);
    }

    const Sequence::SeqAuxInfo *auxInfo = Sequence::getAuxInfo(mapper.getSeqType());
    const unsigned char *num2outputnum = (useDinucMapping && auxInfo != NULL)
        ? auxInfo->num2outputnum
        : NULL;
    return decodeMappedSequenceToRna(mapper, subMat, num2outputnum);
}

static inline void convertTsToUs(std::string &seq) {
    for (size_t i = 0; i < seq.size(); ++i) {
        if (seq[i] == 'T') seq[i] = 'U';
        else if (seq[i] == 't') seq[i] = 'u';
    }
}

// Encode an aligned ACGU/N/'-' row into the byte format MsaFilter expects.
// Residues map to 0..3 (so they pass `< NAA` and count toward coverage/identity);
// 'N' maps to NAA (any-residue sentinel); '-' maps to GAP. Encoding doesn't use
// subMat.aa2num — MsaFilter's identity comparisons are pure byte-equality, and
// the substitution matrix is only consulted for qsc, which we leave at the
// default disabled value when --filter-msa is the rMSA-style cov+id+Ndiff path.
static inline unsigned char encodeAlignedNuc(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'U': case 'T': return 3;
        case '-': return MultipleAlignment::GAP;
        default:  return MultipleAlignment::NAA;  // N or anything else → ANY
    }
}

static std::string buildStockholmText(const std::string &id,
                                      const std::vector<AlnSeq> &seqs,
                                      const std::string &ssCons) {
    std::ostringstream out;
    out << "# STOCKHOLM 1.0\n";
    out << "#=GF ID " << (id.empty() ? "mmseqs_model" : id) << "\n";
    for (size_t i = 0; i < seqs.size(); ++i) {
        std::string rnaAln = seqs[i].aln;
        convertTsToUs(rnaAln);
        out << seqs[i].id << " " << rnaAln << "\n";
    }
    // Force every alignment column to be a CM match column via --hand: RF line
    // of all 'x' makes consensus column c map 1:1 to query position c-1, so
    // cmsearch traces are directly query-coord and result_t records stay
    // compatible with result2dnamsa / convertalis.
    if (!seqs.empty() && !seqs[0].aln.empty()) {
        out << "#=GC RF " << std::string(seqs[0].aln.size(), 'x') << "\n";
    }
    if (!ssCons.empty()) {
        out << "#=GC SS_cons " << ssCons << "\n";
    }
    out << "//\n";
    return out.str();
}

} // namespace

int cmbuild(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    // Defaults for the optional MSA row filter (off; rMSA-style cov+id+Ndiff
    // when the user opts in via --filter-msa 1). qsc disabled because our
    // substitution matrix is dinucleotide so qsc scores would be in the wrong
    // space; cov + filterMaxSeqId + Ndiff are pure byte comparisons and unaffected.
    par.filterMsa = 0;
    par.qsc = -50.0f;
    par.covMSAThr = 0.5f;
    par.filterMaxSeqId = 0.99f;
    par.Ndiff = 128;
    par.parseParameters(argc, argv, command, true, 0, 0);

#ifndef HAVE_INFERNAL_BRIDGE
    Debug(Debug::ERROR) << "cmbuild requires Infernal bridge support\n";
    return EXIT_FAILURE;
#else

    // Fork Infernal workers BEFORE loading DBs. Each worker stays ~small
    // (few MB) instead of CoW-ing a post-mmap fat parent, which serializes
    // forks under the kernel mm lock.
    InfernalBridge::WorkerPool *workerPool = InfernalBridge::startWorkerPool(par.threads);
    if (workerPool == NULL) {
        Debug(Debug::ERROR) << "cmbuild: failed to start Infernal worker pool\n";
        return EXIT_FAILURE;
    }

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    qDbr.open(DBReader<unsigned int>::NOSORT);

    const bool sameDatabase = (par.db1 == par.db2);
    DBReader<unsigned int> *tDbr = &qDbr;
    if (!sameDatabase) {
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                          DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter cmWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads,
                      par.compressed, Parameters::DBTYPE_GENERIC_DB);
    cmWriter.open();

    Debug(Debug::INFO) << "Query database size: " << qDbr.getSize() << " type: " << qDbr.getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << tDbr->getDbTypeName() << "\n";
    {
        std::ostringstream msaEvalStream;
        msaEvalStream << std::scientific << par.cmliteMsaEvalThr;
        Debug(Debug::INFO) << "cmbuild seed E-value threshold: " << msaEvalStream.str() << "\n";
    }

    Debug::Progress progress(resultReader.getSize());

    // Stockholm rows are filtered: drop any target with non-gap coverage
    // below this threshold. Infernal's cm_parsetree_Doctor() can fail when
    // many short fragments dominate the column statistics.
    const float minColCoverage = 0.30f;
    // Stricter filter for the alifold input: noisy rows degrade the
    // covariation signal, so feed only well-covered rows to alifold
    // (default 0.70, override with RIBOSEEK_ALIFOLD_MINCOV).
    float alifoldMinCov = 0.70f;
    if (const char *envCov = std::getenv("RIBOSEEK_ALIFOLD_MINCOV")) {
        float v = std::atof(envCov);
        if (v > 0.0f && v <= 1.0f) alifoldMinCov = v;
    }

    // Optional rMSA/hhfilter-style row filter on the cmbuild input MSA. Knobs
    // map 1:1 to MMseqs MsaFilter: --cov-msa, --filter-max-seqid, --diff (Ndiff),
    // --qid, --qsc, --filter-min-enable. Default is off.
    const bool doMsaFilter = (par.filterMsa != 0);
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0f);
    const unsigned int targetExt = DBReader<unsigned int>::getExtendedDbtype(tDbr->getDbtype());
    const bool targetGpuDb = (targetExt & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    const bool decodeTargetDinuc = ((targetExt & Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0) || targetGpuDb;
    const int targetSeqType = effectiveDecodeSeqType(tDbr->getDbtype(), decodeTargetDinuc);
    std::vector<std::string> qid_str_vec = Util::split(par.qid, ",");
    std::vector<int> qid_vec;
    qid_vec.reserve(qid_str_vec.size());
    for (size_t i = 0; i < qid_str_vec.size(); ++i) {
        float qidf = static_cast<float>(std::atof(qid_str_vec[i].c_str()));
        qid_vec.push_back(static_cast<int>(qidf * 100));
    }
    if (qid_vec.empty()) qid_vec.push_back(0);
    std::sort(qid_vec.begin(), qid_vec.end());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);
        std::set<unsigned int> seenKeys;

        // Per-thread MsaFilter scratch. Allocated once, reused across queries.
        // maxSeqLen and maxSetSize grow inside MsaFilter as needed (increaseSetSize),
        // so seeding with conservative values is fine.
        MsaFilter msaFilter(8192, 1024, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
        Sequence targetMapper(tDbr->getMaxSeqLen(), targetSeqType, &subMat, 0, false, false);
        std::vector<unsigned char> encBuf;
        std::vector<const char *> rowPtrs;

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            alnResults.clear();
            seenKeys.clear();

            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr.getId(queryKey);
            if (queryId == UINT_MAX) {
                Debug(Debug::WARNING) << "cmbuild: invalid query " << queryKey << "\n";
                continue;
            }

            // Render the query row directly from raw nucleotide DB chars.
            // The MSA stored in MultipleAlignment::MSAResult uses the
            // dinucleotide alphabet (24-symbol), so num2aa[]-based decoding
            // collapses most positions to 'N' — fatal for alifold/Infernal.
            // Mirror result2dnamsa.cpp instead: emit raw chars from the DB
            // and walk each backtrace in query coordinates.
            const char *querySeq = qDbr.getData(queryId, thread_idx);
            const size_t qLen = qDbr.getSeqLen(queryId);
            std::string queryRow;
            queryRow.reserve(qLen);
            for (size_t pos = 0; pos < qLen; ++pos) {
                char c = normalizeBase(querySeq[pos]);
                if (c != 'A' && c != 'C' && c != 'G' && c != 'U') c = 'N';
                queryRow.push_back(c);
            }

            Matcher::readAlignmentResults(alnResults, resultReader.getData(id, thread_idx), false);

            std::vector<AlnSeq> stoSeqs;
            std::vector<std::string> alifoldRows;
            stoSeqs.reserve(alnResults.size() + 1);
            alifoldRows.reserve(alnResults.size() + 1);

            AlnSeq queryAs;
            queryAs.id = "query_" + std::to_string(queryKey);
            queryAs.aln = queryRow;
            stoSeqs.push_back(queryAs);
            alifoldRows.push_back(queryRow);

            for (size_t i = 0; i < alnResults.size(); i++) {
                Matcher::result_t res = alnResults[i];
                if (res.backtrace.empty()) continue;
                if (res.eval > par.cmliteMsaEvalThr) continue;
                const size_t targetId = tDbr->getId(res.dbKey);
                if (targetId == UINT_MAX) continue;
                // Infernal requires unique sequence names; skip duplicates.
                if (seenKeys.find(res.dbKey) != seenKeys.end()) continue;
                seenKeys.insert(res.dbKey);

                std::string decodedTargetSeq = decodeDbSequence(*tDbr, targetId, thread_idx,
                                                                targetMapper, decodeTargetDinuc,
                                                                targetGpuDb, subMat);
                const char *targetSeq = decodedTargetSeq.c_str();

                // Fold reverse-strand alignments onto the forward query frame
                // (matches result2dnamsa.cpp). Strand handling lives here so
                // splitstrand-produced reverse hits land in the same column
                // space as forward hits.
                bool queryIsReversed = (res.qStartPos > res.qEndPos);
                bool targetIsReversed = (res.dbStartPos > res.dbEndPos);
                bool isReverseStrand = false;
                if (queryIsReversed && targetIsReversed) {
                    std::swap(res.dbStartPos, res.dbEndPos);
                    std::reverse(res.backtrace.begin(), res.backtrace.end());
                } else if (queryIsReversed && !targetIsReversed) {
                    isReverseStrand = true;
                    std::swap(res.dbStartPos, res.dbEndPos);
                    std::reverse(res.backtrace.begin(), res.backtrace.end());
                } else if (!queryIsReversed && targetIsReversed) {
                    isReverseStrand = true;
                }

                int qStartPos = std::min(res.qStartPos, res.qEndPos);
                int qEndPos = std::max(res.qStartPos, res.qEndPos);

                std::string row;
                row.reserve(qLen);
                for (int pos = 0; pos < qStartPos; ++pos) row.push_back('-');
                size_t residueCount = 0;
                unsigned int seqPos = 0;
                for (size_t bp = 0; bp < res.backtrace.size(); ++bp) {
                    char rawChar = isReverseStrand
                        ? Orf::complement(targetSeq[res.dbStartPos - seqPos])
                        : targetSeq[res.dbStartPos + seqPos];
                    switch (res.backtrace[bp]) {
                        case 'M': {
                            char c = normalizeBase(rawChar);
                            if (c != 'A' && c != 'C' && c != 'G' && c != 'U') c = 'N';
                            row.push_back(c);
                            residueCount++;
                            seqPos++;
                            break;
                        }
                        case 'I':
                            row.push_back('-');
                            break;
                        case 'D':
                            seqPos++;
                            break;
                    }
                }
                for (int pos = qEndPos + 1; pos < (int)res.qLen; ++pos) row.push_back('-');

                // Stockholm requires equal-length rows; res.qLen should match
                // qLen but pad/truncate defensively in case a stale result
                // disagrees with the current query DB.
                if (row.size() < qLen) row.append(qLen - row.size(), '-');
                else if (row.size() > qLen) row.resize(qLen);

                const float cov = qLen > 0
                    ? static_cast<float>(residueCount) / static_cast<float>(qLen)
                    : 0.0f;
                AlnSeq as;
                as.id = "t" + std::to_string(res.dbKey);
                as.aln = std::move(row);
                if (cov >= minColCoverage) {
                    if (cov >= alifoldMinCov) alifoldRows.push_back(as.aln);
                    stoSeqs.push_back(std::move(as));
                } else if (cov >= alifoldMinCov) {
                    alifoldRows.push_back(std::move(as.aln));
                }
            }

            // Consensus SS via alifold over high-coverage rows; falls back to
            // single-seq fold on row 0 if only the query passed the filter.
            std::string ssCons;
            double mfe = 0.0;
            bool ssOk = rnaFoldAlifoldDotBracket(alifoldRows, ssCons, &mfe);
            (void)mfe;
            if (std::getenv("RIBOSEEK_DEBUG_SSCONS") != nullptr) {
                Debug(Debug::WARNING) << "[ssCons] q=" << queryKey
                    << " alifoldRows=" << alifoldRows.size()
                    << " ok=" << (ssOk ? 1 : 0)
                    << " len=" << ssCons.size()
                    << " ss=[" << ssCons << "]\n";
            }

            // Optional MsaFilter pass: rMSA/hhfilter-style row pruning before
            // cmbuild. Encodes ACGU/N/'-' rows into MsaFilter's byte alphabet
            // (residues < NAA, N=NAA, '-'=GAP) so identity comparisons land in
            // pure byte-equality and qsc/qid stay disabled (we set qsc=-50 and
            // qid="0" by default since the loaded subMat is dinucleotide).
            if (doMsaFilter && stoSeqs.size() > 1) {
                const int N = static_cast<int>(stoSeqs.size());
                const int L = static_cast<int>(qLen);
                encBuf.assign(static_cast<size_t>(N) * static_cast<size_t>(L), MultipleAlignment::GAP);
                rowPtrs.assign(N, nullptr);
                for (int i = 0; i < N; ++i) {
                    unsigned char *dst = encBuf.data() + static_cast<size_t>(i) * static_cast<size_t>(L);
                    const std::string &aln = stoSeqs[i].aln;
                    const int copyLen = std::min(L, static_cast<int>(aln.size()));
                    for (int p = 0; p < copyLen; ++p) dst[p] = encodeAlignedNuc(aln[p]);
                    rowPtrs[i] = reinterpret_cast<const char *>(dst);
                }
                msaFilter.filter(N, L,
                                 static_cast<int>(par.covMSAThr * 100),
                                 qid_vec, par.qsc,
                                 static_cast<int>(par.filterMaxSeqId * 100),
                                 par.Ndiff, par.filterMinEnable,
                                 rowPtrs.data(), false);
                bool *keptRaw = new bool[N];
                msaFilter.getKept(keptRaw, N);
                std::vector<AlnSeq> keptSeqs;
                keptSeqs.reserve(N);
                for (int i = 0; i < N; ++i) {
                    if (keptRaw[i]) keptSeqs.push_back(std::move(stoSeqs[i]));
                }
                delete[] keptRaw;
                stoSeqs = std::move(keptSeqs);
            }

            std::string stoText = buildStockholmText(
                "query_" + std::to_string(queryKey), stoSeqs, ssCons);

            if (const char *dumpPath = std::getenv("RIBOSEEK_DUMP_STOCKHOLM")) {
                std::string outPath = std::string(dumpPath) + "/query_" + std::to_string(queryKey) + ".sto";
                FILE *fp = std::fopen(outPath.c_str(), "w");
                if (fp != nullptr) {
                    std::fwrite(stoText.data(), 1, stoText.size(), fp);
                    std::fclose(fp);
                }
            }

            // Call Infernal bridge to build CM. Bridge forks per call so OMP workers
            // can build in parallel despite Infernal's process-wide global state.
            std::string cmText;
            std::string infernalErr;
            bool success = InfernalBridge::buildCmFromStockholmText(workerPool, stoText, cmText, infernalErr, par.calibrateCm);

            if (success) {
                cmWriter.writeData(cmText.c_str(), cmText.size(), queryKey, thread_idx);
            } else {
                Debug(Debug::WARNING) << "cmbuild: Infernal failed for query " << queryKey
                                      << ": " << infernalErr << "\n";
            }
        }
    }

    cmWriter.close(true);
    resultReader.close();
    if (!sameDatabase) {
        tDbr->close();
        delete tDbr;
    }
    qDbr.close();

    InfernalBridge::stopWorkerPool(workerPool);
    return EXIT_SUCCESS;
#endif // HAVE_INFERNAL_BRIDGE
}
