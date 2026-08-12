#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "CmBuildScan.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "RnaFoldingBridge.h"
#include "Matcher.h"
#include "Orf.h"
#include "MsaFilter.h"
#include "MultipleAlignment.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "DinucleotideMapping.h"

#include <set>
#include <unordered_map>
#ifdef RIBOSEEK_CMBUILD_DUMP_STOCKHOLM
#include <fstream>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

extern bool cmbuildFromAlignment(
    const std::string &name,
    const std::vector<std::string> &sqnames,
    const std::vector<std::string> &aseqs,
    const std::string &ssCons,
    std::string &cmText, std::string &err,
    double ereOverride, double symfracOpt, bool noss,
    bool forScan, bool useInside, bool useLocal,
    CM_t **ret_cm
);
extern void cmInfernalGlobalInit();

extern int result2profile(int argc, const char **argv, const Command& command);
int generatecm(int argc, const char **argv, const Command &command) {
    return result2profile(argc, argv, command);
}

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
        if (c != "A"[0] && c != "C"[0] && c != "G"[0] && c != "U"[0]) {
            c = "N"[0];
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
        if (seq[i] == "T"[0]) seq[i] = "U"[0];
        else if (seq[i] == "t"[0]) seq[i] = "u"[0];
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

static std::string trimCopy(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) {
        ++b;
    }
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) {
        --e;
    }
    return s.substr(b, e - b);
}

static std::vector<std::string> splitCommaList(const std::string &s) {
    std::vector<std::string> out;
    size_t pos = 0;
    while (pos <= s.size()) {
        const size_t next = s.find(',', pos);
        const std::string part = trimCopy(s.substr(pos, next - pos));
        if (!part.empty()) {
            out.push_back(part);
        }
        if (next == std::string::npos) {
            break;
        }
        pos = next + 1;
    }
    if (out.empty() && !trimCopy(s).empty()) {
        out.push_back(trimCopy(s));
    }
    return out;
}

static void unlinkIfExists(const std::string &path) {
    if (!path.empty() && FileUtil::fileExists(path.c_str())) {
        unlink(path.c_str());
    }
}

static void removeDbFiles(const std::string &db) {
    const std::vector<std::string> dataFiles = FileUtil::findDatafiles(db.c_str());
    for (size_t i = 0; i < dataFiles.size(); ++i) {
        unlinkIfExists(dataFiles[i]);
    }
    unlinkIfExists(db);
    unlinkIfExists(db + ".index");
    unlinkIfExists(db + ".dbtype");
    unlinkIfExists(db + ".lookup");
    const std::vector<std::string> headerDataFiles = FileUtil::findDatafiles((db + "_h").c_str());
    for (size_t i = 0; i < headerDataFiles.size(); ++i) {
        unlinkIfExists(headerDataFiles[i]);
    }
    unlinkIfExists(db + "_h");
    unlinkIfExists(db + "_h.index");
    unlinkIfExists(db + "_h.dbtype");
    unlinkIfExists(db + "_h.lookup");
}

static bool mergeCmbuildInputs(const std::string &queryDb,
                               const std::vector<std::string> &targetDbs,
                               const std::vector<std::string> &resultDbs,
                               const std::string &mergedTargetDb,
                               const std::string &mergedResultDb,
                               int threads,
                               std::string &error) {
    if (targetDbs.empty() || targetDbs.size() != resultDbs.size()) {
        error = "target/result DB list sizes do not match";
        return false;
    }

    DBReader<unsigned int> qDbr(queryDb.c_str(), (queryDb + ".index").c_str(),
                                std::max(1, threads),
                                DBReader<unsigned int>::USE_INDEX);
    qDbr.open(DBReader<unsigned int>::NOSORT);

    std::vector<DBReader<unsigned int> *> tReaders(targetDbs.size(), NULL);
    std::vector<DBReader<unsigned int> *> rReaders(resultDbs.size(), NULL);
    std::vector<std::unordered_map<unsigned int, unsigned int> > keyMaps(targetDbs.size());
    int mergedTargetDbtype = -1;
    bool qClosed = false;
    auto cleanupReaders = [&]() {
        for (size_t i = 0; i < tReaders.size(); ++i) {
            if (tReaders[i] != NULL) {
                tReaders[i]->close();
                delete tReaders[i];
                tReaders[i] = NULL;
            }
            if (rReaders[i] != NULL) {
                rReaders[i]->close();
                delete rReaders[i];
                rReaders[i] = NULL;
            }
        }
        if (!qClosed) {
            qDbr.close();
            qClosed = true;
        }
    };
#define CMBUILD_MERGE_FAIL(msg) do { error = (msg); cleanupReaders(); return false; } while (0)

        for (size_t i = 0; i < targetDbs.size(); ++i) {
            tReaders[i] = new DBReader<unsigned int>(targetDbs[i].c_str(),
                                                     (targetDbs[i] + ".index").c_str(),
                                                     std::max(1, threads),
                                                     DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
            tReaders[i]->open(DBReader<unsigned int>::NOSORT);
            rReaders[i] = new DBReader<unsigned int>(resultDbs[i].c_str(),
                                                     (resultDbs[i] + ".index").c_str(),
                                                     std::max(1, threads),
                                                     DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
            rReaders[i]->open(DBReader<unsigned int>::NOSORT);
            if (mergedTargetDbtype < 0) {
                mergedTargetDbtype = tReaders[i]->getDbtype();
            } else if (tReaders[i]->getDbtype() != mergedTargetDbtype) {
                error = "target DB types do not match across inputs";
                CMBUILD_MERGE_FAIL(error);
            }
        }

        FILE *targetIndexFile = fopen((mergedTargetDb + ".index").c_str(), "w");
        if (targetIndexFile == NULL) {
            error = "failed to create merged target index";
            CMBUILD_MERGE_FAIL(error);
        }
        unsigned int nextTargetKey = 0;
        size_t nextShardId = 0;
        size_t cumulativeDataSize = 0;
        for (size_t ti = 0; ti < tReaders.size(); ++ti) {
            DBReader<unsigned int> &tDbr = *tReaders[ti];
            std::unordered_map<unsigned int, unsigned int> &map = keyMaps[ti];
            map.reserve(tDbr.getSize() * 2 + 1);
            const std::vector<std::string> dataFiles = FileUtil::findDatafiles(targetDbs[ti].c_str());
            if (dataFiles.empty()) {
                error = "no datafiles found for target DB " + targetDbs[ti];
                CMBUILD_MERGE_FAIL(error);
            }
            size_t thisDbDataSize = 0;
            for (size_t fi = 0; fi < dataFiles.size(); ++fi) {
                const std::string realDataFile = FileUtil::getRealPathFromSymLink(dataFiles[fi]);
                struct stat st;
                if (stat(realDataFile.c_str(), &st) != 0) {
                    error = "failed to stat target data file " + realDataFile;
                    CMBUILD_MERGE_FAIL(error);
                }
                thisDbDataSize += static_cast<size_t>(st.st_size);
                const std::string dstShard = mergedTargetDb + "." + SSTR(nextShardId++);
                if (symlink(realDataFile.c_str(), dstShard.c_str()) != 0) {
                    error = "failed to symlink target data file " + realDataFile + " -> " + dstShard;
                    CMBUILD_MERGE_FAIL(error);
                }
            }
            for (size_t id = 0; id < tDbr.getSize(); ++id) {
                const unsigned int oldKey = tDbr.getDbKey(id);
                const size_t mergedOffset = cumulativeDataSize + tDbr.getOffset(id);
                const size_t len = tDbr.getEntryLen(id);
                fprintf(targetIndexFile, "%u\t%zu\t%zu\n", nextTargetKey, mergedOffset, len);
                map[oldKey] = nextTargetKey;
                ++nextTargetKey;
            }
            cumulativeDataSize += thisDbDataSize;
        }
        if (fclose(targetIndexFile) != 0) {
            error = "failed to close merged target index";
            CMBUILD_MERGE_FAIL(error);
        }
        {
            FILE *dbtypeFile = fopen((mergedTargetDb + ".dbtype").c_str(), "wb");
            if (dbtypeFile == NULL) {
                error = "failed to create merged target dbtype";
                CMBUILD_MERGE_FAIL(error);
            }
            if (fwrite(&mergedTargetDbtype, sizeof(int), 1, dbtypeFile) != 1) {
                fclose(dbtypeFile);
                error = "failed to write merged target dbtype";
                CMBUILD_MERGE_FAIL(error);
            }
            if (fclose(dbtypeFile) != 0) {
                error = "failed to close merged target dbtype";
                CMBUILD_MERGE_FAIL(error);
            }
        }

        FILE *lookupFilePtr = fopen((mergedTargetDb + ".lookup").c_str(), "w");
        if (lookupFilePtr == NULL) {
            error = "failed to create merged target lookup";
            CMBUILD_MERGE_FAIL(error);
        }
        DBWriter headerWriter((mergedTargetDb + "_h").c_str(),
                              (mergedTargetDb + "_h.index").c_str(),
                              1, 0, Parameters::DBTYPE_GENERIC_DB);
        headerWriter.open();
        unsigned int nextFileNumberOffset = 0;
        for (size_t ti = 0; ti < tReaders.size(); ++ti) {
            unsigned int localMaxFileNumber = 0;
            const std::string srcLookup = targetDbs[ti] + ".lookup";
            if (FileUtil::fileExists(srcLookup.c_str())) {
                DBReader<unsigned int> lookupReader(targetDbs[ti].c_str(),
                                                    (targetDbs[ti] + ".index").c_str(),
                                                    1,
                                                    DBReader<unsigned int>::USE_LOOKUP);
                lookupReader.open(DBReader<unsigned int>::NOSORT);
                DBReader<unsigned int>::LookupEntry *lookup = lookupReader.getLookup();
                for (size_t li = 0; li < lookupReader.getLookupSize(); ++li) {
                    const unsigned int oldKey = lookup[li].id;
                    std::unordered_map<unsigned int, unsigned int>::const_iterator it = keyMaps[ti].find(oldKey);
                    if (it == keyMaps[ti].end()) {
                        continue;
                    }
                    const unsigned int mergedKey = it->second;
                    const unsigned int mergedFileNumber = nextFileNumberOffset + lookup[li].fileNumber;
                    if (lookup[li].fileNumber > localMaxFileNumber) {
                        localMaxFileNumber = lookup[li].fileNumber;
                    }
                    fprintf(lookupFilePtr, "%u\t%s\t%u\n",
                            mergedKey,
                            lookup[li].entryName.c_str(),
                            mergedFileNumber);
		    const std::string headerName = lookup[li].entryName + "\n";
                    headerWriter.writeData(headerName.c_str(),
                                           headerName.size(),
                                           mergedKey,
                                           0);
                }
                lookupReader.close();
                nextFileNumberOffset += (localMaxFileNumber + 1);
            } else {
                DBReader<unsigned int> &tDbr = *tReaders[ti];
                for (size_t id = 0; id < tDbr.getSize(); ++id) {
                    const unsigned int oldKey = tDbr.getDbKey(id);
                    std::unordered_map<unsigned int, unsigned int>::const_iterator it = keyMaps[ti].find(oldKey);
                    if (it == keyMaps[ti].end()) {
                        continue;
                    }
                    const unsigned int mergedKey = it->second;
                    fprintf(lookupFilePtr, "%u\t%u\t%u\n", mergedKey, mergedKey, nextFileNumberOffset);
		    const std::string mergedName = SSTR(mergedKey) + "\n";
                    headerWriter.writeData(mergedName.c_str(), mergedName.size(), mergedKey, 0);
                }
                nextFileNumberOffset += 1;
            }
        }
        if (fclose(lookupFilePtr) != 0) {
            error = "failed to close merged target lookup";
            CMBUILD_MERGE_FAIL(error);
        }
        headerWriter.close(true);

        DBWriter resultWriter(mergedResultDb.c_str(), (mergedResultDb + ".index").c_str(),
                              1, 0, Parameters::DBTYPE_ALIGNMENT_RES);
        resultWriter.open();
        for (size_t qi = 0; qi < qDbr.getSize(); ++qi) {
            const unsigned int queryKey = qDbr.getDbKey(qi);
            std::string out;
            for (size_t ri = 0; ri < rReaders.size(); ++ri) {
                DBReader<unsigned int> &rDbr = *rReaders[ri];
                const size_t rid = rDbr.getId(queryKey);
                if (rid == UINT_MAX) {
                    continue;
                }
                char *block = rDbr.getData(rid, 0);
                const size_t blockLen = rDbr.getEntryLen(rid);
                size_t pos = 0;
                const size_t usable = (blockLen == 0 ? 0 : blockLen - 1);
                while (pos < usable) {
                    const size_t lineStart = pos;
                    while (pos < usable && block[pos] != '\n') {
                        ++pos;
                    }
                    const size_t lineEnd = pos;
                    if (pos < usable && block[pos] == '\n') {
                        ++pos;
                    }
                    if (lineEnd <= lineStart) {
                        continue;
                    }
                    size_t tabPos = lineStart;
                    while (tabPos < lineEnd && block[tabPos] != '\t') {
                        ++tabPos;
                    }
                    if (tabPos <= lineStart || tabPos >= lineEnd) {
                        continue;
                    }
                    const unsigned int oldKey = static_cast<unsigned int>(
                        std::strtoul(block + lineStart, NULL, 10));
                    std::unordered_map<unsigned int, unsigned int>::const_iterator it = keyMaps[ri].find(oldKey);
                    if (it == keyMaps[ri].end()) {
                        continue;
                    }
                    out.append(SSTR(it->second));
                    out.append(block + tabPos, lineEnd - tabPos);
                    out.push_back('\n');
                }
            }
            resultWriter.writeData(out.c_str(), out.size(), queryKey, 0, true, true);
        }
        resultWriter.close(true);
    cleanupReaders();
#undef CMBUILD_MERGE_FAIL
    return true;
}

bool buildQueryCmText(LocalParameters &par, CmBuildCtx &ctx, size_t id,
                      unsigned int thread_idx, std::string &cmText, std::string &err,
                      bool forScan, bool useInside, bool useLocal,
                      CM_t **ret_cm) {
    DBReader<unsigned int> &qDbr = *ctx.qDbr;
    DBReader<unsigned int> *tDbr = ctx.tDbr;
    DBReader<unsigned int> &resultReader = *ctx.resultReader;
    SubstitutionMatrix &subMat = *ctx.subMat;
    MsaFilter &msaFilter = *ctx.msaFilter;
    Sequence &targetMapper = *ctx.targetMapper;
    const std::vector<int> &qid_vec = ctx.qid_vec;
    const float minColCoverage = ctx.minColCoverage;
    const bool doMsaFilter = ctx.doMsaFilter;
    const bool decodeTargetDinuc = ctx.decodeTargetDinuc;
    const bool targetGpuDb = ctx.targetGpuDb;
    std::vector<Matcher::result_t> alnResults; alnResults.reserve(300);
    std::set<unsigned int> seenKeys;
    std::vector<unsigned char> encBuf;
    std::vector<const char *> rowPtrs;
            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr.getId(queryKey);
            if (queryId == UINT_MAX) {
                err = "invalid query " + std::to_string(queryKey);
                return false;
            }

            // Render the query row directly from raw nucleotide DB chars.
            // The MSA stored in MultipleAlignment::MSAResult uses the
            // dinucleotide alphabet (24-symbol), so num2aa[]-based decoding
            // collapses most positions to 'N' — fatal for folding/Infernal.
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
            stoSeqs.reserve(alnResults.size() + 1);

            AlnSeq queryAs;
            queryAs.id = "query_" + std::to_string(queryKey);
            queryAs.aln = queryRow;
            stoSeqs.push_back(queryAs);

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
                    stoSeqs.push_back(std::move(as));
                }
            }

            // SS_cons from a single-sequence MFE fold of the query row.
            std::string ssCons;
            double mfe = 0.0;
            bool ssOk = rnaFoldMfeDotBracket(queryRow, ssCons, &mfe);
            (void)mfe;
            if (std::getenv("RIBOSEEK_DEBUG_SSCONS") != nullptr) {
                Debug(Debug::WARNING) << "[ssCons] q=" << queryKey
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

            const std::string modelName = "query_" + std::to_string(queryKey);
            std::string buildErr;
            std::vector<std::string> sqnames;
            std::vector<std::string> aseqs;
            sqnames.reserve(stoSeqs.size());
            aseqs.reserve(stoSeqs.size());
            for (const AlnSeq &s : stoSeqs) {
                sqnames.push_back(s.id);
                std::string row = s.aln;
                convertTsToUs(row);
                aseqs.push_back(std::move(row));
            }
#ifdef RIBOSEEK_CMBUILD_DUMP_STOCKHOLM
            // Serializes seed MSA + SS_cons fed to Infernal cmbuild
            std::ofstream sto(std::string(RIBOSEEK_CMBUILD_DUMP_STOCKHOLM) + "/" + modelName + ".sto");
            if (sto) {
                sto << "# STOCKHOLM 1.0\n#=GF ID " << modelName << "\n";
                for (size_t i = 0; i < sqnames.size(); ++i)
                    sto << sqnames[i] << " " << aseqs[i] << "\n";
                if (!aseqs.empty() && !aseqs[0].empty())
                    sto << "#=GC RF " << std::string(aseqs[0].size(), 'x') << "\n";
                if (!ssCons.empty())
                    sto << "#=GC SS_cons " << ssCons << "\n";
                sto << "//\n";
            }
#endif
            bool success = cmbuildFromAlignment(modelName, sqnames, aseqs, ssCons, cmText, buildErr,
                                                par.cmbuildEre, par.cmbuildSymfrac, par.cmbuildNoss != 0,
                                                forScan, useInside, useLocal, ret_cm);
    if (!success) { err = "build failed for query " + std::to_string(queryKey) + ": " + buildErr; return false; }
    return true;
}

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

    std::string rawTargetDbArg;
    std::string rawResultDbArg;
    std::vector<std::string> parseArgStorage;
    std::vector<const char *> parseArgv;
    if (argc > 2) {
        rawTargetDbArg = argv[1];
        rawResultDbArg = argv[2];
        if (rawTargetDbArg.find(',') != std::string::npos
            || rawResultDbArg.find(',') != std::string::npos) {
            parseArgStorage.reserve(static_cast<size_t>(argc));
            parseArgv.reserve(static_cast<size_t>(argc));
            for (int i = 0; i < argc; ++i) {
                parseArgStorage.push_back(argv[i]);
            }
            const std::vector<std::string> targetParts = splitCommaList(rawTargetDbArg);
            const std::vector<std::string> resultParts = splitCommaList(rawResultDbArg);
            if (!targetParts.empty()) {
                parseArgStorage[1] = targetParts[0];
            }
            if (!resultParts.empty()) {
                parseArgStorage[2] = resultParts[0];
            }
            for (size_t i = 0; i < parseArgStorage.size(); ++i) {
                parseArgv.push_back(parseArgStorage[i].c_str());
            }
        }
    }
    par.parseParameters(argc,
                        parseArgv.empty() ? argv : parseArgv.data(),
                        command,
                        true,
                        0,
                        0);

    // initialize Infernal logsum tables once
    cmInfernalGlobalInit();

    const std::vector<std::string> targetDbList = splitCommaList(rawTargetDbArg.empty() ? par.db2 : rawTargetDbArg);
    const std::vector<std::string> resultDbList = splitCommaList(rawResultDbArg.empty() ? par.db3 : rawResultDbArg);
    if (targetDbList.size() != resultDbList.size()) {
        Debug(Debug::ERROR) << "cmbuild: targetDB/resultDB list sizes differ ("
                            << targetDbList.size() << " vs " << resultDbList.size() << ")\n";
        return EXIT_FAILURE;
    }
    const bool useMergedInputs = (targetDbList.size() > 1 || resultDbList.size() > 1);
    if (useMergedInputs) {
        const std::string mergedTargetDb = par.db4 + "_target_merged";
        const std::string mergedResultDb = par.db4 + "_result_merged";
        removeDbFiles(mergedTargetDb);
        removeDbFiles(mergedResultDb);
        std::string mergeErr;
        Debug(Debug::INFO) << "cmbuild: merging " << targetDbList.size()
                           << " target/result DB pairs next to output CM DB\n";
        if (!mergeCmbuildInputs(par.db1, targetDbList, resultDbList,
                                mergedTargetDb, mergedResultDb,
                                par.threads, mergeErr)) {
            Debug(Debug::ERROR) << "cmbuild: failed to merge multi-DB input: " << mergeErr << "\n";
            return EXIT_FAILURE;
        }
        par.db2 = mergedTargetDb;
        par.db2Index = mergedTargetDb + ".index";
        par.db3 = mergedResultDb;
        par.db3Index = mergedResultDb + ".index";
        Debug(Debug::INFO) << "cmbuild: merged target DB: " << par.db2 << "\n";
        Debug(Debug::INFO) << "cmbuild: merged result DB: " << par.db3 << "\n";
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

    // Stockholm rows are filtered: drop any target with non-gap coverage below
    // --cmbuild-min-row-cov. Infernal's cm_parsetree_Doctor() can fail when many
    // short fragments dominate the column statistics.
    const float minColCoverage = par.cmbuildMinRowCov;

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
        MsaFilter msaFilter(8192, 1024, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
        Sequence targetMapper(tDbr->getMaxSeqLen(), targetSeqType, &subMat, 0, false, false);
        CmBuildCtx ctx;
        ctx.qDbr = &qDbr; ctx.tDbr = tDbr; ctx.resultReader = &resultReader;
        ctx.subMat = &subMat; ctx.msaFilter = &msaFilter; ctx.targetMapper = &targetMapper;
        ctx.qid_vec = qid_vec; ctx.minColCoverage = minColCoverage;
        ctx.doMsaFilter = doMsaFilter; ctx.decodeTargetDinuc = decodeTargetDinuc; ctx.targetSeqType = targetSeqType;
        ctx.targetGpuDb = targetGpuDb;

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            std::string cmText, buildErr;
            if (buildQueryCmText(par, ctx, id, thread_idx, cmText, buildErr)) {
                cmWriter.writeData(cmText.c_str(), cmText.size(), resultReader.getDbKey(id), thread_idx);
            } else if (!buildErr.empty()) {
                Debug(Debug::WARNING) << "cmbuild: " << buildErr << "\n";
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

    return EXIT_SUCCESS;
}
