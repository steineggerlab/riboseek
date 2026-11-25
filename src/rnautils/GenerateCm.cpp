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
#include "infernal/riboseek_cmcalibrate.h"

#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iterator>
#include <cstdlib>
#include <dirent.h>
#include <unistd.h>
#include <sys/wait.h>

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

// Decode only [start, start + len) of the target: a hit on a genomic target
// would otherwise read and decode the whole target for every alignment row.
static std::string decodeDbSequence(DBReader<unsigned int> &dbr,
                                    size_t id,
                                    unsigned int threadIdx,
                                    Sequence &mapper,
                                    bool useDinucMapping,
                                    bool isGpuDb,
                                    const SubstitutionMatrix &subMat,
                                    size_t start,
                                    size_t len) {
    if (isGpuDb) {
        const unsigned char *data =
            reinterpret_cast<const unsigned char *>(dbr.getDataUncompressed(id)) + start;
        mapper.mapSequence(id, dbr.getDbKey(id), std::make_pair(data, static_cast<unsigned int>(len)));
    } else {
        const char *data = dbr.getDataRange(id, threadIdx, start, len);
        mapper.mapSequence(id, dbr.getDbKey(id), data, static_cast<unsigned int>(len));
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

// ---------------------------------------------------------------------------
// R-scape CaCoFold structure (--cmbuild-ss-mode 1)
// ---------------------------------------------------------------------------

static std::string shellQuote(const std::string &s) {
    std::string out = "'";
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\'') out += "'\\''";
        else out.push_back(s[i]);
    }
    out.push_back('\'');
    return out;
}

// Resolve an executable name through PATH (or check a path containing '/').
static bool findExecutable(const std::string &name, std::string &resolved) {
    if (name.empty()) return false;
    if (name.find('/') != std::string::npos) {
        if (access(name.c_str(), X_OK) == 0) { resolved = name; return true; }
        return false;
    }
    const char *path = std::getenv("PATH");
    if (path == NULL) return false;
    const std::vector<std::string> dirs = Util::split(path, ":");
    for (size_t i = 0; i < dirs.size(); ++i) {
        const std::string cand = (dirs[i].empty() ? std::string(".") : dirs[i]) + "/" + name;
        if (access(cand.c_str(), X_OK) == 0) { resolved = cand; return true; }
    }
    return false;
}

// R-scape writes a flat output directory; remove its files and the directory.
static void removeFlatDir(const std::string &dir) {
    DIR *d = opendir(dir.c_str());
    if (d == NULL) return;
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        const std::string n = e->d_name;
        if (n == "." || n == "..") continue;
        unlink((dir + "/" + n).c_str());
    }
    closedir(d);
    rmdir(dir.c_str());
}

// Concatenate the "#=GC <tag>" lines of a (possibly interleaved) Stockholm file.
static bool readStockholmGc(const std::string &path, const std::string &tag, std::string &out) {
    out.clear();
    std::ifstream in(path.c_str());
    if (!in) return false;
    const std::string prefix = "#=GC " + tag;
    std::string line;
    bool found = false;
    while (std::getline(in, line)) {
        if (line.compare(0, prefix.size(), prefix) != 0) continue;
        if (line.size() > prefix.size() && !std::isspace(static_cast<unsigned char>(line[prefix.size()]))) {
            continue;  // e.g. SS_cons_1 / SS_cons_pk_1
        }
        out += trimCopy(line.substr(prefix.size()));
        found = true;
    }
    return found;
}

static std::string tailOfFile(const std::string &path, size_t maxBytes) {
    std::ifstream in(path.c_str());
    if (!in) return std::string();
    std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (all.size() > maxBytes) all = all.substr(all.size() - maxBytes);
    return trimCopy(all);
}

// Run R-scape on the alignment (rows as given to cmbuild) and return the
// CaCoFold SS_cons (all alignment columns; pseudoknot letters are dropped
// later by Infernal's HandModelmaker, as in cmbuild).
static bool rscapeCacofoldSs(const LocalParameters &par, const CmBuildCtx &ctx,
                             const std::string &modelName,
                             const std::vector<AlnSeq> &rows,
                             std::string &ssCons, std::string &err) {
    ssCons.clear();
    if (rows.empty() || rows[0].aln.empty()) { err = "empty alignment"; return false; }
    const size_t alen = rows[0].aln.size();
    const std::string dir = ctx.rscapeWorkDir + "/" + modelName;
    removeFlatDir(dir);
    if (!FileUtil::makeDir(dir.c_str())) { err = "cannot create " + dir; return false; }
    const std::string stoPath = dir + "/input.sto";
    {
        std::ofstream sto(stoPath.c_str());
        if (!sto) { err = "cannot write " + stoPath; return false; }
        sto << "# STOCKHOLM 1.0\n#=GF ID " << modelName << "\n";
        // Exact duplicate aligned rows (common in raw search results) would
        // distort the tree R-scape builds for its null; write each distinct
        // row once. Row 0 (query/reference) is always written first.
        std::unordered_set<std::string> seenRows;
        size_t nWritten = 0;
        for (size_t i = 0; i < rows.size(); ++i) {
            if (!seenRows.insert(rows[i].aln).second) continue;
            sto << rows[i].id << " " << rows[i].aln << "\n";
            ++nWritten;
        }
        if (std::getenv("RIBOSEEK_DEBUG_SSCONS") != nullptr) {
            Debug(Debug::WARNING) << "[rscape] " << modelName << " rows=" << rows.size()
                                  << " distinct=" << nWritten << "\n";
        }
        sto << "#=GC RF " << std::string(alen, 'x') << "\n//\n";
        if (!sto) { err = "cannot write " + stoPath; return false; }
    }
    const std::string logPath = dir + "/rscape.log";
    const std::string cmd = shellQuote(ctx.rscapeExe) + " " + par.rscapeArgs
                          + " --outdir " + shellQuote(dir) + " --outname model "
                          + shellQuote(stoPath) + " > " + shellQuote(logPath) + " 2>&1";
    const int rc = std::system(cmd.c_str());
    const std::string outSto = dir + "/model.cacofold.sto";
    if (rc != 0 || !readStockholmGc(outSto, "SS_cons", ssCons)) {
        err = "R-scape failed (exit " + SSTR(WIFEXITED(rc) ? WEXITSTATUS(rc) : rc) + "), files kept in "
            + dir + ": " + tailOfFile(logPath, 300);
        ssCons.clear();
        return false;
    }
    if (ssCons.size() != alen) {
        err = "R-scape SS_cons has " + SSTR(ssCons.size()) + " columns, alignment has " + SSTR(alen)
            + "; files kept in " + dir;
        ssCons.clear();
        return false;
    }
    removeFlatDir(dir);
    return true;
}

// Nested base pairs in a WUSS string. Pseudoknot letters (Aa, Bb, ...) are not
// counted: Infernal's HandModelmaker removes them before building.
static size_t countNestedPairs(const std::string &ss) {
    size_t n = 0;
    for (size_t i = 0; i < ss.size(); ++i) {
        const char c = ss[i];
        if (c == '<' || c == '(' || c == '[' || c == '{') ++n;
    }
    return n;
}

// Shared tail of both inputs (result DB and MSA DB): optional row filter,
// SS_cons, then Infernal cmbuild. Row 0 is the reference (query) row.
static bool finishCmBuild(LocalParameters &par, CmBuildCtx &ctx, unsigned int key,
                          std::vector<AlnSeq> &stoSeqs, const std::string &givenSsCons,
                          std::string &cmText, std::string &err,
                          bool forScan, bool useInside, bool useLocal, CM_t **ret_cm) {
    if (stoSeqs.empty() || stoSeqs[0].aln.empty()) {
        err = "empty alignment for key " + std::to_string(key);
        return false;
    }
    const std::string modelName = "query_" + std::to_string(key);

    // Optional MsaFilter pass: rMSA/hhfilter-style row pruning before
    // cmbuild. Encodes ACGU/N/'-' rows into MsaFilter's byte alphabet
    // (residues < NAA, N=NAA, '-'=GAP) so identity comparisons land in
    // pure byte-equality and qsc/qid stay disabled (we set qsc=-50 and
    // qid="0" by default since the loaded subMat is dinucleotide).
    if (ctx.doMsaFilter && stoSeqs.size() > 1) {
        std::vector<unsigned char> encBuf;
        std::vector<const char *> rowPtrs;
        const int N = static_cast<int>(stoSeqs.size());
        const int L = static_cast<int>(stoSeqs[0].aln.size());
        encBuf.assign(static_cast<size_t>(N) * static_cast<size_t>(L), MultipleAlignment::GAP);
        rowPtrs.assign(N, nullptr);
        for (int i = 0; i < N; ++i) {
            unsigned char *dst = encBuf.data() + static_cast<size_t>(i) * static_cast<size_t>(L);
            const std::string &aln = stoSeqs[i].aln;
            const int copyLen = std::min(L, static_cast<int>(aln.size()));
            for (int p = 0; p < copyLen; ++p) dst[p] = encodeAlignedNuc(aln[p]);
            rowPtrs[i] = reinterpret_cast<const char *>(dst);
        }
        ctx.msaFilter->filter(N, L,
                              static_cast<int>(par.covMSAThr * 100),
                              ctx.qid_vec, par.qsc,
                              static_cast<int>(par.filterMaxSeqId * 100),
                              par.Ndiff, par.filterMinEnable,
                              rowPtrs.data(), false);
        bool *keptRaw = new bool[N];
        ctx.msaFilter->getKept(keptRaw, N);
        std::vector<AlnSeq> keptSeqs;
        keptSeqs.reserve(N);
        for (int i = 0; i < N; ++i) {
            if (keptRaw[i]) keptSeqs.push_back(std::move(stoSeqs[i]));
        }
        delete[] keptRaw;
        stoSeqs = std::move(keptSeqs);
        if (stoSeqs.empty()) {
            err = "MSA filter removed all rows for key " + std::to_string(key);
            return false;
        }
    }

    // SS_cons: R-scape CaCoFold (mode 1), else a given SS_cons (MSA DB
    // entry), else a single-sequence MFE fold of row 0.
    std::string ssCons;
    bool ssOk = false;
    bool querySsNoss = false;
    if (par.cmbuildNoss == 0) {
        if (par.cmbuildSsMode == 1) {
            std::string rsErr;
            ssOk = rscapeCacofoldSs(par, ctx, modelName, stoSeqs, ssCons, rsErr);
            if (ssOk && countNestedPairs(ssCons) == 0) {
                // CaCoFold found no nested base pair: build this query like cmbuild --noss
                Debug(Debug::INFO) << "cmbuild: " << modelName
                                   << ": CaCoFold structure has no base pairs; building it without structure (--noss)\n";
                querySsNoss = true;
            }
            if (!ssOk) {
                Debug(Debug::WARNING) << "cmbuild: " << modelName << ": " << rsErr
                                      << "; using the "
                                      << (givenSsCons.empty() ? "MFE fold of the first row" : "input SS_cons")
                                      << " instead\n";
            }
        }
        if (!ssOk && !givenSsCons.empty()) {
            ssCons = givenSsCons;
            ssOk = true;
        }
        if (!ssOk) {
            double mfe = 0.0;
            ssOk = rnaFoldMfeDotBracket(stoSeqs[0].aln, ssCons, &mfe);
        }
    }
    if (std::getenv("RIBOSEEK_DEBUG_SSCONS") != nullptr) {
        Debug(Debug::WARNING) << "[ssCons] q=" << key
            << " ok=" << (ssOk ? 1 : 0)
            << " len=" << ssCons.size()
            << " ss=[" << ssCons << "]\n";
    }

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
                                        par.cmbuildEre, par.cmbuildSymfrac, par.cmbuildNoss != 0 || querySsNoss,
                                        forScan, useInside, useLocal, ret_cm);
    if (!success) { err = "build failed for query " + std::to_string(key) + ": " + buildErr; return false; }
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
    Sequence &targetMapper = *ctx.targetMapper;
    const float minColCoverage = ctx.minColCoverage;
    const bool decodeTargetDinuc = ctx.decodeTargetDinuc;
    const bool targetGpuDb = ctx.targetGpuDb;
    std::vector<Matcher::result_t> alnResults; alnResults.reserve(300);
    std::set<unsigned int> seenKeys;
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

                const int tLo = std::min(res.dbStartPos, res.dbEndPos);
                const int tHi = std::max(res.dbStartPos, res.dbEndPos);
                if (tLo < 0) continue;
                // aligned span plus one base each side, as in result2msa, so the
                // dinucleotide codes at the span edges decode as in the full target
                const size_t ctxStart = (tLo > 0) ? static_cast<size_t>(tLo - 1) : 0;
                const size_t ctxEnd = std::min(static_cast<size_t>(tDbr->getSeqLen(targetId)),
                                               static_cast<size_t>(tHi) + 2);
                if (ctxEnd <= ctxStart) continue;
                std::string decodedTargetSeq = decodeDbSequence(*tDbr, targetId, thread_idx,
                                                                targetMapper, decodeTargetDinuc,
                                                                targetGpuDb, subMat,
                                                                ctxStart, ctxEnd - ctxStart);
                // targetSeq[p] is target position ctxStart + p
                const char *targetSeq = decodedTargetSeq.c_str();
                const int ctxOffset = static_cast<int>(ctxStart);

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
                    switch (res.backtrace[bp]) {
                        case 'M': {
                            const int p = isReverseStrand
                                ? res.dbStartPos - ctxOffset - static_cast<int>(seqPos)
                                : res.dbStartPos - ctxOffset + static_cast<int>(seqPos);
                            const char rawChar = isReverseStrand ? Orf::complement(targetSeq[p]) : targetSeq[p];
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

            return finishCmBuild(par, ctx, queryKey, stoSeqs, std::string(), cmText, err,
                                 forScan, useInside, useLocal, ret_cm);
}

// ---------------------------------------------------------------------------
// MSA DB input: riboseek cmbuild msaDB outputCMDB
// ---------------------------------------------------------------------------

static inline char normalizeMsaChar(char c) {
    if (c == '-' || c == '.' || c == '_' || c == '~') return '-';
    c = normalizeBase(c);
    if (c != 'A' && c != 'C' && c != 'G' && c != 'U') c = 'N';
    return c;
}

// Parse one MSA DB entry: Stockholm (convertmsa, result2msa --msa-format-mode 4),
// aligned FASTA (modes 2/3) or A3M (modes 5/6; lowercase insert columns are
// dropped). Row 0 is the reference row. Returns the entry's SS_cons if present.
static bool parseMsaEntry(const char *data, size_t len, std::vector<AlnSeq> &rows,
                          std::string &ssCons, std::string &err) {
    rows.clear();
    ssCons.clear();
    std::vector<std::string> names;
    std::vector<std::string> seqs;
    std::unordered_map<std::string, size_t> stoIndex;
    bool isStockholm = false;
    bool sawFirst = false;
    size_t pos = 0;
    while (pos < len && data[pos] != '\0') {
        size_t eol = pos;
        while (eol < len && data[eol] != '\n' && data[eol] != '\0') ++eol;
        std::string line(data + pos, eol - pos);
        pos = (eol < len && data[eol] == '\n') ? eol + 1 : eol;
        if (!line.empty() && line[line.size() - 1] == '\r') line.resize(line.size() - 1);
        const std::string t = trimCopy(line);
        if (t.empty()) continue;
        if (!sawFirst) {
            sawFirst = true;
            if (t.compare(0, 11, "# STOCKHOLM") == 0) { isStockholm = true; continue; }
        }
        if (isStockholm) {
            if (t == "//") break;
            if (t.compare(0, 12, "#=GC SS_cons") == 0
                && (t.size() == 12 || std::isspace(static_cast<unsigned char>(t[12])))) {
                ssCons += trimCopy(t.substr(12));
                continue;
            }
            if (t[0] == '#') continue;
            size_t sp = 0;
            while (sp < t.size() && !std::isspace(static_cast<unsigned char>(t[sp]))) ++sp;
            const std::string name = t.substr(0, sp);
            std::string seq;
            for (size_t i = sp; i < t.size(); ++i) {
                if (!std::isspace(static_cast<unsigned char>(t[i]))) seq.push_back(t[i]);
            }
            std::unordered_map<std::string, size_t>::const_iterator it = stoIndex.find(name);
            if (it == stoIndex.end()) {
                stoIndex[name] = names.size();
                names.push_back(name);
                seqs.push_back(seq);
            } else {
                seqs[it->second] += seq;
            }
        } else {
            if (t[0] == '#') continue;  // A3M comment / consensus header lines
            if (t[0] == '>') {
                size_t e = 1;
                while (e < t.size() && !std::isspace(static_cast<unsigned char>(t[e]))) ++e;
                names.push_back(t.substr(1, e - 1));
                seqs.push_back(std::string());
            } else {
                if (seqs.empty()) { err = "sequence data before the first header"; return false; }
                for (size_t i = 0; i < t.size(); ++i) {
                    if (!std::isspace(static_cast<unsigned char>(t[i]))) seqs.back().push_back(t[i]);
                }
            }
        }
    }
    if (seqs.empty()) { err = "no sequences"; return false; }

    // Aligned FASTA: every row has the same length, keep all columns.
    // A3M: lengths differ; drop lowercase insert states and '.'.
    bool sameLen = true;
    for (size_t i = 1; i < seqs.size(); ++i) {
        if (seqs[i].size() != seqs[0].size()) { sameLen = false; break; }
    }
    if (!isStockholm && !sameLen) {
        for (size_t i = 0; i < seqs.size(); ++i) {
            std::string kept;
            kept.reserve(seqs[i].size());
            for (size_t j = 0; j < seqs[i].size(); ++j) {
                const char c = seqs[i][j];
                if (std::islower(static_cast<unsigned char>(c)) || c == '.') continue;
                kept.push_back(c);
            }
            seqs[i].swap(kept);
        }
        sameLen = true;
        for (size_t i = 1; i < seqs.size(); ++i) {
            if (seqs[i].size() != seqs[0].size()) { sameLen = false; break; }
        }
    }
    if (!sameLen) { err = "rows have different aligned lengths"; return false; }
    if (seqs[0].empty()) { err = "empty alignment"; return false; }
    if (!ssCons.empty() && ssCons.size() != seqs[0].size()) {
        err = "SS_cons length " + SSTR(ssCons.size()) + " differs from alignment length " + SSTR(seqs[0].size());
        return false;
    }

    // Unique, whitespace-free names (Infernal and R-scape need them).
    std::set<std::string> used;
    rows.reserve(seqs.size());
    for (size_t i = 0; i < seqs.size(); ++i) {
        std::string name;
        for (size_t j = 0; j < names[i].size(); ++j) {
            const unsigned char c = static_cast<unsigned char>(names[i][j]);
            name.push_back(std::isgraph(c) ? static_cast<char>(c) : '_');
        }
        if (name.empty()) name = "seq" + SSTR(i);
        if (used.count(name) != 0) {
            size_t k = 2;
            while (used.count(name + "_" + SSTR(k)) != 0) ++k;
            name += "_" + SSTR(k);
        }
        used.insert(name);
        AlnSeq a;
        a.id = name;
        a.aln.resize(seqs[i].size());
        for (size_t j = 0; j < seqs[i].size(); ++j) a.aln[j] = normalizeMsaChar(seqs[i][j]);
        rows.push_back(std::move(a));
    }
    return true;
}

// Resolve R-scape and create its scratch directory next to the output DB.
static bool setupRscape(const LocalParameters &par, const std::string &outDb,
                        std::string &exe, std::string &workDir) {
    if (par.cmbuildSsMode != 1 || par.cmbuildNoss != 0) return true;
    if (!findExecutable(par.rscapeBin, exe)) {
        Debug(Debug::ERROR) << "cmbuild: --cmbuild-ss-mode 1 needs R-scape; \"" << par.rscapeBin
                            << "\" was not found (set --rscape-bin)\n";
        return false;
    }
    workDir = outDb + "_rscape";
    if (!FileUtil::directoryExists(workDir.c_str()) && !FileUtil::makeDir(workDir.c_str())) {
        Debug(Debug::ERROR) << "cmbuild: cannot create R-scape directory " << workDir << "\n";
        return false;
    }
    Debug(Debug::INFO) << "cmbuild: SS_cons from R-scape (" << exe << " " << par.rscapeArgs << ")\n";
    return true;
}

static void calibrateIfRequested(const LocalParameters &par, std::string &cmText, unsigned int key) {
    if (par.cmCalibrate == 0) return;
    // Infernal cmcalibrate in memory; on failure keep the uncalibrated CM (cmsearch then reports E = 1)
    char *calText = NULL;
    size_t calLen = 0;
    char calErr[128];
    if (riboseek_cmcalibrate(cmText.c_str(), cmText.size(), par.cmCalibrateLen,
                             par.cmCalibrate == 2 ? 1 : 0, &calText, &calLen, calErr) == 0) {
        cmText.assign(calText, calLen);
        free(calText);
    } else {
        Debug(Debug::WARNING) << "cmbuild: cmcalibrate failed for query key " << key << ": " << calErr << "\n";
    }
}

static int cmbuildFromMsaDb(LocalParameters &par) {
    DBReader<unsigned int> msaReader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                     DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    msaReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    if (Parameters::isEqualDbtype(msaReader.getDbtype(), Parameters::DBTYPE_MSA_DB) == false) {
        Debug(Debug::WARNING) << "cmbuild: " << par.db1 << " is a " << msaReader.getDbTypeName()
                              << " DB, reading its entries as MSAs\n";
    }
    if (msaReader.getSize() > 0) {
        // binary cA3M (result2msa --msa-format-mode 0/1) is not text
        const char *d = msaReader.getData(0, 0);
        if (d != NULL && msaReader.getEntryLen(0) > 1 && d[0] != '#' && d[0] != '>') {
            Debug(Debug::ERROR) << "cmbuild: MSA entries must be Stockholm, aligned FASTA or A3M text "
                                   "(result2msa --msa-format-mode 2-6); binary cA3M is not supported\n";
            msaReader.close();
            return EXIT_FAILURE;
        }
    }

    std::string rscapeExe, rscapeWorkDir;
    if (!setupRscape(par, par.db2, rscapeExe, rscapeWorkDir)) {
        msaReader.close();
        return EXIT_FAILURE;
    }

    DBWriter cmWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                      par.compressed, Parameters::DBTYPE_GENERIC_DB);
    cmWriter.open();
    Debug(Debug::INFO) << "MSA database size: " << msaReader.getSize() << "\n";
    if (par.cmCalibrate != 0) {
        Debug(Debug::INFO) << "cmbuild: calibrating CM E-values (cmcalibrate -L " << par.cmCalibrateLen
                           << (par.cmCalibrate == 2 ? " --random" : "") << "); expect ~10-90 CPU-min per CM\n";
    }

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0f);
    std::vector<int> qid_vec;
    {
        std::vector<std::string> qid_str_vec = Util::split(par.qid, ",");
        for (size_t i = 0; i < qid_str_vec.size(); ++i) {
            qid_vec.push_back(static_cast<int>(static_cast<float>(std::atof(qid_str_vec[i].c_str())) * 100));
        }
        if (qid_vec.empty()) qid_vec.push_back(0);
        std::sort(qid_vec.begin(), qid_vec.end());
    }
    Debug::Progress progress(msaReader.getSize());
    size_t nFailed = 0;

#pragma omp parallel reduction(+:nFailed)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        MsaFilter msaFilter(8192, 1024, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
        CmBuildCtx ctx;
        ctx.subMat = &subMat; ctx.msaFilter = &msaFilter;
        ctx.qid_vec = qid_vec; ctx.minColCoverage = par.cmbuildMinRowCov;
        ctx.doMsaFilter = (par.filterMsa != 0);
        ctx.rscapeExe = rscapeExe; ctx.rscapeWorkDir = rscapeWorkDir;
        std::vector<AlnSeq> rows;

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < msaReader.getSize(); id++) {
            progress.updateProgress();
            const unsigned int key = msaReader.getDbKey(id);
            const char *data = msaReader.getData(id, thread_idx);
            const size_t len = msaReader.getEntryLen(id);
            std::string givenSs, err;
            if (!parseMsaEntry(data, len, rows, givenSs, err)) {
                Debug(Debug::WARNING) << "cmbuild: MSA key " << key << ": " << err << "\n";
                nFailed++;
                continue;
            }
            // drop non-reference rows with coverage below --cmbuild-min-row-cov
            const size_t alen = rows[0].aln.size();
            std::vector<AlnSeq> kept;
            kept.reserve(rows.size());
            for (size_t i = 0; i < rows.size(); ++i) {
                size_t res = 0;
                for (size_t j = 0; j < alen; ++j) res += (rows[i].aln[j] != '-');
                if (i == 0 || static_cast<float>(res) / static_cast<float>(alen) >= ctx.minColCoverage) {
                    kept.push_back(std::move(rows[i]));
                }
            }
            std::string cmText;
            if (finishCmBuild(par, ctx, key, kept, givenSs, cmText, err, false, true, false, NULL)) {
                calibrateIfRequested(par, cmText, key);
                cmWriter.writeData(cmText.c_str(), cmText.size(), key, thread_idx);
            } else {
                Debug(Debug::WARNING) << "cmbuild: " << err << "\n";
                nFailed++;
            }
        }
    }
    const size_t nMsa = msaReader.getSize();
    cmWriter.close(true);
    msaReader.close();
    if (!rscapeWorkDir.empty()) rmdir(rscapeWorkDir.c_str());  // stays if failed runs were kept
    if (nFailed > 0) {
        Debug(Debug::WARNING) << "cmbuild: " << nFailed << " of " << nMsa
                              << " MSAs produced no CM\n";
    }
    return EXIT_SUCCESS;
}

// Positional (non-option) arguments in argv; bool flags take no value
// unless the next argument is 0/1/true/false.
static size_t countPositionalArgs(int argc, const char **argv, const std::vector<MMseqsParameter *> &params) {
    size_t n = 0;
    for (int i = 0; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.size() > 1 && a[0] == '-') {
            bool isBool = false;
            for (size_t p = 0; p < params.size(); ++p) {
                if (a == params[p]->name) { isBool = (params[p]->type == typeid(bool)); break; }
            }
            if (!isBool) {
                ++i;
            } else if (i + 1 < argc) {
                const std::string v = argv[i + 1];
                if (v == "0" || v == "1" || v == "true" || v == "false") ++i;
            }
            continue;
        }
        ++n;
    }
    return n;
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

    // riboseek cmbuild msaDB outputCMDB
    if (countPositionalArgs(argc, argv, *command.params) == 2) {
        Command msaCommand = command;
        msaCommand.databases.clear();
        msaCommand.databases.push_back({"msaDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL});
        msaCommand.databases.push_back({"outputCMDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb});
        par.parseParameters(argc, argv, msaCommand, true, 0, 0);
        cmInfernalGlobalInit();
        return cmbuildFromMsaDb(par);
    }

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

    std::string rscapeExe, rscapeWorkDir;
    if (!setupRscape(par, par.db4, rscapeExe, rscapeWorkDir)) {
        return EXIT_FAILURE;
    }

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

    if (par.cmCalibrate != 0) {
        Debug(Debug::INFO) << "cmbuild: calibrating CM E-values (cmcalibrate -L " << par.cmCalibrateLen
                           << (par.cmCalibrate == 2 ? " --random" : "") << "); expect ~10-90 CPU-min per CM\n";
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
        ctx.rscapeExe = rscapeExe; ctx.rscapeWorkDir = rscapeWorkDir;

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            std::string cmText, buildErr;
            if (buildQueryCmText(par, ctx, id, thread_idx, cmText, buildErr)) {
                calibrateIfRequested(par, cmText, resultReader.getDbKey(id));
                cmWriter.writeData(cmText.c_str(), cmText.size(), resultReader.getDbKey(id), thread_idx);
            } else if (!buildErr.empty()) {
                Debug(Debug::WARNING) << "cmbuild: " << buildErr << "\n";
            }
        }
    }

    cmWriter.close(true);
    if (!rscapeWorkDir.empty()) rmdir(rscapeWorkDir.c_str());  // stays if failed runs were kept
    resultReader.close();
    if (!sameDatabase) {
        tDbr->close();
        delete tDbr;
    }
    qDbr.close();

    return EXIT_SUCCESS;
}
