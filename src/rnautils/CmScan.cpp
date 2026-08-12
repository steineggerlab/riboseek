#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"
#include "DBReader.h"
#include "FileUtil.h"
#include "Matcher.h"
#include "Util.h"
#include "NucleotideMatrix.h"
#include "CmBuildScan.h"
#include "Util.h"
#include "Sequence.h"

#include <fstream>
#include <unordered_set>

#ifdef OPENMP
#include <omp.h>
#endif

extern "C" {
#include "infernal.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"
}

extern void cmInfernalGlobalInit();

struct FastaSeq {
    unsigned int key = 0;
    std::string id;
    std::string seq;
};

inline std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])) != 0) ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) != 0) --e;
    return s.substr(b, e - b);
}

inline char normalizeBase(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (c == 'T') return 'U';
    return c;
}

char complementBase(char c) {
    switch (normalizeBase(c)) {
        case 'A': return 'U';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'U': return 'A';
        default:  return 'N';
    }
}

std::string reverseComplement(const std::string &seq) {
    std::string rc(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); ++i) {
        rc[seq.size() - 1 - i] = complementBase(seq[i]);
    }
    return rc;
}

// cigar encoding
std::string rleTraceOps(const std::string &ops) {
    if (ops.empty()) return "NA";
    std::string out;
    out.reserve(ops.size() * 2);
    char cur = ops[0];
    int run = 1;
    for (size_t i = 1; i < ops.size(); ++i) {
        if (ops[i] == cur) { ++run; }
        else { out += std::to_string(run); out.push_back(cur); cur = ops[i]; run = 1; }
    }
    out += std::to_string(run);
    out.push_back(cur);
    return out;
}

// cigar decoding
std::string reverseRleTraceOps(const std::string &rle) {
    struct Run { int n; char op; };
    std::vector<Run> runs;
    runs.reserve(rle.size() / 2 + 1);
    size_t pos = 0;
    while (pos < rle.size()) {
        if (!std::isdigit(static_cast<unsigned char>(rle[pos]))) return rle;
        int n = 0;
        while (pos < rle.size() && std::isdigit(static_cast<unsigned char>(rle[pos]))) {
            n = n * 10 + (rle[pos] - '0'); ++pos;
        }
        if (pos >= rle.size() || n <= 0) return rle;
        runs.push_back(Run{n, rle[pos++]});
    }
    std::string out;
    out.reserve(rle.size());
    for (std::vector<Run>::const_reverse_iterator it = runs.rbegin(); it != runs.rend(); ++it) {
        if (!out.empty() && out.back() == it->op) {
            size_t digitEnd = out.size() - 1;
            size_t digitStart = digitEnd;
            while (digitStart > 0 && std::isdigit(static_cast<unsigned char>(out[digitStart - 1]))) --digitStart;
            int prev = std::atoi(out.c_str() + digitStart);
            out.erase(digitStart);
            out += std::to_string(prev + it->n);
            out.push_back(it->op);
        } else {
            out += std::to_string(it->n);
            out.push_back(it->op);
        }
    }
    return out;
}

int effectiveDecodeSeqType(int dbtype, bool useDinucMapping) {
    if (!useDinucMapping) return dbtype;
    unsigned int ext = DBReader<unsigned int>::getExtendedDbtype(dbtype);
    ext |= Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE)) {
        return DBReader<unsigned int>::setExtendedDbtype(dbtype, ext);
    }
    return DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS, ext);
}

std::string decodeMappedSequenceToRna(const Sequence &seq, const BaseMatrix &subMat,
                                      const unsigned char *num2outputnum) {
    std::string out;
    out.reserve(static_cast<size_t>(seq.L));
    for (int i = 0; i < seq.L; ++i) {
        unsigned char code = seq.numSequence[i];
        if (num2outputnum != NULL) code = num2outputnum[code];
        char c = normalizeBase(subMat.num2aa[code]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'U') c = 'N';
        out.push_back(c);
    }
    return out;
}

FastaSeq decodeOneSequence(DBReader<unsigned int> &dbr, size_t id, BaseMatrix &subMat,
                           Sequence &seqObj, bool useDinucMapping, bool isGpuDb,
                           unsigned int thread_idx) {
    FastaSeq cur;
    cur.key = dbr.getDbKey(id);
    if (dbr.getLookupSize() > 0) {
        const size_t lid = dbr.getLookupIdByKey(cur.key);
        cur.id = dbr.getLookupEntryName(lid);
    } else {
        cur.id = std::to_string(cur.key);
    }
    const size_t seqLen = dbr.getSeqLen(id);
    if (isGpuDb) {
        const unsigned char *data = reinterpret_cast<const unsigned char *>(dbr.getDataUncompressed(id));
        seqObj.mapSequence(id, cur.key, std::make_pair(data, seqLen));
    } else {
        const char *data = dbr.getData(id, thread_idx);
        seqObj.mapSequence(id, cur.key, data, seqLen);
    }
    const Sequence::SeqAuxInfo *auxInfo = Sequence::getAuxInfo(seqObj.getSeqType());
    const unsigned char *num2outputnum = (useDinucMapping && auxInfo != NULL) ? auxInfo->num2outputnum : NULL;
    cur.seq = decodeMappedSequenceToRna(seqObj, subMat, num2outputnum);
    return cur;
}

bool hasDbIndex(const std::string &path) {
    return FileUtil::fileExists((path + ".index").c_str());
}

bool looksLikeInfernalCm(const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in.good()) return false;
    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        return t.rfind("INFERNAL", 0) == 0 || t == "CM" || t.rfind("NAME", 0) == 0;
    }
    return false;
}

struct InfModel {
    // shared read-only across threads (immutable)
    ESL_ALPHABET *abc = NULL;
    // configured template; per-thread clones own their scratch   
    CM_t *cm  = NULL;
    int  clen = 0;
    bool valid = false;
    // scan/alignment modes Inside vs CYK
    bool useInside = true;
    // flags for tcm->align_opts
    int  alignOpts = CM_ALIGN_NONBANDED | CM_ALIGN_CYK;
    // HMM-banded? (alidisplay tau)
    bool alignBanded = false;
    // minimum window bit score
    float scanCutoff = 0.0f;
};

// Load an INFERNAL1/a CM from in-memory text and configure it for scannin non-banded alignment
InfModel loadModel(
    const std::string &cmtext,
    bool useInside,
    int alignOpts,
    bool alignBanded,
    bool useLocal,
    float scanCutoff) {
    InfModel m;
    m.useInside = useInside;
    m.alignOpts = alignOpts;
    m.alignBanded = alignBanded;
    m.scanCutoff = scanCutoff;
    std::vector<char> buf(cmtext.begin(), cmtext.end());
    buf.push_back('\0');
    char errbuf[eslERRBUFSIZE];
    errbuf[0] = '\0';

    CM_FILE *cmfp = NULL;
    if (cm_file_OpenBuffer(buf.data(), static_cast<int>(cmtext.size()), FALSE, &cmfp) != eslOK) {
        Debug(Debug::ERROR) << "cmscan: cm_file_OpenBuffer failed\n";
        return m;
    }
    // null so close is safe
    cmfp->hfp = NULL;
    if (cm_file_Read(cmfp, FALSE, &m.abc, &m.cm) != eslOK || m.cm == NULL) {
        Debug(Debug::ERROR) << "cmscan: cm_file_Read failed: " << cmfp->errbuf << "\n";
        cm_file_Close(cmfp);
        return m;
    }
    cm_file_Close(cmfp);
    m.cm->config_opts |= CM_CONFIG_SCANMX | CM_CONFIG_NONBANDEDMX;
    // glocal by default; local adds CM local begins/ends + matching CP9 (HMM-banded aln) local config
    if (useLocal) m.cm->config_opts |= CM_CONFIG_LOCAL | CM_CONFIG_HMMLOCAL | CM_CONFIG_HMMEL;
    if (useInside) m.cm->search_opts |= CM_SEARCH_INSIDE;  // else CYK scan (flag unset)
    if (cm_Configure(m.cm, errbuf, -1) != eslOK) {
        Debug(Debug::ERROR) << "cmscan: cm_Configure failed: " << errbuf << "\n";
        FreeCM(m.cm); m.cm = NULL;
        esl_alphabet_Destroy(m.abc); m.abc = NULL;
        return m;
    }
    m.clen = m.cm->clen;
    m.valid = true;
    return m;
}

void freeModel(InfModel &m) {
    if (m.cm != NULL) { FreeCM(m.cm); m.cm = NULL; }
    if (m.abc != NULL) { esl_alphabet_Destroy(m.abc); m.abc = NULL; }
    m.valid = false;
}

// Per-region hit
struct InfHit {
    bool valid = false;
    // Inside score, null3-corrected
    double score = -std::numeric_limits<double>::infinity();
    // window 1-based target coord of first MATCH residue
    int start1 = 0;
    // window 1-based target coord of last MATCH residue
    int end1 = 0;
    // 0-based first consensus column (cfrom_emit-1)
    int qStart = -1;
    // 0-based last consensus column (cto_emit-1)
    int qEnd = -1;
    // I=target-consume, D=query-consume, over [firstM..lastM]
    std::string cigar = "NA";
    int cigarAlnLen = 0;
    float seqId = -1.0f;
    unsigned int dbKey = 0;
    unsigned int dbLen = 0;
};

// scan one strand-oriented window with Infernal
InfHit scanRegionInfernal(const InfModel &model, CM_t *tcm, const std::string &regionSeq) {
    InfHit h;
    const int Lwin = static_cast<int>(regionSeq.size());
    if (Lwin <= 0) return h;

    char errbuf[eslERRBUFSIZE];
    errbuf[0] = '\0';

    ESL_DSQ *dsq = NULL;
    if (esl_abc_CreateDsq(model.abc, regionSeq.c_str(), &dsq) != eslOK || dsq == NULL) {
        return h;
    }

    CM_TOPHITS *th = cm_tophits_Create();
    float sc = 0.0f;
    int scanStatus;
    if (model.useInside) {
        scanStatus = FastIInsideScan(
            tcm, errbuf, tcm->smx, SMX_QDB1_TIGHT, dsq,
            1, Lwin, model.scanCutoff, th, TRUE /*do_null3*/, 0.0f, NULL, NULL, NULL, &sc
        );
    } else {
        scanStatus = FastCYKScan(
            tcm, errbuf, tcm->smx, SMX_QDB1_TIGHT, dsq,
            1, Lwin, model.scanCutoff, th, TRUE /*do_null3*/, 0.0f, NULL, NULL, NULL, &sc
        );
    }
    if (scanStatus != eslOK || th->N == 0) {
        cm_tophits_Destroy(th);
        free(dsq);
        return h;
    }

    // best hit by score
    CM_HIT *best = NULL;
    for (uint64_t i = 0; i < th->N; ++i) {
        CM_HIT *hh = &th->unsrt[i];
        if (best == NULL || hh->score > best->score) best = hh;
    }
    if (best == NULL || best->start < 1 || best->stop > Lwin || best->stop < best->start) {
        cm_tophits_Destroy(th);
        free(dsq);
        return h;
    }
    const double hitScore = best->score;
    const int64_t hStart = best->start;
    const int64_t hStop  = best->stop;
    const char    hMode  = (char) best->mode;

    // non-banded CYK align of the best hit
    ESL_SQ *winSq  = esl_sq_CreateDigitalFrom(model.abc, "win", dsq, Lwin, NULL, NULL, NULL);
    ESL_SQ *sq2aln = esl_sq_CreateDigitalFrom(model.abc, "hit", dsq + hStart - 1,
                                              hStop - hStart + 1, NULL, NULL, NULL);
    tcm->align_opts = model.alignOpts;
    ESL_STOPWATCH *watch = esl_stopwatch_Create();
    esl_stopwatch_Start(watch);
    CM_ALNDATA *adata = NULL;
    const float mxsize = 8192.0f;
    const int alnStatus = DispatchSqAlignment(tcm, errbuf, sq2aln, -1, mxsize, hMode,
                                              PLI_PASS_STD_ANY, FALSE /*cp9b invalid*/,
                                              NULL, NULL, NULL, &adata);
    esl_stopwatch_Stop(watch);
    if (alnStatus != eslOK || adata == NULL) {
        if (adata != NULL) cm_alndata_Destroy(adata, FALSE);
        esl_sq_Destroy(sq2aln); esl_sq_Destroy(winSq); esl_stopwatch_Destroy(watch);
        cm_tophits_Destroy(th); free(dsq);
        Debug(Debug::WARNING) << "cmscan: alignment failed, skipping hit"
                              << (errbuf[0] != '\0' ? ": " : "") << errbuf << "\n";
        return h;
    }

    CM_ALIDISPLAY *ad = NULL;
    const int adStatus = cm_alidisplay_Create(tcm, errbuf, adata, winSq, hStart,
                                              model.alignBanded ? tcm->tau : -1.0, watch->elapsed, &ad);
    if (adStatus != eslOK || ad == NULL) {
        if (ad != NULL) cm_alidisplay_Destroy(ad);
        cm_alndata_Destroy(adata, FALSE);
        esl_sq_Destroy(sq2aln); esl_sq_Destroy(winSq); esl_stopwatch_Destroy(watch);
        cm_tophits_Destroy(th); free(dsq);
        Debug(Debug::WARNING) << "cmscan: alidisplay failed, skipping hit"
                              << (errbuf[0] != '\0' ? ": " : "") << errbuf << "\n";
        return h;
    }

    const int N = ad->N;
    std::vector<char> ops;
    // window 1-based coord of the residue consumed (M/I); -1 for D
    std::vector<long> targetCoord;
    ops.reserve(N);
    targetCoord.reserve(N);
    long curTarget = static_cast<long>(ad->sqfrom) - 1;
    int idCount = 0, idDenom = 0;
    for (int k = 0; k < N; ++k) {
        const char mc = ad->model[k];
        const char ac = ad->aseq[k];
        const bool modelCons = (std::isalpha(static_cast<unsigned char>(mc)) != 0);
        const bool aseqRes   = (std::isalpha(static_cast<unsigned char>(ac)) != 0);
        char op = '\0';
        if (modelCons && aseqRes) op = 'M';
        else if (!modelCons && aseqRes) op = 'I';   // target residue, no consensus
        else if (modelCons && !aseqRes) op = 'D';   // consensus, no target residue
        if (op == '\0') continue;
        long coord = -1;
        if (op == 'M' || op == 'I') { ++curTarget; coord = curTarget; }
        ops.push_back(op);
        targetCoord.push_back(coord);
    }
    // trim to the first/last MATCH column so the reported span starts/ends on M
    int firstM = -1, lastM = -1;
    for (size_t k = 0; k < ops.size(); ++k) {
        if (ops[k] == 'M') { if (firstM < 0) firstM = static_cast<int>(k); lastM = static_cast<int>(k); }
    }
    if (firstM >= 0) {
        std::string internal;
        internal.reserve(lastM - firstM + 1);
        for (int k = firstM; k <= lastM; ++k) internal.push_back(ops[k]);
        // identity/seqId: re-walk the alidisplay columns
        {
            int opIdx = -1;
            for (int k = 0; k < N; ++k) {
                const char mc = ad->model[k];
                const char ac = ad->aseq[k];
                const bool modelCons = (std::isalpha(static_cast<unsigned char>(mc)) != 0);
                const bool aseqRes   = (std::isalpha(static_cast<unsigned char>(ac)) != 0);
                char op = '\0';
                if (modelCons && aseqRes) op = 'M';
                else if (!modelCons && aseqRes) op = 'I';
                else if (modelCons && !aseqRes) op = 'D';
                if (op == '\0') continue;
                ++opIdx;
                if (opIdx < firstM || opIdx > lastM) continue;
                ++idDenom;
                if (op == 'M') {
                    const char a = normalizeBase(mc);
                    const char b = normalizeBase(ac);
                    if (a == b && a != 'N') ++idCount;
                }
            }
        }
        h.cigar = rleTraceOps(internal);
        h.cigarAlnLen = static_cast<int>(internal.size());
        h.start1 = static_cast<int>(targetCoord[firstM]);
        h.end1 = static_cast<int>(targetCoord[lastM]);
        // Query coords must match the cigar
        const int qConsOffset = (ad->cfrom_emit > 0) ? (ad->cfrom_emit - 1) : 0;
        int qConsBefore = 0;
        for (int k = 0; k < firstM; ++k) {
            if (ops[k] == 'M' || ops[k] == 'D') {
                ++qConsBefore;
            }
        }
        int qConsSpan = 0;
        for (int k = firstM; k <= lastM; ++k) {
            if (ops[k] == 'M' || ops[k] == 'D') {
                ++qConsSpan;
            }
        }
        h.qStart = qConsOffset + qConsBefore;
        h.qEnd = h.qStart + qConsSpan - 1;
        h.seqId = (idDenom > 0) ? static_cast<float>(idCount) / static_cast<float>(idDenom) : -1.0f;
        h.score = hitScore;
        h.valid = true;
    }

    cm_alidisplay_Destroy(ad);
    cm_alndata_Destroy(adata, FALSE);
    esl_sq_Destroy(sq2aln);
    esl_sq_Destroy(winSq);
    esl_stopwatch_Destroy(watch);
    cm_tophits_Destroy(th);
    free(dsq);
    return h;
}

// Coordinate remap window-local to absolute target coords
void finalizeHit(InfHit &h, int strand, unsigned int tKey, unsigned int fullLen,
                 int offset, int regionLen) {
    h.dbKey = tKey;
    h.dbLen = fullLen;
    if (strand > 0) {
        if (offset > 0) { h.start1 += offset; h.end1 += offset; }
    } else {
        const int fwdStart = regionLen - h.start1 + 1 + offset;
        const int fwdEnd = regionLen - h.end1 + 1 + offset;
        h.start1 = fwdStart;
        h.end1 = fwdEnd;
    }
}

// Full per-query scan pipeline: parse cmText -> loadModel -> gather candidates ->
// scan -> write alignments. Shared by cmscan (cmText from CM DB) and cmbuildscan
// (cmText built in-memory). Identical quantization: same cmText -> same cm_file parse.
// Wrap an already-built, already-scan-configured CM (from the direct build path) as
// an InfModel, no ASCII parse. abc stays alive via cm->abc; freeModel frees both.
static InfModel makeModelFromCM(CM_t *cm, bool useInside, int alignOpts, bool alignBanded) {
    InfModel m;
    m.cm = cm; m.abc = (cm != NULL) ? const_cast<ESL_ALPHABET*>(cm->abc) : NULL; m.clen = (cm != NULL) ? cm->clen : 0;
    m.useInside = useInside; m.alignOpts = alignOpts; m.alignBanded = alignBanded;
    m.valid = (cm != NULL);
    return m;
}

static void processQueryModel(LocalParameters &par, const InfModel &model,
                              unsigned int queryKey,
                              DBReader<unsigned int> &seqDbr, DBReader<unsigned int> &resultReader,
                              DBWriter &resultWriter, BaseMatrix &subMat, NucleotideMatrix &nucMat,
                              int seqDecodeType, bool seqGpuDb, bool decodeSeqDinuc, size_t nThreads) {
    (void) nucMat;
        const unsigned int modelLen = static_cast<unsigned int>(std::max(0, model.clen));
        const int W = model.cm->W;   // CM max hit length (Infernal band setup)
        // --cm-region sizes the CM alignment window: the model is aligned to the
        // prefilter hit region [dbStart,dbEnd] widened by (regionFlank * W) nt on
        // EACH side. 1.0 (default) = pad by a full W (a hit can't reach further);
        // <1 tightens the window (faster, may clip edge hits); <=0 = no window
        // (align the whole target).
        const float regionFlank = par.cmRegionFlanking;

        // Gather prefilter candidate regions for this query, first-per-target.
        struct Cand {
            unsigned int tKey;
            int dbStart = 0, dbEnd = 0, qStart = 0, qEnd = 0, qLen = 0;
            bool hasRegionCoord = false;
            bool prefilterIsRev = false;
        };
        std::vector<Cand> cands;
        std::unordered_set<unsigned int> seen;
        const size_t rid = resultReader.getId(queryKey);
        if (rid != UINT_MAX) {
            char *data = resultReader.getData(rid, 0);
            while (data != NULL && *data != '\0') {
                while (*data == ' ' || *data == '\t') { ++data; }
                const char *lineStart = data;
                char *endptr = NULL;
                const unsigned long k = std::strtoul(data, &endptr, 10);
                if (endptr == data || k > static_cast<unsigned long>(UINT_MAX)) {
                    data = Util::skipLine(data);
                    continue;
                }
                const unsigned int tKey = static_cast<unsigned int>(k);
                if (seen.insert(tKey).second == false) {
                    data = Util::skipLine(data);
                    continue;
                }
                Cand cand;
                cand.tKey = tKey;
                if (*endptr == '\t') {
                    const char *p = endptr;
                    int col = 1;
                    const char *colStart[12] = {NULL};
                    colStart[0] = lineStart;
                    while (*p != '\n' && *p != '\0' && col < 11) {
                        if (*p == '\t') { colStart[col] = p + 1; col++; }
                        p++;
                    }
                    if (col >= 10 && colStart[4] && colStart[5] && colStart[7] && colStart[8]) {
                        cand.dbStart = Util::fast_atoi<int>(colStart[7]);
                        cand.dbEnd   = Util::fast_atoi<int>(colStart[8]);
                        cand.qStart  = Util::fast_atoi<int>(colStart[4]);
                        cand.qEnd    = Util::fast_atoi<int>(colStart[5]);
                        cand.qLen    = (colStart[6] != NULL) ? Util::fast_atoi<int>(colStart[6]) : 0;
                        const bool dbRev = (cand.dbStart > cand.dbEnd);
                        const bool qRev  = (cand.qStart > cand.qEnd);
                        cand.prefilterIsRev = (dbRev != qRev);
                        if (dbRev) std::swap(cand.dbStart, cand.dbEnd);
                        if (qRev)  std::swap(cand.qStart, cand.qEnd);
                        cand.hasRegionCoord = true;
                    }
                }
                data = Util::skipLine(data);
                cands.push_back(cand);
            }
        }

        std::vector<InfHit> hits;
#pragma omp parallel num_threads(static_cast<int>(nThreads))
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence seqObjLocal(seqDbr.getMaxSeqLen(), seqDecodeType, &nucMat, 0, false, false);
            std::vector<InfHit> localHits;

            // per-thread CM clone: owns its own smx + align matrices
            char cerr[eslERRBUFSIZE];
            cerr[0] = '\0';
            CM_t *tcm = NULL;
            const int cloneStatus = cm_Clone(model.cm, cerr, &tcm);
            if (cloneStatus != eslOK || tcm == NULL) {
#pragma omp critical
                Debug(Debug::ERROR) << "cmscan: cm_Clone failed: " << cerr << "\n";
            } else {
#pragma omp for schedule(dynamic, 1) nowait
                for (long ci = 0; ci < static_cast<long>(cands.size()); ++ci) {
                    const Cand &cand = cands[static_cast<size_t>(ci)];
                    const size_t tId = seqDbr.getId(cand.tKey);
                    if (tId == UINT_MAX) continue;
                    FastaSeq fs = decodeOneSequence(seqDbr, tId, subMat, seqObjLocal,
                                                    decodeSeqDinuc, seqGpuDb, thread_idx);
                    const int seqLen = static_cast<int>(fs.seq.size());
                    if (seqLen <= 0) continue;

                    // Rescore a WINDOW around the prefilter envelope, padded by W.
                    // Bounds DP by model size, not the target length.
                    std::string regionBuf;
                    const std::string *regionSeqPtr = &fs.seq;
                    int offset = 0;
                    if (cand.hasRegionCoord && regionFlank > 0.0f) {
                        const int basePad = (W > 0) ? W : std::max(1, model.clen);
                        const int pad = std::max(1, static_cast<int>(regionFlank * basePad));
                        const int regStart = std::max(0, cand.dbStart - pad);
                        const int regEnd = std::min(seqLen, cand.dbEnd + 1 + pad);
                        if (regEnd <= regStart) continue;
                        if (regStart != 0 || regEnd != seqLen) {
                            regionBuf = fs.seq.substr(static_cast<size_t>(regStart),
                                                      static_cast<size_t>(regEnd - regStart));
                            regionSeqPtr = &regionBuf;
                        }
                        offset = regStart;
                    }
                    const std::string &regionSeq = *regionSeqPtr;
                    const int regionLen = static_cast<int>(regionSeq.size());

                    const bool scanFwd = !cand.hasRegionCoord || !cand.prefilterIsRev;
                    const bool scanRev = !cand.hasRegionCoord ||  cand.prefilterIsRev;
                    if (scanFwd) {
                        InfHit h = scanRegionInfernal(model, tcm, regionSeq);
                        if (h.valid) {
                            finalizeHit(h, +1, cand.tKey, static_cast<unsigned int>(seqLen), offset, regionLen);
                            localHits.push_back(h);
                        }
                    }
                    if (scanRev) {
                        const std::string revSeq = reverseComplement(regionSeq);
                        InfHit h = scanRegionInfernal(model, tcm, revSeq);
                        if (h.valid) {
                            finalizeHit(h, -1, cand.tKey, static_cast<unsigned int>(seqLen), offset, regionLen);
                            localHits.push_back(h);
                        }
                    }
                }
            }
            if (tcm != NULL) FreeCM(tcm);
#pragma omp critical
            {
                hits.insert(hits.end(), localHits.begin(), localHits.end());
            }
        }

        std::sort(hits.begin(), hits.end(), [](const InfHit &a, const InfHit &b) {
            if (a.dbKey != b.dbKey) return a.dbKey < b.dbKey;
            if (a.start1 != b.start1) return a.start1 < b.start1;
            return a.end1 < b.end1;
        });

        char buffer[1024 + 32768 * 4];
        resultWriter.writeStart(0);
        for (size_t i = 0; i < hits.size(); ++i) {
            const InfHit &h = hits[i];
            int dbStartOut = std::max(0, h.start1 - 1);
            int dbEndOut = std::max(0, h.end1 - 1);
            const unsigned int qLen = (modelLen > 0)
                ? modelLen
                : static_cast<unsigned int>(std::max(1, dbEndOut - dbStartOut + 1));
            int qStartOut = (h.qStart >= 0) ? h.qStart : 0;
            int qEndOut = (h.qEnd >= 0) ? h.qEnd : static_cast<int>(qLen) - 1;
            if (qEndOut < qStartOut) qEndOut = qStartOut;
            const unsigned int alnLen = (h.cigarAlnLen > 0)
                ? static_cast<unsigned int>(h.cigarAlnLen)
                : static_cast<unsigned int>(std::max(1, std::abs(dbEndOut - dbStartOut) + 1));
            const unsigned int qSpan = static_cast<unsigned int>(qEndOut - qStartOut + 1);
            const unsigned int dbSpan = static_cast<unsigned int>(std::max(1, std::abs(dbEndOut - dbStartOut) + 1));
            const float qcov = (qLen > 0) ? static_cast<float>(qSpan) / static_cast<float>(qLen) : 0.0f;
            const float dbcov = (h.dbLen > 0) ? static_cast<float>(dbSpan) / static_cast<float>(h.dbLen) : 0.0f;
            const int bitScore = static_cast<int>(std::lrint(h.score));
            const double evalue = 1.0; // E-values out of scope
            const bool hasBacktrace = (!h.cigar.empty() && h.cigar != "NA");
            // Swap I/D (I=target-consume, D=query-consume) to MMseqs2 convention
            std::string emitCigar;
            if (hasBacktrace) {
                emitCigar.reserve(h.cigar.size());
                for (char c : h.cigar) {
                    if (c == 'I') emitCigar.push_back('D');
                    else if (c == 'D') emitCigar.push_back('I');
                    else emitCigar.push_back(c);
                }
            }
            if (dbStartOut > dbEndOut && qStartOut <= qEndOut) {
                std::swap(qStartOut, qEndOut);
                std::swap(dbStartOut, dbEndOut);
                if (hasBacktrace) emitCigar = reverseRleTraceOps(emitCigar);
            }
            float seqIdVal = h.seqId;
            if (seqIdVal < 0.0f) {
                const unsigned int bitScorePos = static_cast<unsigned int>(std::max(0, bitScore));
                seqIdVal = Matcher::estimateSeqIdByScorePerCol(
                    static_cast<uint16_t>(std::min(bitScorePos, 65535u)),
                    std::max(1u, alnLen), std::max(1u, alnLen));
            }
            Matcher::result_t res(h.dbKey, bitScore,
                                  std::min(1.0f, std::max(0.0f, qcov)),
                                  std::min(1.0f, std::max(0.0f, dbcov)),
                                  std::min(1.0f, std::max(0.0f, seqIdVal)),
                                  evalue, alnLen, qStartOut, qEndOut, qLen,
                                  dbStartOut, dbEndOut, h.dbLen,
                                  hasBacktrace ? emitCigar : std::string());
            const size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
            resultWriter.writeAdd(buffer, len, 0);
        }
        resultWriter.writeEnd(queryKey, 0);
}

int cmscan(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    cmInfernalGlobalInit();

    const size_t nThreads = (par.threads > 0) ? static_cast<size_t>(par.threads) : 1u;

    struct QueryModelRef { unsigned int key; size_t dbIdx; };
    std::vector<QueryModelRef> queryRefs;
    DBReader<unsigned int> *cmReader = NULL;
    std::string singleCmText;
    bool haveSingle = false;
    if (hasDbIndex(par.db1)) {
        cmReader = new DBReader<unsigned int>(par.db1.c_str(), (par.db1 + ".index").c_str(),
                                              static_cast<int>(nThreads),
                                              DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        cmReader->open(DBReader<unsigned int>::NOSORT);
        queryRefs.reserve(cmReader->getSize());
        for (size_t i = 0; i < cmReader->getSize(); ++i) {
            queryRefs.push_back(QueryModelRef{cmReader->getDbKey(i), i});
        }
        Debug(Debug::INFO) << "cmscan: CM database with " << queryRefs.size() << " models\n";
    } else {
        if (!looksLikeInfernalCm(par.db1)) {
            Debug(Debug::ERROR) << "cmscan requires an Infernal CM (run cmbuild first): " << par.db1 << "\n";
            return EXIT_FAILURE;
        }
        std::ifstream in(par.db1.c_str());
        std::stringstream ss;
        ss << in.rdbuf();
        singleCmText = ss.str();
        haveSingle = true;
        queryRefs.push_back(QueryModelRef{0, SIZE_MAX});
    }

    DBReader<unsigned int> seqDbr(par.db2.c_str(), (par.db2 + ".index").c_str(),
                                  static_cast<int>(nThreads),
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX
                                      | DBReader<unsigned int>::USE_LOOKUP);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    const unsigned int seqExt = DBReader<unsigned int>::getExtendedDbtype(seqDbr.getDbtype());
    const bool seqGpuDb = (seqExt & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    const bool decodeSeqDinuc = ((seqExt & Parameters::DBTYPE_EXTENDED_DINUCLEOTIDE) != 0) || seqGpuDb;
    const int seqDecodeType = effectiveDecodeSeqType(seqDbr.getDbtype(), decodeSeqDinuc);
    NucleotideMatrix nucMat(Parameters::getInstance().scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    BaseMatrix &subMat = static_cast<BaseMatrix &>(nucMat);

    DBReader<unsigned int> resultReader(par.db3.c_str(), (par.db3 + ".index").c_str(),
                                        static_cast<int>(nThreads),
                                        DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(),
                          static_cast<unsigned int>(nThreads), par.compressed,
                          Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    for (size_t qi = 0; qi < queryRefs.size(); ++qi) {
        const QueryModelRef &ref = queryRefs[qi];
        std::string cmText;
        if (ref.dbIdx != SIZE_MAX && cmReader != NULL) {
            const char *raw = cmReader->getData(ref.dbIdx, 0);
            size_t len = cmReader->getEntryLen(ref.dbIdx);
            cmText.assign(raw, len);
            const size_t nul = cmText.find('\0');
            if (nul != std::string::npos) cmText.resize(nul);
        } else if (haveSingle) {
            cmText = singleCmText;
        } else {
            continue;
        }

        const bool useInside = (par.cmMode == 1);
        int alignOpts = par.cmAlign ? CM_ALIGN_OPTACC : CM_ALIGN_CYK;
        if (par.cmAlignBanded == 0) alignOpts |= CM_ALIGN_NONBANDED;
        InfModel model = loadModel(
            cmText, useInside, alignOpts,
            par.cmAlignBanded != 0, par.cmLocal != 0,
            par.cmScanCutoff
        );
        if (!model.valid) {
            Debug(Debug::WARNING) << "cmscan: could not load model key " << ref.key << ", skipping\n";
            continue;
        }
        processQueryModel(par, model, ref.key, seqDbr, resultReader, resultWriter,
                          subMat, nucMat, seqDecodeType, seqGpuDb, decodeSeqDinuc, nThreads);
        freeModel(model);
    }

    resultWriter.close();
    resultReader.close();
    seqDbr.close();
    if (cmReader != NULL) { cmReader->close(); delete cmReader; }
    return EXIT_SUCCESS;
}

// compat
int cmsearch(int argc, const char **argv, const Command &command) {
    return cmscan(argc, argv, command);
}
