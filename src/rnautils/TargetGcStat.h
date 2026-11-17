#ifndef RIBOSEEK_TARGET_GC_STAT_H
#define RIBOSEEK_TARGET_GC_STAT_H

#include "DBReader.h"
#include "Parameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "EvalueModel.h"
#include "RnaEvalueComputation.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

// A fixed-size sufficient statistic for the target database's GC composition:
// residues and sequence count at each of 256 GC levels. 4 KB regardless of
// database size, and enough to derive any bin count, any edge policy and any
// R0 afterwards.
//
// Why not store the finished 4-bin histogram: both the bin edges and
// R0_TARGET_LEN are model decisions, not database facts. R0 moves whenever the
// nulls are refit -- it has already caused a silent 4% error once -- so a
// database file carrying a baked histogram would quietly disagree with a later
// model with nothing to detect it. The statistic below is invariant to both.
struct TargetGcStat {
    // enum, not static const int: these are ODR-used (array bounds, comparisons
    // passed through templates) and an in-class static const would need an
    // out-of-line definition before C++17.
    // Version 2 adds dinucleotide counts per GC level. They fed the profile model's
    // lambda_KA input (S1f); the current model (S1n, 2026-10-20) does not use them,
    // but they stay in the format so existing files and tools keep working.
    // Version 1 files still load, with hasDinuc == false.
    enum { NLEVEL = 256, NSTATE = 16, FORMAT_VERSION = 2 };

    double residues[NLEVEL];     // residues whose sequence has this GC level
    double counts[NLEVEL];       // sequences at this GC level
    double dinuc[NLEVEL][NSTATE];// overlapping dinucleotides (profile column state order) at this level
    double totalResidues;        // ACGTU residues over the whole database
    size_t nSeq;                 // sequences contributing
    bool   sampled;              // true if built from a sample, not a full scan
    bool   hasDinuc;             // dinuc[][] filled (built here, or loaded from a v2 file)

    TargetGcStat() : totalResidues(0.0), nSeq(0), sampled(false), hasDinuc(true) {
        std::fill(residues, residues + NLEVEL, 0.0);
        std::fill(counts, counts + NLEVEL, 0.0);
        std::fill(&dinuc[0][0], &dinuc[0][0] + NLEVEL * NSTATE, 0.0);
    }

    static int level(double gc) {
        int l = (int)(gc * (NLEVEL - 1) + 0.5);
        return l < 0 ? 0 : (l >= NLEVEL ? NLEVEL - 1 : l);
    }
    static double levelGc(int l) { return (double)l / (double)(NLEVEL - 1); }

    void add(double gc, double nres, const double *dn = NULL) {
        const int l = level(gc);
        residues[l] += nres; counts[l] += 1.0; totalResidues += nres; nSeq++;
        if (dn != NULL) for (int a = 0; a < NSTATE; a++) dinuc[l][a] += dn[a];
    }

    // Overlapping dinucleotide counts of one raw-nucleotide sequence, in profile
    // column state order (RnaEvalueComputation::dinucState). A pair touching any
    // non-ACGTU character is skipped, as rawscores/repro/ka_lambda.target_background does.
    static void seqDinuc(const char *data, size_t len, double *dn) {
        std::fill(dn, dn + NSTATE, 0.0);
        int prev = -1;
        for (size_t j = 0; j < len; j++) {
            int n;
            switch (data[j]) {
                case 'A': case 'a': n = 0; break;
                case 'C': case 'c': n = 1; break;
                case 'G': case 'g': n = 2; break;
                case 'T': case 't': case 'U': case 'u': n = 3; break;
                default: n = -1; break;
            }
            if (prev >= 0 && n >= 0) dn[RnaEvalueComputation::dinucState(prev, n)] += 1.0;
            prev = n;
        }
    }

    // GC of one raw-nucleotide sequence. The target DB stores RAW characters --
    // despite the name, target_db_dinuc is not converted on disk; Sequence::
    // mapSequence does the dinucleotide mapping lazily when a sequence is read.
    // Reading these bytes through the dinucleotide tables instead would score
    // 'A' as the state whose dinucleotide is GC and report 0.56 for a 0.28 DB.
    // N and anything else count toward neither numerator nor denominator.
    static bool seqGc(const char *data, size_t len, double &gc, double &n) {
        size_t g = 0, m = 0;
        for (size_t j = 0; j < len; j++) {
            switch (data[j]) {
                case 'G': case 'g': case 'C': case 'c': g++; m++; break;
                case 'A': case 'a': case 'T': case 't':
                case 'U': case 'u':                     m++; break;
                default: break;
            }
        }
        if (m == 0) return false;
        gc = (double)g / (double)m; n = (double)m; return true;
    }

    // A GPU (padded) target DB, as written by `makepaddedseqdb`, does not hold raw
    // nucleotides: each position stores the dinucleotide code of the pair starting
    // there, and DBReader::getData (getUnpadded) hands it back as that code's matrix
    // LETTER (A C D E F G H I K L M N P Q R S T V W Y for codes 0-19; lowercase when
    // soft-masked; pairs assigned in DinucleotideMapping.cpp, T V W Y being the
    // last-base sentinels). Counting those letters as nucleotides gave GC 0.48-0.85
    // for a 0.34-0.63 DB and shifted every GPU E-value. Map each letter to the first
    // base of its pair so seqGc/seqDinuc see what they would for the unpadded DB
    // (exact on all 96,701 benchmark targets, 2026-10-15). Raw numeric codes (0-24,
    // +32 masked) are accepted too, in case a caller reads the data uncompressed.
    static bool isPadded(DBReader<unsigned int> *r) {
        return (DBReader<unsigned int>::getExtendedDbtype(r->getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    }
    static const char *decodePadded(const char *data, size_t len, std::string &buf) {
        // first base of the pair for codes 0-19 (letters A C D E F G H I K L M N P Q R S T V W Y)
        static const char firstBase[20] = {'G','A','C','U','C','A','U','G','C','A',
                                           'G','U','G','A','C','U','A','C','G','U'};
        static const signed char letterCode[26] = {
            /*A*/ 0, /*B*/-1, /*C*/ 1, /*D*/ 2, /*E*/ 3, /*F*/ 4, /*G*/ 5, /*H*/ 6, /*I*/ 7,
            /*J*/-1, /*K*/ 8, /*L*/ 9, /*M*/10, /*N*/11, /*O*/-1, /*P*/12, /*Q*/13, /*R*/14,
            /*S*/15, /*T*/16, /*U*/-1, /*V*/17, /*W*/18, /*X*/-1, /*Y*/19, /*Z*/-1 };
        buf.resize(len);
        for (size_t j = 0; j < len; j++) {
            unsigned char c = (unsigned char)data[j];
            int code = -1;
            if (c >= 'a' && c <= 'z') code = letterCode[c - 'a'];
            else if (c >= 'A' && c <= 'Z') code = letterCode[c - 'A'];
            else if (c < 64) code = (c >= 32) ? c - 32 : c;
            buf[j] = (code >= 0 && code < 20) ? firstBase[code] : 'N';
        }
        return buf.data();
    }

    // One raw-nucleotide sequence, exactly as seqGc + seqDinuc + add would count it, but in a
    // single table-driven, branch-free pass: createdb calls this on its one reader thread for
    // every sequence, where the two branchy loops above cost ~2.8 s per GB (2026-10-15). Counts are
    // integers, so the result is identical (checked byte-identical against rnagcstat).
    void addSequence(const char *data, size_t len) {
        if (data == NULL || len == 0) return;
        // Branch-free: code 4 marks a non-ACGTU byte; it is counted in its own slots and
        // dropped below, which is what skipping it (and any pair touching it) does.
        static const unsigned char *code = nucCodeTable();
        size_t nc[5] = {0}, pc[25] = {0};
        unsigned int prev = 4;
        for (size_t j = 0; j < len; j++) {
            const unsigned int n = code[(unsigned char)data[j]];
            nc[n]++; pc[prev * 5 + n]++;
            prev = n;
        }
        const size_t m = nc[0] + nc[1] + nc[2] + nc[3];
        if (m == 0) return;
        double dn[NSTATE];
        for (int x = 0; x < 4; x++)
            for (int y = 0; y < 4; y++) dn[RnaEvalueComputation::dinucState(x, y)] = (double)pc[x * 5 + y];
        add((double)(nc[1] + nc[2]) / (double)m, (double)m, dn);
    }
    // A C G T/U (either case) -> 0 1 2 3, anything else -> 4, as seqGc/seqDinuc classify them.
    struct NucCode {
        unsigned char t[256];
        NucCode() {
            for (int i = 0; i < 256; i++) t[i] = 4;
            t['A'] = t['a'] = 0; t['C'] = t['c'] = 1; t['G'] = t['g'] = 2;
            t['T'] = t['t'] = t['U'] = t['u'] = 3;
        }
    };
    static const unsigned char *nucCodeTable() {
        static const NucCode tbl;   // thread-safe one-time init (C++11 function-local static)
        return tbl.t;
    }

    // Full scan, sequential in key order so readahead works, parallel over
    // sequences. This is the expensive path; it is what `rnagcstat` pays once.
    static TargetGcStat fromFullScan(DBReader<unsigned int> *r, int threads) {
        TargetGcStat out;
        const bool padded = isPadded(r);
        std::vector<TargetGcStat> part((size_t)std::max(threads, 1));
#pragma omp parallel num_threads(threads)
        {
            unsigned int thread_idx = 0;
            unsigned int thread_count = 1;
#ifdef OPENMP
            thread_idx = (unsigned int)omp_get_thread_num();
            thread_count = (unsigned int)omp_get_num_threads();
#endif
            TargetGcStat &acc = part[thread_idx];
            std::string buf;
            // Balance residues, not sequence counts: a padded DB is sorted by length, so a
            // static split by count handed the last thread every long sequence. A fixed
            // range per thread also keeps the merged sums identical from run to run.
            size_t start = 0, count = 0;
            r->decomposeDomainByAminoAcid(thread_idx, thread_count, &start, &count);
            for (size_t id = start; id < start + count; id++) {
                const char *data = r->getData(id, thread_idx);
                const size_t len = r->getSeqLen(id);
                if (data == NULL || len == 0) continue;
                if (padded) data = decodePadded(data, len, buf);
                double gc, n, dn[NSTATE];
                if (seqGc(data, len, gc, n)) { seqDinuc(data, len, dn); acc.add(gc, n, dn); }
            }
        }
        for (size_t t = 0; t < part.size(); t++) out.merge(part[t]);
        out.sampled = false;
        return out;
    }

    // Strided sample, for when no precomputed statistic exists. Kept only as a
    // fallback: on a multi-TB database this is 200k scattered reads, each one
    // amplified to a 128 KB readahead, and -- because a sample that large is
    // certain to contain both GC tails -- it is also the case where min/max bin
    // edges collapse, which is why deriveHistogram uses quantiles instead.
    // threadIdx: the DBReader buffer to read through (padded and compressed DBs decode into a
    // per-thread buffer, so a caller running alongside other reader threads passes its own).
    static TargetGcStat fromSample(DBReader<unsigned int> *r, size_t maxScan, int threadIdx = 0) {
        TargetGcStat out;
        const size_t nSeq = r->getSize();
        const size_t step = (maxScan > 0 && nSeq > maxScan) ? (nSeq / maxScan) : 1;
        const bool padded = isPadded(r);
        std::string buf;
        for (size_t id = 0; id < nSeq; id += step) {
            const char *data = r->getData(id, threadIdx);
            const size_t len = r->getSeqLen(id);
            if (data == NULL || len == 0) continue;
            if (padded) data = decodePadded(data, len, buf);
            double gc, n, dn[NSTATE];
            if (seqGc(data, len, gc, n)) { seqDinuc(data, len, dn); out.add(gc, n, dn); }
        }
        out.sampled = (step > 1);
        return out;
    }

    void merge(const TargetGcStat &o) {
        for (int i = 0; i < NLEVEL; i++) {
            residues[i] += o.residues[i]; counts[i] += o.counts[i];
            for (int a = 0; a < NSTATE; a++) dinuc[i][a] += o.dinuc[i][a];
        }
        hasDinuc = hasDinuc && o.hasDinuc;
        totalResidues += o.totalResidues; nSeq += o.nSeq;
    }

    bool empty() const { return nSeq == 0 || totalResidues <= 0.0; }

    // The statistic for target database `db` (read through `r`, whose residue total is dbRes), as
    // rnaalign and the fragment realigner both need it (2026-10-20). Prefers a precomputed file
    // (createdb / makepaddedseqdb / rnagcstat write <db>.gcstat; with a createindex index the target
    // is read from <db>.idx, so <db minus .idx>.gcstat is tried as well), and samples the database
    // otherwise -- 200,000 scattered reads, each amplified to a 128 KB kernel readahead, which on a
    // multi-TB database is about 25 GB of disk traffic per search. A file whose residue total is
    // off by more than 2x belongs to another database and is ignored: it would skew every E-value.
    // The file describes the ORIGINAL target while `r` may be the chunked one, so only the residue
    // total is comparable, and loosely, since splitsequence chunks can overlap.
    // Returns an empty statistic when the database has no A/C/G/U residues; callers must stop.
    static TargetGcStat forTarget(const std::string &db, DBReader<unsigned int> *r, double dbRes,
                                  bool &precomputed, int threadIdx = 0) {
        std::string statPath = pathFor(db);
        const std::string idxSuffix = ".idx";
        if (FileUtil::fileExists(statPath.c_str()) == false && db.size() > idxSuffix.size()
            && db.compare(db.size() - idxSuffix.size(), idxSuffix.size(), idxSuffix) == 0) {
            statPath = pathFor(db.substr(0, db.size() - idxSuffix.size()));
        }
        TargetGcStat st;
        precomputed = st.load(statPath);
        if (precomputed) {
            const double ratio = (dbRes > 0.0) ? st.totalResidues / dbRes : 0.0;
            if (ratio < 0.5 || ratio > 2.0) {
                Debug(Debug::WARNING) << statPath << " holds " << st.totalResidues
                                      << " residues but the target database has "
                                      << dbRes << "; ignoring it and sampling instead. "
                                      << "Delete the file and re-run `riboseek rnagcstat`.\n";
                precomputed = false;
            }
        }
        if (precomputed == false) {
            const size_t MAX_SCAN = 200000;
            st = fromSample(r, MAX_SCAN, threadIdx);
        }
        return st;
    }

    // ---- the derivation both paths share ---------------------------------
    //
    // Bin edges are RESIDUE QUANTILES, not the observed min/max. Min/max edges
    // are outlier-fragile and get worse as the database grows: a 200k sample of
    // a 3 TB database contains a GC<0.20 and a GC>0.75 entry with probability
    // 1.000, which stretches the range to 0.05-0.95 and drops 92% of sequences
    // into a single bin -- the one-bin configuration measured at 3.08x spread
    // against 1.57x for four real bins.
    //
    // The bin's representative GC is the RESIDUE-weighted mean, consistent with
    // the weight being residues. Taking an unweighted mean over sequences would
    // let a bin holding one 10 kb sequence at GC 0.50 and a hundred 100 nt ones
    // at 0.42 report 0.421 while its residues are dominated by the 0.50 entry.
    //
    // dbResTotal is the database's true residue count, used to scale a sampled
    // statistic up to full size; pass 0 to use the statistic's own total.
    std::vector<std::pair<double, double> >
    deriveHistogram(int nbin, double r0, double dbResTotal = 0.0,
                    std::vector<std::vector<double> > *binBg = NULL) const {
        if (binBg != NULL) binBg->clear();
        std::vector<std::pair<double, double> > hist;
        if (empty() || nbin < 1 || r0 <= 0.0) return hist;

        // residue-quantile edges over the occupied levels
        std::vector<int> edge; edge.push_back(0);
        const double per = totalResidues / (double)nbin;
        double run = 0.0; int want = 1;
        for (int l = 0; l < NLEVEL && want < nbin; l++) {
            run += residues[l];
            while (want < nbin && run >= per * (double)want) { edge.push_back(l + 1); want++; }
        }
        while ((int)edge.size() < nbin) edge.push_back(NLEVEL);   // degenerate: all one GC
        edge.push_back(NLEVEL);

        const double scale = (dbResTotal > 0.0) ? dbResTotal / totalResidues : 1.0;
        for (int b = 0; b < nbin; b++) {
            const int lo = edge[b], hi = edge[b + 1];
            if (hi <= lo) continue;
            double res = 0.0, wsum = 0.0;
            std::vector<double> bg(NSTATE, 0.0);
            for (int l = lo; l < hi; l++) {
                res += residues[l]; wsum += levelGc(l) * residues[l];
                for (int a = 0; a < NSTATE; a++) bg[a] += dinuc[l][a];
            }
            if (res <= 0.0) continue;
            hist.push_back(std::make_pair(wsum / res, res * scale / r0));
            if (binBg != NULL) {
                double t = 0.0;
                for (int a = 0; a < NSTATE; a++) t += bg[a];
                for (int a = 0; a < NSTATE; a++) bg[a] = (t > 0.0) ? bg[a] / t : 0.0;
                binBg->push_back(bg);
            }
        }
        return hist;
    }

    // ---- persistence ------------------------------------------------------
    // Small versioned text file. 4 KB is nothing next to the database, and a
    // readable format means a mismatch can be diagnosed by eye.
    bool save(const std::string &path) const {
        const std::string tmp = path + ".tmp";
        FILE *f = fopen(tmp.c_str(), "w");
        if (f == NULL) return false;
        fprintf(f, "RIBOSEEK_GCSTAT %d\n", FORMAT_VERSION);
        fprintf(f, "nseq %zu\n", nSeq);
        fprintf(f, "residues %.17g\n", totalResidues);
        fprintf(f, "sampled %d\n", sampled ? 1 : 0);
        fprintf(f, "levels %d\n", NLEVEL);
        for (int l = 0; l < NLEVEL; l++) {
            if (residues[l] <= 0.0 && counts[l] <= 0.0) continue;
            fprintf(f, "%d %.17g %.17g", l, residues[l], counts[l]);
            for (int a = 0; a < NSTATE; a++) fprintf(f, " %.17g", dinuc[l][a]);
            fprintf(f, "\n");
        }
        if (fclose(f) != 0) { FileUtil::remove(tmp.c_str()); return false; }
        // atomic publish, so a concurrent reader never sees a half-written file
        if (rename(tmp.c_str(), path.c_str()) != 0) { FileUtil::remove(tmp.c_str()); return false; }
        return true;
    }

    bool load(const std::string &path) {
        FILE *f = fopen(path.c_str(), "r");
        if (f == NULL) return false;
        int ver = 0, nlev = 0, samp = 0;
        char tag[64];
        if (fscanf(f, "%63s %d", tag, &ver) != 2 || std::string(tag) != "RIBOSEEK_GCSTAT"
            || ver < 1 || ver > FORMAT_VERSION) { fclose(f); return false; }
        unsigned long long ns = 0;
        if (fscanf(f, " nseq %llu", &ns) != 1) { fclose(f); return false; }
        if (fscanf(f, " residues %lf", &totalResidues) != 1) { fclose(f); return false; }
        if (fscanf(f, " sampled %d", &samp) != 1) { fclose(f); return false; }
        if (fscanf(f, " levels %d", &nlev) != 1 || nlev != NLEVEL) { fclose(f); return false; }
        nSeq = (size_t)ns; sampled = (samp != 0);
        std::fill(residues, residues + NLEVEL, 0.0);
        std::fill(counts, counts + NLEVEL, 0.0);
        std::fill(&dinuc[0][0], &dinuc[0][0] + NLEVEL * NSTATE, 0.0);
        hasDinuc = (ver >= 2);
        int l; double r, c;
        while (fscanf(f, "%d %lf %lf", &l, &r, &c) == 3) {
            double dn[NSTATE];
            if (hasDinuc) {
                for (int a = 0; a < NSTATE; a++) {
                    if (fscanf(f, "%lf", &dn[a]) != 1) { fclose(f); return false; }
                }
            }
            if (l >= 0 && l < NLEVEL) {
                residues[l] = r; counts[l] = c;
                if (hasDinuc) for (int a = 0; a < NSTATE; a++) dinuc[l][a] = dn[a];
            }
        }
        fclose(f);
        return empty() == false;
    }

    static std::string pathFor(const std::string &dbName) { return dbName + ".gcstat"; }
};

#endif
