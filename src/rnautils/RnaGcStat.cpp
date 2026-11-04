#include "LocalParameters.h"
#include "LocalCommandDeclarations.h"
#include "Debug.h"
#include "DBReader.h"
#include "FileUtil.h"
#include "TargetGcStat.h"
#include "EvalueModel.h"
#include "CommandDeclarations.h"
#include "CreatedbObserver.h"
#include "MakepaddedseqdbObserver.h"

#include <algorithm>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

// Precompute the target database's GC statistic, once, so that searches do not
// have to sample it. Writes <targetDB>.gcstat (about 4 KB).
//
// Without it, rnaalign falls back to a strided sample capped at 200,000
// sequences. On a multi-TB database that sample is 200,000 scattered reads,
// each amplified to a 128 KB kernel readahead -- about 25 GB pulled off disk to
// look at 200 MB of it -- and it is paid again on every search. It is also
// statistically worse the larger the database gets: a sample that size is
// certain to contain both GC tails, so min/max bin edges would span the whole
// axis. (deriveHistogram uses residue quantiles, which is robust to that, but
// a full scan is still the better input.)
int rnagcstat(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::NOSORT);

    const std::string out = TargetGcStat::pathFor(par.db1);
    if (FileUtil::fileExists(out.c_str())) {
        TargetGcStat old;
        if (old.load(out) && old.hasDinuc) {
            Debug(Debug::INFO) << out << " already exists (" << old.nSeq
                               << " sequences). Delete it to rebuild.\n";
            reader.close();
            return EXIT_SUCCESS;
        }
        Debug(Debug::WARNING) << out << " exists but could not be read or predates the dinucleotide "
                              << "counts the profile E-value model needs; rebuilding.\n";
    }

    Debug(Debug::INFO) << "Scanning " << reader.getSize() << " sequences for GC composition ...\n";
    TargetGcStat st = TargetGcStat::fromFullScan(&reader, par.threads);
    if (st.empty()) {
        Debug(Debug::ERROR) << "No usable nucleotide residues in " << par.db1 << ".\n";
        reader.close();
        return EXIT_FAILURE;
    }

    if (st.save(out) == false) {
        Debug(Debug::ERROR) << "Could not write " << out << ".\n";
        reader.close();
        return EXIT_FAILURE;
    }

    // Report the histogram this will actually produce, so a bad database is
    // visible here rather than as a silently skewed E-value later.
    std::vector<std::pair<double, double> > h =
        st.deriveHistogram(4, EvalueModel::R0_TARGET_LEN);
    Debug(Debug::INFO) << "Wrote " << out << ": " << st.nSeq << " sequences, "
                       << st.totalResidues << " residues.\n";
    Debug(Debug::INFO) << "GC histogram (weights are residues / "
                       << EvalueModel::R0_TARGET_LEN << "):";
    for (size_t b = 0; b < h.size(); b++) {
        Debug(Debug::INFO) << " [" << h[b].first << " x " << h[b].second << "]";
    }
    Debug(Debug::INFO) << "\n";

    reader.close();
    return EXIT_SUCCESS;
}

// `riboseek createdb` = MMseqs createdb, gathering the target GC statistic in the same pass
// (the CreatedbObserver hook sees each sequence as it is written) and writing <db>.gcstat for a
// nucleotide database, so it is ready for calibrated E-values without a separate step or a
// second read of the database (2026-10-15). Always rewritten: createdb has just replaced the
// database, so an existing .gcstat describes the old contents. A failure here only warns -- the
// search falls back to sampling -- and never fails the createdb itself.
// `riboseek base:createdb` still runs the plain MMseqs module.
namespace {
struct GcStatObserver : public CreatedbObserver {
    TargetGcStat st;
    void reset() { st = TargetGcStat(); }
    void addSequence(const char *seq, size_t len) { st.addSequence(seq, len); }
};
}

int riboseekCreatedb(int argc, const char **argv, const Command &command) {
    GcStatObserver obs;
    createdbObserver = &obs;
    const int ret = createdb(argc, argv, command);
    createdbObserver = NULL;
    if (ret != EXIT_SUCCESS) {
        return ret;
    }
    Parameters &par = Parameters::getInstance();
    const std::string db = par.filenames.back();
    const std::string out = TargetGcStat::pathFor(db);
    if (FileUtil::fileExists(out.c_str())) {
        FileUtil::remove(out.c_str());
    }
    const int dbtype = FileUtil::parseDbType(db.c_str());
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_NUCLEOTIDES) == false) {
        return EXIT_SUCCESS;
    }
    if (par.gpu) {
        // createdb --gpu pads single-nucleotide codes, not the dinucleotide codes riboseek's
        // GPU search reads (it scored ROC1 0.02 on the benchmark, 2026-10-15); the observer is
        // not fed in that mode, so there is nothing to write.
        Debug(Debug::WARNING) << "No GC statistic for a createdb --gpu database: riboseek's GPU search needs "
                              << "`riboseek createdb in.fa db` then `riboseek makepaddedseqdb db db_pad`.\n";
        return EXIT_SUCCESS;
    }
    obs.st.sampled = false;
    if (obs.st.empty()) {
        Debug(Debug::WARNING) << db << " has no A/C/G/U residues, so it has no GC statistic: a search against it "
                              << "cannot compute E-values and will stop.\n";
        return EXIT_SUCCESS;
    }
    if (obs.st.save(out) == false) {
        Debug(Debug::WARNING) << "Could not write the GC statistic " << out
                              << "; searches will sample the database instead. "
                              << "Run `riboseek rnagcstat " << db << "` to retry.\n";
        return EXIT_SUCCESS;
    }
    Debug(Debug::INFO) << "Wrote " << out << " (" << obs.st.nSeq << " sequences) for E-value calibration\n";
    return EXIT_SUCCESS;
}

// `riboseek makepaddedseqdb` = MMseqs makepaddedseqdb, then hand the padded (GPU) DB the
// source DB's GC statistic, so GPU searches load it instead of sampling the target on
// every run (2026-10-15). Padding re-orders and re-encodes the sequences but keeps every
// base, so the source's statistic is exact for the padded DB. A source without one (made
// by base:createdb, or before createdb wrote it) gets it gathered while it is padded --
// the MakepaddedseqdbObserver hook sees each raw source entry in the padding loop, so the
// database is read once (2026-11-04). The padded DB's own statistic is always rewritten so
// it can never describe an earlier source.
namespace {
struct PaddedGcStatObserver : public MakepaddedseqdbObserver {
    std::vector<TargetGcStat> part;
    bool active = false;
    void begin(unsigned int threads) {
        // Decided here, after makepaddedseqdb has parsed its parameters: gather only when
        // the source has no usable statistic to copy and holds raw nucleotides. A padded
        // source hands back dinucleotide letters, which addSequence would misread; it is
        // left to the full scan below, which decodes them.
        Parameters &par = Parameters::getInstance();
        const std::string src = TargetGcStat::pathFor(par.db1);
        TargetGcStat old;
        if (FileUtil::fileExists(src.c_str()) && old.load(src) && old.hasDinuc) {
            return;
        }
        const int dbtype = FileUtil::parseDbType(par.db1.c_str());
        if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_NUCLEOTIDES) == false
            || (DBReader<unsigned int>::getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_GPU) != 0) {
            return;
        }
        unsigned int slots = std::max(threads, 1u);
#ifdef OPENMP
        slots = std::max(slots, (unsigned int)omp_get_max_threads());
#endif
        part.assign(slots, TargetGcStat());
        active = true;
    }
    void addSequence(unsigned int threadIdx, const char *seq, size_t len) {
        if (active) {
            part[threadIdx].addSequence(seq, len);
        }
    }
};
}

int riboseekMakepaddedseqdb(int argc, const char **argv, const Command &command) {
    PaddedGcStatObserver obs;
    makepaddedseqdbObserver = &obs;
    const int ret = makepaddedseqdb(argc, argv, command);
    makepaddedseqdbObserver = NULL;
    if (ret != EXIT_SUCCESS) {
        return ret;
    }
    Parameters &par = Parameters::getInstance();
    const std::string src = TargetGcStat::pathFor(par.db1);
    const std::string out = TargetGcStat::pathFor(par.db2);
    if (FileUtil::fileExists(out.c_str())) {
        FileUtil::remove(out.c_str());
    }
    TargetGcStat st;
    if (obs.active) {
        // per-thread parts merged in thread order: the ranges are fixed per thread count,
        // so the statistic is the same from run to run
        for (size_t t = 0; t < obs.part.size(); t++) {
            st.merge(obs.part[t]);
        }
        st.sampled = false;
    } else {
        TargetGcStat old;
        if (FileUtil::fileExists(src.c_str()) && old.load(src) && old.hasDinuc) {
            FileUtil::copyFile(src, out);
            Debug(Debug::INFO) << "Copied " << src << " to " << out << " for E-value calibration\n";
            return EXIT_SUCCESS;
        }
        if (Parameters::isEqualDbtype(FileUtil::parseDbType(par.db1.c_str()), Parameters::DBTYPE_NUCLEOTIDES) == false) {
            return EXIT_SUCCESS;
        }
        DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                      DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        reader.open(DBReader<unsigned int>::NOSORT);
        st = TargetGcStat::fromFullScan(&reader, par.threads);
        reader.close();
    }
    if (st.empty()) {
        Debug(Debug::WARNING) << par.db1 << " has no A/C/G/U residues, so it has no GC statistic: a search against "
                              << par.db2 << " cannot compute E-values and will stop.\n";
        return EXIT_SUCCESS;
    }
    if (st.save(out) == false) {
        Debug(Debug::WARNING) << "Could not write the GC statistic " << out
                              << "; GPU searches will sample the database instead.\n";
        return EXIT_SUCCESS;
    }
    Debug(Debug::INFO) << "Wrote " << out << " (" << st.nSeq << " sequences, from " << par.db1
                       << ") for E-value calibration\n";
    return EXIT_SUCCESS;
}
