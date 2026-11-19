#include "LocalParameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <climits>
#include <string>

#ifdef OPENMP
#include <omp.h>
#endif

// Materialize the dinucleotide encoding that DinucleotideMapping.cpp otherwise applies
// on the fly: every position i becomes the letter of the pair (seq[i], seq[i+1]), the
// last position pairs with an unknown nucleotide. The output is a plain amino acid
// database over the 25 letter dinucleotide alphabet of dinuc.out, so all modules that
// score raw characters (rescorediagonal, align2clust, ...) work unchanged.
// --strand (default forward, keys preserved): 0 encodes the reverse complement under the same keys, 2 writes both
// strands, the forward one under key 2k and the reverse complement under 2k+1 (riboseek cluster's both-strand queries).
static char complementNucleotide(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': case 'U': return 'A';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': case 'u': return 'a';
        case 'R': return 'Y'; case 'Y': return 'R'; case 'K': return 'M'; case 'M': return 'K';
        case 'B': return 'V'; case 'V': return 'B'; case 'D': return 'H'; case 'H': return 'D';
        case 'r': return 'y'; case 'y': return 'r'; case 'k': return 'm'; case 'm': return 'k';
        case 'b': return 'v'; case 'v': return 'b'; case 'd': return 'h'; case 'h': return 'd';
        default: return c;
    }
}

static void encodeDinucleotides(const char *data, size_t seqLen, const unsigned char *aa2num, const char *num2aa,
                                std::string &buffer) {
    buffer.clear();
    buffer.reserve(seqLen + 1);
    for (size_t pos = 0; pos + 1 < seqLen; ++pos) {
        unsigned short pair = static_cast<unsigned short>(
                (static_cast<unsigned char>(data[pos]) << 8) | static_cast<unsigned char>(data[pos + 1]));
        buffer.push_back(num2aa[aa2num[pair]]);
    }
    unsigned short lastPair = static_cast<unsigned short>(
            (static_cast<unsigned char>(data[seqLen - 1]) << 8) | static_cast<unsigned char>('X'));
    buffer.push_back(num2aa[aa2num[lastPair]]);
    buffer.push_back('\n');
}

int dinucdb(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    if (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES) == false) {
        Debug(Debug::ERROR) << "Input " << par.db1 << " is not a nucleotide sequence database\n";
        reader.close();
        return EXIT_FAILURE;
    }

    // The setup callback registered by registerDinucleotideMapping() fills the 16 bit
    // (first << 8 | second) entries of aa2num while the matrix is constructed.
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0f);
    const unsigned char *aa2num = subMat.aa2num;
    const char *num2aa = subMat.num2aa;

    const int strand = par.PARAM_STRAND.wasSet ? par.strand : 1;
    if (strand == 2 && reader.getLastKey() > (UINT_MAX - 1) / 2) {
        Debug(Debug::ERROR) << "--strand 2 writes keys 2k and 2k+1, the largest key of " << par.db1 << " is too large\n";
        reader.close();
        return EXIT_FAILURE;
    }

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed,
                    Parameters::DBTYPE_AMINO_ACIDS);
    writer.open();

    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        std::string reverse;
        buffer.reserve(4096);

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char *data = reader.getData(i, thread_idx);
            size_t seqLen = reader.getSeqLen(i);
            const unsigned int forwardKey = (strand == 2) ? 2 * key : key;
            const unsigned int reverseKey = (strand == 2) ? 2 * key + 1 : key;
            if (seqLen == 0) {
                if (strand != 0) {
                    writer.writeData("\n", 1, forwardKey, thread_idx);
                }
                if (strand != 1) {
                    writer.writeData("\n", 1, reverseKey, thread_idx);
                }
                continue;
            }

            if (strand != 0) {
                encodeDinucleotides(data, seqLen, aa2num, num2aa, buffer);
                writer.writeData(buffer.c_str(), buffer.size(), forwardKey, thread_idx);
            }
            if (strand != 1) {
                reverse.resize(seqLen);
                for (size_t pos = 0; pos < seqLen; ++pos) {
                    reverse[pos] = complementNucleotide(data[seqLen - 1 - pos]);
                }
                encodeDinucleotides(reverse.c_str(), seqLen, aa2num, num2aa, buffer);
                writer.writeData(buffer.c_str(), buffer.size(), reverseKey, thread_idx);
            }
        }
    }
    writer.close(true);
    reader.close();

    // keys are preserved, so headers/lookup/source of the input stay valid (not for --strand 2: keys 2k, 2k+1)
    if (strand != 2) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);
    }

    return EXIT_SUCCESS;
}
