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
        buffer.reserve(4096);

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char *data = reader.getData(i, thread_idx);
            size_t seqLen = reader.getSeqLen(i);
            if (seqLen == 0) {
                writer.writeData("\n", 1, key, thread_idx);
                continue;
            }

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

            writer.writeData(buffer.c_str(), buffer.size(), key, thread_idx);
        }
    }
    writer.close(true);
    reader.close();

    // keys are preserved, so headers/lookup/source of the input stay valid
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);

    return EXIT_SUCCESS;
}
