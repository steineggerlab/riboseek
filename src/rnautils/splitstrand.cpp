#include "Debug.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Orf.h"
#include "Sequence.h"

#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

int splitstrand(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    const int inputDbtype = reader.getDbtype();
    const bool inputIsSequence = Parameters::isEqualDbtype(inputDbtype, Parameters::DBTYPE_AMINO_ACIDS)
                                 || Parameters::isEqualDbtype(inputDbtype, Parameters::DBTYPE_NUCLEOTIDES);

    // Keep the base type: NUCLEOTIDES for sequence input, HMM_PROFILE for profile input.
    // Data on disk stays as nucleotide characters (ACGT).
    // The DINUCLEOTIDE flag triggers on-the-fly dinucleotide conversion when sequences are read.
    int outputDbtype = inputIsSequence ? Parameters::DBTYPE_NUCLEOTIDES : Parameters::DBTYPE_HMM_PROFILE;
    unsigned int extFlags = DBReader<unsigned int>::getExtendedDbtype(inputDbtype);
    // Only add DINUCLEOTIDE flag for sequence input — profile input is already in
    // dinucleotide space and should use standard mapProfile (matching old system).
    if (inputIsSequence) {
        extFlags |= LocalParameters::DBTYPE_EXTENDED_DINUCLEOTIDE;
    }
    extFlags |= LocalParameters::DBTYPE_EXTENDED_STRAND_SPLIT;
    outputDbtype = DBReader<unsigned int>::setExtendedDbtype(outputDbtype, extFlags);

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, outputDbtype);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        size_t querySize = 0;
        size_t queryFrom = 0;
        reader.decomposeDomainByAminoAcid(thread_idx, par.threads, &queryFrom, &querySize);
        if (querySize == 0) {
            queryFrom = 0;
        }

        char buffer[1024];
        size_t bufferLen;

        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t dataLen = reader.getEntryLen(i);  // includes null terminator
            size_t seqLen = reader.getSeqLen(i);

            // Forward strand (will get even key after renumbering)
            if (par.strand == 1 || par.strand == 2) {
                sequenceWriter.writeData(data, dataLen - 1, key, thread_idx);
                bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t>(0), seqLen - 1, 0, 0);
                headerWriter.writeData(buffer, bufferLen, key, thread_idx);
            }

            // Reverse complement strand (will get odd key after renumbering)
            if (par.strand == 0 || par.strand == 2) {
                // For sequence data, store the SAME forward-strand characters.
                // dinucBuildProfileFromSequence handles reversal of dinucleotide
                // position order + revcomp column swaps (matching extractqueryprofiles).
                // For profile data, reverse the profile entries.
                if (Parameters::isEqualDbtype(inputDbtype, Parameters::DBTYPE_HMM_PROFILE)) {
                    // Profile: reverse the order of PROFILE_READIN_SIZE-sized entries
                    std::string revData;
                    revData.reserve(seqLen * Sequence::PROFILE_READIN_SIZE);
                    for (int pos = (int)seqLen - 1; pos >= 0; pos--) {
                        revData.append(data + pos * Sequence::PROFILE_READIN_SIZE, Sequence::PROFILE_READIN_SIZE);
                    }
                    sequenceWriter.writeData(revData.c_str(), revData.length(), key, thread_idx);
                } else {
                    // Sequence: write the same forward-strand data for reverse entry.
                    // The dinucleotide encoding produces forward-strand numSequence,
                    // which dinucBuildProfileFromSequence then reverses to get the
                    // correct reversed position order (matching old extractqueryprofiles).
                    sequenceWriter.writeData(data, dataLen - 1, key, thread_idx);
                }
                bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 1, static_cast<size_t>(0), 0, 0);
                headerWriter.writeData(buffer, bufferLen, key, thread_idx);
            }
        }
    }
    headerWriter.close(true);
    sequenceWriter.close(true);
    reader.close();

    // Renumber so forward gets even keys (0, 2, 4, ...) and reverse gets odd keys (1, 3, 5, ...)
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                DBWriter::createRenumberedDB(par.hdr2, par.hdr2Index, "", "");
            }

#pragma omp task
            {
                DBWriter::createRenumberedDB(par.db2, par.db2Index, par.createLookup ? par.db1 : "", par.createLookup ? par.db1Index : "");
            }
        }
    }
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SOURCE);

    return EXIT_SUCCESS;
}
