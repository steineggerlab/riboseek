#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Util.h"
#include "Orf.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Prefiltering.h"
#include "MathUtil.h"
#include "PSSMCalculator.h"
#include "MultipleAlignment.h"

#include <unistd.h>
#include <climits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

void toBuffer(const char* pssm, const unsigned char* sequence, const unsigned char* consensus, const float* neffM, size_t seqLen, std::string& result) {
    size_t currPos = 0;
    for (size_t pos = 0; pos < seqLen; pos++) {
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            result.push_back(pssm[currPos + aa]);
        }
        result.push_back(static_cast<unsigned char>(sequence[pos]));
        result.push_back(static_cast<unsigned char>(consensus[pos]));
        result.push_back(static_cast<unsigned char>(MathUtil::convertNeffToChar(neffM[pos])));
        currPos += Sequence::PROFILE_READIN_SIZE;
    }
}

int extractqueryprofiles(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    const int inputDbtype = reader.getDbtype();

    unsigned int maxSeqLength = reader.getMaxSeqLen();
    // for SIMD memory alignment
    maxSeqLength = (maxSeqLength) / (VECSIZE_INT * 4) + 2;
    maxSeqLength *= (VECSIZE_INT * 4);

    int outputDbtype = Parameters::DBTYPE_HMM_PROFILE;
    if (par.pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC) {
        outputDbtype = DBReader<unsigned int>::setExtendedDbtype(outputDbtype, Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS);
    }
    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, outputDbtype);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);

    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        Sequence seq(maxSeqLength + 1, inputDbtype, &subMat, 0, false, par.compBiasCorrection != 0);
        size_t querySize = 0;
        size_t queryFrom = 0;
        reader.decomposeDomainByAminoAcid(thread_idx, par.threads, &queryFrom, &querySize);
        if (querySize == 0) {
            queryFrom = 0;
        }

        char* profile_char = new char[(maxSeqLength + 1) * Sequence::PROFILE_READIN_SIZE];
        float* neffM = new float[maxSeqLength + 1];
        // Initialize with 1.0
        for (size_t i = 0; i < maxSeqLength + 1; i++) {
            neffM[i] = 1.0f;
        }
        float* compositionBias = new float[maxSeqLength + 1];

        char buffer[1024];
        std::string result;
        result.reserve((par.maxSeqLen + 1) * Sequence::PROFILE_READIN_SIZE * sizeof(char));
        size_t bufferLen;

        // Currently only think about single sequence input, so maxSetSize is 1
        PSSMCalculator calculator(
            &subMat, maxSeqLength + 1, 1, par.pcmode, par.pca, par.pcb
        );

        char* data_cpy = (char*)malloc(sizeof(char) * maxSeqLength * Sequence::PROFILE_READIN_SIZE);
        
        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t seqLen = reader.getSeqLen(i);

            seq.mapSequence(i, key, data, seqLen);

            if (Parameters::isEqualDbtype(inputDbtype, Parameters::DBTYPE_HMM_PROFILE)) {
                // Copy the data
                memcpy(data_cpy, data, sizeof(char) * seqLen * Sequence::PROFILE_READIN_SIZE);
                if (par.strand == 1 || par.strand == 2) {
                    // Nothing to do for the forward strand
                    toBuffer(data_cpy, seq.numSequence, seq.numConsensusSequence, seq.neffM, seqLen, result);
                    sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);
                    result.clear();

                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t>(0), seqLen - 1, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                }

                if (par.strand == 0 || par.strand == 2) {
                    // Do the same thing with the reversed profile
                    // Reverse the numSequence, numConsensusSequence, seq.neffM, data_cpy
                    std::reverse(seq.numSequence, seq.numSequence + seqLen);
                    std::reverse(seq.numConsensusSequence, seq.numConsensusSequence + seqLen);
                    std::reverse(seq.neffM, seq.neffM + seqLen);
                    char tmpPssm[Sequence::PROFILE_READIN_SIZE];
                    int i_curr = 0;
                    int j_curr = (seqLen - 1) * Sequence::PROFILE_READIN_SIZE;
                    for (size_t pos = 0; pos < seqLen/2; pos++) {
                        memcpy(&tmpPssm[0], data_cpy + i_curr, Sequence::PROFILE_READIN_SIZE * sizeof(char));
                        memcpy(data_cpy + i_curr, data_cpy + j_curr, Sequence::PROFILE_READIN_SIZE * sizeof(char));
                        memcpy(data_cpy + j_curr, &tmpPssm[0], Sequence::PROFILE_READIN_SIZE * sizeof(char));
                        i_curr += Sequence::PROFILE_READIN_SIZE;
                        j_curr -= Sequence::PROFILE_READIN_SIZE;
                    }
                    size_t currPos = 0, l = 0;
                    // for (size_t pos = 0; pos < seqLen; pos++) {
                    while (l < seqLen && l < maxSeqLength) {
                        // Swap following position pairs: (1,15), (2,6), (4,12), (5,7), (8,9), (10,11), => canonical dinucleotides
                        //                                (16,23), (17,22), (18,21), (19,20) => non-canonical dinucleotides
                        std::swap(data_cpy[currPos + 1], data_cpy[currPos + 15]);
                        std::swap(data_cpy[currPos + 2], data_cpy[currPos + 6]);
                        std::swap(data_cpy[currPos + 4], data_cpy[currPos + 12]);
                        std::swap(data_cpy[currPos + 5], data_cpy[currPos + 7]);
                        std::swap(data_cpy[currPos + 8], data_cpy[currPos + 9]);
                        std::swap(data_cpy[currPos + 10], data_cpy[currPos + 11]);
                        std::swap(data_cpy[currPos + 16], data_cpy[currPos + 23]);
                        std::swap(data_cpy[currPos + 17], data_cpy[currPos + 22]);
                        std::swap(data_cpy[currPos + 18], data_cpy[currPos + 21]);
                        std::swap(data_cpy[currPos + 19], data_cpy[currPos + 20]);
                        currPos += Sequence::PROFILE_READIN_SIZE;
                        l++;
                    }
                    toBuffer(data_cpy, seq.numSequence, seq.numConsensusSequence, seq.neffM, seqLen, result);
                    sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);
                    
                    bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 1, static_cast<size_t>(0), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    result.clear();
                }
            } else {
                // Get local composition bias
                if (par.compBiasCorrection == true) {
                    SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, seq.numSequence, seqLen, compositionBias, par.compBiasCorrectionScale);
                } else {
                    memset(compositionBias, 0, sizeof(float) * (seqLen + 1));
                }
                // Make a basic profile
                for (size_t pos = 0; pos < seqLen; pos++) {
                    unsigned int aaIdx = seq.numSequence[pos];
                    for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                        profile_char[pos * Sequence::PROFILE_READIN_SIZE + aa_num] = subMat.subMatrix[aaIdx][aa_num] + static_cast<char>((compositionBias[pos]) > 0.0 ?
                                                    compositionBias[pos]+0.5 : compositionBias[pos]-0.5);
                    }
                }
                // seq.numConsensusSequence same as seq.numSequence for single sequence
                memcpy(seq.numConsensusSequence, seq.numSequence, sizeof(unsigned char) * seqLen);
                if (par.strand == 1 || par.strand == 2) {
                    // seq.neffM is 1.0 for single sequence
                    toBuffer(profile_char, seq.numSequence, seq.numConsensusSequence, neffM, seqLen, result);
                    sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);

                    bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t>(0), seqLen - 1, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    result.clear();
                }
                
                if (par.strand == 0 || par.strand == 2) {
                    // Now do the same thing for the reversed profile
                    // Reverse the seq.numSequence and profile_char
                    std::reverse(seq.numSequence, seq.numSequence + seqLen);
                    std::reverse(seq.numConsensusSequence, seq.numConsensusSequence + seqLen);
                    char tmpPssm[Sequence::PROFILE_READIN_SIZE];
                    int i_curr = 0;
                    int j_curr = (seqLen - 1) * Sequence::PROFILE_READIN_SIZE;
                    for (size_t pos = 0; pos < seqLen/2; pos++) {
                        memcpy(&tmpPssm[0], profile_char + i_curr, Sequence::PROFILE_READIN_SIZE * sizeof(char));
                        memcpy(profile_char + i_curr, profile_char + j_curr, Sequence::PROFILE_READIN_SIZE * sizeof(char));
                        memcpy(profile_char + j_curr, &tmpPssm[0], Sequence::PROFILE_READIN_SIZE * sizeof(char));
                        i_curr += Sequence::PROFILE_READIN_SIZE;
                        j_curr -= Sequence::PROFILE_READIN_SIZE;
                    }
                    for (size_t pos = 0; pos < seqLen; pos++) {
                        // Swap following position pairs: (1,15), (2,6), (4,12), (5,7), (8,9), (10,11), => canonical dinucleotides
                        //                                (16,23), (17,22), (18,21), (19,20) => non-canonical dinucleotides
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 1], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 15]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 2], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 6]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 4], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 12]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 5], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 7]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 8], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 9]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 10], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 11]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 16], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 23]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 17], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 22]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 18], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 21]);
                        std::swap(profile_char[pos * Sequence::PROFILE_READIN_SIZE + 19], profile_char[pos * Sequence::PROFILE_READIN_SIZE + 20]);
                    }
                    toBuffer(profile_char, seq.numSequence, seq.numConsensusSequence, neffM, seqLen, result);
                    sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);

                    bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 1, static_cast<size_t>(0), 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                    result.clear();
                }
            }
        }
        delete[] profile_char;
        delete[] neffM;
        delete[] compositionBias;
        free(data_cpy);
    }
    headerWriter.close(true);
    sequenceWriter.close(true);
    reader.close();


    // make identifiers stable
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

