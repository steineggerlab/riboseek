#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Util.h"
#include "itoa.h"

#include "Orf.h"

#include <unistd.h>
#include <climits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

int splitsequence(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.maxSeqLen = 10000;
    par.sequenceOverlap = 300;
    par.parseParameters(argc, argv, command, true, 0, 0);
    int mode = DBReader<unsigned int>::USE_INDEX;
    if (par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_HARD) {
        mode |= DBReader<unsigned int>::USE_DATA;
    }
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, mode);
    reader.open(DBReader<unsigned int>::NOSORT);
    bool sizeLarger = false;
    for (size_t i = 0; i < reader.getSize(); i++) {
        sizeLarger |= (reader.getSeqLen(i) > par.maxSeqLen);
    }

    // if no sequence needs to be splitted
    if (sizeLarger == false) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_DB);
        reader.close();
        return EXIT_SUCCESS;
    }

    DBReader<unsigned int> headerReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::NOSORT);

    if (par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && par.compressed == true) {
        Debug(Debug::WARNING) << "Sequence split mode (--sequence-split-mode 0) and compressed (--compressed 1) can not be combined.\nTurn compressed to 0";
        par.compressed = 0;
    }

    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    size_t sequenceOverlap = par.sequenceOverlap;

    // makepaddedseqdb pads every sequence to a multiple of ALIGN and asserts that each
    // entry starts on an ALIGN-byte boundary; the GPU (marv) kernels rely on it for
    // vectorised loads. Splitting such a DB writes entries at offset + split*step, so
    // the step must preserve that alignment or the kernel faults with
    // "CUDA error: misaligned address". Round the step DOWN to a multiple of ALIGN,
    // which only ever widens the requested overlap, never narrows it.
    const size_t PADDED_DB_ALIGN = 4;  // keep in sync with ALIGN in makepaddedseqdb.cpp
#ifdef RIBOSEEK
    const char PADDED_DB_PAD_SYMBOL = 24;
#else
    const char PADDED_DB_PAD_SYMBOL = 20;
#endif
    const bool isPaddedDb =
        (DBReader<unsigned int>::getExtendedDbtype(reader.getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    size_t splitStep = par.maxSeqLen - sequenceOverlap;
    if (isPaddedDb) {
        size_t alignedStep = splitStep & ~(PADDED_DB_ALIGN - 1);
        if (alignedStep == 0) {
            Debug(Debug::ERROR) << "--max-seq-len minus --sequence-overlap must be at least "
                                << PADDED_DB_ALIGN << " for a GPU (padded) database.\n";
            EXIT(EXIT_FAILURE);
        }
        if (alignedStep != splitStep) {
            Debug(Debug::INFO) << "GPU database: reducing split step from " << splitStep << " to "
                               << alignedStep << " to keep entries " << PADDED_DB_ALIGN
                               << "-byte aligned (effective overlap "
                               << (par.maxSeqLen - alignedStep) << ").\n";
        }
        splitStep = alignedStep;
    }
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

        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data=NULL;
            if (par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_HARD) {
                // getData() decodes a padded (GPU) database back to letters via
                // DBReader::getUnpadded(); a hard split must copy the raw encoded bytes
                // instead, or the result is no longer a valid GPU database.
                data = isPaddedDb ? reader.getDataUncompressed(i)
                                  : reader.getData(i, thread_idx);
            }
            size_t seqLen = reader.getSeqLen(i);
            char* header = headerReader.getData(i, thread_idx);
            size_t headerLen = headerReader.getEntryLen(i) - 1;
            Orf::SequenceLocation loc;
            loc.id = UINT_MAX;
            loc.strand = Orf::STRAND_PLUS;
            size_t from = 0;
            unsigned int dbKey = key;
            if (par.headerSplitMode == 0) {
                loc = Orf::parseOrfHeader(header);
                if (loc.id != UINT_MAX) {
                    from = (loc.strand == Orf::STRAND_MINUS) ? loc.to : loc.from;
                    dbKey = loc.id;
                }
            }
            size_t splitCnt = (size_t) ceilf(static_cast<float>(seqLen) / static_cast<float>(splitStep));

            for (size_t split = 0; split < splitCnt; split++) {
                size_t startPos = split * splitStep;
                size_t len = std::min(par.maxSeqLen, seqLen - startPos);
                if (par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
                    // +2 to emulate the \n\0
                    sequenceWriter.writeIndexEntry(key, reader.getOffset(i) + startPos, len+2, thread_idx);
                } else if (isPaddedDb) {
                    // A hard split writes entries back to back, which would destroy the
                    // ALIGN-byte start alignment makepaddedseqdb guarantees and the GPU
                    // kernels rely on. Mirror makepaddedseqdb: pad the written payload up
                    // to ALIGN, but record only the logical length (len + \n\0) in the
                    // index so the padding is never read as sequence.
                    // makepaddedseqdb stores sequence bytes padded up to ALIGN with the
                    // pad symbol and NO newline/null; the "+2" in the index length is a
                    // convention DBReader subtracts back off. Reproduce that exactly.
                    const size_t padding = (len % PADDED_DB_ALIGN == 0)
                                         ? 0 : PADDED_DB_ALIGN - (len % PADDED_DB_ALIGN);
                    std::string chunk;
                    chunk.reserve(len + padding);
                    chunk.append(data + startPos, len);
                    chunk.append(padding, PADDED_DB_PAD_SYMBOL);
                    sequenceWriter.writeData(chunk.c_str(), chunk.size(), key, thread_idx, false, false);
                    sequenceWriter.writeIndexEntry(key, sequenceWriter.getStart(thread_idx), len + 2, thread_idx);
                } else {
                    sequenceWriter.writeStart(thread_idx);
                    sequenceWriter.writeAdd(data + startPos, len, thread_idx);
                    char newLine = '\n';
                    sequenceWriter.writeAdd(&newLine, 1, thread_idx);
                    sequenceWriter.writeEnd(key, thread_idx, true);
                }

                if (par.headerSplitMode == 0) {
                    size_t fromPos = from + startPos;
                    size_t toPos = (from + startPos) + (len - 1);
                    if (loc.id != UINT_MAX && loc.strand == Orf::STRAND_MINUS) {
                        fromPos = (seqLen - 1) - (from + startPos);
                        toPos   = fromPos - std::min(fromPos, len);
                    }

                    size_t bufferLen = Orf::writeOrfHeader(buffer, dbKey, fromPos, toPos, 0, 0);
                    headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                } else {
                    headerWriter.writeData(header, headerLen, key, thread_idx);
                }
            }
        }
    }
    headerWriter.close(true);
    sequenceWriter.close(true);
    headerReader.close();
    reader.close();
    if (par.sequenceSplitMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::DATA);
    }
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

