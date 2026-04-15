#include "LocalCommandDeclarations.h"
#include "CommandDeclarations.h"
#include "Debug.h"
#include "DBReader.h"
#include "MMseqsMPI.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "RNAFoldBridge.h"
#include "Matcher.h"
#include "MultipleAlignment.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#ifdef HAVE_INFERNAL_BRIDGE
#include "infernal/InfernalBridge.h"
#endif

#ifdef OPENMP
#include <omp.h>
#endif

#include <cctype>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Native generatecm command:
// maps DB -> DB covariance-model generation to MMseqs profile construction.
int generatecm(int argc, const char **argv, const Command &command) {
    return result2profile(argc, argv, command);
}

namespace {

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

static std::string buildStockholmText(const std::string &id,
                                      const std::vector<AlnSeq> &seqs,
                                      const std::string &ssCons) {
    std::ostringstream out;
    out << "# STOCKHOLM 1.0\n";
    out << "#=GF ID " << (id.empty() ? "mmseqs_model" : id) << "\n";
    for (size_t i = 0; i < seqs.size(); ++i) {
        out << seqs[i].id << " " << seqs[i].aln << "\n";
    }
    if (!ssCons.empty()) {
        out << "#=GC SS_cons " << ssCons << "\n";
    }
    out << "//\n";
    return out.str();
}

} // namespace

int cmbuild(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

#ifndef HAVE_INFERNAL_BRIDGE
    Debug(Debug::ERROR) << "cmbuild requires Infernal bridge support\n";
    return EXIT_FAILURE;
#else

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

    const unsigned int maxSequenceLength = std::max(tDbr->getMaxSeqLen(), qDbr.getMaxSeqLen());
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);

    DBWriter cmWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads,
                      par.compressed, Parameters::DBTYPE_GENERIC_DB);
    cmWriter.open();

    Debug(Debug::INFO) << "Query database size: " << qDbr.getSize() << " type: " << qDbr.getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << tDbr->getDbTypeName() << "\n";

    Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        Sequence centerSequence(maxSequenceLength, qDbr.getDbtype(), &subMat, 0, false, false);
        Sequence edgeSequence(maxSequenceLength, tDbr->getDbtype(), &subMat, 0, false, false);
        MultipleAlignment aligner(maxSequenceLength, &subMat);

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);
        std::vector<Matcher::result_t> filteredResults;
        filteredResults.reserve(300);
        std::vector<std::vector<unsigned char>> seqSet;
        seqSet.reserve(300);
        std::vector<unsigned int> seqKeys;
        seqKeys.reserve(300);
        std::set<unsigned int> seenKeys;

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            alnResults.clear();
            filteredResults.clear();
            seqSet.clear();
            seqKeys.clear();
            seenKeys.clear();

            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr.getId(queryKey);
            if (queryId == UINT_MAX) {
                Debug(Debug::WARNING) << "cmbuild: invalid query " << queryKey << "\n";
                continue;
            }

            centerSequence.mapSequence(queryId, queryKey,
                                       qDbr.getData(queryId, thread_idx),
                                       qDbr.getSeqLen(queryId));

            // Parse alignment results and map target sequences
            Matcher::readAlignmentResults(alnResults, resultReader.getData(id, thread_idx), false);
            for (size_t i = 0; i < alnResults.size(); i++) {
                Matcher::result_t &res = alnResults[i];
                if (res.backtrace.empty()) {
                    continue;
                }
                const size_t targetId = tDbr->getId(res.dbKey);
                if (targetId == UINT_MAX) {
                    continue;
                }
                // Skip duplicate target keys (Infernal requires unique names)
                if (seenKeys.find(res.dbKey) != seenKeys.end()) {
                    continue;
                }
                seenKeys.insert(res.dbKey);
                edgeSequence.mapSequence(targetId, res.dbKey,
                                         tDbr->getData(targetId, thread_idx),
                                         tDbr->getSeqLen(targetId));
                seqSet.emplace_back(edgeSequence.numSequence,
                                    edgeSequence.numSequence + edgeSequence.L);
                seqKeys.emplace_back(res.dbKey);
                filteredResults.push_back(res);
            }

            if (filteredResults.empty()) {
                continue;
            }

            // Compute MSA using MultipleAlignment (same as result2msa)
            MultipleAlignment::MSAResult msaRes = aligner.computeMSA(
                &centerSequence, seqSet, filteredResults, true);

            // Build consensus sequence for ViennaRNA SS prediction
            std::string consensus;
            consensus.reserve(msaRes.centerLength);
            for (size_t pos = 0; pos < msaRes.centerLength; pos++) {
                char aa = msaRes.msaSequence[0][pos];
                char c = (aa < MultipleAlignment::GAP) ? subMat.num2aa[(int)aa] : 'N';
                c = normalizeBase(c);
                consensus.push_back(c);
            }

            // Predict secondary structure with ViennaRNA
            std::string ssCons;
            double mfe = 0.0;
            rnaFoldPredictDotBracket(consensus, ssCons, &mfe);

            // Build Stockholm text from MSA, filtering very gappy sequences
            // Infernal's cm_parsetree_Doctor() can fail with many short fragments
            const float minColCoverage = 0.30f;
            std::vector<AlnSeq> stoSeqs;
            stoSeqs.reserve(msaRes.setSize);
            for (size_t i = 0; i < msaRes.setSize; i++) {
                AlnSeq as;
                if (i == 0) {
                    as.id = "query_" + std::to_string(queryKey);
                } else {
                    as.id = "t" + std::to_string(seqKeys[i - 1]);
                }
                as.aln.reserve(msaRes.centerLength);
                size_t residueCount = 0;
                for (size_t pos = 0; pos < msaRes.centerLength; pos++) {
                    char aa = msaRes.msaSequence[i][pos];
                    if (aa < MultipleAlignment::GAP) {
                        char c = normalizeBase(subMat.num2aa[(int)aa]);
                        if (c == 'A' || c == 'C' || c == 'G' || c == 'U') {
                            as.aln.push_back(c);
                        } else {
                            as.aln.push_back('N');
                        }
                        residueCount++;
                    } else {
                        as.aln.push_back('-');
                    }
                }
                // Always keep the query (i == 0); filter gappy targets
                if (i == 0 || (msaRes.centerLength > 0 &&
                    static_cast<float>(residueCount) / static_cast<float>(msaRes.centerLength) >= minColCoverage)) {
                    stoSeqs.push_back(as);
                }
            }

            std::string stoText = buildStockholmText(
                "query_" + std::to_string(queryKey), stoSeqs, ssCons);

            // Call Infernal bridge to build CM (not thread-safe, needs critical section)
            std::string cmText;
            std::string infernalErr;
            bool success = false;
#pragma omp critical(infernal_bridge)
            {
                success = InfernalBridge::buildCmFromStockholmText(stoText, cmText, infernalErr);
            }

            MultipleAlignment::deleteMSA(&msaRes);

            if (success) {
                cmWriter.writeData(cmText.c_str(), cmText.size(), queryKey, thread_idx);
            } else {
                Debug(Debug::WARNING) << "cmbuild: Infernal failed for query " << queryKey
                                      << ": " << infernalErr << "\n";
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
#endif // HAVE_INFERNAL_BRIDGE
}
