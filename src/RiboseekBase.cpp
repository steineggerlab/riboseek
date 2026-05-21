#include "LocalParameters.h"
#include "LocalCommandDeclarations.h"
#include "DownloadDatabase.h"
#include "Prefiltering.h"

LocalParameters& localPar = LocalParameters::getLocalInstance();

void updateValidation() {}
void (*validatorUpdate)(void) = updateValidation;

std::vector<KmerThreshold> externalThreshold = {};
std::vector<DatabaseDownload> externalDownloads = {};

std::vector<Command> riboseekCommands = {
    {"splitstrand", splitstrand, &localPar.splitstrand, COMMAND_SEQUENCE,
            "Split database into forward and reverse complement strand entries",
            NULL,
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryDB> <o:strandDB>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"strandDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
    {"rnaalign", rnaalign, &localPar.rnaalign, COMMAND_ALIGNMENT,
            "RNA-aware gapped local alignment with dinucleotide scoring",
            NULL,
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
    {"rnasearch", rnasearch, &localPar.searchworkflow, COMMAND_MAIN,
            "RNA homology search with dinucleotide profiles",
            "# Search RNA query against RNA target database\n"
            "riboseek rnasearch queryDB targetDB resultDB tmp\n\n"
            "# Iterative profile search\n"
            "riboseek rnasearch queryDB targetDB resultDB tmp --num-iterations 3\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryDB> <i:targetDB> <o:resultDB> <tmpDir>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                      {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
    {"lolalign", lolalign, &localPar.rnaalign, COMMAND_ALIGNMENT,
            "LoL-style iterative local RNA profile/structure realignment",
            "riboseek lolalign queryDB targetDB resultDB alignmentDB --cm-region 1.5 -e inf --lolalign-msa-eval 1e-3\n"
            "Environment overrides:\n"
            "  MMSEQS_LOLALIGN_STRUCTURE_BACKEND=auto|rnafold|linearfold|canonical\n"
            "  MMSEQS_LOLALIGN_ITERATIONS=N\n"
            "  MMSEQS_LOLALIGN_REL_WEIGHT=F\n"
            "  MMSEQS_LOLALIGN_GUIDE_WEIGHT=F\n",
            "OpenAI Codex",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}}},
    {"lolcmsearch", lolcmsearch, &localPar.rnaalign, COMMAND_ALIGNMENT,
            "Full CM CYK refinement with LoLalign-derived state bands",
            "riboseek lolcmsearch queryDB targetDB resultDB alignmentDB --cm-region 1.5 -e inf --lolalign-msa-eval 1e-3\n"
            "targetDB and resultDB may also be comma-separated lists of equal length.\n"
            "Environment overrides:\n"
            "  MMSEQS_LOLCMSEARCH_BAND_PAD=N\n"
            "  MMSEQS_LOLCMSEARCH_KEEP_TMP=1\n"
            "  MMSEQS_LOLCMSEARCH_SKIP_BUILD=1\n"
            "  MMSEQS_LOLCMSEARCH_SKIP_ROUGH=1\n"
            "Argument:\n"
            "  --lolcmsearch-tmpdir DIR\n",
            "OpenAI Codex",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}}},
    {"cmsearch", cmsearch, &localPar.align, COMMAND_ALIGNMENT,
            "CM search with in-tree CYK/Inside dynamic programming",
            "riboseek cmsearch queryCMDB targetDB resultDB alignmentDB\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryCMDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
            CITATION_MMSEQS2, {{"queryCMDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
    {"cmscan", cmscan, &localPar.align, COMMAND_HIDDEN,
            "Compatibility alias for cmsearch",
            "riboseek cmscan queryCMDB targetDB resultDB alignmentDB\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryCMDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
            CITATION_MMSEQS2, {{"queryCMDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
    {"generatecm", generatecm, &localPar.result2profile, COMMAND_PROFILE,
            "Generate CM-like profile DB from query/target/result DBs",
            "riboseek generatecm queryDB targetDB resultDB queryCM\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:queryCM>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"queryCM", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::profileDb }}},
    {"cmbuild", cmbuild, &localPar.cmbuild, COMMAND_PROFILE,
            "Build CM from aligned Stockholm/FASTA input",
            "riboseek cmbuild queryDB targetDB resultDB outputCMDB\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:outputCMDB>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"outputCMDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb }}},
};
