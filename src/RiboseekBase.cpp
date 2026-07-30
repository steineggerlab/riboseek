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
    {"search", riboseekSearch, &localPar.searchworkflow, COMMAND_MAIN,
            "RNA homology search with dinucleotide profiles",
            "# Search RNA query against RNA target database\n"
            "riboseek search queryDB targetDB resultDB tmp\n\n"
            "# Iterative profile search\n"
            "riboseek search queryDB targetDB resultDB tmp --num-iterations 3\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryDB> <i:targetDB> <o:resultDB> <tmpDir>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                      {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
    {"cmsearch", cmsearch, &localPar.cmscan, COMMAND_ALIGNMENT,
            "CM search with in-tree CYK/Inside dynamic programming",
            "riboseek cmsearch queryCMDB targetDB resultDB alignmentDB\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryCMDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
            CITATION_MMSEQS2, {{"queryCMDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
    {"cmscan", cmscan, &localPar.cmscan, COMMAND_HIDDEN,
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
            "riboseek cmbuild queryDB targetDB resultDB outputCMDB\n"
            "targetDB and resultDB may be comma-separated lists of equal length.\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:queryDB> <i:targetDB> <i:resultDB> <o:outputCMDB>",
            CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                      {"outputCMDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb }}},
};
