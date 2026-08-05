#include "LocalParameters.h"
#include "LocalCommandDeclarations.h"
#include "DownloadDatabase.h"
#include "Prefiltering.h"

LocalParameters& localPar = LocalParameters::getLocalInstance();

void updateValidation() {}
void (*validatorUpdate)(void) = updateValidation;

// Dinucleotide k-mer thresholds for the prefilter. dinuc.out has a nearly flat diagonal
// (~19-24 per dinucleotide at the seeding bit factor) instead of the wide BLOSUM62 spread,
// so MMseqs2's built-in thresholds land in the wrong range: below -s ~2 they exceed the
// self-score of every dinucleotide k-mer and the index stays completely empty. base and
// sensPerStep are anchored so that -s 1 still indexes every k-mer (measured saturation:
// 95 / 116 / 135 for k = 5 / 6 / 7) while -s 7.5 keeps the threshold MMseqs2 would use.
std::vector<KmerThreshold> externalThreshold = {
    {Parameters::DBTYPE_AMINO_ACIDS, 5,  99.6f, 4.60f},
    {Parameters::DBTYPE_AMINO_ACIDS, 6, 119.0f, 3.03f},
    {Parameters::DBTYPE_AMINO_ACIDS, 7, 140.1f, 5.08f},
};
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
    {"cluster", riboseekCluster, &localPar.clusterworkflow, COMMAND_MAIN,
            "Cascaded RNA clustering in dinucleotide space",
            "# Cascaded clustering of a riboseek createdb database\n"
            "riboseek cluster sequenceDB clusterDB tmp\n\n"
            "# Sequence identity and coverage are measured on dinucleotides:\n"
            "# a single nucleotide mismatch breaks the two dinucleotides it is part of,\n"
            "# so id_nucleotide is roughly (1 + id_dinucleotide) / 2\n"
            "riboseek cluster sequenceDB clusterDB tmp --min-seq-id 0.8\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:sequenceDB> <o:clusterDB> <tmpDir>",
            CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                      {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
    {"linclust", riboseekLinclust, &localPar.linclustworkflow, COMMAND_MAIN,
            "Linear-time RNA clustering in dinucleotide space",
            "# Linear-time clustering of a riboseek createdb database\n"
            "riboseek linclust sequenceDB clusterDB tmp\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:sequenceDB> <o:clusterDB> <tmpDir>",
            CITATION_MMSEQS2|CITATION_LINCLUST, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                      {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                      {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
    {"dinucdb", dinucdb, &localPar.dinucdb, COMMAND_SEQUENCE,
            "Encode a nucleotide database in the dinucleotide alphabet",
            "# Amino acid typed database over the 25 dinucleotide letters of dinuc.out.\n"
            "# Keys are preserved, so results carry over to the nucleotide database.\n"
            "riboseek dinucdb sequenceDB dinucDB\n",
            "Martin Steinegger <martin.steinegger@snu.ac.kr>",
            "<i:sequenceDB> <o:dinucSeqDB>",
            CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                      {"dinucSeqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::aaDb }}},
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
