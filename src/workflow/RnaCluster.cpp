#include "LocalParameters.h"
#include "LocalCommandDeclarations.h"
#include "BaseMatrix.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include "rnaclust.sh.h"

#include <algorithm>
#include <cassert>

// Dinucleotide k-mers are indexed positionally in base (alphabetSize - 1) = 24, so
// 24^k has to fit into a size_t: k = 14 spans 15 nucleotides and is the maximum,
// k = 15 overflows.
static const int DINUC_KMER_SIZE = 14;

// Mirrors setAutomaticThreshold() of MMseqs2's cluster workflow, but with a floor: the
// k-mer score thresholds behind -s are calibrated on BLOSUM62 statistics, and on the
// dinucleotide matrix the protein setting for high identity (-s 1) makes the prefilter
// miss even identical sequence pairs.
static float dinucSensitivity(float seqIdThr) {
    float sens;
    if (seqIdThr <= 0.3f) {
        sens = 6.0f;
    } else if (seqIdThr > 0.8f) {
        sens = 1.0f;
    } else {
        sens = 1.0f + 10.0f * (0.7f - seqIdThr);
    }
    return std::max(4.0f, sens);
}

static bool contains(const std::vector<MMseqsParameter*> &list, const MMseqsParameter *par) {
    for (size_t i = 0; i < list.size(); i++) {
        if (list[i]->uniqid == par->uniqid) {
            return true;
        }
    }
    return false;
}

// Parameters riboseek resolves itself for dinucleotide space. An explicit user setting
// always wins, everything else keeps the MMseqs2 protein defaults.
static std::string dinucParameterString(LocalParameters &par, bool linear,
                                        const std::vector<MMseqsParameter*> &workflow) {
    std::vector<MMseqsParameter*> owned;

    // full 25 letter dinucleotide alphabet: its letters are concrete base pairs, so
    // the mutual-information alphabet reduction MMseqs2 uses for proteins does not apply
    if (par.PARAM_ALPH_SIZE.wasSet == false) {
        par.alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(25, 5));
    }
    owned.push_back(&par.PARAM_ALPH_SIZE);

    // dinuc.out is scaled like a protein matrix (match ~+5.6 bits), the gap costs match
    // the ones riboseek search uses with it
    if (par.PARAM_GAP_OPEN.wasSet == false) {
        par.gapOpen = MultiParam<NuclAA<int>>(NuclAA<int>(23, 5));
    }
    owned.push_back(&par.PARAM_GAP_OPEN);
    if (par.PARAM_GAP_EXTEND.wasSet == false) {
        par.gapExtend = MultiParam<NuclAA<int>>(NuclAA<int>(1, 2));
    }
    owned.push_back(&par.PARAM_GAP_EXTEND);

    // rRNA/ncRNA sequences are long, so seed as densely as the nucleotide path does
    // instead of the ~20 k-mers per sequence the protein path picks
    if (par.PARAM_KMER_PER_SEQ.wasSet == false) {
        par.kmersPerSequence = 60;
    }
    owned.push_back(&par.PARAM_KMER_PER_SEQ);
    if (par.PARAM_KMER_PER_SEQ_SCALE.wasSet == false) {
        par.kmersPerSequenceScale = MultiParam<NuclAA<float>>(NuclAA<float>(0.2f, 0.2f));
    }
    owned.push_back(&par.PARAM_KMER_PER_SEQ_SCALE);

    if (linear == false) {
        if (par.PARAM_S.wasSet == false) {
            par.sensitivity = dinucSensitivity(par.seqIdThr);
        }
        owned.push_back(&par.PARAM_S);
        if (par.clusterVersion == Parameters::CLUSTER_VERSION2) {
            Debug(Debug::WARNING) << "--cluster-version 2 derives the prefilter sensitivity from --min-seq-id itself, "
                                  << "which is calibrated for proteins and under-clusters in dinucleotide space.\n";
        }
    } else {
        // k is only used by kmermatcher here; the cascaded workflow needs it unset so
        // that prefilter can pick its own (much smaller) k-mer size
        if (par.PARAM_K.wasSet == false) {
            par.kmerSize = DINUC_KMER_SIZE;
        }
        if (par.kmerSize > DINUC_KMER_SIZE && par.alphabetSize.values.aminoacid() >= 21) {
            Debug(Debug::WARNING) << "k-mer size " << par.kmerSize << " overflows the dinucleotide k-mer index, using "
                                  << DINUC_KMER_SIZE << " instead\n";
            par.kmerSize = DINUC_KMER_SIZE;
        }
        owned.push_back(&par.PARAM_K);
    }

    // forward everything the caller set explicitly, it overrides the values above
    std::vector<MMseqsParameter*> forward;
    for (size_t i = 0; i < workflow.size(); i++) {
        if (workflow[i]->wasSet == false || contains(owned, workflow[i])) {
            continue;
        }
        forward.push_back(workflow[i]);
    }

    return par.createParameterString(owned) + par.createParameterString(forward);
}

static int rnaClusterWorkflow(int argc, const char **argv, const Command &command, bool linear) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    std::vector<MMseqsParameter*> &workflow = linear ? par.linclustworkflow : par.clusterworkflow;
    par.parseParameters(argc, argv, command, false, 0, 0);

    const int dbType = FileUtil::parseDbType(par.db1.c_str());
    // nucleotide input gets encoded into the dinucleotide alphabet first; an amino acid
    // input is already dinucleotide encoded (cascaded clustering calls linclust on the
    // converted database)
    const bool convert = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);
    if (convert == false && Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_AMINO_ACIDS) == false) {
        Debug(Debug::ERROR) << "Input " << par.db1 << " is not a sequence database\n";
        return EXIT_FAILURE;
    }

    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, workflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    const std::string clustPar = dinucParameterString(par, linear, workflow);
    if (convert) {
        Debug(Debug::INFO) << "Clustering in dinucleotide space ("
                           << BaseMatrix::unserializeName(par.scoringMatrixFile.values.aminoacid().c_str())
                           << ", 25 letters). Sequence identity is measured on dinucleotides.\n";
    }

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("CONVERT", convert ? "TRUE" : NULL);
    cmd.addVariable("DINUCDB_PAR", par.createParameterString(par.dinucdb).c_str());
    cmd.addVariable("CLUST_MODULE", linear ? "linclustbase" : "clusterbase");
    cmd.addVariable("CLUST_PAR", clustPar.c_str());

    std::string program = tmpDir + "/rnaclust.sh";
    FileUtil::writeFile(program, rnaclust_sh, rnaclust_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}

int riboseekCluster(int argc, const char **argv, const Command &command) {
    return rnaClusterWorkflow(argc, argv, command, false);
}

int riboseekLinclust(int argc, const char **argv, const Command &command) {
    return rnaClusterWorkflow(argc, argv, command, true);
}
