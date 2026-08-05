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

    // -s needs no override: the dinucleotide k-mer thresholds registered in
    // externalThreshold (RiboseekBase.cpp) make the MMseqs2 sensitivity automagic work on
    // the dinuc.out score scale.
    if (linear) {
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

    // MMseqs2's MultiParam parser silently turns a partial "--alph-size aa:13" into
    // aa:0,nucl:0, which makes the k-mer alphabet reduction corrupt the heap. Catch it
    // here instead of letting kmermatcher die.
    const int aaAlphabetSize = par.alphabetSize.values.aminoacid();
    if (par.PARAM_ALPH_SIZE.wasSet && (aaAlphabetSize < 2 || aaAlphabetSize > 25)) {
        Debug(Debug::ERROR) << "--alph-size aa:" << aaAlphabetSize << " is out of range for the 25 letter "
                            << "dinucleotide alphabet.\nPass both values, e.g. --alph-size aa:13,nucl:5.\n";
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
