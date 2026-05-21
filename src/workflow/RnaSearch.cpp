#include "DBReader.h"
#include "CommandCaller.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"
#include "LocalParameters.h"

#include "blastdigp.sh.h"

#include <iomanip>
#include <climits>
#include <cassert>

static void setRnaSearchDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 7.5;
    p->evalThr = 0.001;
    p->evalProfile = 0.1;
    p->searchType = Parameters::SEARCH_TYPE_NUCLEOTIDES;
    // Match old mmseqs pipeline defaults for dinucleotide search
    p->gapOpen = MultiParam<NuclAA<int>>(NuclAA<int>(23, 23));
    p->gapExtend = MultiParam<NuclAA<int>>(NuclAA<int>(1, 1));
    p->pca = MultiParam<PseudoCounts>(PseudoCounts(1.1, 1.1));
    p->pcb = MultiParam<PseudoCounts>(PseudoCounts(1.8, 1.8));
    p->maskProfile = 0;
    if (p->PARAM_MAX_SEQ_LEN.wasSet == false) {
        p->maxSeqLen = 10000;
    }
}

int rnasearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setRnaSearchDefaults(&par);
    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MIN_SEQ_ID.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, false, 0,
                        MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);

    std::string indexStr = PrefilteringIndexReader::searchForIndex(par.db2);
    const bool targetGpuDb =
        (DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db2.c_str()))
         & Parameters::DBTYPE_EXTENDED_GPU) != 0;
    const bool gpuWorkflow = (par.gpu == 1) || targetGpuDb;

    // RNA always uses both strands; mark as set so result2profile inherits --strand 2
    par.strand = 2;
    par.PARAM_STRAND.wasSet = true;

    // validate and set parameters for iterative search
    if (par.numIterations > 1) {
        par.addBacktrace = true;
    }

    par.printParameters(command.cmd, argc, argv, par.searchworkflow);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.searchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("VERB_COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    cmd.addVariable("GPU", targetGpuDb ? "TRUE" : NULL);

    // RNA always uses rnaalign
    cmd.addVariable("ALIGN_MODULE", "rnaalign");

    // GPU can only use the ungapped prefilter
    if (gpuWorkflow && par.PARAM_PREF_MODE.wasSet == false) {
        if (par.numIterations > 1
            || par.alignmentMode != Parameters::ALIGNMENT_MODE_SCORE_ONLY
            || par.altAlignment > 0
            || par.scoreBias != 0.0
            || par.realign == true
            || par.addBacktrace == true
            ) {
            par.prefMode = Parameters::PREF_MODE_UNGAPPED;
        } else {
            par.prefMode = Parameters::PREF_MODE_UNGAPPED_AND_GAPPED;
        }
    }

    switch (par.prefMode) {
        case Parameters::PREF_MODE_KMER:
            cmd.addVariable("PREFMODE", "KMER");
            break;
        case Parameters::PREF_MODE_UNGAPPED:
            cmd.addVariable("PREFMODE", "UNGAPPED");
            break;
        case Parameters::PREF_MODE_UNGAPPED_AND_GAPPED:
            cmd.addVariable("PREFMODE", "UNGAPPED_AND_GAPPED");
            break;
        case Parameters::PREF_MODE_EXHAUSTIVE:
            cmd.addVariable("PREFMODE", "EXHAUSTIVE");
            break;
    }

    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    std::string targetDB = (indexStr == "") ? par.db2.c_str() : indexStr.c_str();
    par.filenames[1] = targetDB;

    // Always use blastdigp.sh — it handles both single and multi-iteration,
    // and includes splitstrand, splitsequence, and offsetalignment
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    cmd.addVariable("SPLITSTRAND", "TRUE");
    cmd.addVariable("SPLITSEQUENCE_PAR", par.createParameterString(par.splitsequence).c_str());
    // Match fork's behavior: single-iter splits target (via blastdi.sh);
    // multi-iter skips target split (blastdigp.sh has the block commented out).
    if (indexStr == "" && par.numIterations == 1) {
        cmd.addVariable("NEEDTARGETSPLIT", "TRUE");
    }
    cmd.addVariable("NEEDQUERYSPLIT", "TRUE");
    cmd.addVariable("SPLIT_STRAND_PAR",
                    par.createParameterString(par.splitstrand).c_str());
    cmd.addVariable("OFFSETALIGNMENT_PAR",
                    par.createParameterString(par.offsetalignment).c_str());

    double originalEval = par.evalThr;
    if (par.numIterations > 1) {
        par.evalThr = (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    }

    for (int i = 0; i < par.numIterations; i++) {
        // Match old MMseqs2 nucl-nucl iterative search: realign disabled across
        // all iterations (Search.cpp lines 511-518 in the forked MMseqs2).
	if (i == 0) {
            par.realign = true;
	} else {
            par.realign = false;
	}

        if (i == (par.numIterations - 1)) {
            par.evalThr = originalEval;
        }

        if (par.prefMode == Parameters::PREF_MODE_KMER) {
            cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(),
                            par.createParameterString(par.prefilter).c_str());
        } else if (par.prefMode == Parameters::PREF_MODE_UNGAPPED) {
            cmd.addVariable(std::string("UNGAPPEDPREFILTER_PAR_" + SSTR(i)).c_str(),
                            par.createParameterString(par.ungappedprefilter).c_str());
        }

        cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(),
                        par.createParameterString(par.align).c_str());
        cmd.addVariable(std::string("PROFILE_PAR_" + SSTR(i)).c_str(),
                        par.createParameterString(par.result2profile).c_str());
    }

    FileUtil::writeFile(tmpDir + "/blastdigp.sh", blastdigp_sh, blastdigp_sh_len);
    std::string program = std::string(tmpDir + "/blastdigp.sh");

    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
