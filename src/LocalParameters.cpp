#include "LocalParameters.h"
#include "dinuc.out.h"

#include <cfloat>

LocalParameters::LocalParameters() : Parameters(),
    PARAM_CM_REGION(PARAM_CM_REGION_ID, "--cm-region", "CM region flanking",
        "Extract target subregion around prefilter hit for CM alignment.\n"
        "Value is flanking fraction of prefilter alignment length (qend-qstart+1) added on each side (0.0 = disabled, use full target)",
        typeid(float), (void *) &cmRegionFlanking,
        "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CM_MODE(PARAM_CM_MODE_ID, "--cm-mode", "CM mode",
        "CM alignment mode: 0 = CYK, 1 = inside",
        typeid(int), (void *) &cmMode,
        "^[0-1]$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_DB_SIZE(PARAM_DB_SIZE_ID, "--db-size", "Database size",
        "Effective database size (0: use actual size)",
        typeid(size_t), (void *) &dbSize,
        "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CALIBRATE_CM(PARAM_CALIBRATE_CM_ID, "--calibrate-cm", "Calibrate CM",
        "Run Infernal cmcalibrate after cmbuild to fit E-value exp-tail parameters.\n"
        "Slow (~1-2s per query). Off by default; CmScan uses a native E-value fallback.",
        typeid(bool), (void *) &calibrateCm,
        "", MMseqsParameter::COMMAND_MISC),
    PARAM_CMLITE_MSA_EVAL(200001, "--cmlite-msa-eval", "CmLite MSA E-value",
        "Include only hits with <= this E-value when building the CmLite seed MSA/profile.\n"
        "All hits from resultDB are still realigned afterward; this only filters the profile-building subset.",
        typeid(double), (void *) &cmliteMsaEvalThr,
        "^(([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|([iI][nN][fF]))$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_LOLALIGN_MSA_EVAL(200002, "--lolalign-msa-eval", "LoLalign MSA E-value",
        "Include only hits with <= this E-value when building the LoLalign seed MSA/profile.\n"
        "All hits from resultDB are still realigned afterward; this only filters the profile-building subset.",
        typeid(double), (void *) &lolalignMsaEvalThr,
        "^(([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|([iI][nN][fF]))$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_LOLCMSEARCH_TMPDIR(200003, "--lolcmsearch-tmpdir", "LoLCMSearch tmp dir",
        "Directory to use for lolcmsearch intermediate databases.\n"
        "If set, this directory is created if needed and is not auto-deleted.",
        typeid(std::string), (void *) &lolcmsearchTmpDir,
        "", MMseqsParameter::COMMAND_ALIGN)
{
    cmRegionFlanking = 0.0f;
    cmMode = 0;
    dbSize = 0;
    calibrateCm = false;
    cmliteMsaEvalThr = DBL_MAX;
    lolalignMsaEvalThr = DBL_MAX;
    lolcmsearchTmpDir.clear();

    // Register dinuc.out as compiled-in matrix and set as default
    scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("dinuc.out", "dinuc.out"));
    seedScoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("dinuc.out", "dinuc.out"));
    substitutionMatrices.emplace_back("dinuc.out", dinuc_out, dinuc_out_len);

    // Match reference defaults for dinucleotide search
    alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(25, 5));
    maskMode = 0;

    align.push_back(&PARAM_CM_REGION);

    // rnaalign inherits the standard align parameters plus --db-size
    rnaalign = align;
    rnaalign.push_back(&PARAM_DB_SIZE);
    rnaalign.push_back(&PARAM_CMLITE_MSA_EVAL);
    rnaalign.push_back(&PARAM_LOLALIGN_MSA_EVAL);
    rnaalign.push_back(&PARAM_LOLCMSEARCH_TMPDIR);

    splitstrand.push_back(&PARAM_STRAND);
    splitstrand.push_back(&PARAM_THREADS);
    splitstrand.push_back(&PARAM_COMPRESSED);
    splitstrand.push_back(&PARAM_V);

    cmbuild.push_back(&PARAM_CALIBRATE_CM);
    // hhfilter/rMSA-style row filter on the cmbuild input MSA (off by default)
    cmbuild.push_back(&PARAM_FILTER_MSA);
    cmbuild.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    cmbuild.push_back(&PARAM_FILTER_QID);
    cmbuild.push_back(&PARAM_FILTER_QSC);
    cmbuild.push_back(&PARAM_FILTER_COV);
    cmbuild.push_back(&PARAM_FILTER_NDIFF);
    cmbuild.push_back(&PARAM_FILTER_MIN_ENABLE);
    cmbuild.push_back(&PARAM_CMLITE_MSA_EVAL);
    cmbuild.push_back(&PARAM_THREADS);
    cmbuild.push_back(&PARAM_COMPRESSED);
    cmbuild.push_back(&PARAM_V);

    // result2profile needs --strand so the RNA-corrected E-value gets the
    // both-strands doubling (matches MMseqs2 RNA fork behavior)
    result2profile.push_back(&PARAM_STRAND);
}
