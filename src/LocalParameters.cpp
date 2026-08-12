#include "LocalParameters.h"
#include "dinuc.out.h"

#include <cfloat>

LocalParameters::LocalParameters() : Parameters(),
    PARAM_CM_REGION(PARAM_CM_REGION_ID, "--cm-region", "CM alignment window",
        "CM alignment window: flanking pad on each side of the prefilter region, as a fraction of the\n"
        "model max hit length W.\n 1.0: full W (default, safest); smaller tightens the window (faster on\n"
        "long/genomic targets, may clip hits near the prefilter edge)\n 0: full target (no windowing)",
        typeid(float), (void *) &cmRegionFlanking,
        "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CM_SCAN_CUTOFF(PARAM_CM_SCAN_CUTOFF_ID, "--cm-scan-cutoff", "CM scan score cutoff",
        "Minimum bit score a window must reach in the CM scan to be considered for alignment",
        typeid(float), (void *) &cmScanCutoff,
        "^-?[0-9]*(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CM_MODE(PARAM_CM_MODE_ID, "--cm-mode", "CM scan scoring",
        "CM detection scoring: 0 = CYK (max, faster), 1 = Inside (sum, more sensitive)",
        typeid(int), (void *) &cmMode,
        "^[0-1]$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CM_ALIGN(PARAM_CM_ALIGN_ID, "--cm-align", "CM alignment method",
        "CM hit alignment:\n0: CYK (max-likelihood parse)\n1: Optimal Accuracy (posterior decoding, higher quality)",
        typeid(int), (void *) &cmAlign,
        "^[0-1]$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CM_ALIGN_BANDED(PARAM_CM_ALIGN_BANDED_ID, "--cm-align-banded", "CM alignment banding",
        "CM alignment banding:\n0: nonbanded (exact, fails >8GB on large models)\n1: HMM-banded",
        typeid(int), (void *) &cmAlignBanded,
        "^[0-1]$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CM_LOCAL(PARAM_CM_LOCAL_ID, "--cm-local", "CM alignment mode",
        "CM configuration:\n0: glocal (whole model must align; best for full-length rRNA)\n"
        "1: local (partial model matches via local begins/ends; better for fragmentary/divergent hits)",
        typeid(int), (void *) &cmLocal,
        "^[0-1]$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_DB_SIZE(PARAM_DB_SIZE_ID, "--db-size", "Database size",
        "Effective database size (0: use actual size)",
        typeid(size_t), (void *) &dbSize,
        "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CMLITE_MSA_EVAL(PARAM_CMLITE_MSA_EVAL_ID, "--cmlite-msa-eval", "CmLite MSA E-value",
        "Include only hits with <= this E-value when building the cmbuild seed CM.\n"
        "All hits from resultDB are still searched afterward; this only filters the CM-building subset.",
        typeid(double), (void *) &cmliteMsaEvalThr,
        "^([0-9eE+.-]+|[iI][nN][fF])$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CMBUILD_ERE(PARAM_CMBUILD_ERE_ID, "--cmbuild-ere", "cmbuild target relative entropy",
        "cmbuild: target mean match-state relative entropy in bits.\n"
        "Lower: more permissive/sensitive\nhigher: more specific\n<0: Infernal default (0.59 struct / 0.38 no-ss).",
        typeid(double), (void *) &cmbuildEre,
        "^-?[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
    PARAM_CMBUILD_SYMFRAC(PARAM_CMBUILD_SYMFRAC_ID, "--cmbuild-symfrac", "cmbuild match-column fraction",
        "cmbuild: fraction of non-gap residues for an MSA column to become a consensus/match column.\n<0: keep ALL columns (query-centric)\n[0,1]: occupancy-based.",
        typeid(double), (void *) &cmbuildSymfrac,
        "^-?[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_PROFILE),
    PARAM_CMBUILD_NOSS(PARAM_CMBUILD_NOSS_ID, "--cmbuild-noss", "cmbuild without secondary structure",
        "cmbuild:\n0: use SS_cons base pairs\n1: build a sequence-only CM (ignore/require-no structure).",
        typeid(int), (void *) &cmbuildNoss,
        "^[0-1]$", MMseqsParameter::COMMAND_PROFILE)
{
    cmRegionFlanking = 1.0f;  // 1.0 => pad by full W (prior hardcoded behavior)
    cmScanCutoff = 0.0f;
    cmMode = 1;
    cmAlign = 0;
    cmAlignBanded = 1;
    cmLocal = 0;        // config: glocal (current shipped default)
    dbSize = 0;
    cmliteMsaEvalThr = DBL_MAX;
    cmbuildEre = -1.0;      // <0 => Infernal default target rel-entropy
    cmbuildSymfrac = -1.0;  // <0 => keep all columns (current default)
    cmbuildNoss = 0;        // use SS_cons (current default)

    // Register dinuc.out as compiled-in matrix and set as default
    scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("dinuc.out", "dinuc.out"));
    seedScoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("dinuc.out", "dinuc.out"));
    substitutionMatrices.emplace_back("dinuc.out", dinuc_out, dinuc_out_len);

    // Match reference defaults for dinucleotide search
    alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(25, 5));
    maskMode = 0;

    cmscan = align;
    cmscan.push_back(&PARAM_CM_REGION);
    cmscan.push_back(&PARAM_CM_SCAN_CUTOFF);
    cmscan.push_back(&PARAM_CM_MODE);
    cmscan.push_back(&PARAM_CM_ALIGN);
    cmscan.push_back(&PARAM_CM_ALIGN_BANDED);
    cmscan.push_back(&PARAM_CM_LOCAL);

    rnaalign = align;
    rnaalign.push_back(&PARAM_DB_SIZE);

    splitstrand.push_back(&PARAM_STRAND);
    splitstrand.push_back(&PARAM_THREADS);
    splitstrand.push_back(&PARAM_COMPRESSED);
    splitstrand.push_back(&PARAM_V);

    // hhfilter/rMSA-style row filter on the cmbuild input MSA (off by default)
    cmbuild.push_back(&PARAM_FILTER_MSA);
    cmbuild.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    cmbuild.push_back(&PARAM_FILTER_QID);
    cmbuild.push_back(&PARAM_FILTER_QSC);
    cmbuild.push_back(&PARAM_FILTER_COV);
    cmbuild.push_back(&PARAM_FILTER_NDIFF);
    cmbuild.push_back(&PARAM_FILTER_MIN_ENABLE);
    cmbuild.push_back(&PARAM_CMLITE_MSA_EVAL);
    cmbuild.push_back(&PARAM_CMBUILD_ERE);
    cmbuild.push_back(&PARAM_CMBUILD_SYMFRAC);
    cmbuild.push_back(&PARAM_CMBUILD_NOSS);
    cmbuild.push_back(&PARAM_THREADS);
    cmbuild.push_back(&PARAM_COMPRESSED);
    cmbuild.push_back(&PARAM_V);

    // result2profile needs --strand so the RNA-corrected E-value gets the
    // both-strands doubling (matches MMseqs2 RNA fork behavior)
    result2profile.push_back(&PARAM_STRAND);
}
