#include "LocalParameters.h"
#include "dinuc.out.h"

LocalParameters::LocalParameters() : Parameters(),
    PARAM_CM_REGION(PARAM_CM_REGION_ID, "--cm-region", "CM region flanking",
        "Extract target subregion around prefilter hit for CM alignment.\n"
        "Value is flanking fraction of CLEN added on each side (0.0 = disabled, use full target)",
        typeid(float), (void *) &cmRegionFlanking,
        "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CM_MODE(PARAM_CM_MODE_ID, "--cm-mode", "CM mode",
        "CM alignment mode: 0 = CYK, 1 = inside",
        typeid(int), (void *) &cmMode,
        "^[0-1]$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_DB_SIZE(PARAM_DB_SIZE_ID, "--db-size", "Database size",
        "Effective database size (0: use actual size)",
        typeid(size_t), (void *) &dbSize,
        "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN)
{
    cmRegionFlanking = 0.0f;
    cmMode = 0;
    dbSize = 0;

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

    splitstrand.push_back(&PARAM_STRAND);
    splitstrand.push_back(&PARAM_THREADS);
    splitstrand.push_back(&PARAM_COMPRESSED);
    splitstrand.push_back(&PARAM_V);

    // result2profile needs --strand so the RNA-corrected E-value gets the
    // both-strands doubling (matches MMseqs2 RNA fork behavior)
    result2profile.push_back(&PARAM_STRAND);
}
