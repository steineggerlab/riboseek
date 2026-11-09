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
        "^[0-1]$", MMseqsParameter::COMMAND_PROFILE),
    PARAM_LOCAL_MIN(PARAM_LOCAL_MIN_ID, "--local-min", "Whole-target length",
        "envalign: targets of at most this many residues are aligned whole, longer ones inside envelopes",
        typeid(int), (void *) &localMin, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_BAND(PARAM_CHAIN_BAND_ID, "--chain-band", "Chain diagonal band",
        "envalign: window hits whose diagonals differ by at most this much join one chain (indel tolerance)",
        typeid(int), (void *) &chainBand, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_GAP(PARAM_CHAIN_GAP_ID, "--chain-gap", "Chain query gap",
        "envalign: largest query gap between two window hits of one chain",
        typeid(int), (void *) &chainGap, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_ENV_MARGIN(PARAM_ENV_MARGIN_ID, "--env-margin", "Envelope margin",
        "envalign: residues added on each side of a chain before alignment; doubled while the alignment reaches a "
        "clipped envelope edge",
        typeid(int), (void *) &envMargin, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_ENV_MAX_ITER(PARAM_ENV_MAX_ITER_ID, "--env-max-iter", "Envelope widenings",
        "envalign: most margin doublings per envelope",
        typeid(int), (void *) &envMaxIter, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_WORD(PARAM_CHAIN_WORD_ID, "--chain-word", "Chain word length",
        "chainalign: exact word hits of this many letters seed the ungapped HSPs (1-12)",
        typeid(int), (void *) &chainWord, "^([1-9]|1[0-2])$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_XDROP(PARAM_CHAIN_XDROP_ID, "--chain-xdrop", "Chain ungapped X-drop",
        "chainalign: an HSP extension stops when its raw score falls this far below its best",
        typeid(int), (void *) &chainXdrop, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_MIN_HSP(PARAM_CHAIN_MIN_HSP_ID, "--chain-min-hsp", "Chain HSP score",
        "chainalign: HSPs below this raw score are dropped",
        typeid(int), (void *) &chainMinHsp, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_BRIDGE(PARAM_CHAIN_BRIDGE_ID, "--chain-bridge", "Chain bridge length",
        "chainalign: longest space (letters, either sequence) between two chained HSPs",
        typeid(int), (void *) &chainBridge, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_MAX_HSP(PARAM_CHAIN_MAX_HSP_ID, "--chain-max-hsp", "Chain HSPs per pair",
        "chainalign: the best this many HSPs of a pair are chained",
        typeid(int), (void *) &chainMaxHsp, "^[1-9][0-9]*$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_EXTEND(PARAM_CHAIN_EXTEND_ID, "--chain-extend", "Chain end extension X-drop",
        "chainalign: gapped X-drop extension (raw score) from both chain ends, as xdropalign; 0: no extension",
        typeid(int), (void *) &chainExtend, "^[0-9]+$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_CHAIN_EXTEND_EVAL(PARAM_CHAIN_EXTEND_EVAL_ID, "--chain-extend-eval", "Chain extension E-value",
        "chainalign: only chains with an E-value at most this before extension are extended",
        typeid(double), (void *) &chainExtendEval, "^([0-9eE+.-]+|[iI][nN][fF])$", MMseqsParameter::COMMAND_ALIGN),
    PARAM_KMER_ALPH(PARAM_KMER_ALPH_ID, "--kmer-alph", "Reduced k-mer alphabet",
        "redkmermatcher: k-mer alphabet built from the nucleotides of the dinucleotide db\n"
        "n4: nucleotides (k nt per k-mer)\nd16: dinucleotides (k + 1 nt)\nd12: dinucleotides, 12 classes (k + 1 nt)\n"
        "d8: dinucleotides, 8 classes (k + 1 nt)",
        typeid(std::string), (void *) &kmerAlph, "^(n4|d16|d12|d8)$", MMseqsParameter::COMMAND_CLUSTLINEAR)
{
    cmRegionFlanking = 1.0f;  // 1.0 => pad by full W (prior hardcoded behavior)
    cmMode = 1;
    cmAlign = 0;
    cmAlignBanded = 1;
    cmLocal = 0;        // config: glocal (current shipped default)
    dbSize = 0;
    cmliteMsaEvalThr = DBL_MAX;
    cmbuildEre = -1.0;      // <0 => Infernal default target rel-entropy
    cmbuildSymfrac = -1.0;  // <0 => keep all columns (current default)
    cmbuildNoss = 0;        // use SS_cons (current default)
    localMin = 5000;
    chainBand = 100;
    chainGap = 1000;
    envMargin = 200;
    envMaxIter = 8;
    chainWord = 8;
    chainXdrop = 30;
    chainMinHsp = 0;
    chainBridge = 500;
    chainMaxHsp = 500;
    chainExtend = 0;
    chainExtendEval = 1.0;
    kmerAlph = "n4";

    // Register dinuc.out as compiled-in matrix and set as default
    scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("dinuc.out", "dinuc.out"));
    seedScoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("dinuc.out", "dinuc.out"));
    substitutionMatrices.emplace_back("dinuc.out", dinuc_out, dinuc_out_len);

    // Match reference defaults for dinucleotide search
    alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(25, 5));
    maskMode = 0;

    cmscan = align;
    cmscan.push_back(&PARAM_CM_REGION);
    cmscan.push_back(&PARAM_CM_MODE);
    cmscan.push_back(&PARAM_CM_ALIGN);
    cmscan.push_back(&PARAM_CM_ALIGN_BANDED);
    cmscan.push_back(&PARAM_CM_LOCAL);

    rnaalign = align;

    envalign = align;
    envalign.push_back(&PARAM_LOCAL_MIN);
    envalign.push_back(&PARAM_CHAIN_BAND);
    envalign.push_back(&PARAM_CHAIN_GAP);
    envalign.push_back(&PARAM_ENV_MARGIN);
    envalign.push_back(&PARAM_ENV_MAX_ITER);

    chainalign = align;
    chainalign.push_back(&PARAM_CHAIN_WORD);
    chainalign.push_back(&PARAM_CHAIN_XDROP);
    chainalign.push_back(&PARAM_CHAIN_MIN_HSP);
    chainalign.push_back(&PARAM_CHAIN_BRIDGE);
    chainalign.push_back(&PARAM_CHAIN_MAX_HSP);
    chainalign.push_back(&PARAM_CHAIN_EXTEND);
    chainalign.push_back(&PARAM_CHAIN_EXTEND_EVAL);

    redkmermatcher = kmermatcher;
    redkmermatcher.push_back(&PARAM_KMER_ALPH);
    // linclust --kmer-alph: both kmermatcher levels run redkmermatcher (RnaCluster.cpp)
    linclustworkflow.push_back(&PARAM_KMER_ALPH);
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

    dinucdb.push_back(&PARAM_SUB_MAT);
    dinucdb.push_back(&PARAM_THREADS);
    dinucdb.push_back(&PARAM_COMPRESSED);
    dinucdb.push_back(&PARAM_V);

    // result2profile needs --strand so the RNA-corrected E-value gets the
    // both-strands doubling (matches MMseqs2 RNA fork behavior)
    result2profile.push_back(&PARAM_STRAND);
}
