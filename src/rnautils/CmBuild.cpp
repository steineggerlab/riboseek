extern "C" {
#include "infernal.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_vectorops.h"
#include "esl_bitfield.h"
}

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <vector>

// One-time global init of logsum tables
void cmInfernalGlobalInit() {
    static std::once_flag once;
    std::call_once(once, []() {
        // SIMD dispatch (static inline; harmless to call)
        impl_Init();
        // integer Inside logsum table
        init_ilogsum();
        // float logsum table
        FLogsumInit();
         // hmmer p7 float table
        p7_FLogsumInit();
    });
}

void flattenInsertEmissions(CM_t *cm) {
    esl_vec_FNorm(cm->null, cm->abc->K);
    for (int v = 0; v < cm->M; v++) {
        if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
            esl_vec_FSet(cm->e[v], cm->abc->K * cm->abc->K, 0.0f);
            esl_vec_FCopy(cm->null, cm->abc->K, cm->e[v]);
        }
    }
}

double targetRelent(int clen, int nbps) {
    const double esigma = 45.0;
    const double reTarget = (nbps > 0) ? 0.59 : 0.38;
    double etarget = (esigma - eslCONST_LOG2R *
                      std::log(2.0 / ((double) clen * (double) (clen + 1)))) / (double) clen;
    return (etarget > reTarget) ? etarget : reTarget;
}

// Round a probability the same way as ASCII CM
static inline float asciiRoundProb(float p, float nullLog, float nullMul) {
    if (p == 0.0f) return 0.0f;
    double sc = std::log((double) p / (double) nullLog) * 1.44269504;
    double r  = std::round(sc * 1000.0) / 1000.0;
    return (float) (std::exp(r / 1.44269504) * (double) nullMul);
}

// Quantize to the ASCII CM-file precision
static void quantizeCM(CM_t *cm) {
    const int K = cm->abc->K;
    std::vector<float> nullOrig(cm->null, cm->null + K);
    for (int x = 0; x < K; x++) {
        cm->null[x] = asciiRoundProb(cm->null[x], 1.0f / (float) K, 1.0f / (float) K);
    }
    for (int v = 0; v < cm->M; v++) {
        if (cm->sttype[v] != B_st)
            for (int x = 0; x < cm->cnum[v]; x++) {
                cm->t[v][x] = asciiRoundProb(cm->t[v][x], 1.0f, 1.0f);
            }
        if (cm->sttype[v] == MP_st) {
            for (int x = 0; x < K; x++) {
                for (int y = 0; y < K; y++) {
                   cm->e[v][x*K+y] = asciiRoundProb(cm->e[v][x*K+y], nullOrig[x]*nullOrig[y], cm->null[x]*cm->null[y]);
                }
            }
        } else if (cm->sttype[v]==ML_st || cm->sttype[v]==MR_st || cm->sttype[v]==IL_st || cm->sttype[v]==IR_st) {
            for (int x = 0; x < K; x++) {
                cm->e[v][x] = asciiRoundProb(cm->e[v][x], nullOrig[x], cm->null[x]);
            }
        }
    }
}

bool cmbuildFromAlignment(const std::string &name,
                                  const std::vector<std::string> &sqnames,
                                  const std::vector<std::string> &aseqs,
                                  const std::string &ssCons,
                                  std::string &cmText, std::string &err,
                                  double ereOverride, double symfracOpt, bool noss,
                                  bool forScan, bool useInside, bool useLocal,
                                  CM_t **ret_cm) {
    if (aseqs.empty()) { err = "no sequences"; return false; }
    if (sqnames.size() != aseqs.size()) { err = "name/row count mismatch"; return false; }
    const int nseq = (int) aseqs.size();
    const int alen = (int) aseqs[0].size();
    if (alen == 0) { err = "empty alignment rows"; return false; }
    for (const std::string &r : aseqs)
        if ((int) r.size() != alen) { err = "ragged alignment rows"; return false; }
    if (!noss) {
        if (ssCons.empty()) { err = "no SS_cons annotation"; return false; }
        if ((int) ssCons.size() != alen) { err = "SS_cons length mismatch"; return false; }
    }

    cmInfernalGlobalInit();

    char errbuf[eslERRBUFSIZE];
    errbuf[0] = '\0';

    ESL_ALPHABET *abc = esl_alphabet_Create(eslRNA);
    ESL_MSA *msa = esl_msa_CreateDigital(abc, nseq, alen);
    for (int i = 0; i < nseq; i++) {
        esl_abc_Digitize(abc, aseqs[i].c_str(), msa->ax[i]);
        msa->wgt[i] = 1.0;
        esl_strdup(sqnames[i].c_str(), -1, &msa->sqname[i]);
    }

    // --hand: all-'x' RF => every alignment column is a consensus/match column,
    // so clen == alen == query length and cmscan traces map 1:1 to query coords.
    msa->rf = (char *) malloc((size_t) alen + 1);
    memset(msa->rf, 'x', (size_t) alen);
    msa->rf[alen] = '\0';
    msa->ss_cons = (char *) malloc((size_t) alen + 1);
    // no base pairs -> sequence-only CM
    if (noss) {
        memset(msa->ss_cons, '.', (size_t) alen);
    } else {
        memcpy(msa->ss_cons, ssCons.c_str(), (size_t) alen);
    }
    msa->ss_cons[alen] = '\0';
    esl_msa_SetName(msa, name.empty() ? "query" : name.c_str(), -1);

    uint32_t checksum = 0;
    esl_msa_Checksum(msa, &checksum);

    // PB (Henikoff position-based) relative weights over the RF consensus cols
    ESL_MSAWEIGHT_CFG *wcfg = esl_msaweight_cfg_Create();
    wcfg->ignore_rf = FALSE;   // use the RF line we set
    esl_msaweight_PB_adv(wcfg, msa, NULL);
    esl_msaweight_cfg_Destroy(wcfg);

    // overwrite each fragment's terminal gaps with the missing symbol for parsetree truncation
    ESL_BITFIELD *fragassign = NULL;
    esl_msa_MarkFragments(msa, 0.5f, &fragassign);
    const ESL_DSQ missing = (ESL_DSQ) esl_abc_XGetMissing(abc);
    for (int i = 0; i < nseq; i++) {
        if (!esl_bitfield_IsSet(fragassign, i)) continue;
        for (int apos = 1; apos <= alen; apos++) {
            if (esl_abc_XIsResidue(abc, msa->ax[i][apos])) break;
            msa->ax[i][apos] = missing;
        }
        for (int apos = alen; apos >= 1; apos--) {
            if (esl_abc_XIsResidue(abc, msa->ax[i][apos])) break;
            msa->ax[i][apos] = missing;
        }
    }
    esl_bitfield_Destroy(fragassign);

    // guide tree + expanded CM
    CM_t *cm = NULL;
    Parsetree_t *mtr = NULL;
    const int   useRf = (symfracOpt < 0.0) ? TRUE : FALSE;  // <0: all columns via RF; else occupancy
    const float symfr = (symfracOpt < 0.0) ? 0.5f : (float) symfracOpt;
    if (HandModelmaker(msa, errbuf, useRf, FALSE, FALSE, symfr, &cm, &mtr) != eslOK) {
        err = std::string("HandModelmaker: ") + errbuf;
        esl_msa_Destroy(msa); esl_alphabet_Destroy(abc);
        return false;
    }

    float *null = NULL;
    DefaultNullModel(abc, &null);
    CMSetNullModel(cm, null);

    // rebalance
    CM_t *balanced = NULL;
    if (CMRebalance(cm, errbuf, &balanced) != eslOK) {
        err = std::string("CMRebalance: ") + errbuf;
        FreeParsetree(mtr); FreeCM(cm); free(null);
        esl_msa_Destroy(msa); esl_alphabet_Destroy(abc);
        return false;
    }
    FreeCM(cm);
    cm = balanced;

    int nbps = 0;
    for (int nd = 0; nd < cm->nodes; nd++) {
        if (cm->ndtype[nd] == MATP_nd) nbps++;
    }
    Prior_t *pri = Prior_Default(nbps == 0 ? TRUE : FALSE);

    // count transitions and joint emissions
    int *used_el = (int *) malloc(((size_t) alen + 1) * sizeof(int));
    esl_vec_ISet(used_el, alen + 1, FALSE);
    std::vector<Parsetree_t *> tr(nseq, NULL);
    for (int i = 0; i < nseq; i++) {
        if (Transmogrify(cm, errbuf, mtr, msa->ax[i], used_el, alen, &tr[i]) != eslOK) {
            err = std::string("Transmogrify: ") + errbuf;
            for (int k = 0; k < i; k++) FreeParsetree(tr[k]);
            free(used_el); Prior_Destroy(pri); FreeParsetree(mtr); FreeCM(cm); free(null);
            esl_msa_Destroy(msa); esl_alphabet_Destroy(abc);
            return false;
        }
        ParsetreeCountExceptTruncatedMPs(cm, tr[i], msa->ax[i], msa->wgt[i]);
    }

    // truncated-MP marginals against a stable double copy of the MP counts
    std::vector<double *> dbl_e(cm->M, NULL);
    const int K = cm->abc->K;
    for (int v = 0; v < cm->M; v++) {
        if (cm->sttype[v] == MP_st) {
            dbl_e[v] = (double *) malloc((size_t) K * K * sizeof(double));
            for (int a = 0; a < K * K; a++) dbl_e[v][a] = (double) cm->e[v][a];
        }
    }
    for (int i = 0; i < nseq; i++) {
        ParsetreeCountOnlyTruncatedMPs(cm, tr[i], msa->ax[i], msa->wgt[i], dbl_e.data(), pri);
    }
    for (int v = 0; v < cm->M; v++) free(dbl_e[v]);

    // flanking ROOT_IL/IR inserts are not learned from the alignment.
    cm_zero_flanking_insert_counts(cm, errbuf);

    // detach-check dual inserts
    cm_find_and_detach_dual_inserts(cm, TRUE, FALSE);

    // 1emit map, EL self-transition score, sequence counts
    cm->emap = CreateEmitMap(cm);
    cm->el_selfsc = sreLOG2(0.94);
    cm->nseq = nseq;
    cm->eff_nseq = (float) nseq;

    // null2/null3 omega
    cm->null2_omega = DEFAULT_NULL2_OMEGA;
    cm->null3_omega = DEFAULT_NULL3_OMEGA;

    // name.
    cm_SetName(cm, msa->name);

    // entropy weighting
    {
        const int clen = cm->clen;
        const double etarget = (ereOverride >= 0.0) ? ereOverride : targetRelent(clen, nbps);
        double hmm_re = 0.0, neff = 0.0;
        cm_EntropyWeight(cm, pri, etarget, 0.1, (double) cm->nseq, FALSE, &hmm_re, &neff);
        cm->eff_nseq = (float) neff;
        cm_Rescale(cm, (float) (neff / (double) nseq));
    }

    // priors, detach, flatten inserts, renormalize
    PriorifyCM(cm, pri);
    cm_find_and_detach_dual_inserts(cm, FALSE, TRUE);
    flattenInsertEmissions(cm);
    CMRenormalize(cm);
    if (forScan) {
        // round to CM-file precision + renorm, exactly as cm_file_Read
        quantizeCM(cm);
        CMRenormalize(cm);
    }

    // QDB bands + W + logodds + consensus
    cm->beta_W = 1e-7;
    cm->qdbinfo->beta1 = 1e-7;
    cm->qdbinfo->beta2 = 1e-15;
    // global mode: no CM_CONFIG_LOCAL
    cm->config_opts |= CM_CONFIG_QDB;
    if (forScan) {
        cm->config_opts |= CM_CONFIG_SCANMX | CM_CONFIG_NONBANDEDMX;
        if (useLocal)  {
            cm->config_opts |= CM_CONFIG_LOCAL | CM_CONFIG_HMMLOCAL | CM_CONFIG_HMMEL;
        }
        if (useInside) {
            cm->search_opts |= CM_SEARCH_INSIDE;
        }
    }
    if (cm_Configure(cm, errbuf, -1) != eslOK) {
        err = std::string("cm_Configure: ") + errbuf;
        for (int i = 0; i < nseq; i++) FreeParsetree(tr[i]);
        free(used_el); Prior_Destroy(pri); FreeParsetree(mtr); FreeCM(cm); free(null);
        esl_msa_Destroy(msa); esl_alphabet_Destroy(abc);
        return false;
    }

    // consensus for CM_ALIDISPLAY
    cm_SetConsensus(cm, cm->cmcons, NULL);

    if (forScan) {
        // caller owns cm; abc stays alive via cm->abc
        *ret_cm = cm;
        for (int i = 0; i < nseq; i++) {
            FreeParsetree(tr[i]);
        }
        free(used_el);
        Prior_Destroy(pri);
        FreeParsetree(mtr);
        free(null);
        esl_msa_Destroy(msa);
        return true;
    }

    // Disable unused p7 HMM filter with very few iteration
    if (cm->mlp7 != NULL) {
        // not thread safe?
        static std::mutex p7Mutex;
        std::lock_guard<std::mutex> lock(p7Mutex);
        double gfmu = 0.0, gflambda = 0.0;
        if (cm_p7_Calibrate(cm->mlp7, errbuf, 200, 200, 100, ESL_MAX(100, 2 * cm->clen),
                            10, 10, 10, 10, 0.055, 0.065, &gfmu, &gflambda) == eslOK) {
            cm_p7_hmm_SetConsensus(cm->mlp7);
            cm_SetFilterHMM(cm, cm->mlp7, gfmu, gflambda);
        }
    }

    // checksum, validate, write
    cm->checksum = checksum;
    cm->flags |= CMH_CHKSUM;
    if (cm_Validate(cm, 1e-4, errbuf) != eslOK) {
        err = std::string("cm_Validate: ") + errbuf;
        for (int i = 0; i < nseq; i++) FreeParsetree(tr[i]);
        free(used_el); Prior_Destroy(pri); FreeParsetree(mtr); FreeCM(cm); free(null);
        esl_msa_Destroy(msa); esl_alphabet_Destroy(abc);
        return false;
    }

    // serialize to an in-memory buffer
    bool ok = true;
    {
        char *buf = NULL;
        size_t sz = 0;
        FILE *fp = open_memstream(&buf, &sz);
        if (fp == NULL) {
            err = "open_memstream failed";
            ok = false;
        } else {
            if (cm_file_WriteASCII(fp, -1, cm) != eslOK) {
                err = "cm_file_WriteASCII failed";
                ok = false;
            }
            fclose(fp);
            if (ok) {
                cmText.assign(buf, sz);
            }
            free(buf);
        }
    }

    // Cleanup
    for (int i = 0; i < nseq; i++) {
        FreeParsetree(tr[i]);
    }
    free(used_el);
    Prior_Destroy(pri);
    FreeParsetree(mtr);
    FreeCM(cm);
    free(null);
    esl_msa_Destroy(msa);
    esl_alphabet_Destroy(abc);
    return ok;
}
