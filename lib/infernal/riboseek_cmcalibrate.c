/* In-memory cmcalibrate for riboseek, see riboseek_cmcalibrate.h.
 *
 * Includes the vendored cmcalibrate.c (its main() renamed away) so the
 * calibration reuses Infernal's own static helpers unchanged: initialize_cm,
 * generate_sequences, process_search_workunit and fit_histogram. The driver
 * below is serial_master() for a single CM without the progress output,
 * the --merge/--part/--forecast branches and the file rewrite.
 */
#define main riboseek_cmcalibrate_unused_main
#include "cmcalibrate.c"
#undef main

#include <pthread.h>
#include "riboseek_cmcalibrate.h"

/* CM parsing and cm_Configure() are run under a lock; the searches and fits are not. */
static pthread_mutex_t rs_cal_mutex = PTHREAD_MUTEX_INITIALIZER;

static int
rs_read_cm(const char *text, size_t len, ESL_ALPHABET **abc, CM_t **ret_cm, char *errbuf)
{
  int      status;
  CM_FILE *cmfp = NULL;
  char    *buf  = NULL;

  *ret_cm = NULL;
  ESL_ALLOC(buf, len + 1);
  memcpy(buf, text, len);
  buf[len] = '\0';
  if ((status = cm_file_OpenBuffer(buf, (int) len, FALSE, &cmfp)) != eslOK) {
    free(buf);
    ESL_FAIL(status, errbuf, "could not open CM text");
  }
  cmfp->hfp = NULL; /* cm_file_OpenBuffer() leaves it unset; cm_file_Close() frees it */
  /* the filter HMM is not needed for calibration; the output keeps the input text */
  if ((status = cm_file_Read(cmfp, FALSE, abc, ret_cm)) != eslOK || *ret_cm == NULL) {
    snprintf(errbuf, eslERRBUFSIZE, "could not parse CM: %s", cmfp->errbuf);
    cm_file_Close(cmfp);
    free(buf);
    return (status == eslOK) ? eslEFORMAT : status;
  }
  cm_file_Close(cmfp);
  free(buf);
  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "out of memory");
}

/* serial_loop() without the progress bar */
static int
rs_search_block(WORKER_INFO *info, char *errbuf, ESL_SQ_BLOCK *sq_block)
{
  int status;
  int i, h;
  CM_TOPHITS *th = NULL;

  info->nhits = 0;
  if (info->scA != NULL) { free(info->scA); info->scA = NULL; }
  for (i = 0; i < sq_block->count; i++) {
    if ((status = process_search_workunit(info->cm, errbuf, sq_block->list[i].dsq, sq_block->list[i].L, info->cutoff, &th)) != eslOK) return status;
    if (th->N > 0) {
      ESL_REALLOC(info->scA, sizeof(float) * (info->nhits + th->N));
      for (h = 0; h < th->N; h++) info->scA[info->nhits + h] = th->unsrt[h].score;
      info->nhits += th->N;
    }
    cm_tophits_Destroy(th);
    th = NULL;
  }
  return eslOK;

 ERROR:
  if (th != NULL) cm_tophits_Destroy(th);
  ESL_FAIL(eslEMEM, errbuf, "out of memory");
}

int
riboseek_cmcalibrate(const char *cmtext, size_t len, double L_mb, int use_random,
                     char **ret_text, size_t *ret_len, char *out_errbuf)
{
  int            status   = eslOK;
  char           errbuf[eslERRBUFSIZE];
  ESL_GETOPTS   *go       = NULL;
  struct cfg_s   cfg;
  CM_t          *cm       = NULL;
  CM_t          *nc_cm    = NULL;
  ESL_SQ_BLOCK  *sq_block = NULL;
  WORKER_INFO    info;
  int            exp_mode;
  int            i;
  double         mu, lambda;
  int            nrandhits;
  float          tailp;
  ExpInfo_t     *expA[EXP_NMODES];
  char           Lbuf[32];
  char           ecm[4][160];
  int            ecmlen = 0;
  const char    *ins;
  char          *out = NULL;
  int            locked = FALSE;
  static const int   order[4] = { EXP_CM_LC, EXP_CM_GC, EXP_CM_LI, EXP_CM_GI };
  static const char *tags[4]  = { "ECMLC", "ECMGC", "ECMLI", "ECMGI" };

  errbuf[0] = '\0';
  *ret_text = NULL;
  *ret_len  = 0;
  for (i = 0; i < EXP_NMODES; i++) expA[i] = NULL;
  memset(&cfg, 0, sizeof(cfg));
  memset(&info, 0, sizeof(info));

  /* options exactly as on the command line, so all defaults are Infernal's */
  snprintf(Lbuf, sizeof(Lbuf), "%g", L_mb);
  {
    char *argv[6];
    int   argc = 0;
    argv[argc++] = "cmcalibrate";
    argv[argc++] = "-L";
    argv[argc++] = Lbuf;
    if (use_random) argv[argc++] = "--random";
    argv[argc++] = "in-memory.cm";
    argv[argc]   = NULL;
    if ((go = esl_getopts_Create(options)) == NULL)          { status = eslEMEM; snprintf(errbuf, sizeof(errbuf), "esl_getopts_Create failed"); goto CLEANUP; }
    if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
        esl_opt_VerifyConfig(go)               != eslOK)     { status = eslEINVAL; snprintf(errbuf, sizeof(errbuf), "bad cmcalibrate options: %s", go->errbuf); goto CLEANUP; }
  }

  cfg.L         = EXPTAIL_CHUNKLEN;
  cfg.N         = (int) (((esl_opt_GetReal(go, "-L") * 1000000.) / (float) cfg.L) + 0.5);
  cfg.sc_cutoff = -eslINFINITY;
  cfg.r         = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  cfg.r_est     = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  if (cfg.r == NULL || cfg.r_est == NULL) { status = eslEMEM; snprintf(errbuf, sizeof(errbuf), "failed to create RNG"); goto CLEANUP; }

  pthread_mutex_lock(&rs_cal_mutex); locked = TRUE;
  if ((status = rs_read_cm(cmtext, len, &cfg.abc, &cm, errbuf))  != eslOK) goto CLEANUP;
  if ((status = cm_Clone(cm, errbuf, &nc_cm))                     != eslOK) goto CLEANUP;
  if ((status = initialize_cm(go, &cfg, errbuf, cm, FALSE))      != eslOK) goto CLEANUP;
  pthread_mutex_unlock(&rs_cal_mutex); locked = FALSE;

  if ((status = CreateGenomicHMM(cfg.abc, errbuf, &cfg.ghmm_sA, &cfg.ghmm_tAA, &cfg.ghmm_eAA, &cfg.ghmm_nstates)) != eslOK) goto CLEANUP;
  if ((status = set_dnull(&cfg, cm, errbuf))                       != eslOK) goto CLEANUP;
  if ((status = generate_sequences(go, &cfg, errbuf, cm, &sq_block)) != eslOK) goto CLEANUP;
  for (i = 0; i < EXP_NMODES; i++) {
    if ((expA[i] = CreateExpInfo()) == NULL) { status = eslEMEM; snprintf(errbuf, sizeof(errbuf), "out of memory"); goto CLEANUP; }
  }

  info.cutoff = cfg.sc_cutoff;
  info.ipart  = -1;
  info.npart  = -1;
  for (exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
    /* glocal modes first; switch to the local configuration once, as serial_master() does */
    if (exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
      FreeCM(cm);
      cm = nc_cm;
      nc_cm = NULL;
      pthread_mutex_lock(&rs_cal_mutex); locked = TRUE;
      status = initialize_cm(go, &cfg, errbuf, cm, TRUE);
      pthread_mutex_unlock(&rs_cal_mutex); locked = FALSE;
      if (status != eslOK) goto CLEANUP;
    }
    if (ExpModeIsInside(exp_mode)) cm->search_opts |=  CM_SEARCH_INSIDE;
    else                           cm->search_opts &= ~CM_SEARCH_INSIDE;
    info.cm = cm;
    if ((status = rs_search_block(&info, errbuf, sq_block)) != eslOK) goto CLEANUP;
    if ((status = fit_histogram(go, &cfg, errbuf, info.scA, info.nhits, exp_mode, &mu, &lambda, &nrandhits, &tailp)) != eslOK) goto CLEANUP;
    SetExpInfo(expA[exp_mode], lambda, mu, (double) (cfg.L * cfg.N), nrandhits, tailp);
  }

  /* ECM lines formatted as cm_file_WriteASCII() writes them */
  for (i = 0; i < 4; i++) {
    const ExpInfo_t *e = expA[order[i]];
    snprintf(ecm[i], sizeof(ecm[i]), "%s    %.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n",
             tags[i], e->lambda, e->mu_extrap, e->mu_orig, (long) (e->dbsize + 0.5), e->nrandhits, e->tailp);
    ecmlen += (int) strlen(ecm[i]);
  }
  /* insert before the "CM" line that ends the header */
  ins = NULL;
  {
    const char *p = cmtext, *end = cmtext + len;
    while (p < end) {
      const char *nl = memchr(p, '\n', (size_t) (end - p));
      size_t      ll = (nl != NULL) ? (size_t) (nl - p) : (size_t) (end - p);
      if (ll == 2 && p[0] == 'C' && p[1] == 'M') { ins = p; break; }
      if (ll >= 3 && strncmp(p, "ECM", 3) == 0) { status = eslEINVAL; snprintf(errbuf, sizeof(errbuf), "CM is already calibrated"); goto CLEANUP; }
      if (nl == NULL) break;
      p = nl + 1;
    }
  }
  if (ins == NULL) { status = eslEFORMAT; snprintf(errbuf, sizeof(errbuf), "no CM line in CM text"); goto CLEANUP; }
  ESL_ALLOC(out, len + (size_t) ecmlen + 1);
  {
    size_t head = (size_t) (ins - cmtext), off = head;
    memcpy(out, cmtext, head);
    for (i = 0; i < 4; i++) { size_t l = strlen(ecm[i]); memcpy(out + off, ecm[i], l); off += l; }
    memcpy(out + off, ins, len - head);
    off += len - head;
    out[off] = '\0';
    *ret_text = out;
    *ret_len  = off;
    out = NULL;
  }
  status = eslOK;

 CLEANUP:
 ERROR:
  if (status == eslEMEM && errbuf[0] == '\0') snprintf(errbuf, sizeof(errbuf), "out of memory");
  if (locked) pthread_mutex_unlock(&rs_cal_mutex);
  if (out      != NULL) free(out);
  if (info.scA != NULL) free(info.scA);
  if (sq_block != NULL) esl_sq_DestroyBlock(sq_block);
  for (i = 0; i < EXP_NMODES; i++) if (expA[i] != NULL) free(expA[i]);
  if (cm       != NULL) FreeCM(cm);
  if (nc_cm    != NULL) FreeCM(nc_cm);
  if (cfg.dnull   != NULL) free(cfg.dnull);
  if (cfg.ghmm_sA != NULL) free(cfg.ghmm_sA);
  if (cfg.ghmm_tAA != NULL) { for (i = 0; i < cfg.ghmm_nstates; i++) if (cfg.ghmm_tAA[i] != NULL) free(cfg.ghmm_tAA[i]); free(cfg.ghmm_tAA); }
  if (cfg.ghmm_eAA != NULL) { for (i = 0; i < cfg.ghmm_nstates; i++) if (cfg.ghmm_eAA[i] != NULL) free(cfg.ghmm_eAA[i]); free(cfg.ghmm_eAA); }
  if (cfg.r     != NULL) esl_randomness_Destroy(cfg.r);
  if (cfg.r_est != NULL) esl_randomness_Destroy(cfg.r_est);
  if (cfg.abc   != NULL) esl_alphabet_Destroy(cfg.abc);
  if (go        != NULL) esl_getopts_Destroy(go);
  if (status != eslOK && out_errbuf != NULL) snprintf(out_errbuf, 128, "%s", errbuf[0] ? errbuf : "cmcalibrate failed");
  return (status == eslOK) ? 0 : 1;
}
