#ifndef RIBOSEEK_CMCALIBRATE_H
#define RIBOSEEK_CMCALIBRATE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* In-memory, reentrant Infernal cmcalibrate for one CM (ASCII INFERNAL1/a text).
 *
 * Runs the same exponential-tail calibration as `cmcalibrate` (serial path,
 * default options, --seed 181) and returns the input text with the four ECM
 * lines inserted before the "CM" line, i.e. what `cmcalibrate` writes minus
 * its COM line. Nothing is printed and no files are touched, so several
 * threads may calibrate different CMs at once.
 *
 *   L_mb        random sequence length in Mb (cmcalibrate -L, default 1.6)
 *   use_random  0: Infernal default background (genome-like HMM),
 *               1: iid sequence from the CM's null model (cmcalibrate --random)
 *   ret_text    malloc'ed calibrated CM text, caller frees
 *   errbuf      at least 128 bytes; set on failure
 *
 * Returns 0 on success.
 */
int riboseek_cmcalibrate(const char *cmtext, size_t len, double L_mb, int use_random,
                         char **ret_text, size_t *ret_len, char *errbuf);

#ifdef __cplusplus
}
#endif

#endif
