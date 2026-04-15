/*
 * CmCalibrateCapture.c
 *
 * Provides the capture-to-memory shim for embedded infernal cmcalibrate.
 * CMakeLists.txt renames cmcalibrate's main() to infernal_cmcalibrate_main_impl().
 * This file wraps that with infernal_cmcalibrate_main(), which after a successful
 * run reads the rewritten CM file (always the last positional argv) into a heap
 * buffer so that InfernalBridge can retrieve it without an extra disk read.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Defined by cmcalibrate.c after preprocessor rename */
int infernal_cmcalibrate_main_impl(int argc, char **argv);

static int    g_capture_enabled = 0;
static char  *g_capture_buf     = NULL;
static size_t g_capture_len     = 0;

void infernal_cmcalibrate_set_capture_to_memory(int enabled) {
    g_capture_enabled = enabled;
    if (!enabled) {
        free(g_capture_buf);
        g_capture_buf = NULL;
        g_capture_len = 0;
    }
}

const char *infernal_cmcalibrate_get_captured_cm_text(size_t *opt_len) {
    if (opt_len) { *opt_len = g_capture_len; }
    return g_capture_buf;
}

int infernal_cmcalibrate_main(int argc, char **argv) {
    int ret = infernal_cmcalibrate_main_impl(argc, argv);
    if (g_capture_enabled && ret == 0 && argc > 0) {
        /* cmcalibrate rewrites argv[argc-1] (the CM file path) in-place */
        const char *cmpath = argv[argc - 1];
        FILE *f = fopen(cmpath, "rb");
        if (f) {
            fseek(f, 0, SEEK_END);
            long sz = ftell(f);
            rewind(f);
            if (sz > 0) {
                char *buf = (char *)realloc(g_capture_buf, (size_t)(sz) + 1);
                if (buf) {
                    g_capture_buf = buf;
                    g_capture_len = (size_t)fread(g_capture_buf, 1, (size_t)sz, f);
                    g_capture_buf[g_capture_len] = '\0';
                }
            }
            fclose(f);
        }
    }
    return ret;
}
