#ifndef P7_CONFIGH_INCLUDED
#define P7_CONFIGH_INCLUDED

/* 1. Compile-time constants controlling computational behavior / output.
 * (p7_RAMLIMIT and p7_NCPU are omitted: the first is unreferenced by any source
 * we compile, the second is only used by HMMER's own program mains, which we
 * don't build.) */
#ifndef p7_ALILENGTH
#define p7_ALILENGTH       50
#endif

/* 2. Empirically tuned default parameters. */
#define p7_ETARGET_AMINO  0.59 /* bits,  from the work of Steve Johnson. */
#define p7_ETARGET_DNA    0.62 /* bits,  from the work of Travis Wheeler and Robert Hubley. */
#define p7_ETARGET_OTHER  1.0  /* bits */

#define p7_SEQDBENV          "BLASTDB"
#define p7_HMMDBENV          "PFAMDB"

/* 3. Fundamental parameters. */
#define p7_MAXABET    20      /* maximum size of alphabet (4 or 20)              */
#define p7_MAXCODE    29      /* maximum degenerate alphabet size (18 or 29)     */
#define p7_MAX_SC_TXTLEN   11
#define p7_MAXDCHLET  20      /* maximum # Dirichlet components in mixture prior */

/* 4. Version info. */
#define HMMER_VERSION "3.4"
#define HMMER_DATE "Aug 2023"
#define HMMER_COPYRIGHT "Copyright (C) 2023 Howard Hughes Medical Institute."
#define HMMER_LICENSE "Freely distributed under the BSD open source license."
#define HMMER_URL "http://hmmer.org/"

/* Choice of optimized implementation. Must match Easel's esl_config.h and the
 * impl_* sources CMake compiles. Selected from the compiler's target macros. */
#if defined(__aarch64__) || defined(__ARM_NEON) || defined(_M_ARM64)
#  define eslENABLE_NEON 1
#elif defined(__SSE2__) || defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#  define eslENABLE_SSE 1
#  define HAVE_FLUSH_ZERO_MODE 1
#else
#  error "riboseek: unsupported target architecture for vendored HMMER (need NEON or SSE2)"
#endif

/* System headers (present on our Linux/macOS targets) */
#define HAVE_NETINET_IN_H 1
#define HAVE_SYS_PARAM_H 1

/* HMMER_THREADS DISABLED */
/* #undef HMMER_THREADS */

#endif /*P7_CONFIGH_INCLUDED*/
