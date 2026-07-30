#ifndef eslCONFIG_INCLUDED
#define eslCONFIG_INCLUDED

/* Version info. */
#define EASEL_VERSION "0.49"
#define EASEL_DATE "Aug 2023"
#define EASEL_COPYRIGHT "Copyright (C) 2023 Howard Hughes Medical Institute."
#define EASEL_LICENSE "Freely distributed under the BSD open source license."
#define EASEL_URL "http://bioeasel.org/"

/* Control of debugging instrumentation */
#define eslDEBUGLEVEL 0

#if defined(__aarch64__) || defined(__ARM_NEON) || defined(_M_ARM64)
#  define eslENABLE_NEON 1
#  if defined(__aarch64__) || defined(_M_ARM64)
#    define eslHAVE_NEON_AARCH64 1
#  endif
#elif defined(__SSE2__) || defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#  define eslENABLE_SSE 1
/* _MM_SET_FLUSH_ZERO_MODE, always avail w/ SSE */
#  define HAVE_FLUSH_ZERO_MODE 1
#else
#  error "unsupported target architecture for Easel (need NEON or SSE2)"
#endif

/* HAVE_PTHREAD DISABLED */
/* #undef HAVE_PTHREAD */

/* HAVE_GZIP DISABLED */
/* #undef HAVE_GZIP */

/* Headers */
#if defined(__has_include)
#  if __has_include(<endian.h>)
/* glibc; absent on macOS/BSD (they use machine/endian.h) */
#    define HAVE_ENDIAN_H 1
#  endif
#elif defined(__linux__)
#  define HAVE_ENDIAN_H 1
#endif
#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_UNISTD_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_STRINGS_H 1
#define HAVE_NETINET_IN_H 1

#define HAVE_SYS_PARAM_H 1

/* Compiler characteristics */
#define HAVE_FUNC_ATTRIBUTE_NORETURN 1
#define HAVE_FUNC_ATTRIBUTE_FORMAT 1

/* Functions (all present on Linux/macOS targets) */
#define HAVE_GETCWD 1
#define HAVE_GETPID 1
#define HAVE_POPEN 1
#define HAVE_STRCASECMP 1
#define HAVE_STRSEP 1
#define HAVE_SYSCONF 1
#define HAVE_TIMES 1
#define HAVE_FSEEKO 1

/* Function behavior */
#define eslSTOPWATCH_HIGHRES

#endif /*eslCONFIG_INCLUDED*/
