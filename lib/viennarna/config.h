/* Minimal config.h for ViennaRNA built within MMseqs2 */
#ifndef VRNA_MMSEQS_CONFIG_H
#define VRNA_MMSEQS_CONFIG_H

/* Version */
#define VRNA_VERSION "2.7.2"
#define VRNA_VERSION_MAJOR 2
#define VRNA_VERSION_MINOR 7
#define VRNA_VERSION_PATCH 2

/* We always have math.h */
#define HAVE_MATH_H 1

/* Check for erand48 (POSIX) */
#if defined(__unix__) || defined(__APPLE__)
#define HAVE_ERAND48 1
#endif

/* No SVM support */
/* #undef VRNA_WITH_SVM */

/* No GSL support */
/* #undef VRNA_WITH_GSL */

/* Naview layout support */
#define VRNA_WITH_NAVIEW_LAYOUT 1

/* JSON support */
#define VRNA_WITH_JSON_SUPPORT 1

/* No colored TTY output in library mode */
#define VRNA_WITHOUT_TTY_COLORS 1

#endif
