/*****************************************************************/
/*   Implementation-dependent definitions                        */
/*                                                               */
/*   Copyright (c) William J. Welch 1991--99.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*   Version: 2001.01.08                                         */
/*****************************************************************/

/* Uncomment the following for UNIX and/or SUN environments. */
#define UNIX_DEFINED
/*
#define SUN_DEFINED
*/

/* Last year of licence: 100 for Year 2000, etc. */
#define YEAR_LIC    116

#ifndef DBL_EPSILON
     #include <float.h>
#endif

#ifndef UINT_MAX
     #include <limits.h>
#endif

#ifdef NULL
     #undef NULL
#endif
#define NULL             0L

typedef double      real;
#define EPSILON     DBL_EPSILON    /* Machine precision. */
#define REAL_MAX    DBL_MAX
#define REAL_MIN    DBL_MIN
#define LN_MAX      (700.0)        /* Approx. ln(REAL_MAX) */
#define LN_MIN      (-700.0)       /* Approx. ln(REAL_MIN) */

#define DEF_IN_DIR   "."  /* Default directory for input files. */
#define DEF_OUT_DIR  "."  /* Default directory for output files. */

#define INPUT_COLS  256   /* Maximum input line width (in characters). */
#define OUTPUT_COLS  80   /* Line width (in characters) for output. */
#define MAXTOK      256   /* Maximum length for filenames, */
                          /* converted numbers, etc.       */

#define PRECISION     6   /* Default precision for %g, etc. */

#define LOG_FILE      "logfile.out"

/* Directory separator. */
#ifdef UNIX_DEFINED
     #define DIR_SEP      "/"
#else
     /* DOS */
     #define DIR_SEP      "\\" /* Directory separator. */
#endif

/* SIZE_T_MAX depends on the definition of size_t: */
/* UINT_MAX for NeXT gcc and DOS Microsoft C;      */
/* INT_MAX for Sun gcc.                            */
#ifdef SUN_DEFINED
     #define SIZE_T_MAX  INT_MAX
#else
     #define SIZE_T_MAX  UINT_MAX
#endif

/* Defects in Sun gcc environments. */
#ifdef SUN_DEFINED
     #define difftime(time2, time1)   ((double) (time2) - (time1))
     #define strtoul(a, b, c) (unsigned long) strtol(a, b, c)
     double strtod(const char *s, char **endp);
#endif
