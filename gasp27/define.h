/*****************************************************************/
/*   Copyright (c) William J. Welch 1990--99.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*   Version: 1999.07.01                                         */
/*****************************************************************/

/* The following names, text, etc. are gathered here to */
/* facilitate consistency and possible translation into */
/* foreign languages.                                   */

/* Change ABSOLUTE to ABSOLUTE_TXT, etc. */

#define ABSOLUTE         "Absolute"
#define ANALYZE          "Analyze"
#define ARCSIN_STR       "Arcsin"
#define ARCTAN_STR       "Arctan"
#define AVERAGE          "Average"
#define CAND_GROUP       "CandidateGroup"
#define CASES            "Cases"
#define CONTINUOUS_STR   "Continuous"
#define DISCRETE_STR     "Discrete"
#define DISTRIBUTION     "Distribution"
#define FIXED_STR        "Fixed"
#define GRID_STR         "Grid"
#define GROUP            "Group"
#define ILLEGAL_COND_TXT "Illegal condition"
#define INFINITY_TXT     "Infinity"
#define INCLUSIVE        "Inclusive"
#define LOGUNIFORM_STR   "LogUniform"
#define LOG_STR          "Log"
#define MAX              "Max"
#define MIN              "Min"
#define NO_STR           "No"
#define NONE_STR         "None"
#define NORMAL_STR       "Normal"
#define NOT_AVAIL        "NA"
#define NUM_LEVELS       "NumberLevels"
#define NUM_CATS         "NumberCategories"
#define RELATIVE         "Relative"
#define SECONDS          "Seconds"
#define SUM              "Sum"
#define STEP             "Step"
#define SUPPORT          "Support"
#define TERM             "Term"
#define TRANSFORMATION   "Transformation"
#define UNIFORM_STR      "Uniform"
#define UNITS            "Units"
#define VARIABLE         "Variable"
#define WEIGHT           "Weight"
#define YES_STR          "Yes"

#define NUMERIC_ERR_TXT  "Numerical instability.\n"

/* Not-available for numeric types.      */
/* 1999.07.01: Made smaller than maxima, */
/* to distinguish NA from maximum.       */ 
#define NA_INT      (INT_MAX - 1)
#define NA_SIZE_T   (SIZE_T_MAX - 1)
#define NA_REAL     (REAL_MAX / 10)

/* Error numbers: */

#define OK               0
#define FILE_ERR         -10
#define INCOMPAT_ERR     -15
#define INPUT_ERR        -20
#define LOGIC_ERR        -25
#define MEMORY_ERR       -30
#define NONUNIQ_ERR      -35
#define NUMERIC_ERR      -40
#define RANGE_ERR        -45

#define ARG1_ERR         -101
#define ARG2_ERR         -102
#define ARG3_ERR         -103
#define ARG4_ERR         -104
#define ARG5_ERR         -105

#define ALL_DONE         -1000

#define INDEX_ERR        SIZE_T_MAX
#define INDEX_OK         (SIZE_T_MAX - 1)

/* Mnemonics for special characters: */
 
#define EOL              '\n'
#define TAB              '\t'
#define BLANK            ' '
#define FORMFEED         '\f'
 
 
/* Macros for math: */

#define HALFPI           1.570796327
#define MAX_COND_NUM     1.0e10

#ifndef is_even
#define is_even(x)  (((x) % 2 == 0) ? 1 : 0)
#endif

#ifndef is_odd
#define is_odd(x)   (((x) % 2 == 1) ? 1 : 0)
#endif

#ifndef max
#define max(a, b)   ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b)   ((a) < (b) ? (a) : (b))
#endif

#ifndef sign
#define sign(a, b)  ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

#ifndef sqr
#define sqr(x)      ((x) * (x))
#endif


/* typedefs: */

typedef unsigned long    ulong;
typedef char *           string;
typedef int              boolean;

#define YES              1
#define NO               0
#define TRUE             YES
#define FALSE            NO

/* Data types: */

#define INTEGER 0
#define REAL    1
#define SIZE_T  2
#define STRING  3

#define NUM_DATA_TYPES   4

#define MIXED   5


/* Data classes: */

#define MATRIX  0
#define SCALAR  1

#define UNKNOWN -1

