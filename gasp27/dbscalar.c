/*****************************************************************/
/*   ROUTINES TO MANAGE SCALAR OBJECTS IN THE SYSTEM DATABASE    */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--2011.                  */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"
#include "kriging.h"
#include "optdes.h"
#include "alex.h"

extern Matrix  XDescrip;
extern string  Prompt;

static List    *ScalTemplates = NULL;

/* The scalars are actually stored as vectors of length 1. */
static vector       *Scal    = NULL;
static size_t       nScalars = 0;


/* Table of int scalars with defaults */
/* (illegal value = no default).      */

int       Seed = 100;

static struct IntStruct
{
     string    Name;
     int       Min;
     int       Max;
     int      *i;
}
IntScalar[] =
{
     {RAN_NUM_SEED,      1,   30000,    &Seed}
};

#define NUM_INTS    (sizeof(IntScalar) / sizeof(struct IntStruct))


/* Table of real scalars with defaults */
/* (illegal value = no default).       */
real AlphaMax            =  1.0;   /* Min p is 1. */
real AlphaMin            =  0.0;
real CensLimit           = -1.0;
real CoverDist           = -1.0;
real CritLogLikeDiff     =  1.0;
real InterPerc           =  5.0;
real Lambda              =  1.0;
real LogLikeTol          =  0.00001;
real MainPerc            =  2.5;
real Metric              =  1.0;
     /* Dispersion parameter for the GCD criterion. */
real Phi                 = -1.0;
     /* Tolerances consistent with optdes.h. */
     /* Should change design/mle algorithms. */
real SPVarPropMin        = 0.0;
real SPVarPropMax        = 1.0;
real TolAbs              =  EPSILON;
real TolRel              = 1.0e-4;
/* 1996.04.09: Changed from 10.0   */
/* 1996.04.10: Changed from 1000.0 */
/* 1999.06.23: Changed from 100.0  */
/* 1999.06.24: Changed from 1000.0 */
real ThetaStandMin       = 0.0;
real ThetaStandMax       = REAL_MAX;
real yCritical           = 0.0;

static struct RealStruct
{
     string    Name;
     real      Min;
     real      Max;
     real      *r;
}
RealScalar[] =
{
     {ALPHA "." MAX,      0.0, 1.999999, &AlphaMax       },
     {ALPHA "." MIN,      0.0, 1.999999, &AlphaMin       },
     {CENSORING_LIMIT,    0.0, REAL_MAX, &CensLimit      },
     {COVER_DIST,         0.0, REAL_MAX, &CoverDist      },
     {CRIT_LOG_LIKE_DIFF, 0.0, REAL_MAX, &CritLogLikeDiff},
     {INTER_EFF_PERC,     0.0,    100.0, &InterPerc      },
     {LAMBDA,             1.0, REAL_MAX, &Lambda         },
     {LOG_LIKE_TOL,       0.0, REAL_MAX, &LogLikeTol     },
     {MAIN_EFF_PERC,      0.0,    100.0, &MainPerc       },
     {METRIC,             1.0, REAL_MAX, &Metric         },
     {"Phi",              0.0, REAL_MAX, &Phi            },
     {SP_VAR_PROP "." MAX,0.0,      1.0, &SPVarPropMax   },
     {SP_VAR_PROP "." MIN,0.0,      1.0, &SPVarPropMin   },
     {THETA "." STANDARDIZED "." MAX,
                          0.0, REAL_MAX, &ThetaStandMax  },
     {THETA "." STANDARDIZED "." MIN,
                          0.0, REAL_MAX, &ThetaStandMin  },
     {TOL "." ABSOLUTE,   0.0, REAL_MAX, &TolAbs         },
     {TOL "." RELATIVE,   0.0, REAL_MAX, &TolRel         },
     {"yCritical",       -REAL_MAX, REAL_MAX, &yCritical }
};

#define NUM_REALS   (sizeof(RealScalar) / sizeof(struct RealStruct))


/* Table of size_t scalars with defaults */
/* (illegal value = no default).         */

size_t    derivMin       = 0;      /* Matern correlation derivatives */
size_t    derivMax       = 3;      /* Codes infinity! */
size_t    k              = 0;      /* Replace! */
size_t    kf             = 0;      /* Replace! */
size_t    ProjDimMax     = 2;
size_t    ProjDimMin     = 2;
size_t    nProtected     = 0;
size_t    n              = 0;
size_t    nRefit         = 1;
size_t    s              = 0;      /* Replace! */
size_t    Tries          = 1;
size_t    nXVars         = 0;

static struct Size_tStruct
{
     string    Name;
     size_t    Min;
     size_t    Max;
     size_t    *z;
}
Size_tScalar[] =
{
     {"Derivatives.Min", 0,                 3,    &derivMin      },
     {"Derivatives.Max", 0,                 3,    &derivMax      },
     {"k",               1,        SIZE_T_MAX,    &k             },
     {"kf",              1,        SIZE_T_MAX,    &kf            },
     {PROJ_DIM "." MAX,  1,        SIZE_T_MAX,    &ProjDimMax    },
     {PROJ_DIM "." MIN,  1,        SIZE_T_MAX,    &ProjDimMin    },
     {PROTECTED_RUNS,    0,        SIZE_T_MAX,    &nProtected    },
     {REFIT_RUNS,        1,        SIZE_T_MAX,    &nRefit        },
     {RUNS,              1,        SIZE_T_MAX,    &n             },
     {"s",               1,        SIZE_T_MAX,    &s             },
     {TRIES,             1,        SIZE_T_MAX,    &Tries         },
     {N_X_VARS,          1,        SIZE_T_MAX,    &nXVars        }
};

#define NUM_SIZE_TS (sizeof(Size_tScalar) / sizeof(struct Size_tStruct))


/* Table of string (sometimes converted to case index) scalars.  */
/* Default of zero means that *first* name is default.           */
/* INDEX_ERR = no default.                                       */
/* Why are some of these Num and some Size_t? */

size_t CorFamNum              = 0;
size_t CritNum                = INDEX_ERR;
size_t DesAlgNum              = INDEX_ERR;
size_t GenPredCoefsSize_t     = 0;
size_t InDirSize_t            = INDEX_ERR;
size_t LifeDist               = INDEX_ERR;
size_t LikeNum                = 0;
size_t LinkNum                = 0;
size_t ModCompCritNum         = 0;
size_t NormalizedRangesSize_t = INDEX_ERR;
size_t RanErrSize_t           = INDEX_ERR;
size_t RespFuncSize_t         = INDEX_ERR;
size_t SeqCritNum             = INDEX_ERR;
size_t OutDirSize_t           = INDEX_ERR;
size_t VarFnNum               = 0;

boolean RanErr           = NO;
boolean GenPredCoefs     = NO;
boolean NormalizedRanges = NO;

string  RespFunc         = NULL;
string  InDir            = DEF_IN_DIR;
string  OutDir           = DEF_OUT_DIR;
string  SeqCrit          = NULL;

static string DesAlgName[]         = DES_ALG_NAMES;
static string CorFamName[]         = COR_FAM_NAMES;
static string LifeDistName[]       = {"Exponential", "Weibull"};
static string LikeName[]           = LIKE_NAMES;
static string LinkName[]           = LINK_FN_NAMES;
static string ModCompCritName[]    = MOD_COMP_CRIT_NAMES;
static string NoYes[]              = {NO_STR, YES_STR};
static string SeqCritName[]        = {"Minimize", "Discriminate"};
static string VarFnName[]          = VAR_FN_NAMES;

static struct StrStruct
{
     string    Name;
     size_t    nLegalStrs;
     string    *LegalStr;
     size_t    *StrNum;
}
StrScalar[] =
{
     {COR_FAM,           NumStr(CorFamName),      CorFamName,
                                                  &CorFamNum          },
     {DESIGN_CRIT,       0,                       NULL,
                                                  &CritNum            },
     {DESIGN_ALG,        NumStr(DesAlgName),      DesAlgName,
                                                  &DesAlgNum          },
     {GEN_PRED_COEF,     2,                       NoYes,
                                                  &GenPredCoefsSize_t },
     {IN_DIR,            0,                       NULL,
                                                  &InDirSize_t        },
     {LIFE_DIST,         NumStr(LifeDistName),    LifeDistName,
                                                  &LifeDist           },
     {"Likelihood",      NumStr(LikeName),        LikeName,
                                                  &LikeNum            },
     {"LinkFunction",    NumStr(LinkName),        LinkName,
                                                  &LinkNum            },
     {MOD_COMP_CRIT,     NumStr(ModCompCritName), ModCompCritName,
                                                  &ModCompCritNum     },
     {NORMALIZED_RANGES, 2,                       NoYes,
                                                  &NormalizedRangesSize_t},
     {OUT_DIR,           0,                       NULL,
                                                  &OutDirSize_t       },
     {RESP_FUNC,         0,                       NULL,
                                                  &RespFuncSize_t     },
     {SEQ_CRIT,          NumStr(SeqCritName),     SeqCritName,
                                                  &SeqCritNum         },
     {RAN_ERR,           2,                       NoYes,
                                                  &RanErrSize_t       },
     {"VarianceFunction",NumStr(VarFnName),       VarFnName,
                                                  &VarFnNum   }
};

#define NUM_STRS    (sizeof(StrScalar) / sizeof(struct StrStruct))


/******************************+++********************************/
void DbScalarInit(void)
/*****************************************************************/
/*   Purpose:  Initialize the scalar database.                   */
/*                                                               */
/*   Version:  1995 May 11                                       */
/*****************************************************************/
{
     size_t    i, nLegalStrs;
     string    Name;
     string    *LegalStr;
     vector    *v;

     /* Allocate int scalars and their templates. */
     for (i = 0; i < NUM_INTS; i++)
     {
          Scal = (vector *) AllocGeneric(nScalars + 1,
                    sizeof(vector), Scal);
          v = &Scal[nScalars++];
          VecPutName(v, IntScalar[i].Name);
          VecPutType(v, INTEGER);
          VecPutLength(v, 1);
          VecPutInts(v, IntScalar[i].i);

          ScalTemplates = TemplIntAlloc(IntScalar[i].Name,
                    IntScalar[i].Min, IntScalar[i].Max,
                    ScalTemplates);
     }

     /* Allocate real scalars and their templates. */
     for (i = 0; i < NUM_REALS; i++)
     {
          Scal = (vector *) AllocGeneric(nScalars + 1,
                    sizeof(vector), Scal);
          v = &Scal[nScalars++];
          VecPutName(v, RealScalar[i].Name);
          VecPutType(v, REAL);
          VecPutLength(v, 1);
          VecPutReals(v, RealScalar[i].r);

          ScalTemplates = TemplRealAlloc(RealScalar[i].Name,
                    RealScalar[i].Min, RealScalar[i].Max,
                    ScalTemplates);
     }

     /* Allocate size_t scalars and their templates. */
     for (i = 0; i < NUM_SIZE_TS; i++)
     {
          Scal = (vector *) AllocGeneric(nScalars + 1,
                    sizeof(vector), Scal);
          v = &Scal[nScalars++];
          VecPutName(v, Size_tScalar[i].Name);
          VecPutType(v, SIZE_T);
          VecPutLength(v, 1);
          VecPutSize_ts(v, Size_tScalar[i].z);

          ScalTemplates = TemplSize_tAlloc(Size_tScalar[i].Name,
                    Size_tScalar[i].Min, Size_tScalar[i].Max,
                    ScalTemplates);
     }

     /* Allocate string scalars and their templates. */
     for (i = 0; i < NUM_STRS; i++)
     {
          if (stricmp(StrScalar[i].Name, DESIGN_CRIT) == 0)
          {
               nLegalStrs = 0;
               LegalStr  = NULL;
               while ( (Name = CritName(nLegalStrs)) != NULL)
               {
                    LegalStr = AllocStr(nLegalStrs + 1, LegalStr);
                    LegalStr[nLegalStrs++] = Name;
               }
               StrScalar[i].nLegalStrs = nLegalStrs;
               StrScalar[i].LegalStr  = LegalStr;
          }

          Scal = (vector *) AllocGeneric(nScalars + 1,
                    sizeof(vector), Scal);
          v = &Scal[nScalars++];
          VecPutName(v, StrScalar[i].Name);
          VecPutType(v, STRING);
          VecPutLength(v, 1);

          /* Set default string. */
          VecPutStrs(v, AllocStr(1, NULL));
          if (stricmp(StrScalar[i].Name, IN_DIR) == 0)
               VecPutStr(v, 0, DEF_IN_DIR);
          else if (stricmp(StrScalar[i].Name, OUT_DIR) == 0)
               VecPutStr(v, 0, DEF_OUT_DIR);
          else if (*(StrScalar[i].StrNum) != INDEX_ERR)
               /* Use LegalStr[0] as the default (may be NULL). */
               VecPutStr(v, 0, StrScalar[i].LegalStr[0]);
          else
               /* No default. */
               VecPutStr(v, 0, NOT_AVAIL);

          /* Naughty: use size_t to store StrNum (pointer). */
          VecPutSize_ts(v, StrScalar[i].StrNum);

          ScalTemplates = TemplStrAlloc(StrScalar[i].Name,
                    StrScalar[i].LegalStr, StrScalar[i].nLegalStrs,
                    ScalTemplates);
     }

     return;
}

/******************************+++********************************/
size_t DbScalIndex(const string Name, boolean ErrorMessage)
/*****************************************************************/
/*   Purpose:  Return the Scal index for Name.                   */
/*                                                               */
/*   Returns:  INDEX_ERR if not found;                           */
/*             OK        otherwise.                              */
/*                                                               */
/*   Version:  1994 September 6                                  */
/*****************************************************************/
{
     size_t    j;

     for (j = 0; j < nScalars; j++)
          if (stricmp(Name, VecName(&Scal[j])) == 0)
               return j;

     /* Not found. */
     if (ErrorMessage)
          Error(DB_SCALAR, Name);

     return INDEX_ERR;
}

/******************************+++********************************/
int DbScalCheck(size_t Index)
/*****************************************************************/
/*   Purpose:  Check whether scalar Index in Scal is legal and   */
/*             compatible.                                       */
/*                                                               */
/*   Returns:  INPUT_ERR, INCOMPAT_ERR, or OK.                   */
/*                                                               */
/*   96.01.17: CodeBug replaced by CodeCheck.                    */
/*   96.04.05: Unassigned string gives INPUT_ERR.                */
/*                                                               */
/*   Version:  1996.04.05                                        */
/*****************************************************************/
{
     int       ErrNum;
     string    Name;
     template  *T;
     vector    *ScalIndex;

     ScalIndex = &Scal[Index];

     Name = VecName(ScalIndex);

     T = TemplPtr(Name, ScalTemplates);
     CodeCheck(T != NULL);

     if (VecTempl(ScalIndex, T) != INDEX_OK)
          ErrNum = INPUT_ERR;
     else if (VecType(ScalIndex) == STRING &&
               stricmp(VecStr(ScalIndex, 0), NOT_AVAIL) == 0)
          ErrNum = INPUT_ERR;
     else
          ErrNum = DbScalCompat(Index);

     return ErrNum;
}

/******************************+++********************************/
string DbScalValue(size_t Index)
/*****************************************************************/
/*   Purpose:  Return the value of the scalar Index in Scal as a */
/*             string.                                           */
/*                                                               */
/*   Version:  1995 March 2                                      */
/*****************************************************************/
{
     string Value;

     VecToStr(&Scal[Index], &Value);

     return Value;
}

/******************************+++********************************/
int DbScalParse(size_t Index, const string Token)
/*****************************************************************/
/*   Purpose:  Parse Token and, if legal, put it in the scalar   */
/*             database.                                         */
/*                                                               */
/*   Returns:  INPUT_ERR    if Token is illegal;                 */
/*             OK           otherwise.                           */
/*                                                               */
/*   96.01.17: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996.04.15                                        */
/*****************************************************************/
{
     int       ErrNum;
     template  *T;
     vector    Old;
     vector    *ScalIndex;

     ScalIndex = &Scal[Index];

     T = TemplPtr(VecName(ScalIndex), ScalTemplates);
     CodeCheck(T != NULL);

     Old = *ScalIndex;

     if (VecFromStr(&Token, ScalIndex) != INDEX_OK ||
               VecTempl(ScalIndex, T) != INDEX_OK)
     {
          ErrNum = INPUT_ERR;
          TemplError(T);
          *ScalIndex = Old;
     }
     else
     {
          ErrNum = OK;

          if (stricmp(VecName(ScalIndex), RAN_NUM_SEED) == 0)
               RandInit(Seed, Seed, Seed);

          if (VecType(ScalIndex) == STRING)
          {
               VecPutSize_t(ScalIndex, 0,
                         TemplStrIndex(VecStr(ScalIndex, 0), T));

               if (stricmp(VecName(ScalIndex), RAN_ERR) == 0)
                     RanErr = (boolean) VecSize_t(ScalIndex, 0);
               else if (stricmp(VecName(ScalIndex), GEN_PRED_COEF)
                         == 0)
                    GenPredCoefs = (boolean) VecSize_t(ScalIndex, 0);
               else if (stricmp(VecName(ScalIndex), NORMALIZED_RANGES)
                         == 0)
                    NormalizedRanges = (boolean) VecSize_t(ScalIndex, 0);
               else if (stricmp(VecName(ScalIndex), RESP_FUNC) == 0)
                    RespFunc = VecStr(ScalIndex, 0);
               else if (stricmp(VecName(ScalIndex), IN_DIR) == 0)
                    InDir = VecStr(ScalIndex, 0);
               else if (stricmp(VecName(ScalIndex), OUT_DIR) == 0)
                    OutDir = VecStr(ScalIndex, 0);
               else if (stricmp(VecName(ScalIndex), SEQ_CRIT) == 0)
                    SeqCrit = VecStr(ScalIndex, 0);
          }
     }

     return ErrNum;
}

/******************************+++********************************/
int DbScalCompat(size_t Index)
/*****************************************************************/
/*   Purpose:  Check Scal[Index] for compatibility with other    */
/*             database objects.                                 */
/*                                                               */
/*   Returns:  INCOMPAT_ERR or OK.                               */
/*                                                               */
/*   Version:  1995 March 10                                     */
/*****************************************************************/
{
     int       ErrNum;
     size_t    ExtOffset, IndexMin, IndexMax, z;
     string    Name, NameCopy;
     vector    *ScalMax, *ScalMin;

     ErrNum = OK;

     Name = VecName(&Scal[Index]);

     /* Look for ".Min" or ".Max". */
     ExtOffset = strlen(Name) - 4;
     if (stricmp(Name + ExtOffset, "." MIN) == 0 ||
               stricmp(Name + ExtOffset, "." MAX) == 0)
     {
          NameCopy = StrDup(Name);
          strcpy(NameCopy + ExtOffset, "." MIN);
          IndexMin = DbScalIndex(NameCopy, YES);
          strcpy(NameCopy + ExtOffset, "." MAX);
          IndexMax = DbScalIndex(NameCopy, YES);
          AllocFree(NameCopy);
          ScalMin = &Scal[IndexMin];
          ScalMax = &Scal[IndexMax];

          if ((VecType(ScalMin) == REAL &&
                    VecReal(ScalMin, 0) > VecReal(ScalMax, 0)) ||
                    (VecType(ScalMin) == SIZE_T &&
                    VecSize_t(ScalMin, 0) > VecSize_t(ScalMax, 0)))
          {
               Incompatibility(DB_MIN_MAX, VecName(ScalMin),
                         VecName(ScalMax));
               ErrNum = INCOMPAT_ERR;
          }

     }

     if (stricmp(Name, PROTECTED_RUNS) == 0 && n > 0 &&
               nProtected >= n)
     {
          Incompatibility(DB_PROT_RUNS);
          ErrNum = INCOMPAT_ERR;
     }
     else if (stricmp(Name, PROJ_DIM "." MIN) == 0 ||
               stricmp(Name, PROJ_DIM "." MAX) == 0)
     {
          z = VecSize_t(&Scal[Index], 0);
          if (!MatEmpty(&XDescrip) && z > MatNumRows(&XDescrip))
          {
               Incompatibility(DB_PROJ_XVARS, Name);
               ErrNum = INCOMPAT_ERR;
          }
     }

     return ErrNum;
}
