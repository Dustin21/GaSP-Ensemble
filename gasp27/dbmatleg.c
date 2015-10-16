/*****************************************************************/
/*   ROUTINES TO CHECK THE LEGALITY OF A DATABASE MATRIX         */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--2009.                  */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"
#include "alex.h"

extern boolean      DesignJob;

extern LinModel     RegMod;
extern LinModel     SPMod;

List                *ColTemplates = NULL;

/* Table for real-column templates. */
static struct RealColStruct
{
     string    Name;
     real      Min;
     real      Max;
}
RealCol[] =
{
     {ALPHA,             0.0,           1.99      },
     {AVERAGE,           -REAL_MAX,     REAL_MAX  },
     {ANOVA_TOTAL_PERC,  0.0,           REAL_MAX  },
     {BETA,              -REAL_MAX,     REAL_MAX  },
     {COEF,              -REAL_MAX,     REAL_MAX  },
     {COND_NUM,          -REAL_MAX,     REAL_MAX  },
     {CV_MAX_ERR,        -REAL_MAX,     REAL_MAX  },
     {CV_ROOT_MSE,       0.0,           REAL_MAX  },
     /* DbMatColChk checks integer values, too. */
     {"Derivatives",     0.0,           3.0       }, 
     {ERR_VAR,           0.0,           REAL_MAX  },
     {LOG_LIKE,          -REAL_MAX,     REAL_MAX  },
     {MAX,               -REAL_MAX,     REAL_MAX  },
     {MAX_ERR,           -REAL_MAX,     REAL_MAX  },
     {MIN,               -REAL_MAX,     REAL_MAX  },
     {PRED,              -REAL_MAX,     REAL_MAX  },
     {ROOT_MSE,          0.0,           REAL_MAX  },
     {STD_ERR,           0.0,           REAL_MAX  },
     {SP_VAR,            0.0,           REAL_MAX  },
     {THETA,             0.0,           REAL_MAX  },
     {WEIGHT,            0.0,           REAL_MAX  },
     {"x_i",             -REAL_MAX,     REAL_MAX  },
     {"x_j",             -REAL_MAX,     REAL_MAX  },
     {"y",               -REAL_MAX,     REAL_MAX  }
};

#define NUM_REAL_COLS    (sizeof(RealCol) / sizeof(struct RealColStruct))

/* Table for size_t-column templates. */
static struct Size_tColStruct
{
     string    Name;
     size_t    Min;
     size_t    Max;
}
Size_tCol[] =
{
     {CAND_GROUP,        0,   SIZE_T_MAX     },
     {CASES,             0,   SIZE_T_MAX     },
     {NUM_CATS,          0,   SIZE_T_MAX     },
     {NUM_LEVELS,        0,   SIZE_T_MAX     },
};

#define NUM_SIZE_T_COLS  (sizeof(Size_tCol)       \
                         / sizeof(struct Size_tColStruct))

/* NUM_LEVELS must be at least 2 for a GRID variable. */
/* This is checked in RegExtract.                     */

static string DistribName[]   = DISTRIB_NAMES;
static string NoYes[]         = {NO_STR, YES_STR};
static string SuppName[]      = SUPPORT_NAMES;
static string TransName[]     = TRANSFORM_NAMES;

/* Table for string-column templates. */
static struct StrColStruct
{
     string    Name;
     size_t    nLegalStrs;
     string    *LegalStr;
}
StrCol[] =
{
     {ANALYZE,           2,                   NoYes           },
     {CASE_MAX_ERR,      0,                   NULL            },
     {CASE_CV_MAX_ERR,   0,                   NULL            },
     {DISTRIBUTION,      NumStr(DistribName), DistribName     },
     {INCLUSIVE,         2,                   NoYes           },
     {SUPPORT,           NumStr(SuppName),    SuppName        },
     {TERM,              0,                   NULL            },
     {TRANSFORMATION,    NumStr(TransName),   TransName       },
     {UNITS,             0,                   NULL            },
     {VARIABLE,          0,                   NULL            }
};

#define NUM_STR_COLS  (sizeof(StrCol) / sizeof(struct StrColStruct))

/*******************************+++*******************************/
void DbColTemplInit(void)
/*****************************************************************/
/*   Purpose:  Initialize column templates.                      */
/*                                                               */
/*   Version:  1994 December 13                                  */
/*****************************************************************/
{
     size_t    i;

     /* Allocate templates for real columns. */
     for (i = 0; i < NUM_REAL_COLS; i++)
          ColTemplates = TemplRealAlloc(RealCol[i].Name,
                    RealCol[i].Min, RealCol[i].Max,
                    ColTemplates);

     /* Allocate templates for size_t columns. */
     for (i = 0; i < NUM_SIZE_T_COLS; i++)
          ColTemplates = TemplSize_tAlloc(Size_tCol[i].Name,
                    Size_tCol[i].Min, Size_tCol[i].Max,
                    ColTemplates);

     /* Allocate templates for string columns. */
     for (i = 0; i < NUM_STR_COLS; i++)
          ColTemplates = TemplStrAlloc(StrCol[i].Name,
                    StrCol[i].LegalStr, StrCol[i].nLegalStrs,
                    ColTemplates);

     return;
}

/*******************************+++*******************************/
int DbMatLegal(DbMatrix *D)
/*****************************************************************/
/*   Purpose:  Check that D is legal.                            */
/*                                                               */
/*   Returns:  INPUT_ERR    if the matrix is illegal;            */
/*             OK           otherwise.                           */
/*                                                               */
/*   Comment:  May change the matrix, e.g., by setting up a      */
/*             default.                                          */
/*                                                               */
/*   Version:  1995 October 12                                   */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    *M;
     string    Name;
     string    *Term;

     M    = D->M;
     Name = D->Name;

     /* Check column names and legality of contents. */
     if ( (ErrNum = DbMatColChk(D)) != OK)
          return ErrNum;

     /* Further specific checks: */

     else if (stricmp(Name, REG_MOD) == 0)
     {
          Term = (MatNumRows(M) > 0) ?
                    MatStrColFind(M, TERM, YES) : NULL;
          ErrNum = ModParse1(MatNumRows(M), Term, REG_MOD, &RegMod);
     }

     else if (stricmp(Name, SP_MOD) == 0)
     {
          Term = (MatNumRows(M) > 0) ?
                    MatStrColFind(M, TERM, YES) : NULL;
          ErrNum = ModParse1(MatNumRows(M), Term, SP_MOD, &SPMod);
     }

     /*
     else if (stricmp(Name, T_MAT) == 0 &&
               MatNumRows(M) > MatNumCols(M))
          ErrNum = INPUT_ERR;
     */

     else if (stricmp(Name, X_COR) == 0)
          ErrNum = DbMatCor(M);

     else if (!DesignJob &&
               (stricmp(Name, X_MAT) == 0 ||
               stricmp(Name, Y_MAT) == 0) &&
               (MatNumRows(M) == 0 || MatNumCols(M) == 0) )
          ErrNum = INPUT_ERR;

     return ErrNum;
}

/*******************************+++*******************************/
int DbMatColChk(DbMatrix *D)
/*****************************************************************/
/* Purpose:  Check each column independently.                    */
/*                                                               */
/* Returns:  OK or INPUT_ERR.                                    */
/*                                                               */
/* 1996.01.17: CodeBug replaced by CodeCheck.                    */
/* 1996.01.22: IllegalType replaced by CodeBug.                  */
/* 1996.03.28: Column names checked against optional columns     */
/*             both with and without extent (previously only     */
/*             without).                                         */
/* 2009.05.08: Derivatives (Matern) checked                      */                        
/*****************************************************************/
{
     boolean   InvalidName;
     char      *DotPtr;
     int       ErrNum;
     Matrix    *M;
     real      *col;
     size_t    BadIndex, i, ii, j, nCompCols, nOptCols, nRows;
     string    Name, NameNoExt;
     string    *CompCol, *OptCol, *s, *t;
     template  *T;

     ErrNum = OK;

     M = D->M;
     nRows = MatNumRows(M);

     /* Check for duplicate column names. */
     s = MatColNames(M);
     for (j = 1; j < MatNumCols(M); j++)
     {
          if (StrIndex(s[j], s, j) != INDEX_ERR)
          {
               Error(DB_DUP_COL, s[j]);
               ErrNum = INPUT_ERR;
          }
     }

     /* Check compulsory columns are present. */
     nCompCols = D->nCompCols;
     CompCol   = D->CompCol;
     for (j = 0; j < nCompCols; j++)
          if (StrIndex(CompCol[j], MatColNames(M), MatNumCols(M))
                    == INDEX_ERR)
          {
               Error(DB_COMPULSORY, D->Name, CompCol[j]);
               ErrNum = INPUT_ERR;
          }

     nOptCols = D->nOptCols;
     OptCol   = D->OptCol;
     if (nCompCols > 0 || nOptCols > 0)
     {
          /* Check all columns are compulsory or optional. */
          for (j = 0; j < MatNumCols(M); j++)
          {
               Name = MatColName(M, j);

               /* If column name includes extension, remove it. */
               NameNoExt = StrDup(Name);
               if ( (DotPtr = strrchr(NameNoExt, '.')) != NULL)
                    *DotPtr = NULL;

               /* Check for column that is neither compulsory */
               /* nor optional.                               */
               InvalidName = (StrIndex(Name, CompCol, nCompCols)
                         == INDEX_ERR &&
                         StrIndex(Name, OptCol, nOptCols)
                         == INDEX_ERR &&
                         StrIndex(NameNoExt, OptCol, nOptCols)
                         == INDEX_ERR);

               if (InvalidName)
               {
                    Error(DB_COL, Name);
                    ErrNum = INPUT_ERR;
               }
               else
               {
                    T = TemplPtr(NameNoExt, ColTemplates);
                    CodeCheck(T != NULL);
               }

               AllocFree(NameNoExt);

               if (InvalidName)
                    continue;

               /* Convert column to appropriate type. */
               if (MatColType(M, j) == STRING &&
                         TemplType(T) != STRING)
                    BadIndex = MatColConvert(j, TemplType(T), M);
               else
                    BadIndex = INDEX_OK;

               /* Check column contents. */
               if (BadIndex == INDEX_OK)
                    switch (MatColType(M, j))
                    {
                         case REAL:
                              BadIndex = TemplRealCheck(T,
                                        MatCol(M, j), nRows);
                              break;

                         case SIZE_T:
                              BadIndex = TemplSize_tCheck(T,
                                        MatSize_tCol(M, j), nRows);
                              break;

                         case STRING:
                              BadIndex = TemplStrCheck(T,
                                        MatStrCol(M, j), nRows);
                              break;

                         default:
                              CodeBug("Illegal type");
                    }

               if (BadIndex != INDEX_OK)
               {
                    TemplError(T);
                    Output(DB_MAT_ERR_AT, Name,
                              MatRowName(M, BadIndex));
                    ErrNum = INPUT_ERR;
               }

               if (stricmp(Name, "Derivatives") == 0)
               {
                    /* Must be integer valued. */
                    col = MatCol(M, j);
                    for (i = 1; i < nRows && ErrNum == OK; i++)
                         if (col[i]  != 0.0 && col[i] != 1.0 &&
                                   col[i]  != 2.0 && col[i] != 3.0)
                         {
                              Error("Derivatives must be integer valued.");
                              ErrNum = INPUT_ERR;
                         }
               }

               if (stricmp(Name, TERM) == 0)
               {
                    /* Duplicates not allowed. */
                    s = MatStrCol(M, j);
                    for (i = 1; i < nRows && ErrNum == OK; i++)
                         if (StrIndex(s[i], s, i) != INDEX_ERR)
                         {
                              Error(DB_DUP_TERM, s[i]);
                              ErrNum = INPUT_ERR;
                         }
               }

               if (stricmp(Name, VARIABLE) == 0)
               {
                    /* Duplicate VARIABLE-TRANSFORMATION not allowed. */
                    s = MatStrCol(M, j);
                    t = MatStrColFind(M, TRANSFORMATION, NO);

                    for (i = 1; i < nRows && ErrNum == OK; i++)
                         for (ii = 0; ii < i; ii++)
                              if (stricmp(s[ii], s[i]) == 0 &&
                                        (t == NULL ||
                                        stricmp(t[ii], t[i]) == 0))
                              {
                                   Error(DB_DUP_VAR, s[i]);
                                   ErrNum = INPUT_ERR;
                                   break;
                              }
               }
          }
     }

     return ErrNum;
}

/*******************************+++*******************************/
int DbMatCor(const Matrix *M)
/*****************************************************************/
/*   Purpose:  Check legality of a correlation matrix.           */
/*                                                               */
/*   Returns:  INPUT_ERR    if illegal;                          */
/*             OK           otherwise.                           */
/*                                                               */
/*   Version:  1994 September 5                                  */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    R;
     real      e;
     size_t    BadIndex, i, j, nRows;

     if ( (nRows = MatNumRows(M)) != MatNumCols(M))
     {
          Error(DB_SQ);
          ErrNum = INPUT_ERR;
     }

     else if ( (BadIndex = StrVecCmp(MatRowNames(M),
               MatColNames(M), nRows)) != INDEX_OK)
     {
          Error(DB_COR_LABEL, MatRowName(M, BadIndex),
                    MatColName(M, BadIndex));
          ErrNum = INPUT_ERR;
     }

     else
           ErrNum = OK;

     for (i = 0; i < nRows && ErrNum == OK; i++)
     {
          if (MatElem(M, i, i) != 1.0)
          {
               Error(DB_COR_DIAG, i + 1);
               ErrNum = INPUT_ERR;
          }

          for (j = 0; j < i && ErrNum == OK; j++)
               if (fabs(e = MatElem(M, i, j)) > 1.0)
               {
                    Error(DB_COR_ONE, i + 1, j + 1);
                    ErrNum = INPUT_ERR;
               }
               else if (e != MatElem(M, j, i))
               {
                    Error(DB_COR_SYM, i + 1, j + 1, j + 1,
                              i + 1);
                    ErrNum = INPUT_ERR;
               }
     }

     if (nRows > 0 && ErrNum == OK)
     {
          /* Check positive definiteness. */
          MatAlloc(nRows, nRows, UP_TRIANG, &R);
          if (TriCholesky(M, 0, &R) != OK)
               Error(DB_POS_DEF);
          MatFree(&R);
     }

     return ErrNum;
}
