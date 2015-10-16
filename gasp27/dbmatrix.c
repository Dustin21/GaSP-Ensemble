/*****************************************************************/
/*   ROUTINES TO INPUT AND OUTPUT DATABASE MATRICES              */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--6.                     */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"
#include "alex.h"

extern int          ErrorSeverityLevel;
extern List         *ColTemplates;

extern boolean      DesignJob;
extern boolean      RanErr;
extern size_t       nProtected;
extern size_t       n;
extern size_t       nXVars;
extern string       FuncName;
extern string       InDir;
extern string       OutDir;
extern string       Prompt;

LinModel  RegMod;
LinModel  SPMod;

Matrix    ANOVAPerc;
Matrix    Cand;
Matrix    CV;
Matrix    ExpReg;
Matrix    JointEff;
Matrix    KPhi;
Matrix    MainEff;
Matrix    PredReg;
Matrix    PredCoef;
Matrix    PriorSamp;
Matrix    RegModMat;
Matrix    SPModMat;
Matrix    T;
Matrix    X;
Matrix    XCor;
Matrix    XDescrip;
Matrix    XPred;
Matrix    YMat;
Matrix    YDescrip;
Matrix    YPred;
Matrix    YTrueMat;

/* Names of optional and compulsory columns. */

static string KPhiCompCol[]   = {"Pi", "PhiStar", NULL};  /* Fix: */

static string JointCompCol[]  = {VARIABLE ".x_i", VARIABLE ".x_j",
                                   VARIABLE ".y", "x_i", "x_j",
                                   "y", STD_ERR, NULL};
static string MainCompCol[]   = {VARIABLE ".x_i", VARIABLE ".y",
                                   "x_i", "y", STD_ERR, NULL};
static string PredCoefOptCol[]= {COEF, NULL};
static string PredOptCol[]    = {PRED, STD_ERR, NULL};
static string RegModOptCol[]  = {BETA, NULL};
static string SPModOptCol[]   = {ALPHA, THETA, "Derivatives", NULL};
static string TermCompCol[]   = {TERM, NULL};
static string VarCompCol[]    = {VARIABLE, NULL};
static string XDescripOptCol[]= {TRANSFORMATION, UNITS, NUM_CATS, 
          SUPPORT,      SUPPORT      "." FIT, SUPPORT      "." PRED,
          MIN,          MIN          "." FIT, MIN          "." PRED,
          MAX,          MAX          "." FIT, MAX          "." PRED,
          NUM_LEVELS,   NUM_LEVELS   "." FIT, NUM_LEVELS   "." PRED,
          DISTRIBUTION, DISTRIBUTION "." FIT, DISTRIBUTION "." PRED,
          CAND_GROUP,   CAND_GROUP   "." FIT, CAND_GROUP   "." PRED,
          INCLUSIVE, WEIGHT, NULL};
static string YDescripCompCol[]= {VARIABLE, SP_VAR, ERR_VAR, NULL};
static string YDescripOptCol[]= {TRANSFORMATION, UNITS, MIN, MAX, ANALYZE,
          CASES, SP_VAR, ERR_VAR, LOG_LIKE, COND_NUM,
          ANOVA_TOTAL_PERC, AVERAGE, STD_ERR "." AVERAGE,
          CASE_MAX_ERR, MAX_ERR, ROOT_MSE, CASE_CV_MAX_ERR,
          CV_MAX_ERR, CV_ROOT_MSE, NULL};

/* Table of matrices. */
static DbMatrix DbMat[] =
{
     {ANOVA_PERC, &ANOVAPerc,  REAL, ANOVA_PERC_TITLE},
     /* No checking! */
     {      CAND,      &Cand,  REAL,             NULL},
     {    CV_MAT,        &CV,  REAL,         CV_TITLE},
     {   EXP_REG,    &ExpReg, MIXED,             NULL},
     { JOINT_EFF,  &JointEff, MIXED,  JOINT_EFF_TITLE},
     {     K_PHI,      &KPhi,  REAL,             NULL},
     {  MAIN_EFF,   &MainEff, MIXED,   MAIN_EFF_TITLE},
     { PRED_COEF,  &PredCoef,  REAL,  PRED_COEF_TITLE},
     {  PRED_REG,   &PredReg, MIXED,             NULL},
     {PRIOR_SAMP, &PriorSamp,  REAL,             NULL},
     {   REG_MOD, &RegModMat, MIXED,    REG_MOD_TITLE},
     {    SP_MOD,  &SPModMat, MIXED,     SP_MOD_TITLE},
     {     T_MAT,         &T,  REAL,             NULL},
     {     X_MAT,         &X,  REAL,          X_TITLE},
     {     X_COR,      &XCor,  REAL,             NULL},
     { X_DESCRIP,  &XDescrip, MIXED,  X_DESCRIP_TITLE},
     {    X_PRED,     &XPred,  REAL,             NULL},
     {     Y_MAT,      &YMat,  REAL,             NULL},
     { Y_DESCRIP,  &YDescrip, MIXED,  Y_DESCRIP_TITLE},
     {    Y_PRED,     &YPred,  REAL,     Y_PRED_TITLE},
     {    Y_TRUE,  &YTrueMat,  REAL,             NULL}
};

#define NUM_MATS    (sizeof(DbMat) / sizeof(DbMatrix))

/******************************+++********************************/
void DbMatInit(void)
/*****************************************************************/
/*   Purpose:  Initialize the matrix database.                   */
/*                                                               */
/*   96.02.16: D->IsX = YES added for CAND so that CAND is       */
/*             processed.                                        */
/*   96.03.06: D->IsX = YES removed for X_PRED so that Min/Max   */
/*             for fitting not applied.                          */
/*   96.03.07: D->IsX = YES restored for X_PRED.                 */
/*   96.03.07: D->IsFit / D->IsPred added.                       */
/*   96.04.07: RemoveNAs etc. not members of DbMatrix.           */
/*                                                               */
/*   Version:  1996.04.07                                        */
/*****************************************************************/
{
     DbMatrix  *D, *DCrossRef;
     size_t    i;
     string    Title;

     D = DbMatFind(ANOVA_PERC, YES);
     D->IsOutput = YES;

     D = DbMatFind(CAND, YES);
     D->IsX = YES;

     D = DbMatFind(CV_MAT, YES);
     D->OptCol = PredOptCol;
     D->RowLabelsComp1 = DbMatFind(X_MAT, YES);
     D->RowLabelsComp2 = DbMatFind(Y_MAT, YES);
     D->IsOutput = YES;

     D = DbMatFind(JOINT_EFF, YES);
     D->CompCol = JointCompCol;
     D->IsOutput = YES;

     D = DbMatFind(K_PHI, YES);
     D->CompCol = KPhiCompCol;

     D = DbMatFind(MAIN_EFF, YES);
     D->CompCol = MainCompCol;
     D->IsOutput = YES;

     D = DbMatFind(PRED_COEF, YES);
     D->OptCol = PredCoefOptCol;
     D->RowLabelsComp1 = DbMatFind(X_MAT, YES);
     D->RowLabelsComp2 = DbMatFind(Y_MAT, YES);
     D->IsOutput = YES;

     D = DbMatFind(REG_MOD, YES);
     D->CompCol = TermCompCol;
     D->OptCol  = RegModOptCol;

     D = DbMatFind(SP_MOD, YES);
     D->CompCol = TermCompCol;
     D->OptCol  = SPModOptCol;

     D = DbMatFind(T_MAT, YES);
     D->RowLabelsComp1 = DbMatFind(X_MAT, YES);
     D->RowLabelsComp2 = DbMatFind(Y_MAT, YES);

     D = DbMatFind(X_MAT, YES);
     D->IsX = YES;
     D->IsFit = YES;
     if (DesignJob)
          D->IsOutput = YES;

     D = DbMatFind(X_DESCRIP, YES);
     D->CompCol = VarCompCol;
     D->OptCol  = XDescripOptCol;

     D = DbMatFind(X_PRED, YES);
     D->IsX = YES;
     D->IsPred = YES;

     D = DbMatFind(Y_MAT, YES);
     D->IsY = YES;
     D->IsFit = YES;
     D->RowLabelsComp1 = DbMatFind(X_MAT, YES);

     D = DbMatFind(Y_DESCRIP, YES);
     D->CompCol = VarCompCol;
     D->OptCol  = YDescripOptCol;
     D->IsOutput = YES;

     D = DbMatFind(Y_PRED, YES);
     D->OptCol = PredOptCol;
     D->RowLabelsComp1 = DCrossRef = DbMatFind(X_PRED, YES);
     D->IsOutput = YES;

     D = DbMatFind(Y_TRUE, YES);
     D->RowLabelsComp1 = DbMatFind(X_PRED, YES);
     D->RowLabelsComp2 = DbMatFind(Y_PRED, YES);
     D->IsY = YES;
     D->IsPred = YES;

     for (i = 0; i < NUM_MATS; i++)
     {
          D = &DbMat[i];
          D->nCompCols = StrNumberOf(D->CompCol);
          D->nOptCols  = StrNumberOf(D->OptCol);

          MatInit(RECT, D->Type, YES, D->M);

          if (D->Title != NULL)
          {
               Title = StrPaste(2, D->Title, "\n\n");
               MatPutText(D->M, Title);
               AllocFree(Title);
          }

          D->FileName = StrDup("");
     }

     return;
}


/******************************+++********************************/
DbMatrix *DbMatFind(const string Name, boolean HardFail)
/*****************************************************************/
/*   Purpose:  Return the DbMat entry for Name.                  */
/*                                                               */
/*   Returns:  Pointer to DbMat[?] or NULL.                      */
/*                                                               */
/*   96.01.17: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 17                                   */
/*****************************************************************/
{
     size_t    j;

     for (j = 0; j < NUM_MATS; j++)
          if (stricmp(Name, DbMat[j].Name) == 0)
               return &DbMat[j];

     /* Not found. */
     Error(DB_MATRIX, Name);

     CodeCheck(!HardFail);

     return NULL;
}

/******************************+++********************************/
int DbMatCheck(DbMatrix *D, string *FileName)
/*****************************************************************/
/*   Purpose:  Check whether matrix D in DbMat is legal and      */
/*             compatible.                                       */
/*             On return, *FileName contains the input file name.*/
/*                                                               */
/*   Returns:  INPUT_ERR, INCOMPAT_ERR, or OK.                   */
/*                                                               */
/*   Version:  1996.05.20                                        */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    *M;
     string    Name;

     Name = D->Name;
     M    = D->M;

     ErrNum = OK;

     if (MatEmpty(M) && DbMatDefault(D) == NO)
     {
          D->FileName = StrDup(EMPTY_STR);
          ErrNum = INPUT_ERR;
     }
     else
     {
          /* Jobs requiring previous fit. ?? */
          if (stricmp(Name, Y_DESCRIP) == 0 &&
                    !DesignJob && stricmp(FuncName, "Fit") != 0)
               D->CompCol = YDescripCompCol;

          if ( (ErrNum = DbMatLegal(D)) == OK)
               ErrNum = DbMatCompat(D);

          if (stricmp(Name, Y_DESCRIP) == 0)
               D->CompCol = VarCompCol;
     }

     *FileName = StrDup(D->FileName);

     return ErrNum;
}

/******************************+++********************************/
boolean DbMatDefault(DbMatrix *D)
/*****************************************************************/
/*   Purpose:  Assign a default to a matrix.                     */
/*                                                               */
/*   Returns:  YES if a default is assigned;                     */
/*             NO  otherwise.                                    */
/*                                                               */
/*   96.03.22: ANOVA_PERC matrix not set up here.                */
/*   96.04.02: Default X_DESCRIP matrix based on N_X_VARS.       */
/*   96.04.07: X instead of XX.                                  */
/*   96.04.07: nRowsOrig etc. not members of DbMatrix.           */
/*   96.06.09: Y (global) renamed YMat to avoid conflict with y. */
/*                                                               */
/*   Version:  1996.06.09                                        */
/*****************************************************************/
{
     boolean   HasDefault;
     DbMatrix  *DCrossRef;
     int       ColType;
     int       *NewColType;
     Matrix    *M;
     size_t    j, nCols, nRows;
     string    Name;
     string    *xName;

     HasDefault = NO;

     Name = D->Name;
     M    = D->M;

     if (stricmp(Name, CAND) == 0 ||
               stricmp(Name, ANOVA_PERC) == 0 ||
               (DesignJob && stricmp(Name, X_MAT) == 0) ||
               (DesignJob && stricmp(Name, Y_MAT) == 0))
     {
          D->FileName = StrDup(EMPTY_STR);
          HasDefault = YES;
     }

     else if (stricmp(Name, SP_MOD) == 0)
     {
          nXVars = MatNumRows(&XDescrip);
          xName  = MatStrColFind(&XDescrip, VARIABLE, NO);

          /* Default model is all x variables. */
          ColType = STRING;
          /* Using MatAllocate loses the title. */
          MatReAllocate(nXVars, 1, &ColType, M);
          MatPutColName(M, 0, TERM);
          VecStrCopy(xName, nXVars, MatStrCol(M, 0));
          D->FileName = StrDup(DEFAULT_STR);
          HasDefault = YES;
     }

     else if (stricmp(Name, X_COR) == 0)
     {
          D->FileName = StrDup(EMPTY_STR);
          HasDefault = YES;
     }

     else if (stricmp(Name, X_DESCRIP) == 0)
     {
          if (!MatEmpty(&X))
          {
               nXVars = MatNumCols(&X);
               xName = MatColNames(&X);
          }
          else if (DesignJob)
               xName = NULL;

          if (nXVars > 0)
          {
               /* Set up default XDescrip based on default x names. */
               DbMatDefDescrip(nXVars, xName, M);
               D->FileName = StrDup(DEFAULT_STR);
               HasDefault = YES;
          }
     }

     else if (stricmp(Name, X_PRED) == 0)
     {
          D->FileName = StrDup(EMPTY_STR);
          HasDefault = YES;
     }

     else if (stricmp(Name, Y_DESCRIP) == 0)
     {
          if (!MatEmpty(&YMat))
          {
               DbMatDefDescrip(MatNumCols(&YMat), MatColNames(&YMat), M);
               D->FileName = StrDup(DEFAULT_STR);
               HasDefault = YES;
          }
     }

     else if (stricmp(Name, Y_TRUE) == 0)
     {
          D->FileName = StrDup(EMPTY_STR);
          HasDefault = YES;
     }

     else if (D->IsOutput)
     {
          if ( (DCrossRef = D->RowLabelsComp1) != NULL)
               nRows = MatNumRows(DCrossRef->M);
          else
               nRows = 0;

          nCols = D->nCompCols;
          NewColType = AllocInt(nCols, NULL);
          for (j = 0; j < nCols; j++)
               NewColType[j] = TemplType(TemplPtr(D->CompCol[j],
                         ColTemplates));
          MatReAllocate(nRows, nCols, NewColType, M);
          AllocFree(NewColType);

          if (DCrossRef != NULL)
               VecStrCopy(MatRowNames(DCrossRef->M), nRows,
                         MatRowNames(M));

          VecStrCopy(D->CompCol, nCols, MatColNames(M));

          D->FileName = StrDup(OUTPUT_STR);
          HasDefault = YES;
     }

     return HasDefault;
}

/******************************+++********************************/
void DbMatDefDescrip(size_t nXVars, const string *xName,
     Matrix *Descrip)
/*****************************************************************/
/*   Purpose:  Assign a default XDescrip or YDescrip.            */
/*                                                               */
/*   Version:  1996.04.05                                        */
/*****************************************************************/
{
     int       ColType;
     size_t    i;
     string    VarName;

     ColType = STRING;

     /* Using MatAllocate loses the title. */
     MatReAllocate(nXVars, 1, &ColType, Descrip);

     MatPutColName(Descrip, 0, VARIABLE);

     if (xName != NULL)
          VecStrCopy(xName, nXVars, MatStrCol(Descrip, 0));
     else
     {
          for (i = 0; i < nXVars; i++)
          {
               VarName = StrPaste(2, "x", StrFromSize_t(i + 1));
               MatPutStrElem(Descrip, i, 0, VarName);
               AllocFree(VarName);
          }
     }
}

/******************************+++********************************/
int DbMatRead(DbMatrix *D, const string FileName)
/*****************************************************************/
/*   Purpose:  Read a matrix from a file and, if legal, put it   */
/*             in the database.                                  */
/*                                                               */
/*   Returns:  INPUT_ERR    if the matrix is illegal;            */
/*             INCOMPAT_ERR if the matrix is incompatible        */
/*                          (with XDescrip or YDescrip);         */
/*             OK           otherwise.                           */
/*                                                               */
/*   96.04.07: nRowsOrig etc. not members of DbMatrix.           */
/*                                                               */
/*   Version:  1996.04.07                                        */
/*****************************************************************/
{
     FILE      *InpFile;
     int       ErrNum, TempType, Type;
     string    DirFileName;
     Matrix    *M;

     D->FileName = StrReplace(FileName, D->FileName);

     if (stricmp(InDir, DEF_IN_DIR) != 0)
     {
          DirFileName = StrPaste(3, InDir, DIR_SEP, FileName);
          InpFile = FileOpen(DirFileName, "r");
          AllocFree(DirFileName);
     }
     else
          InpFile = FileOpen(FileName, "r");

     if (InpFile == NULL)
          return INPUT_ERR;

     M = D->M;

     MatFree(M);

     /* A mixed-type matrix will be read as a string matrix. */
     /* It will be converted in DbMatColChk.                 */
     Type = MatType(M);
     TempType = (Type == MIXED) ? STRING : Type;
     ErrNum = MatRead(InpFile, TempType, M);

     /* Restore the proper matrix type. */
     MatPutType(M, Type);

     if (ErrNum == OK)
     {
          Output(DB_ROWS_COLS,
                    MatNumRows(M), StrPlural(MatNumRows(M)),
                    MatNumCols(M), StrPlural(MatNumCols(M)));

          ErrNum = DbMatLegal(D);
     }

     if (ErrNum == INPUT_ERR)
          MatFree(M);

     if (ErrNum == OK && stricmp(D->Name, X_DESCRIP) == 0)
          DbNewDescrip(X_DESCRIP);

     else if (ErrNum == OK && stricmp(D->Name, Y_DESCRIP) == 0)
          DbNewDescrip(Y_DESCRIP);

     return ErrNum;
}

/******************************+++********************************/
void DbNewDescrip(const string DescripName)
/*****************************************************************/
/*   Purpose:  Free data matrices for a new Descrip matrix.      */
/*                                                               */
/*   Version:  1995 October 12                                   */
/*****************************************************************/
{
     DbMatrix  *D;
     int       ErrNum, ErrSevSave;
     size_t    i;
     string    Title;

     ErrNum = OK;
     for (i = 0; i < NUM_MATS && ErrNum == OK; i++)
     {
          D = &DbMat[i];
          if ( (stricmp(DescripName, X_DESCRIP) == 0 && D->IsX) ||
                    (stricmp(DescripName, Y_DESCRIP) == 0 && D->IsY))
          {
               if (D->IsProcessed)
               {
                    /* Data matrix had already been processed. */
                    MatFree(D->M);

                    ErrSevSave = ErrorSeverityLevel;
                    ErrorSeverityLevel = SEV_WARNING;
                    Error("New %s: Matrix %s is being emptied.\n",
                              DescripName, D->Name);
                    ErrorSeverityLevel = ErrSevSave;

                    if (D->Title != NULL)
                    {
                         Title = StrPaste(2, D->Title, "\n\n");
                         MatPutText(D->M, Title);
                         AllocFree(Title);
                    }
                    D->IsProcessed = NO;
               }
               /*
               else if (MatNumRows(D->M) > 0)
                    Sets IsProcessed.
                    ErrNum = DbProcessDataMat(DescripName, Descrip, D);
               */
          }
     }
}

/******************************+++********************************/
int DbMatWrite(DbMatrix *D, const string FileName,
          const string BlockingOption)
/*****************************************************************/
/*   Purpose:  Output a matrix to a file.                        */
/*             BlockingOption should be BLOCKED or UNBLOCKED;    */
/*             a NULL value is the same as BLOCKED.              */
/*                                                               */
/*   Returns:  OK        if successful;                          */
/*             INPUT_ERR otherwise.                              */
/*                                                               */
/*   Version:  1995 March 10                                     */
/*****************************************************************/
{
     boolean   Blocked;
     FILE      *OutFile;
     Matrix    *M;
     string    DirFileName;

     if (stricmp(FileName, SCREEN) == 0)
          OutFile = stdout;
     else
     {
          DirFileName = StrPaste(3, OutDir, DIR_SEP, FileName);
          OutFile = FileOpen(DirFileName, "w");
          AllocFree(DirFileName);
          if (OutFile == NULL)
               return INPUT_ERR;
     }

     /* Default is blocked output. */
     if (BlockingOption == NULL || *BlockingOption == NULL ||
               stricmp(BlockingOption, BLOCKED) == 0)
          Blocked = YES;
     else if (stricmp(BlockingOption, UNBLOCKED) == 0)
          Blocked = NO;
     else
     {
          Error(DB_BLOCK);
          return INPUT_ERR;
     }

     M = D->M;

     if (Blocked)
           MatWriteBlock(M, YES, OutFile);
     else
           MatWrite(M, OutFile);

     return OK;
}

