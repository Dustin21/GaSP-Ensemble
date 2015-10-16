/*****************************************************************/
/*   ROUTINES TO CHECK COMPATIBILITY OF DATABASE MATRICES        */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--2009.                  */
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

extern boolean      DesignJob;
extern boolean      RanErr;
extern LinModel     RegMod;
extern LinModel     SPMod;
extern Matrix       ANOVAPerc;
extern Matrix       Cand;
extern Matrix       ExpReg;
extern Matrix       PredReg;
extern Matrix       SPModMat;
extern Matrix       XDescrip;
extern Matrix       X;
extern Matrix       YDescrip;
extern size_t       CorFamNum;
extern size_t       n;
extern size_t       nProtected;

extern string       FuncName;

/*****************************************************************/
int DbMatCompat(DbMatrix *D)
/*****************************************************************/
/* Purpose:    Check compatibility of a matrix with other        */
/*             information in the database.                      */
/*                                                               */
/* Returns:    INCOMPAT_ERR or OK.                               */
/*                                                               */
/* 1996.01.18: If FuncName is "Optimize" do not try to extract   */
/*             correlation parameters from SPModMat or the       */
/*             stochastic-process variance from YDescrip.        */
/* 1996.02.13: RegCandCompat called for EXP_REG or PRED_REG      */
/*             rather than for CAND.                             */
/* 1996.03.26: Default ANOVA_PERC set up here, and ANOVA_PERC    */
/*             takes account of averaging w.r.t. groups.         */
/* 1996.04.15: "Optimize" renamed "SequentialDesign".            */
/* 1996.04.07: DbProcessDataMat call replaced by                 */
/*             DbProcessDescrip.                                 */
/*             XX replaced by X.                                 */
/* 1996.02.20: ANOVAPerc matrix not checked if PredReg is empty. */
/*             RegMod, SPMod, and XCor not checked if XDescrip   */
/*             is empty.                                         */
/*             DbMatCatCols not called if XDescrip is empty.     */
/* 2009.05.07: CorParAlloc replaces PEAlloc (multiple            */
/*             correlation families)                             */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    CorParOneY;
     Matrix    *M;
     size_t    GroupSize, h, i, j, nEffects, nGroups;
     size_t    nXVars, nYVars;
     size_t    *GroupVarIndex, *nCats;
     string    Name;
     string    s;
     string    *RowName, *Term, *xName, *yName;

     M    = D->M;
     Name = D->Name;

     /* Check consistency of row labels. */
     if ( (ErrNum = DbRowLabels(D)) != OK)
          return ErrNum;

     nXVars = MatNumRows(&XDescrip);
     xName  = MatStrColFind(&XDescrip, VARIABLE, NO);
     nCats  = MatSize_tColFind(&XDescrip, NUM_CATS, NO);

     /* Specific checks: */

     if (stricmp(Name, ANOVA_PERC) == 0 && !MatEmpty(&PredReg))
     {
          nXVars = MatNumRows(&PredReg);
          xName  = MatStrColFind(&PredReg, VARIABLE, YES);
          
          GroupVarIndex = AllocSize_t(nXVars, NULL);
          for (nGroups = 0, j = 0; j < nXVars; j++)
          {
               RegGroupIndices(&PredReg, j, GroupVarIndex);
               if (GroupVarIndex[0] == j)
                    nGroups++;
          }

          nEffects = nGroups * (nGroups + 1) / 2;

          if (MatEmpty(&ANOVAPerc))
               MatReAlloc(nEffects, 0, &ANOVAPerc);
          else if (nEffects != MatNumRows(&ANOVAPerc))
               ErrNum = INCOMPAT_ERR;

          RowName = MatRowNames(&ANOVAPerc);

          for (i = 0, j = 0; j < nXVars && ErrNum == OK; j++)
          {
               GroupSize = RegGroupIndices(&PredReg, j,
                         GroupVarIndex);

               if (GroupVarIndex[0] != j)
                    continue;

               if (GroupSize == 1)
                    s = StrDup(xName[j]);
               else
                    s = StrPaste(2, GROUP, StrFromSize_t(
                              RegCandGroup(&PredReg, j)));

               if (RowName[i] == NULL)
                    MatPutRowName(&ANOVAPerc, i, s);
               else if (stricmp(RowName[i], s) != 0)
                    ErrNum = INCOMPAT_ERR;

               AllocFree(s);
               i++;
          }

          for (h = nGroups, i = 0; i < nGroups - 1; i++)
          {
               for (j = i + 1; j < nGroups && ErrNum == OK; j++, h++)
               {
                    s = StrPaste(3, RowName[i], ".", RowName[j]);
                    if (RowName[h] == NULL)
                         MatPutRowName(&ANOVAPerc, h, s);
                    else if (stricmp(RowName[h], s) != 0)
                         ErrNum = INCOMPAT_ERR;
                    AllocFree(s);
               }
          }

          AllocFree(GroupVarIndex);

          if (ErrNum != OK)
               Incompatibility(DB_ANOVA);
     }

     else if (stricmp(Name, REG_MOD) == 0 && !MatEmpty(&XDescrip))
          ErrNum = ModParse2(nXVars, xName, nCats, REG_MOD, &RegMod);

     else if (stricmp(Name, SP_MOD) == 0 && !MatEmpty(&XDescrip))
     {
          ErrNum = ModParse2(nXVars, xName, nCats, SP_MOD, &SPMod);
          if (ErrNum == OK && !DesignJob &&
                    stricmp(FuncName, "Fit") != 0 &&
                    stricmp(FuncName, "SequentialDesign") != 0)
          {
               nYVars = MatNumRows(&YDescrip);
               yName  = MatStrColFind(&YDescrip, VARIABLE, NO);

               /* Try to extract correlation-parameter matrix */
               /* for each response.                          */
               for (i = 0; i < nYVars && ErrNum == OK; i++)
               {
                    Term = MatStrColFind(&SPModMat, TERM, YES);
                    CorParAlloc(CorFamNum, MatNumRows(&SPModMat), Term, &CorParOneY);
                    ErrNum = CorParExtract(&SPModMat, yName[i], YES,
                              &CorParOneY);
                    MatFree(&CorParOneY);
               }
          }
     }

     else if (stricmp(Name, X_COR) == 0 && !MatEmpty(&XDescrip))
     {
          for (j = 0; j < MatNumCols(M); j++)
               if (StrIndex(MatColName(M, j), xName, nXVars)
                         == INDEX_ERR)
               {
                    Error(DB_COR_VAR, MatColName(M, j));
                    return INPUT_ERR;
               }
     }

     else if (stricmp(Name, Y_DESCRIP) == 0)
     {
          if (!DesignJob && stricmp(FuncName, "Fit") != 0 &&
                   stricmp(FuncName, "SequentialDesign") != 0 &&
                   MatColFind(&YDescrip, SP_VAR, NO) == NULL)
          {
               Incompatibility(DB_COMPULSORY, Y_DESCRIP, SP_VAR);
               ErrNum = INCOMPAT_ERR;
          }
          if (!DesignJob && stricmp(FuncName, "Fit") != 0 &&
                    RanErr == YES &&
                    MatColFind(&YDescrip, ERR_VAR, NO) == NULL)
          {
               Incompatibility(DB_COMPULSORY, Y_DESCRIP, ERR_VAR);
               ErrNum = INCOMPAT_ERR;
          }
     }

     if (ErrNum == OK && D->IsX && !D->IsProcessed &&
               !MatEmpty(&XDescrip) && !MatEmpty(M))
          ErrNum = DbProcessDescrip(X_DESCRIP, &XDescrip, D);

     if (ErrNum == OK && D->IsY && !D->IsProcessed &&
               !MatEmpty(&YDescrip) && !MatEmpty(M))
          ErrNum = DbProcessDescrip(Y_DESCRIP, &YDescrip, D);

     if (ErrNum == OK && MatType(M) == REAL && !MatEmpty(&XDescrip))
          ErrNum = DbMatCatCols(nXVars, xName, nCats, M);

     if (ErrNum == OK && DesignJob && stricmp(Name, X_MAT) == 0)
          ErrNum = DbMatX(&X, &ExpReg);

     return ErrNum;
}

/*******************************+++*******************************/
int DbRowLabels(const DbMatrix *D)
/*****************************************************************/
/*   Purpose:  Check consistency of row labels with those of up  */
/*             to two other database matrices.                   */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if an inconsistency is found;        */
/*             OK           otherwise.                           */
/*                                                               */
/*   Version:  1994 December 13                                  */
/*****************************************************************/
{
     int  ErrNum;

     ErrNum = OK;

     if (D->RowLabelsComp1 != NULL)
          ErrNum = DbRowLabelsComp(D, D->RowLabelsComp1);

     if (ErrNum == OK && D->RowLabelsComp2 != NULL)
          ErrNum = DbRowLabelsComp(D, D->RowLabelsComp2);

     return ErrNum;
}

/*******************************+++*******************************/
int DbRowLabelsComp(const DbMatrix *D1, const DbMatrix *D2)
/*****************************************************************/
/*   Purpose:  Compare the row labels of D1 and D2.              */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if an inconsistency is found;        */
/*             OK           otherwise.                           */
/*                                                               */
/*   96.04.07: nRowsOrig etc. not members of DbMatrix.           */
/*                                                               */
/*   Version:  1996.04.07                                        */
/*****************************************************************/
{
     int       ErrNum;
     size_t    i, nRows1, nRows2;
     string    *RowLabels1, *RowLabels2;

     nRows1 = MatNumRows(D1->M);
     nRows2 = MatNumRows(D2->M);

     if (nRows1 == 0 || nRows2 == 0)
          ErrNum = OK;

     else if (nRows1 != nRows2)
     {
          Incompatibility(DB_N_ROWS, D1->Name, D2->Name);
          ErrNum = INCOMPAT_ERR;
     }

     else
     {
          ErrNum = OK;

          RowLabels1 = MatRowNames(D1->M);
          RowLabels2 = MatRowNames(D2->M);

          if ( (i = StrVecCmp(RowLabels1, RowLabels2, nRows1))
                    != INDEX_OK)
          {
               Incompatibility(DB_ROW_NAME, D1->Name, D2->Name,
                                   RowLabels1[i], RowLabels2[i]);
               ErrNum = INCOMPAT_ERR;
          }
     }

     return ErrNum;
}

/*******************************+++*******************************/
int DbMatX(const Matrix *X, Matrix *XReg)
/*****************************************************************/
/* Purpose:    Check that X the number of rows in X is           */
/*             compatible with the protected runs and any fixed  */
/*             variables in XReg.                                */
/*                                                               */
/* Returns:    OK or INCOMPAT_ERR.                               */
/*                                                               */
/* 1996.02.20: Does not check X against XReg if XReg is empty.   */
/*                                                               */
/* Version:    1996.02.20                                        */
/*****************************************************************/
{
     int       ErrNum;
     real      *Col;
     size_t    i, nXVars;

     ErrNum = OK;

     if (MatNumRows(X) < nProtected)
     {
          Incompatibility(DB_X_PROT, nProtected);
          ErrNum = INCOMPAT_ERR;
     }

     if (MatEmpty(XReg))
          return ErrNum;

     nXVars = MatNumRows(XReg);

     for (i = 0; i < nXVars && ErrNum == OK; i++)
     {
          if (RegSupport(XReg, i) == FIXED)
          {
               if (MatNumRows(X) != n)
               {
                    Incompatibility(DB_X_N, n);
                    ErrNum = INCOMPAT_ERR;
               }

               else
               {
                    /* Transformation: won't find, e.g., LogX !! */
                    Col = MatColFind(X, RegVar(XReg, i), YES);
                    RegPutMin(XReg, i, VecMin(Col, n));
                    RegPutMax(XReg, i, VecMax(Col, n));
                    RegPutNumLevels(XReg, i, n);
               }
          }
     }

     return ErrNum;
}

/*****************************************************************/
int DbMatCatCols(size_t nXVars, const string *xName,
          const size_t *nCats, const Matrix *M)
/*****************************************************************/
/*   Purpose:  Check that any categorical columns have legal     */
/*             values.                                           */
/*                                                               */
/*   Returns:  INCOMPAT_ERR or OK.                               */
/*                                                               */
/*   96.01.17: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 17                                   */
/*****************************************************************/
{
     int       ErrNum;
     real      *Col;
     size_t    CatIndex, i, j;

     ErrNum = OK;
     for (j = 0; j < MatNumCols(M) && ErrNum == OK; j++)
     {
          if ( (CatIndex = StrIndex(MatColName(M, j), xName,
                    nXVars)) != INDEX_ERR && nCats != NULL &&
                    nCats[CatIndex] > 0)
          {
               /* Column j is a categorical variable. */
               Col = MatCol(M, j);
               CodeCheck(Col != NULL);

               for (i = 0; i < MatNumRows(M); i++)
                    if (Col[i] < 1.0 || Col[i]
                              > (double) nCats[CatIndex] ||
                              fmod(Col[i], 1.0) != 0.0)
                    {
                         Incompatibility(DB_CAT, xName[CatIndex],
                                   nCats[CatIndex]);
                         ErrNum = INCOMPAT_ERR;
                         break;
                    }
          }
     }

     return ErrNum;
}
