/*****************************************************************/
/*   ROUTINES TO EXECUTE CROSS-VALIDATION COMPUTATIONS           */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--6.                     */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"
#include "kriging.h"
#include "alex.h"

extern boolean      ErrorSave;
extern string       ErrorVar;

extern boolean      RanErr;

extern LinModel     RegMod;
extern LinModel     SPMod;

extern Matrix       CV;
extern Matrix       SPModMat;
extern Matrix       T;
extern Matrix       X;
extern Matrix       YDescrip;

extern real         *y;
extern size_t       CorFamNum;
extern size_t       nCasesXY;
extern size_t       *IndexXY;
extern string       yName;

static string       SummaryStats[] = {VARIABLE, TRANSFORMATION,
                         CASES, CV_ROOT_MSE, CV_MAX_ERR,
                         CASE_CV_MAX_ERR};

/*******************************+++*******************************/
int CrossValidate(void)
/*****************************************************************/
/*   Purpose:  Compute cross validations and put them in the CV  */
/*             matrix.                                           */
/*                                                               */
/*   1996.04.12: Completed removed.                              */
/*   1996.04.04: X and y include NA's.                           */
/*   1996.04.14: KrigModAlloc/KrigModData; DbIndexXY.            */
/*   2009.05.07: Multiple correlation families                   */
/*****************************************************************/
{
     int            ErrNum, ErrReturn;
     KrigingModel   KrigMod;
     real           *CVMaxErr, *CVRootMSE, *ErrVar, *PredCol;
     real           *SE, *SECol, *SPVar, *YHatCV;
     size_t         IndexMaxErr, j;
     string         ColName;
     string         *CaseMaxErr;

     if (MatNumCols(&T) > 0)
     {
          Error("Transformations not allowed.\n");
          return INPUT_ERR;
     }

     ErrVar = MatColFind(&YDescrip, ERR_VAR, NO);
     SPVar  = MatColFind(&YDescrip, SP_VAR, YES);

     /* Add columns to the YDescrip matrix. */
     CVRootMSE  = MatColAdd(CV_ROOT_MSE,        &YDescrip);
     CVMaxErr   = MatColAdd(CV_MAX_ERR,         &YDescrip);
     CaseMaxErr = MatStrColAdd(CASE_CV_MAX_ERR, &YDescrip);

     /* Compute cross-validation predictions for each response. */
     ErrReturn = OK;
     ErrorSave = YES;
     for (j = 0; j < MatNumRows(&YDescrip); j++)
     {
          if (DbIndexXY(j) == 0)
               continue;

          ErrorVar = yName;

          /* Allocations. */
          YHatCV = AllocReal(nCasesXY, NULL);
          SE     = AllocReal(nCasesXY, NULL);

          /* Set up kriging model. */
          KrigModAlloc(nCasesXY, MatNumCols(&X), yName, &T, &RegMod,
                    &SPMod, CorFamNum, RanErr, &KrigMod);
          KrigModData(nCasesXY, IndexXY, &X, y, &KrigMod);

          /* SPModMat contains the correlation parameters. */
          ErrNum = KrigModSetUp(&SPModMat, yName, SPVar[j],
                    (ErrVar != NULL) ? ErrVar[j] : 0.0, &KrigMod);

          if (ErrNum == OK)
               ErrNum = CalcCV(&KrigMod, YHatCV, SE);

          if (ErrNum == OK)
          {
               /* Put cross validations and standard errors */
               /* in CV.                                    */

               ColName = StrPaste(3, PRED, ".", yName);
               PredCol = MatColAdd(ColName, &CV);
               AllocFree(ColName);
               VecCopyIndex(nCasesXY, NULL, YHatCV, IndexXY,
                         PredCol);

               ColName = StrPaste(3, STD_ERR, ".", yName);
               SECol = MatColAdd(ColName, &CV);
               AllocFree(ColName);
               VecCopyIndex(nCasesXY, NULL, SE, IndexXY, SECol);

               /* Compute summary statistics. */
               CVRootMSE[j] = RootMSE(MatNumRows(&CV), PredCol, y,
                         &CVMaxErr[j], &IndexMaxErr);
               if (IndexMaxErr != INDEX_ERR)
                    CaseMaxErr[j] = StrReplace(
                              MatRowName(&X, IndexMaxErr),
                              CaseMaxErr[j]);

          }

          KrigModFree(&KrigMod);

          if (ErrNum != OK)
               ErrReturn = ErrNum;

          AllocFree(YHatCV);
          AllocFree(SE);
     }

     OutputSummary(&YDescrip, NumStr(SummaryStats), SummaryStats);

     return ErrReturn;
}

/*******************************+++*******************************/
int CalcCV(KrigingModel *KrigMod, real *YHatCV, real *SE)
/*****************************************************************/
/* Purpose:    Compute kriging cross-validation predictions and, */
/*             optionally, their standard errors.                */
/*                                                               */
/* Args:       KrigMod   Input: Kriging model without            */
/*                       decompositions.                         */
/*                       Output: Decompositions are garbage.     */
/*             YHatCV    Output: Cross-validation predictions.   */
/*             SE        Output: Standard errors (computed only  */
/*                       if SE != NULL).                         */
/*                                                               */
/* Returns:    OK or an error number.                            */
/*                                                               */
/* Comment:    Calling routine must allocate space for YHatCV    */
/*             and SE.                                           */
/*             Better matrix updating for doing this?            */
/*             KrigMod decompositions are changed.               */
/*             Standard errors include contribution from epsilon */
/*             in predicted observation.                         */
/* 1995.02.21: SigmaSq not recomputed.                           */
/* 1996.04.12: Temporary output showing progress.                */
/*                                                               */
/* Version:    1996.04.12                                        */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    C, FTilde;
     Matrix    *Chol, *F, *Q, *R;
     real      c, s, t;
     real      *Col, *Beta, *f, *r, *RBeta, *ResTilde, *Y, *YTilde;
     size_t    i, ii, j, k, m, n;

     Y    = KrigY(KrigMod);
     F    = KrigF(KrigMod);
     Chol = KrigChol(KrigMod);
     Q    = KrigQ(KrigMod);
     R    = KrigR(KrigMod);

     /* Use workspace in KrigMod. */
     f        = KrigMod->fRow;
     r        = KrigMod->r;
     RBeta    = KrigMod->RBeta;
     Beta     = KrigMod->Beta;
     ResTilde = KrigMod->ResTilde;

     n = MatNumRows(F);
     k = MatNumCols(F);

     if (n == 0)
          return OK;
     else if (n == 1)
     {
          YHatCV[0] = NA_REAL;
          if (SE != NULL)
              SE[0] = NA_REAL;
          return OK;
     }

     MatAlloc(n, n, UP_TRIANG, &C);
     MatAlloc(n, k, RECT, &FTilde);
     YTilde = AllocReal(n, NULL);

     MatPutNumRows(Q, n - 1);

     /* Put correlation matrix in C. */
     KrigCorMat(0, NULL, KrigMod);
     MatCopy(Chol, &C);

     /* Overwrite correlation matrix with Cholesky decomposition. */
     if (TriCholesky(Chol, 0, Chol) != OK)
     {
          Error("Ill-conditioned Cholesky factor.\n");
          ErrNum = NUMERIC_ERR;
     }
     else
          ErrNum = OK;

     /* Compute FTilde and YTilde for all n rows. */
     if (ErrNum == OK)
          ErrNum = KrigSolve(Chol, F, Y, &FTilde, YTilde);

     /* Delete case i and predict Y[i]. */
     for (i = n - 1, ii = 0; ii < n && ErrNum == OK; ii++, i--)
     {
          OutputTemp("Cross validating variable: %s  Run: %d",
                    yName, i + 1);

          /* Permute adjacent columns of Chol until  */
          /* column i is moved to the last column.   */
          for (j = i; j < n - 1; j++)
          {
               TriPerm(j, j + 1, Chol, &c, &s);

               /* Apply the same rotation to YTilde and FTilde. */
               t           =  c * YTilde[j] + s * YTilde[j+1];
               YTilde[j+1] = -s * YTilde[j] + c * YTilde[j+1];
               YTilde[j]   = t;
               for (m = 0; m < k; m++)
               {
                    Col      = MatCol(&FTilde, m);
                    t        =  c * Col[j] + s * Col[j+1];
                    Col[j+1] = -s * Col[j] + c * Col[j+1];
                    Col[j]   = t;
               }
          }

          /* Correlations between case i and the other cases.  */
          /* Note that cases after i are now in reverse order. */
          for (j = 0; j < i; j++)
               r[j] = MatElem(&C, j, i);
          for (j = 0; j < n - 1 - i; j++)
               r[i+j] = MatElem(&C, i, n - 1 - j);

          /* Linear model terms for case i. */
          MatRow(F, i, f);

          /* Pretend we have only n - 1 cases. */
          MatPutNumRows(Chol, n - 1);
          MatPutNumCols(Chol, n - 1);
          MatPutNumRows(&FTilde, n - 1);

          /* Gram-Schmidt QR orthogonalization of FTilde. */
          if (QRLS(&FTilde, YTilde, Q, R, RBeta, ResTilde) != OK)
          {
               Error("Cannot perform QR decomposition.\n");
               ErrNum = NUMERIC_ERR;
          }

          else
          {
               /* Leave-one-out beta's can be obtained as follows. */
               /*
               if (TriBackSolve(R, RBeta, Beta) != OK)
                    Error("Cannot compute regression beta's.\n");
               else
               {
                    for (j = 0; j < k; j++)
                         Output(" %e", Beta[j]);
                    Output("\n");
               }
               */

               if (SE != NULL)
               {
                    /* Standard error required.             */
                    /* KrigMod->SigmaSq is not updated.     */
                    /* RAve = 1.0 for epsilon contribution. */
                    ErrNum = KrigYHatSE(KrigMod, 1.0, f, r,
                              &YHatCV[i], &SE[i]);
               }
               else
                    /* No standard error. */
                    ErrNum = KrigYHatSE(KrigMod, 1.0, f, r,
                              &YHatCV[i], NULL);
          }

          /* Restore sizes of Chol and FTilde. */
          MatPutNumRows(Chol, n);
          MatPutNumCols(Chol, n);
          MatPutNumRows(&FTilde, n);
     }

     OutputTemp("");

     if (ErrNum != OK)
          for (i = 0; i < n; i++)
               YHatCV[i] = SE[i] = NA_REAL;

     MatPutNumRows(Q, n);

     MatFree(&C);
     MatFree(&FTilde);
     AllocFree(YTilde);

     return ErrNum;
}
