/*****************************************************************/
/*   ROUTINES TO EXECUTE PREDICTIONS COMPUTATIONS                */
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
extern boolean      GenPredCoefs;

extern LinModel     RegMod;
extern LinModel     SPMod;

extern Matrix       PredCoef;
extern Matrix       SPModMat;
extern Matrix       T;
extern Matrix       X;
extern Matrix       XPred;
extern Matrix       YPred;
extern Matrix       YDescrip;

extern real         *y;
extern real         *yTrue;
extern size_t       CorFamNum;
extern size_t       nCasesXY;
extern size_t       *IndexXY;
extern string       yName;

static string       SummaryStats[] = {VARIABLE, TRANSFORMATION,
                         CASES, ROOT_MSE, MAX_ERR, CASE_MAX_ERR};

/*******************************+++*******************************/
int Predict(void)
/*****************************************************************/
/*   Purpose:  Compute predictions, their standard errors, and   */
/*             prediction coefficients.                          */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*                                                               */
/*   1996.04.12: Completed replaced by OutputTemp.               */
/*   1096.04.04: X and y include NA's.                           */
/*   1996.04.14: KrigModAlloc/KrigModData; DbIndexXY.            */
/*   2009.05.07: Multiple correlation families                   */
/*****************************************************************/
{
     boolean        NewXs;
     int            ErrNum, ErrReturn;
     KrigingModel   KrigMod;
     real           *ErrVar, *MaxErr, *NewCol, *ResTildeTilde;
     real           *RMSE, *SE, *SPVar, *yHat;
     size_t         IndexMaxErr, j, m;
     string         ColName;
     string         *CaseMaxErr;

     /* Are there new x's at which to predict? */
     m = MatNumRows(&XPred);
     NewXs = (m > 0);

     /* Check that there are predictions to make and/or */
     /* prediction coefficients to generate.            */
     if (!NewXs && !GenPredCoefs)
     {
          Error("%s is empty and %s = %s: nothing to do!\n",
                    X_PRED, GEN_PRED_COEF, NO_STR);
          return INPUT_ERR;
     }

     ErrVar = MatColFind(&YDescrip, ERR_VAR, NO);
     SPVar  = MatColFind(&YDescrip, SP_VAR, YES);

     /* Add a column to the YDescrip matrix. */
     /* MatPutText(&YDescrip, Y_DESCRIP_TITLE); */

     /* Compute predictions/coefficients for each response. */
     ErrReturn = OK;
     ErrorSave = YES;
     for (j = 0; j < MatNumRows(&YDescrip); j++)
     {
          if (DbIndexXY(j) == 0)
               continue;

          OutputTemp("Predicting variable: %s", yName);

          ErrorVar = yName;

          /* Set up kriging model. */
          KrigModAlloc(nCasesXY, MatNumCols(&X), yName, &T, &RegMod,
                    &SPMod, CorFamNum, RanErr, &KrigMod);
          KrigModData(nCasesXY, IndexXY, &X, y, &KrigMod);

          /* SPModMat contains the correlation parameters. */
          ErrNum = KrigModSetUp(&SPModMat, yName, SPVar[j],
                    (ErrVar != NULL) ? ErrVar[j] : 0.0, &KrigMod);

          if (NewXs && ErrNum == OK)
          {
               yHat = AllocReal(m, NULL);
               SE   = AllocReal(m, NULL);

               ErrNum = KrigPredSE(&KrigMod, &XPred, yHat, SE);

               /* Put predictions and standard errors in YPred. */
               if (ErrNum == OK)
               {
                    ColName = StrPaste(3, PRED, ".", yName);
                    NewCol = MatColAdd(ColName, &YPred);
                    AllocFree(ColName);
                    VecCopy(yHat, m, NewCol);

                    ColName = StrPaste(3, STD_ERR, ".", yName);
                    NewCol = MatColAdd(ColName, &YPred);
                    AllocFree(ColName);
                    VecCopy(SE, m, NewCol);

                    if (yTrue != NULL)
                    {
                         /* Add columns to the YDescrip matrix. */
                         RMSE   = MatColAdd(ROOT_MSE, &YDescrip);
                         MaxErr = MatColAdd(MAX_ERR,  &YDescrip);
                         CaseMaxErr = MatStrColAdd(CASE_MAX_ERR,
                                         &YDescrip);

                         /* Compute summary statistics. */
                         RMSE[j] = RootMSE(m, yHat, yTrue, &MaxErr[j],
                                   &IndexMaxErr);
                         if (IndexMaxErr != INDEX_ERR)
                              CaseMaxErr[j] = StrReplace(
                                        MatRowName(&XPred, IndexMaxErr),
                                        CaseMaxErr[j]);
                    }
               }

               AllocFree(yHat);
               AllocFree(SE);
          }

          if (GenPredCoefs && ErrNum == OK)
          {
               ResTildeTilde = AllocReal(nCasesXY, NULL);

               /* Compute prediction coefficients. */
               /* This seems to be unstable!       */
               ErrNum = TriBackSolve(KrigChol(&KrigMod),
                         KrigMod.ResTilde, ResTildeTilde);

               if (ErrNum == OK)
               {
                    ColName = StrPaste(3, COEF, ".", yName);
                    NewCol = MatColAdd(ColName, &PredCoef);
                    AllocFree(ColName);

                    /* NAs in X or Y give zeros not NAs. */
                    VecInit(0.0, MatNumRows(&PredCoef), NewCol);
                    VecCopyIndex(nCasesXY, NULL, ResTildeTilde,
                              IndexXY, NewCol);
               }

               AllocFree(ResTildeTilde);
          }

          KrigModFree(&KrigMod);

          if (ErrNum != OK)
               ErrReturn = ErrNum;
     }

     OutputTemp("");

     OutputSummary(&YDescrip, NumStr(SummaryStats), SummaryStats);

     return ErrReturn;
}

/*******************************+++*******************************/
void OutputSummary(Matrix *Summ, size_t nCols,
          const string *ColName)
/*****************************************************************/
/*   Purpose:  Output summary information.                       */
/*                                                               */
/*   Version:  1995 October 22                                   */
/*****************************************************************/
{
     Matrix    Out;
     size_t    j, jSumm, jOut;

     MatAllocate(MatNumRows(Summ), 0, RECT, MIXED, NULL, YES,
               &Out);
     MatPutText(&Out, "Summary statistics:\n");

     for (j = 0; j < nCols; j++)
     {
          if ( (jSumm = MatColIndex(Summ, ColName[j])) == INDEX_ERR)
               continue;

          jOut = MatColumnAdd(ColName[j], MatColType(Summ, jSumm),
                    &Out);

          MatCopyCol(jSumm, Summ, jOut, &Out);
     }

     /* Do not write case labels. */
     MatWriteBlock(&Out, NO, stdout);
     MatFree(&Out);
     Output("\n");
}


