/*****************************************************************/
/*   ROUTINES TO EXECUTE FIT COMPUTATIONS                        */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--9.                     */
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

extern size_t       nPointers;

extern boolean      ErrorSave;
extern string       ErrorVar;
extern size_t       ErrorTry;

extern boolean      RanErr;

extern LinModel     RegMod;
extern LinModel     SPMod;

extern Matrix       RegModMat;
extern Matrix       SPModMat;
extern Matrix       T;
extern Matrix       X;
extern Matrix       YDescrip;

extern real         CritLogLikeDiff;
extern real         LogLikeTol;
extern real         *y;
extern size_t       CorFamNum;
extern size_t       ModCompCritNum;
extern size_t       Tries;
extern size_t       nCasesXY;
extern size_t       *IndexXY;
extern string       yName;

static string       SummaryStats[] = {VARIABLE, TRANSFORMATION,
                         CASES, LOG_LIKE, CV_ROOT_MSE, COND_NUM};

/*******************************+++*******************************/
int Fit(void)
/*****************************************************************/
/*   Purpose:  Fit the model parameters.                         */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*                                                               */
/*   1996.01.10: KrigGSpacing moved to KrigModAlloc              */
/*   1996.04.05: Completed removed (temporary output in krmle)   */
/*   1996.04.14: KrigModAlloc/KrigModData; DbIndexXY             */
/*   2009.05.07: Multiple correlation families                   */
/*****************************************************************/
{
     int            ErrNum;
     KrigingModel   KrigMod;
     Matrix         CorPar;
     real           NegLogLike;
     real           *Beta, *CondNum, *CVRootMSE, *ErrVar, *NewCol;
     real           *LogLike, *SPVar;
     size_t         j, jj, k;
     string         ColName;
     ulong          TotEvals;
     unsigned       nEvals;

     k = ModDF(&RegMod);

     Beta = AllocReal(k, NULL);

     /* Add columns to the YDescrip matrix. */
     CVRootMSE  = MatColAdd(CV_ROOT_MSE,  &YDescrip);
     LogLike    = MatColAdd(LOG_LIKE,     &YDescrip);
     SPVar      = MatColAdd(SP_VAR,       &YDescrip);
     ErrVar     = MatColAdd(ERR_VAR,      &YDescrip);
     CondNum    = MatColAdd(COND_NUM,     &YDescrip);

     /* CorPar will hold the correlation parameters */
     /* for one response.                           */
     CorParAlloc(CorFamNum, ModDF(&SPMod), ModTermNames(&SPMod), &CorPar);

     Output("%20s%5s%11s%16s\n", "Variable", "Try", "Iteration",
               "LogLikelihood");

     /* Perform a fit for each response. */
     ErrNum = OK;
     TotEvals = 0;
     ErrorSave = YES;
     for (j = 0; j < MatNumRows(&YDescrip); j++)
     {
          if (DbIndexXY(j) == 0)
               continue;

          ErrorVar = yName;

          /* Set up kriging model. */
          KrigModAlloc(nCasesXY, MatNumCols(&X), yName, &T, &RegMod,
                    &SPMod, CorFamNum, RanErr, &KrigMod);
          KrigModData(nCasesXY, IndexXY, &X, y, &KrigMod);

          ErrNum = FitBest(&KrigMod, Tries, Beta, &CorPar,
                    &SPVar[j], &ErrVar[j], &NegLogLike,
                    &CVRootMSE[j], &nEvals, &CondNum[j]);

          TotEvals += nEvals;

          LogLike[j] = -NegLogLike;

          if (k > 0)
          {
               ColName = StrPaste(3, BETA, ".", yName);
               NewCol = MatColAdd(ColName, &RegModMat);
               AllocFree(ColName);
               VecCopy(Beta, k, NewCol);
          }

          for (jj = 0; jj < MatNumCols(&CorPar); jj++)
          {
               ColName = StrPaste(3, MatColName(&CorPar, jj), ".",
                         yName);
               NewCol = MatColAdd(ColName, &SPModMat);
               AllocFree(ColName);
               VecCopy(MatCol(&CorPar, jj), MatNumRows(&CorPar),
                         NewCol);
          }

          KrigModFree(&KrigMod);
     }

     Output("\n");

     OutputSummary(&YDescrip, NumStr(SummaryStats), SummaryStats);

     Output("Evaluations:     %lu\n", TotEvals);

     AllocFree(Beta);
     MatFree(&CorPar);

     return ErrNum;
}

/*******************************+++*******************************/
int FitBest(KrigingModel *KrigMod, size_t Tries, real *Beta,
     Matrix *CorPar, real *SPVar, real *ErrVar, real *NegLogLike,
     real *CVRootMSE, unsigned *nEvals, real *CondNum)
/*****************************************************************/
/* Purpose:    Choose best of several MLE tries.                 */
/*                                                               */
/* Returns:    OK or an error condition.                         */
/*                                                               */
/* 1996.03.07: First try starts from existing model parameters   */
/*             if they are available.                            */
/* 1996.04.05: Completed removed (temporary output in krmle).    */
/* 1999.04.23: Compare models via user-defined criterion.        */
/*                                                               */
/* Version:    1999.04.23                                        */
/*****************************************************************/
{
     boolean   Better;
     int       ErrNum, ErrThisTry;
     Matrix    RegCorPar;
     real      CondNumTry, CVRootMSETry, MaxErr;
     real      NegLogLikeTry;
     real      *YHatCV;
     size_t    IndexMaxErr, j;
     unsigned  nEvalsTry;

     YHatCV = AllocReal(nCasesXY, NULL);

     ErrNum = !OK;
     *CVRootMSE = REAL_MAX;
     *NegLogLike = REAL_MAX;
     *nEvals = 0;
     *CondNum = NA_REAL;
     for (j = 0; j < Tries; j++)
     {
          /* Try number for error matrix. */
          ErrorTry = j + 1;

          MLEStart(KrigMod, &RegCorPar);

          if (j == 0)
          {
               /* First try: If SPModMat contains correlation   */
               /* parameters, then use them as starting values. */
               CorParExtract(&SPModMat, yName, NO,
                         KrigCorPar(KrigMod));

               if (!RanErr && *SPVar != NA_REAL && *ErrVar != NA_REAL)
                    KrigMod->SPVarProp = *SPVar / (*SPVar + *ErrVar);
          }

          ErrThisTry = MLEFit(&RegCorPar, KrigMod, LogLikeTol,
                    CritLogLikeDiff, j + 1, &NegLogLikeTry,
                    &CondNumTry, &nEvalsTry);
          MatFree(&RegCorPar);

          Better = FALSE;
          if (ErrThisTry == OK &&
                    (ErrThisTry = CalcCV(KrigMod, YHatCV, NULL))
                              == OK)
          {
               CVRootMSETry = RootMSE(nCasesXY, YHatCV, KrigY(KrigMod),
                         &MaxErr, &IndexMaxErr);
               switch (ModCompCritNum)
               {
                    case MOD_COMP_CRIT_CV:
                         if (CVRootMSETry < *CVRootMSE)
                              Better = TRUE;
                         break;          

                    case MOD_COMP_CRIT_LIKE:
                         if (NegLogLikeTry < *NegLogLike)
                              Better = TRUE;
                         break;     

                    default:          
                         CodeBug(ILLEGAL_COND_TXT);
               }          
          }
               
          if (ErrThisTry == OK && Better)                 
          {
               /* One good try is sufficient. */
               ErrNum = OK;

               /* Best parameters so far. */

               *CVRootMSE = CVRootMSETry;
               VecCopy(KrigMod->Beta,
                         ModDF(KrigRegMod(KrigMod)), Beta);

               MatCopy(KrigCorPar(KrigMod), CorPar);

               *SPVar  = KrigMod->SigmaSq * KrigMod->SPVarProp;
               *ErrVar = KrigMod->SigmaSq
                         * (1.0 - KrigMod->SPVarProp);
               *NegLogLike = NegLogLikeTry;
               *CondNum    = CondNumTry;
          }
          *nEvals += nEvalsTry;
     }

     AllocFree(YHatCV);

     return ErrNum;
}
