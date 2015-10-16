/*****************************************************************/
/*   ROUTINES FOR COMPUTING PREDICTIONS FROM THE MODEL           */
/*   Y = REGRESSION + STOCHASTIC PROCESS                         */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--94.                    */
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
#include "kriging.h"

/*******************************+++*******************************/
int KrigPredSetUp
(
     const KrigingModel  *KrigMod, /* Kriging model and          */
                                   /* decompositions.            */
     real           *ResTildeTilde /* Output: Inv(C) * GLS       */
                                   /* residuals.                 */
)
/*****************************************************************/
/*   Purpose:  Set up for computing predictions *without*        */
/*             standard errors.                                  */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for the       */
/*             vector ResTildeTilde (length n).                  */
/*             For generating prediction coefficients, this      */
/*             seems to be unstable.                             */
/*                                                               */
/*   Version:  1994 November 18                                  */
/*****************************************************************/
{
       return TriBackSolve(KrigChol(KrigMod), KrigMod->ResTilde,
               ResTildeTilde);
}

/*******************************+++*******************************/
void KrigPred(KrigingModel *KrigMod, const Matrix *XPred,
     const real *ResTildeTilde, real *YHat)
/*****************************************************************/
/* Purpose: Compute kriging predictions *without* standard       */
/*          errors for a matrix of points.                       */
/*                                                               */
/* Comment:  Calling routine must allocate space for YHat.       */
/*                                                               */
/* 1996.04.03: Cases in XPred with NA's generate NA for YHat.    */
/* 2009.05.13: KrigCorVec arguments changed                      */
/*****************************************************************/
{
     LinModel  *RegMod, *SPMod;
     real      *Beta, *fRow, *gRow, *r, *xRow;
     size_t    i, n;

     n = MatNumRows(KrigChol(KrigMod));

     RegMod = KrigRegMod(KrigMod);
     SPMod  = KrigSPMod(KrigMod);

     Beta = KrigMod->Beta;

     /* Use workspace in KrigMod. */
     xRow = KrigMod->xRow;
     fRow = KrigMod->fRow;
     gRow = KrigMod->gRow;
     r    = KrigMod->r;

     /* For each prediction. */
     for (i = 0; i < MatNumRows(XPred); i++)
     {
          MatRow(XPred, i, xRow);

          if (VecHasNA(MatNumCols(XPred), xRow))
               YHat[i] = NA_REAL;
          else
          {
               XToF(RegMod, xRow, fRow);
               XToF(SPMod,  xRow, gRow);

               KrigCorVec(gRow, KrigG(KrigMod), n, 0, NULL, YES,
                    KrigMod, r);

               YHat[i] = DotProd(fRow, Beta, ModDF(RegMod))
                         + DotProd(r, ResTildeTilde, n);
          }
     }

     return;
}

/*******************************+++*******************************/
int KrigPredSE(KrigingModel *KrigMod, const Matrix *XPred,
          real *YHat, real *SE)
/*****************************************************************/
/* Purpose:    Compute kriging predictions *with* standard       */
/*             errors for a matrix of points.                    */
/*                                                               */
/* Args:       KrigMod   Kriging model and decompositions.       */
/*             XPred     Prediction points.                      */
/*             YHat      Output: vector of predictions.          */
/*             SE        Output: vector of standard errors.      */
/*                                                               */
/* Returns:    OK or an error number.                            */
/*                                                               */
/* Comment:    Calling routine must allocate space for YHat and  */
/*             SE.                                               */
/*                                                               */
/* 1996.04.03: Cases in XPred with NA's generate NA for YHat.    */
/* 1996.02.18: Argument RAve in KrigYHatSE is passed             */
/*             KrigMod->SPVarProp instead of 1.0 to predict      */
/*             f(x) beta + Z *without* epsilon.                  */
/* 2009.05.13: KrigCorVec arguments changed                      */
/*****************************************************************/
{
     int       ErrNum;
     LinModel  *RegMod, *SPMod;
     Matrix    *G;
     real      *fRow, *gRow, *r, *xRow;
     size_t    i, m;

     G = KrigG(KrigMod);

     RegMod = KrigRegMod(KrigMod);
     SPMod  = KrigSPMod(KrigMod);

     /* Use workspace in KrigMod. */
     xRow = KrigMod->xRow;
     fRow = KrigMod->fRow;
     gRow = KrigMod->gRow;
     r    = KrigMod->r;

     m = MatNumRows(XPred);

     /* For each prediction. */
     ErrNum = OK;
     for (i = 0; i < m && ErrNum == OK; i++)
     {
          MatRow(XPred, i, xRow);

          if (VecHasNA(MatNumCols(XPred), xRow))
               YHat[i] = SE[i] = NA_REAL;
          else
          {
               XToF(RegMod, xRow, fRow);
               XToF(SPMod,  xRow, gRow);

               KrigCorVec(gRow, G, MatNumRows(G), 0, NULL, YES, KrigMod, r);

               /* Need to transform r with T? */

               /* Compute YHat[i] and SE[i]. */
               ErrNum = KrigYHatSE(KrigMod, KrigMod->SPVarProp,
                         fRow, r, &YHat[i], &SE[i]);
          }
     }

         
     if (ErrNum != OK)
          for (i = 0; i < m; i++)
               YHat[i] = SE[i] = NA_REAL;

     return ErrNum;
}

/*******************************+++*******************************/
int KrigYHatSE(KrigingModel *KrigMod, real RAve, real *f, real *r,
          real *YHat, real *SE)
/*****************************************************************/
/*   Purpose:  Compute a kriging prediction and a standard error */
/*             at a single point.                                */
/*                                                               */
/*   Args:     KrigMod   Kriging model and decompositions.       */
/*             RAve      Average R (1.0 for prediction at a      */
/*                       point; used for s.e. of average).       */
/*             f         Input: vector of linear model terms for */
/*                       the new point.                          */
/*                       Output: used for workspace.             */
/*             r         Input: vector of correlations between   */
/*                       the new point and the design points.    */
/*                       Output: used for workspace.             */
/*             YHat      Output: the prediction.                 */
/*             SE        Output: the standard error (computed    */
/*                       only if SE != NULL).                    */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Version:  1995 October 19                                   */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    *Q;
     real      MSE;
     size_t    j, k, n;

     /* Defaults in case of error. */
     *YHat = NA_REAL;
     if (SE != NULL)
          *SE = NA_REAL;

     Q = KrigQ(KrigMod);
     n = MatNumRows(Q);
     k = MatNumCols(Q);

     if ( (ErrNum = KrigTilde(KrigMod, f, r)) == OK)
     {
          *YHat = DotProd(f, KrigMod->RBeta, k)
                    + DotProd(r, KrigMod->ResTilde, n);

          if (SE != NULL)
          {
               for (j = 0; j < k; j++)
                    f[j] -= DotProd(MatCol(Q, j), r, n);

               MSE = KrigMod->SigmaSq
                         * (RAve - VecSS(r, n) + VecSS(f, k));
               *SE = (MSE > 0.0) ? sqrt(MSE) : 0.0;
                           /* 2006.02.14 (Bela's suggestion, not implemented yet): */
                   /* Replace negative with zero in both terms separately? */
          }
     }

     return ErrNum;
}

/*******************************+++*******************************/
int KrigTilde(const KrigingModel *KrigMod, real *f, real *r)
/*****************************************************************/
/*   Purpose:  Overwrite f with fTilde and r with rTilde.        */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Version:  1995 October 19                                   */
/*****************************************************************/
{
     int ErrNum;

     if ( (ErrNum = TriForSolve(KrigR(KrigMod), f, 0, f)) != OK)
          Error("Ill-conditioned expanded-design matrix.\n");

     else if ( (ErrNum = TriForSolve(KrigChol(KrigMod), r, 0, r))
               != OK)
          Error("Ill-conditioned correlation matrix.\n");

     return ErrNum;
}

