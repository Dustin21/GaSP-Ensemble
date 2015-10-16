/*****************************************************************/
/*   ROUTINES FOR FITTING THE MODEL                              */
/*             Y = REGRESSION + STOCHASTIC PROCESS               */
/*   OPTIMIZING THE CORRELATION PARAMETERS ONE TERM AT A TIME    */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--2011.                  */
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
#include "min.h"
#include "kriging.h"

extern real    SPVarPropMax;
extern real    SPVarPropMin;

extern int     ErrorSeverityLevel;
extern size_t  nPointers;

/* These parameters are used by the continuous-space optimizer. */
/* Note that MAXFUNCS is for a *single* continuous-space        */
/* optimization.                                                */
#define RELTOL      1.0e-10      /* Set small so won't be used. */
#define MAXFUNCS    100

/* These variables are external for communication with the */
/* objective functions, MLELikeObj and MLELikeUpdate, and  */
/* MLELikeScale.                                           */
static KrigingModel *ExtKrigMod;
static Matrix       CPartial;
static size_t       TermIndex;

/* External to avoid repeated allocations in MLEOneTerm. */
static Matrix  RegSub;
static size_t  *Active = NULL;
static real    *CorParRow = NULL;

/* Communicates with MLELike (likelihood calculation). */
static int     OptErr;

/*******************************+++*******************************/
void MLEStart(KrigingModel *KrigMod, Matrix *RegCorPar)
/*****************************************************************/
/*   Purpose:  Put random values of correlation parameters and   */
/*             SPVarProp in KrigMod, and set up optimization     */
/*             region, RegCorPar (which includes SPVarProp).     */
/*                                                               */
/*   Returns:  OK or error condition.                            */
/*                                                               */
/*   1995.03.13: Created?                                        */
/*   2009.05.14: Multiple correlation families                   */
/*   2011.07.06: SPVarProp in [SPVarPropMin, SPVarPropax]        */
/*****************************************************************/
{
     size_t    nPars;

     /* Get starting values of correlation parameters */
     /* and their optimization regions.               */
     CorParStart(KrigCorFam(KrigMod), KrigG(KrigMod),
          KrigCorPar(KrigMod), RegCorPar);

     /* SPVarProp is an extra parameter. If there is no */
     /* random error, the region will fix SPVarProp.    */
     /* This avoids having to continually test RanErr.  */
     nPars = MatNumRows(RegCorPar) + 1;

     /* Add SPVarProp to the optimization region. */
     MatReAlloc(nPars, MatNumCols(RegCorPar), RegCorPar);

     if (KrigRanErr(KrigMod))
     {
          if (SPVarPropMin < SPVarPropMax)
          {
               RegPutMin(    RegCorPar, nPars - 1, SPVarPropMin);
               RegPutMax(    RegCorPar, nPars - 1, SPVarPropMax);
               RegPutSupport(RegCorPar, nPars - 1, CONTINUOUS);
               RegPutDistrib(RegCorPar, nPars - 1, UNIFORM);
               KrigMod->SPVarProp = RegRand(RegCorPar, nPars - 1);
          }
          else
          {
               RegPutSupport(RegCorPar, nPars - 1, FIXED);
               KrigMod->SPVarProp = SPVarPropMin;
          }
     }
     else
     {
          RegPutSupport(RegCorPar, nPars - 1, FIXED);
          KrigMod->SPVarProp = 1.0;
     }

     return;
}

/*******************************+++*******************************/
int MLEFit
(
     Matrix       *RegCorPar, /* Unchanged.                      */
     KrigingModel *KrigMod,
     real         LogLikeTol,
     real         CritLogLikeDiff,
     size_t       Try,
     real         *NegLogLike,/* Output: -log(likelihood).       */
     real         *CondNum,   /* Output: condition number.       */
     unsigned     *TotFuncs   /* Output: likelihood evaluations. */
)
/*****************************************************************/
/* Purpose:    Fit regression and correlation parameters.        */
/*                                                               */
/* Returns:    OK or error condition.                            */
/*                                                               */
/* 1995.04.05: Extrapolation changed to optimizing free          */
/*             parameters                                        */
/* 1996.04.05: Temporary output showing progress.                */
/* 1996.05.27: Active, CorParRow, and RegSub allocated here.     */
/* 1999.07.01: Optimize one term at a time before all terms,     */
/*             to randomize order of terms at beginning.         */
/* 2009.05.14: CorParTest replaces PETest (multiple correlation  */
/*             families)                                         */
/* 2011.08.01: SPVarProp not optimized if support is FIXED       */
/*****************************************************************/
{
     real      AbsTol, CondChol, CondR, SPVarPropSave;
     real      NullNegLogLike, OldNegLogLike;
     real      *CorParVec;
     Matrix    RegSPVarProp;
     Matrix    *Chol, *CorPar, *G;
     size_t    i, Iter, j, kSP, nParsOneTerm, nPars;
     size_t    *Perm, *SupportSave;

     /* Copy to external. */
     ExtKrigMod = KrigMod;

     Chol   = KrigChol(KrigMod);
     CorPar = KrigCorPar(KrigMod);
     G      = KrigG(KrigMod);

     kSP   = MatNumCols(G);
     nPars = MatNumRows(RegCorPar);

     /* Allocations. */
     nParsOneTerm = MatNumCols(CorPar);
     RegAlloc(1, &RegSPVarProp);
     RegAlloc(nParsOneTerm, &RegSub);
     MatAlloc(MatNumRows(Chol), MatNumCols(Chol), UP_TRIANG,
               &CPartial);
     if (kSP > 1)
          Active = AllocSize_t(kSP - 1, NULL);
     CorParRow   = AllocReal(nParsOneTerm, NULL);
     CorParVec   = AllocReal(nPars, NULL);
     Perm        = AllocSize_t(kSP, NULL);
     SupportSave = AllocSize_t(nPars, NULL);

     /* Region for optimizing SPVarProp. */
     MatCopySub(1, MatNumCols(RegCorPar), nPars - 1, 0, RegCorPar,
               0, 0, &RegSPVarProp);

     ErrorSeverityLevel = SEV_WARNING;

     /* Get starting likelihood. */
     KrigCorMat(0, NULL, KrigMod);
     *NegLogLike = MLELike();
     *TotFuncs = 1;

     /* Do until converged. */
     Iter = 1;
     AbsTol = (kSP > 1 || KrigRanErr(KrigMod)) ? 1.0 : LogLikeTol;
     do
     {
          OutputTemp("%20s%5d%11d%16g", KrigYName(KrigMod), Try,
                    Iter, -(*NegLogLike));
          Iter++;

          OldNegLogLike = *NegLogLike;

          /* Each pass through the stochastic-process terms */
          /* will be in a different, random, order.         */
          for (i = 0; i < kSP; i++)
               Perm[i] = i;
          PermRand(kSP, Perm);

          for (i = 0; i < kSP; i++)
          {
               TermIndex = Perm[i];

               /* Optimize correlation parameters */
               /* for term TermIndex.             */
               *TotFuncs += MLEOneTerm(AbsTol, RegCorPar,
                         NegLogLike);

               *TotFuncs += CorParTest(KrigCorFam(KrigMod), RegCorPar, TermIndex,
                         AbsTol, CritLogLikeDiff, CorPar, NegLogLike);
          }

          if (KrigRanErr(KrigMod) && RegSupport(&RegSPVarProp, 0) != FIXED)
          {
               /* Put the unscaled correlation matrix in Chol. */
               SPVarPropSave = KrigMod->SPVarProp;
               KrigMod->SPVarProp = 1.0;
               KrigCorMat(0, NULL, KrigMod);
               KrigMod->SPVarProp = SPVarPropSave;

               /* Then copy Chol to CPartial. */
               /* Change using KrigCorMatC. */
               MatCopy(Chol, &CPartial);

               /* BUG?  Why start with no error variance??? */
               /* What if SPVarProp in some [a, b] ???      */

               /* Negative log likelihood for no error variance. */
               NullNegLogLike = MLELike();
               *TotFuncs += 1;

               /* Optimize SPVarProp. */
               *TotFuncs += MinAnyX(MLELikeScale, AbsTol, RELTOL,
                         MAXFUNCS, &RegSPVarProp, 1, POWELLALG,
                         &KrigMod->SPVarProp, NegLogLike);

               if (NullNegLogLike - *NegLogLike < CritLogLikeDiff)
               {
                    KrigMod->SPVarProp = 1.0;
                    *NegLogLike = NullNegLogLike;
               }
          }

          if (kSP > 1 || KrigRanErr(KrigMod))
          {
               /* Optimize the free parameters          */
               /* (i.e., not on a constraint boundary). */

               MatStack(CorPar, NO, CorParVec);
               CorParVec[nPars-1] = KrigMod->SPVarProp;

               /* Save supports, and change support to fixed */
               /* for parameters on a boundary.              */
               for (i = 0; i < nPars; i++)
               {
                    SupportSave[i] = RegSupport(RegCorPar, i);
                    if (CorParVec[i] == RegMin(RegCorPar, i) ||
                              CorParVec[i] == RegMax(RegCorPar, i))
                         RegPutSupport(RegCorPar, i, FIXED);
               }

               *TotFuncs += MinAnyX(MLELikeObj, AbsTol, RELTOL,
                         MAXFUNCS, RegCorPar, nPars, POWELLALG,
                         CorParVec, NegLogLike);

               MatUnStack(CorParVec, NO, CorPar);
               KrigMod->SPVarProp = CorParVec[nPars-1];
               for (i = 0; i < nPars; i++)
                    RegPutSupport(RegCorPar, i, SupportSave[i]);
          }

          if (OldNegLogLike - *NegLogLike < AbsTol)
               AbsTol /= 4.0;

     } while (AbsTol >= LogLikeTol && (kSP > 1 || KrigRanErr(KrigMod)));

     OutputTemp("");
     Output("%20s%5d%11d%16g\n", KrigYName(KrigMod), Try,
                    Iter, -(*NegLogLike));

     /* Make sure working arrays correspond to optimum,  */
     /* then compute betas, etc.                         */
     ErrorSeverityLevel = SEV_ERROR;
     KrigCorMat(0, NULL, KrigMod);
     *NegLogLike = MLELike();
     if (OptErr != OK)
     {
          Error(NUMERIC_ERR_TXT);
          for (j = 0; j < MatNumCols(CorPar); j++)
               for (i = 0; i < MatNumRows(CorPar); i++)
                    MatPutElem(CorPar, i, j, NA_REAL);
          for (j = 0; j < ModDF(KrigRegMod(KrigMod)); j++)
               KrigMod->Beta[j] = NA_REAL;
          KrigMod->SigmaSq = *NegLogLike = NA_REAL;
     }

     CondChol = TriCond(KrigChol(KrigMod));
     CondR    = TriCond(KrigR(KrigMod));
     *CondNum = max(CondChol, CondR);

     MatFree(&RegSPVarProp);
     MatFree(&RegSub);
     MatFree(&CPartial);

     if (kSP > 1)
          AllocFree(Active);
     AllocFree(CorParRow);
     AllocFree(CorParVec);
     AllocFree(Perm);
     AllocFree(SupportSave);

     return OptErr;
}

/*******************************+++*******************************/
unsigned MLEOneTerm
(
     real AbsTol,             /* Absolute convergence tolerance  */
                              /* on -log(likelihood).            */
     const Matrix *RegCorPar, /* Optimization region for all     */
                              /* parameters.                     */
     real  *NegLogLike        /* Input: Current -log(likelihood);*/
                              /* Output: Optimized value.        */
)
/*****************************************************************/
/*   Purpose:  Optimize correlation parameters for term          */
/*             TermIndex (external), other parameters remaining  */
/*             fixed.                                            */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   96.05.27: Active, CorParRow, and RegSub not allocated here. */
/*                                                               */
/*   Version:  1996.05.27                                        */
/*****************************************************************/
{
     Matrix    *CorPar;
     real      SPVarPropSave;
     size_t    j, kSP, nParsOneTerm;
     unsigned  NumFuncs;

     CorPar = KrigCorPar(ExtKrigMod);

     kSP            = MatNumCols(KrigG(ExtKrigMod));
     nParsOneTerm = MatNumCols(CorPar);

     /* Load row TermIndex of CorPar into row vector CorParRow. */
     MatRow(CorPar, TermIndex, CorParRow);

     /* Load corresponding regions into RegSub. */
     MatCopySub(nParsOneTerm, MatNumCols(RegCorPar),
                nParsOneTerm * TermIndex, 0, RegCorPar, 0, 0,
                &RegSub);

     for(j = 0; j < nParsOneTerm; j++)
          if (RegSupport(&RegSub, j) != FIXED)
               break;
     if (j == nParsOneTerm)
           /* All parameters are fixed. */
           NumFuncs = 0;
     else
     {
          /* Get correlation matrix excluding column TermIndex of G. */
          if (kSP > 1)
          {
               for (j = 0; j < TermIndex; j++)
                    Active[j] = j;
               for (j = TermIndex + 1; j < kSP; j++)
                    Active[j-1] = j;

               /* Scaling will be applied to the correlations */
               /* for TermIndex.                              */
               SPVarPropSave = ExtKrigMod->SPVarProp;
               ExtKrigMod->SPVarProp = 1.0;
               KrigCorMat(kSP - 1, Active, ExtKrigMod);
               ExtKrigMod->SPVarProp = SPVarPropSave;
               MatCopy(KrigChol(ExtKrigMod), &CPartial);
          }

          /* Optimize parameters for term TermIndex. */
          NumFuncs = MinAnyX(MLELikeUpdate, AbsTol, RELTOL, MAXFUNCS,
                    &RegSub, nParsOneTerm, POWELLALG, CorParRow,
                    NegLogLike);

          /* Load CorParRow back into row TermIndex of CorPar. */
          for (j = 0; j < nParsOneTerm; j++)
               MatPutElem(CorPar, TermIndex, j, CorParRow[j]);
     }

     return NumFuncs;
}

/*******************************+++*******************************/
real MLELike()
/*****************************************************************/
/*   Purpose:  Computes the negative log likelihood given by     */
/*             1/2 [log det (C) + n log sigma hat squared].      */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*                                                               */
/*   Comment   Correct correlation matrix must be already loaded */
/*             into ExtKrigMod->Chol.                            */
/*                                                               */
/*   Version:  1995 February 14                                  */
/*****************************************************************/
{
     int       d2;
     real      d1, NegLogLike;
     size_t    n;

     /* Get basic decompositions. */
     if ( (OptErr = KrigDecompose(ExtKrigMod)) != OK)
          /* Have to be careful not to overflow. */
          return sqrt(REAL_MAX);

     /* Chol does not have zeros on diagonal, so d1 > 0.0. */
     TriDet(KrigChol(ExtKrigMod), &d1, &d2);
     n = MatNumRows(KrigChol(ExtKrigMod));
     ExtKrigMod->SigmaSq = VecSS(ExtKrigMod->ResTilde, n) / n;
     NegLogLike = log(d1) + d2 * log(10.0)
               + 0.5 * n * log(ExtKrigMod->SigmaSq);

     /* Put in all constants - especially n. */

     /*
     Output("d1 = %g  d2 = %d  SigmaSq = %g  LogLike = %g\n",
          d1, d2, ExtKrigMod->SigmaSq, -NegLogLike);
     */

     return NegLogLike;
}

/*******************************+++*******************************/
real MLELikeUpdate
(
     real   *CorParRow,  /* Input: Correlation parameters for a  */
                         /* single term (one row of CorPar).     */
     size_t nPars      /* Number of parameters for a single    */
                         /* term.                                */
)
/*****************************************************************/
/*   Purpose:  Called by optimizer to update -log(likelihood)    */
/*             when the correlation parameters change only for a */
/*             single stochastic-process term.                   */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*                                                               */
/*   Version:  1994 September 26                                 */
/*****************************************************************/
{
     Matrix    *Chol, *CorPar;
     real      *Cj, *Cholj;
     size_t    i, j;

     Chol   = KrigChol(ExtKrigMod);
     CorPar = KrigCorPar(ExtKrigMod);

     /* Load row vector CorParRow into row TermIndex of CorPar */
     for (j = 0; j < nPars; j++)
          MatPutElem(CorPar, TermIndex, j, CorParRow[j]);

     /* Put the correlation matrix for the single term */
     /* (scaled for SPVarProp) in Chol.                */
     KrigCorMat(1, &TermIndex, ExtKrigMod);

     if (MatNumRows(CorPar) > 1)
          /* Overwrite Chol with the correlation matrix for */
          /* all variables.                                 */
          for (j = 1; j < MatNumCols(Chol); j++)
          {
               Cholj = MatCol(Chol, j);
               Cj    = MatCol(&CPartial, j);
               for (i = 0; i < j; i++)
                    Cholj[i] *= Cj[i];
          }

     /* Chol will be overwritten by the Cholesky decomposition. */
     return MLELike();
}

/*******************************+++*******************************/
real MLELikeObj
(
     real   *CorParVec,  /* Input: Correlation parameters        */
                         /* arranged as a vector for all terms,  */
                         /* possibly followed by SPVarProp.      */
     size_t nPars      /* Number of parameters.                */
)
/*****************************************************************/
/*   Purpose:  Called by optimizer to compute -log likelihood.   */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*                                                               */
/*   Version:  1994 September 26                                 */
/*****************************************************************/
{
     Matrix    *CorPar;

     CorPar = KrigCorPar(ExtKrigMod);

     /* Load the vector CorParVec into the matrix CorPar. */
     MatUnStack(CorParVec, NO, CorPar);

     /* CorParVec[nPars-1] is the variance proportion. */
     ExtKrigMod->SPVarProp = CorParVec[nPars-1];

     /* Correlation matrix for all terms. */
     KrigCorMat(0, NULL, ExtKrigMod);

     return MLELike();
}

/*******************************+++*******************************/
real MLELikeScale
(
     real   *SPVarProp, /* Scaling proportion for the error      */
                        /* variance.                             */
     size_t nDims       /* Always 1, not used.                   */
)
/*****************************************************************/
/*   Purpose:  Called by optimizer to compute -log likelihood    */
/*             when the off-diagonal elements of CPartial are    */
/*             re-scaled for *SPVarProp.                         */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*                                                               */
/*   Version:  1994 September 26                                 */
/*****************************************************************/
{
     Matrix    *Chol;
     size_t    j;

     Chol = KrigChol(ExtKrigMod);

     /* Copy the unscaled correlation matrix CPartial to Chol. */
     MatCopy(&CPartial, Chol);

     /* Re-scale Chol. */
     if (*SPVarProp < 1.0)
          for (j = 1; j < MatNumCols(Chol); j++)
               VecMultScalar(*SPVarProp, j, MatCol(Chol, j));

     /* Chol will be overwritten by the Cholesky decomposition. */
     return MLELike();

}
