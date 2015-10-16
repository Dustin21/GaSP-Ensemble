/*****************************************************************/
/*   ROUTINES FOR MATERN CORRELATION FUNCTION                    */
/*                                                               */
/*   Copyright (c) William J. Welch 2009.                        */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "min.h"
#include "model.h"
#include "kriging.h"

extern size_t  derivMax;
extern size_t  derivMin;
extern real    ThetaStandMax;
extern real    ThetaStandMin;

#define RELTOL      1.0e-10      /* Set small so won't be used. */
#define MAXFUNCS    100

/*******************************+++*******************************/
void MaternAlloc
(
     size_t         NumTerms,      /* Number of terms.           */
     Matrix         *CorPar        /* Output: allocated and      */
                                   /* labelled correlation-      */
                                   /* parameter matrix.          */
)
/*****************************************************************/
/* Purpose: Allocate correlation matrix and label columns.       */
/*                                                               */
/* 2009.05.08: Created                                           */
/*****************************************************************/
{
     MatAllocate(NumTerms, 2, RECT, REAL, NULL, YES, CorPar);

     MatPutText(CorPar, "Matern-family correlation "
          "parameters.\n\n");

     MatPutColName(CorPar, 0, "Theta");
     MatPutColName(CorPar, 1, "Derivatives");

     return;
}

/*******************************+++*******************************/
void MaternStart
(
     const Matrix *G,    /* Expanded-design matrix for the       */
                         /* stochastic-process model.            */
     Matrix *CorPar,     /* Output: Starting values of the       */
                         /* correlation parameters.              */
     Matrix *CorReg      /* Output: Feasibility region for the   */
                         /* correlation parameters.              */
)
/*****************************************************************/
/* Purpose:    Return starting values for the correlation        */
/*             parameters and their optimization region.         */
/*                                                               */
/* Comment:    CorReg is allocated here.                         */
/*                                                               */
/* 2009.05.14: Created.                                          */
/*****************************************************************/
{
     real      Distinct1, Distinct2, Range, thetaMax, thetaMaxTemp;
     real      *deriv, *GCol, *theta;
     size_t    i, j, n, NumDistinct, nTerms;

     n      = MatNumRows(G);
     nTerms = MatNumCols(G);

     RegAlloc(2 * nTerms, CorReg);

     /* Random starting values, bounds, etc. */
     /* for theta's and deriv's.             */
     theta = MatCol(CorPar, 0);
     deriv = MatCol(CorPar, 1);
     for (i = 0; i < nTerms; i++)
     {
          GCol = MatCol(G, i);

          Range = VecMax(GCol, n) - VecMin(GCol, n);

          /* How many distinct values in column i of G? */
          NumDistinct = 1;
          Distinct1 = GCol[0];
          for (j = 1; j < n; j++)
               if (GCol[j] != Distinct1)
               {
                    if (NumDistinct == 1)
                    {
                         NumDistinct = 2;
                         Distinct2 = GCol[j];
                    }
                    else if (GCol[j] != Distinct2)
                    {
                         NumDistinct = 3;
                         break;
                    }
               }     

          /* Output("NumDistinct = %u.\n", NumDistinct); */

          /* Random starting value and region for theta. */

          RegPutDistrib(CorReg, 2 * i, ARCTAN);
          
          if (NumDistinct == 1)
          {
               /* Inactive term. */
               RegPutSupport(CorReg, 2 * i, FIXED);
               RegPutMin(CorReg, 2 * i, 0.0);
               RegPutMax(CorReg, 2 * i, 0.0);
               theta[i] = 0.0;
          }
          else 
          {
               RegPutSupport(CorReg, 2 * i, CONTINUOUS);
               RegPutMin(CorReg, 2 * i, ThetaStandMin / (Range * Range));

               /* Temporary upper bound for random starting value. */
               thetaMaxTemp = (ThetaStandMin +
                         min(ThetaStandMax - ThetaStandMin, 100)) /
                         (double) nTerms / (Range * Range); 
               RegPutMax(CorReg, 2 * i, thetaMaxTemp);
               
               theta[i] = RegRand(CorReg, 2 * i);

               thetaMax = (ThetaStandMax == REAL_MAX) ?
                         REAL_MAX : ThetaStandMax / (Range * Range);

               RegPutMax(CorReg, 2 * i, thetaMax);
          }

          /* Random starting value and region for deriv. */

          RegPutDistrib(CorReg, 2 * i + 1, UNIFORM);
          RegPutMin(CorReg, 2 * i + 1, derivMin);
          RegPutMax(CorReg, 2 * i + 1, derivMax);
 
          if (derivMin == derivMax || NumDistinct <= 2)
          {
               RegPutSupport(CorReg, 2 * i + 1, FIXED);
               deriv[i] = derivMin;
          }
          else
          {
               RegPutSupport(CorReg, 2 * i + 1, GRID);
               RegPutNumLevels(CorReg, 2 * i + 1,
                    derivMax - derivMin + 1);
               RegPutStep(CorReg, 2 * i + 1, 1.0);
               deriv[i] = RegRand(CorReg, 2 * i + 1);
         }
     }

     return;
}

/*******************************+++*******************************/
void MaternCor(
     const real   *g,        /* A point.                         */
     const Matrix *G,        /* Matrix of points.                */
     size_t       n,         /* The correlations for only the    */
                             /* first n rows of G are computed.  */
     size_t       NumActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     const Matrix *CorPar,   /* Correlation parameters.          */
     real         *Cor       /* Output: correlations.            */
)
/*****************************************************************/
/* Purpose:  Compute correlations between the point g and the    */
/*           points in the first n rows of G.                    */
/*                                                               */
/* 2009.05.08: Created                                           */
/*****************************************************************/
{
     real      *deriv, *theta;
     size_t    i, ii;

     VecInit(1.0, n, Cor);

     theta = MatCol(CorPar, 0);
     deriv = MatCol(CorPar, 1);

     /* This implementation multiplies together the 1-d         */
     /* correlations.  It is inefficient because there is one   */
     /* exp call for every dimension (compare with power exp).  */
     /* Mathematically, we could have just one exp call for all */
     /* dimensions, but I am concerned about underflow/overflow */
     /* for deriv = 1 or 2.  The run-time penalty can be found  */
     /* by running, say, Gaussian as a special cases of power   */
     /* exp and Matern, and comparing.                          */

     if (Active == NULL)
          for (i = 0; i < MatNumCols(G); i++)
               MaternCorOneDim(g[i], MatCol(G, i), n, theta[i],
                         deriv[i], Cor);
     else
          for (ii = 0; ii < NumActive; ii++)
          {
               i = Active[ii];
               MaternCorOneDim(g[i], MatCol(G, i), n, theta[i],
                              deriv[i], Cor);
          }

     return;
}


/*******************************+++*******************************/
void MaternCorOneDim(real h, const real *g, size_t n, real theta,
          real deriv, real *Cor)
/*****************************************************************/
/* Purpose:  Multiply the correlations Cor[0],...,Cor[n-1] by    */
/*           the 1-d Matern correlations from the distances      */
/*            between h and g[0],...,g[n-1].                     */
/*                                                               */
/* 2009.05.08: Created                                           */
/*****************************************************************/
{
     real      wtDist;
     size_t    i;

     if (theta == 0.0)
          return;

     if (deriv == 0.0)
          /* Exponential correlation function. */
          for (i = 0; i < n; i++)
          {
               wtDist = theta * fabs(h - g[i]);
               /* Bug in exp(x) for certain small values of x! */
               Cor[i] *= ((wtDist < EPSILON) ? 1.0 : exp(-wtDist));
          }
     
     else if (deriv == 1.0)
          for (i = 0; i < n; i++)
          {
               wtDist = theta * fabs(h - g[i]);
               Cor[i] *= ((wtDist < EPSILON) ? 1.0 : exp(-wtDist)) * (wtDist + 1.0);
          }

     else if (deriv == 2.0)
          for (i = 0; i < n; i++)
          {
               wtDist = theta * fabs(h - g[i]);
               Cor[i] *= ((wtDist < EPSILON) ? 1.0 : exp(-wtDist)) * (wtDist * wtDist / 3 + wtDist + 1.0);
          }

     else if (deriv == 3.0)
          /* Gaussian correlation function. */
          for (i = 0; i < n; i++)
          {
               wtDist = fabs(h - g[i]);
               wtDist *= theta * wtDist; /* weighted squared distance */
               Cor[i] *= ((wtDist < EPSILON) ? 1.0 : exp(-wtDist));
          }
     
     else
          CodeBug("Illegal deriv in MaternCorOneDim.\n");

     return;
}

/*******************************+++*******************************/
unsigned MaternTest
(
     Matrix *CorReg,     /* Feasibility region for the           */
                         /* correlation parameters.              */
     size_t TermIndex,   /* Index of the tested term.            */
     real   AbsTol,
     real   CritLogLikeDiff,
     Matrix *CorPar,     /* Input:  Correlation parameters;      */
                         /* Output: Row TermIndex may change.    */
     real   *NegLogLike  /* Input: Negative log likelihood;      */
                         /* Output: New value.                   */
)
/*****************************************************************/
/* Purpose:  Test whether deriv = derivMax and/or theta = 0 for  */
/*           a single term.                                      */
/*                                                               */
/* Returns:  Number of function evaluations.                     */
/*                                                               */
/* 2009.05.08: Created                                           */
/*****************************************************************/
{
     boolean   Finished;
     real      derivMax, derivOld, derivNullNegLogLike;
     real      thetaMin, thetaOld, thetaNullNegLogLike;
     real      CorParRow[2];
     real      *deriv, *theta;
     Matrix    CorRegSub;
     unsigned  NumFuncs;

     theta  = CorParRow;
     deriv  = CorParRow + 1;
     *theta = MatElem(CorPar, TermIndex, 0);
     *deriv = MatElem(CorPar, TermIndex, 1);

     RegAlloc(2, &CorRegSub);
     MatCopySub(2, MatNumCols(CorReg), 2 * TermIndex, 0, CorReg,
               0, 0, &CorRegSub);

     thetaMin = RegMin(&CorRegSub, 0);
     derivMax = RegMax(&CorRegSub, 1);

     NumFuncs = 0;

     Finished = (*theta == 0.0 || CritLogLikeDiff == 0.0) ? YES : NO;

     if (!Finished && *theta > thetaMin && thetaMin == 0.0)
     {
          /* Test theta = 0.0. */
          thetaOld = *theta;
          *theta = 0.0;
          thetaNullNegLogLike = MLELikeUpdate(CorParRow, 2);
          NumFuncs++;

          if (thetaNullNegLogLike - *NegLogLike < CritLogLikeDiff)
          {
               *NegLogLike = thetaNullNegLogLike;
               Finished = TRUE;
          }
          else
               *theta = thetaOld;
     }

     if (!Finished && *deriv < derivMax &&
               RegSupport(&CorRegSub, 1) != FIXED)
     {
          /* Test deriv = derivMax, re-optimizing theta. */
          thetaOld = *theta;
          derivOld = *deriv;
          *deriv = derivMax;

          RegPutSupport(&CorRegSub, 1, FIXED);
          derivNullNegLogLike = MLELikeUpdate(CorParRow, 2);
          NumFuncs++;
          NumFuncs += MinAnyX(MLELikeUpdate, AbsTol, RELTOL,
                    MAXFUNCS, &CorRegSub, 2, POWELLALG,
                    CorParRow, &derivNullNegLogLike);
          if (derivNullNegLogLike - *NegLogLike < CritLogLikeDiff)
               *NegLogLike = derivNullNegLogLike;
          else
          {
               *theta = thetaOld;
               *deriv = derivOld;
          }
     }

     if (!Finished && *theta > thetaMin)
     {
          /* Test theta = thetaMin. */
          thetaOld = *theta;
          *theta = thetaMin;
          if (thetaMin > 0.0)
          {
               thetaNullNegLogLike = MLELikeUpdate(CorParRow, 2);
               NumFuncs++;
          }

          if (thetaNullNegLogLike - *NegLogLike < CritLogLikeDiff)
               *NegLogLike = thetaNullNegLogLike;
          else
               *theta = thetaOld;
     }

     if (*theta == 0.0)
          *deriv = derivMax;  /* Just for tidiness. */

     MatPutElem(CorPar, TermIndex, 0, *theta);
     MatPutElem(CorPar, TermIndex, 1, *deriv);

     MatFree(&CorRegSub);

     return NumFuncs;
}

/*******************************+++*******************************/
boolean MaternIsActive
(
     const Matrix *CorPar,  /* Correlation parameters           */
     size_t       TermIndex  /* Index of the term of interest    */

)
/*****************************************************************/
/* Purpose: Is term TermIndex active in the correlation          */
/*             function?                                         */
/*                                                               */
/* 2009.05.14: Created                                           */
/*****************************************************************/
{
     return (MatElem(CorPar, TermIndex, 0) != 0.0);
}

