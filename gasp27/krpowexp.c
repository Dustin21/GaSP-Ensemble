/*****************************************************************/
/*   ROUTINES FOR POWER-EXPONENTIAL CORRELATION FUNCTION         */
/*   THETA * |DISTANCE| ** (2 - ALPHA)                           */
/*                                                               */
/*   Note that distance is raised to 2 - Alpha, not p.           */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--99.                    */
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

extern real    AlphaMax;
extern real    AlphaMin;
extern real    ThetaStandMax;
extern real    ThetaStandMin;

#define RELTOL      1.0e-10      /* Set small so won't be used. */
#define MAXFUNCS    100


/*******************************+++*******************************/
void PEAlloc
(
     size_t         NumTerms,      /* Number of terms.           */
     Matrix         *CorPar        /* Output: allocated and      */
                                   /* labelled correlation-      */
                                   /* parameter matrix.          */
)
/*****************************************************************/
/*   Purpose:  Allocate correlation matrix and label columns.    */
/*                                                               */
/*   1995.01.13: Created                                         */
/*   2009.05.07: Row names not assigned (done in CorMatAlloc)    */
/*****************************************************************/
{
     MatAllocate(NumTerms, 2, RECT, REAL, NULL, YES, CorPar);

     MatPutText(CorPar, "Power-exponential-family correlation "
               "parameters.\n\n");

     MatPutColName(CorPar, 0, "Theta");
     MatPutColName(CorPar, 1, "Alpha");

     return;
}


/*******************************+++*******************************/
void PEStart
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
/* 1995.12.07: Created.                                          */
/* 1999.06.17: Some braces added to avoid ambiguous else.        */
/* 1999.06.23: NumTerms renamed nTerms;                          */
/*             Random starting theta is log uniform on range     */
/*             (instead of uniform on first 0.01 of range).      */
/* 1999.06.26: Distribution for theta is ARCTAN.                 */
/* 1999.07.30: ThetaMinTemp renamed ThetaMaxTemp;                */
/*             Region completely specified for Alpha when        */
/*             AlphaMin == AlphaMax || NumDistinct <= 2.         */
/*****************************************************************/
{
     real      Distinct1, Distinct2, Range, ThetaMax, ThetaMaxTemp;
     real      *Alpha, *GCol, *Theta;
     size_t    i, j, n, NumDistinct, nTerms;

     n      = MatNumRows(G);
     nTerms = MatNumCols(G);

     RegAlloc(2 * nTerms, CorReg);

     /* Random starting values, bounds, etc. */
     /* for Theta's and Alpha's.             */
     Theta = MatCol(CorPar, 0);
     Alpha = MatCol(CorPar, 1);
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

          /* Random starting value and region for Theta. */

          RegPutDistrib(CorReg, 2 * i, ARCTAN);
          
          if (NumDistinct == 1)
          {
               /* Inactive term. */
               RegPutSupport(CorReg, 2 * i, FIXED);
               RegPutMin(CorReg, 2 * i, 0.0);
               RegPutMax(CorReg, 2 * i, 0.0);
               Theta[i] = 0.0;
          }
          else 
          {
               RegPutSupport(CorReg, 2 * i, CONTINUOUS);
               RegPutMin(CorReg, 2 * i, ThetaStandMin / (Range * Range));

               /* Temporary upper bound for random starting value. */
               ThetaMaxTemp = (ThetaStandMin +
                         min(ThetaStandMax - ThetaStandMin, 100)) /
                         (double) nTerms / (Range * Range); 
               RegPutMax(CorReg, 2 * i, ThetaMaxTemp);
               
               Theta[i] = RegRand(CorReg, 2 * i);

               ThetaMax = (ThetaStandMax == REAL_MAX) ?
                         REAL_MAX : ThetaStandMax / (Range * Range);

               RegPutMax(CorReg, 2 * i, ThetaMax);
          }

          /* Random starting value and region for Alpha. */

          RegPutDistrib(CorReg, 2 * i + 1, UNIFORM);
          RegPutMin(CorReg, 2 * i + 1, AlphaMin);
 
          if (AlphaMin == AlphaMax || NumDistinct <= 2)
          {
               RegPutSupport(CorReg, 2 * i + 1, FIXED);
               RegPutMax(CorReg, 2 * i + 1, AlphaMax);
               Alpha[i] = AlphaMin;
          }
          else
          {
               RegPutSupport(CorReg, 2 * i + 1, CONTINUOUS);
               
               /* Temporary upper bound for random starting value. */
               /* Default AlphaMin and AlphaMax => Alpha in [0, 0.5]. */
               RegPutMax(CorReg, 2 * i + 1, 0.5 * (AlphaMin + AlphaMax));
               Alpha[i] = RegRand(CorReg, 2 * i + 1);

               /* Upper bound for optimization. */
               RegPutMax(CorReg, 2 * i + 1, AlphaMax);
          }
     }

     return;
}


/*******************************+++*******************************/
void PECor(
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
/*   Purpose:  Compute correlations between the point g and the  */
/*             points in the first n rows of G.                  */
/*                                                               */
/*   Version:  1994 February 2                                   */
/*****************************************************************/
{
     size_t    i;

     /* Put the distances in Cor. */
     PEDist(g, G, n, NumActive, Active, CorPar, Cor);

     /* Exponentiate the sum. */
     for (i = 0; i < n; i++)
          /* Bug in exp(x) for certain small values of x! */
          Cor[i] = (Cor[i] < EPSILON) ? 1.0 : exp(-Cor[i]);

     return;
}

/*******************************+++*******************************/
void PEDist
(
     const real   *g,        /* A point.                         */
     const Matrix *G,        /* Matrix of points.                */
     size_t       n,         /* The distances for only the first */
                             /* n rows of G are computed.        */
     size_t       NumActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     const Matrix *CorPar,   /* Correlation parameters.          */
     real         *Dist      /* Output: distances.               */
)
/*****************************************************************/
/*   Purpose:  Compute distances from the point g to the points  */
/*             in the first n rows of G.                         */
/*                                                               */
/*   Version:  1994 February 2                                   */
/*****************************************************************/
{
     real      *Alpha, *Theta;
     size_t    i, ii;

     Theta = MatCol(CorPar, 0);
     Alpha = MatCol(CorPar, 1);

     VecInit(0.0, n, Dist);

     if (Active == NULL)
          for (i = 0; i < MatNumCols(G); i++)
               PEDistInc(g[i], MatCol(G, i), n, Theta[i],
                         Alpha[i], Dist);
     else
          for (ii = 0; ii < NumActive; ii++)
          {
               i = Active[ii];
               PEDistInc(g[i], MatCol(G, i), n, Theta[i],
                              Alpha[i], Dist);
          }
}

/*******************************+++*******************************/
void PEDistInc(real h, const real *g, size_t n, real Theta,
          real Alpha, real *Dist)
/*****************************************************************/
/*   Purpose:  Increment the distances Dist[0],...,Dist[n-1] for */
/*             the distances between h and g[0],...,g[n-1].      */
/*                                                               */
/*   Version:  1994 February 2                                   */
/*****************************************************************/
{
     real      diff;
     size_t    i;

     if (Theta == 0.0)
          return;

     if (Alpha == 0.0 && Theta != 1.0)
          for (i = 0; i < n; i++)
          {
               diff = h - g[i];
               Dist[i] += Theta * diff * diff;
          }

     else if (Alpha == 0.0 && Theta == 1.0)
          for (i = 0; i < n; i++)
          {
               diff = h - g[i];
               Dist[i] += diff * diff;
          }

     else if (Alpha == 1.0 && Theta != 1.0)
          for (i = 0; i < n; i++)
               Dist[i] += Theta * fabs(h - g[i]);

     else if (Alpha == 1.0 && Theta == 1.0)
          for (i = 0; i < n; i++)
               Dist[i] += fabs(h - g[i]);

     else
          for (i = 0; i < n; i++)
               Dist[i] += Theta * pow(fabs(h - g[i]), 2.0 - Alpha);

     return;
}

/*******************************+++*******************************/
unsigned PETest
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
/*   Purpose:  Test whether Alpha = 0 and/or Theta = 0 for a     */
/*             single term.                                      */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   96.01.31: Fixed bug where old correlation-parameter values  */
/*             not restored after test.                          */
/*   96.05.25: Fixed bug by setting Finished to TRUE when Theta  */
/*             is 0 on entry.                                    */
/*                                                               */
/*   Version:  1996.05.25                                        */
/*****************************************************************/
{
     boolean   Finished;
     real      AlphaMin, AlphaOld, AlphaNullNegLogLike;
     real      ThetaMin, ThetaOld, ThetaNullNegLogLike;
     real      CorParRow[2];
     real      *Alpha, *Theta;
     Matrix    CorRegSub;
     unsigned  NumFuncs;

     Theta  = CorParRow;
     Alpha  = CorParRow + 1;
     *Theta = MatElem(CorPar, TermIndex, 0);
     *Alpha = MatElem(CorPar, TermIndex, 1);

     RegAlloc(2, &CorRegSub);
     MatCopySub(2, MatNumCols(CorReg), 2 * TermIndex, 0, CorReg,
               0, 0, &CorRegSub);

     ThetaMin = RegMin(&CorRegSub, 0);
     AlphaMin = RegMin(&CorRegSub, 1);

     NumFuncs = 0;

     Finished = (*Theta == 0.0 || CritLogLikeDiff == 0.0) ? YES : NO;

     if (!Finished && *Theta > ThetaMin && ThetaMin == 0.0)
     {
          /* Test Theta = 0.0. */
          ThetaOld = *Theta;
          *Theta = 0.0;
          ThetaNullNegLogLike = MLELikeUpdate(CorParRow, 2);
          NumFuncs++;

          if (ThetaNullNegLogLike - *NegLogLike < CritLogLikeDiff)
          {
               *NegLogLike = ThetaNullNegLogLike;
               Finished = TRUE;
          }
          else
               *Theta = ThetaOld;
     }

     if (!Finished && *Alpha > AlphaMin &&
               RegSupport(&CorRegSub, 1) != FIXED)
     {
          /* Test Alpha = AlphaMin, re-optimizing Theta. */
          ThetaOld = *Theta;
          AlphaOld = *Alpha;
          *Alpha = AlphaMin;

          RegPutSupport(&CorRegSub, 1, FIXED);
          AlphaNullNegLogLike = MLELikeUpdate(CorParRow, 2);
          NumFuncs++;
          NumFuncs += MinAnyX(MLELikeUpdate, AbsTol, RELTOL,
                    MAXFUNCS, &CorRegSub, 2, POWELLALG,
                    CorParRow, &AlphaNullNegLogLike);
          if (AlphaNullNegLogLike - *NegLogLike < CritLogLikeDiff)
               *NegLogLike = AlphaNullNegLogLike;
          else
          {
               *Theta = ThetaOld;
               *Alpha = AlphaOld;
          }
     }

     if (!Finished && *Theta > ThetaMin)
     {
          /* Test Theta = ThetaMin. */
          ThetaOld = *Theta;
          *Theta = ThetaMin;
          if (ThetaMin > 0.0)
          {
               ThetaNullNegLogLike = MLELikeUpdate(CorParRow, 2);
               NumFuncs++;
          }

          if (ThetaNullNegLogLike - *NegLogLike < CritLogLikeDiff)
               *NegLogLike = ThetaNullNegLogLike;
          else
               *Theta = ThetaOld;
     }

     if (*Theta == 0.0)
          *Alpha = AlphaMin;  /* Just for tidiness. */

     MatPutElem(CorPar, TermIndex, 0, *Theta);
     MatPutElem(CorPar, TermIndex, 1, *Alpha);

     MatFree(&CorRegSub);

     return NumFuncs;
}

/*******************************+++*******************************/
boolean PEIsActive
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

