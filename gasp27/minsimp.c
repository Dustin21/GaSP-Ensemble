/*****************************************************************/
/*   NELDER-MEAD SIMPLEX ALGORITHM                               */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "min.h"

/*******************************+++*******************************/
/*                                                               */
/*   unsigned  Simplex(real (*ObjFunc)(real *x, size_t nDims),   */
/*                  real AbsTol, real RelTol, unsigned MaxFuncs, */
/*                  size_t nDims, real **p, real *y)             */
/*                                                               */
/*   Purpose:  Function minimization by the Nelder-Mead simplex  */
/*             method.                                           */
/*                                                               */
/*   Args:     ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxFuncs  Maximum function evaluations.           */
/*             nDims     Number of dimensions.                   */
/*             p         Input:  Simplex of nDims + 1 starting   */
/*                               points;                         */
/*                       Output: Final simplex; the best point   */
/*                               is in the first row.            */
/*             y         Input:  Objective values corresponding  */
/*                               to the points in the simplex;   */
/*                       Output: Final objectives, y[0] is the   */
/*                               best.                           */
/*                                                               */
/*   Return:   Number of function evaluations.                   */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*                                                               */
/*****************************************************************/

/* Parameters for expansions and contractions. */
#define ALPHA 1.0
#define BETA  0.5
#define GAMMA 2.0

real *ptry = NULL;   /* Used by SimpTry for trial point. */

unsigned Simplex(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          size_t nDims, real **p, real *y)
{
     size_t    i, ilo, ihi, inhi, j, NumPts;
     real      *psum = NULL;
     real      *PtrSave = NULL;
     real      ave, range, ysave, ytry;
     unsigned  NumFuncs;

     NumFuncs = 0;
     NumPts   = nDims + 1;

     psum = AllocReal(nDims, NULL);
     ptry = AllocReal(nDims, NULL);

     for (j = 0; j < nDims; j++)
          for (i = 0, psum[j] = 0.0; i < NumPts; i++)
               psum[j] += p[i][j];

     while (1)
     {
          /* Determine which point is highest (worst), */
          /* next highest, and lowest (best).          */
          ilo = 0;
          ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
          for (i = 0; i < NumPts; i++)
          {
               if (y[i] < y[ilo])
                    ilo=i;
               if (y[i] > y[ihi])
               {
                    inhi=ihi;
                    ihi=i;
               }
               else if (y[i] > y[inhi])
                    if (i != ihi)
                         inhi=i;
          }

          /* Stop? */
          range = y[ihi] - y[ilo];
          /* Have to be careful here to avoid overflow when */
          /* y[ilo] = y[ihi] = HUGE_VAL.                    */
          ave = 0.5 * fabs(y[ihi]) + 0.5 * fabs(y[ilo]);
          if (range <= AbsTol || range <= RelTol * ave ||
                    NumFuncs >= MaxFuncs)
               break;

          /* Extrapolate by a factor ALPHA through the face of */
          /* the simplex across from the high point.           */
          ytry = SimpTry(p, y, psum, nDims, ObjFunc, ihi,
                    &NumFuncs, -ALPHA);
          if (ytry <= y[ilo])
               /* Gives an improvement, so try an additional */
               /* extrapolation by a factor GAMMA.           */
               ytry = SimpTry(p, y, psum, nDims, ObjFunc, ihi,
                         &NumFuncs, GAMMA);
          else if (ytry >= y[inhi])
          {
               /* Worse than the second highest, so contract. */
               ysave = y[ihi];
               ytry = SimpTry(p, y, psum, nDims, ObjFunc, ihi,
                         &NumFuncs, BETA);
               if (ytry >= ysave)
               {
                    /* Can't get rid of high point.  */
                    /* Contract around th low point. */
                    for (i = 0; i < NumPts; i++)
                    {
                         if (i != ilo)
                         {
                              for (j = 0; j < nDims; j++)
                              {
                                   psum[j] = 0.5 * (p[i][j]
                                             + p[ilo][j]);
                                   p[i][j] = psum[j];
                               }
                              y[i] = (*ObjFunc)(psum, nDims);
                         }
                    }
                    NumFuncs += nDims;
                    for (j = 0; j < nDims; j++)
                         for (i = 0, psum[j] = 0.0; i < NumPts; i++)
                              psum[j] += p[i][j];
               }
          }
     }

     /* Put best point in first row of simplex, */
     /* and its objective in y[0].              */
     PtrSave = p[0];
     p[0]    = p[ilo];
     p[ilo]  = PtrSave;
     ysave   = y[0];
     y[0]    = y[ilo];
     y[ilo]  = ysave;

     AllocFree(psum);
     AllocFree(ptry);
     return(NumFuncs);
}

real SimpTry(real **p, real *y, real *psum, size_t nDims,
          real (*ObjFunc)(real *x, size_t nDims), size_t ihi,
          unsigned *NumFuncs, real fac)

/* Extrapolates by a factor fac across the face of the simplex */
/* from the high point, tries it, and replaces the high point  */
/* if the new point is better.                                 */

{
     size_t    j;
     real      fac1, fac2, ytry;

     fac1 = (1.0 - fac) / nDims;
     fac2 = fac1 - fac;
     for (j = 0; j < nDims; j++)
          ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
     ytry = (*ObjFunc)(ptry, nDims);
     ++(*NumFuncs);

     if (ytry < y[ihi])
     {
          /* Improvement, so replace the highest point. */
          y[ihi] = ytry;
          for (j = 0; j < nDims; j++)
          {
               psum[j] += ptry[j] - p[ihi][j];
               p[ihi][j] = ptry[j];
          }
     }

     return(ytry);
}
