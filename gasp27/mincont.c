/*****************************************************************/
/*   MINIMIZE A FUNCTION OF CONTINUOUS VARIABLES (SUBJECT TO     */
/*   SIMPLE LOWER AND UPPER BOUNDS)                              */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--94.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "min.h"

/* These external variables communicate with ObjFuncUncon. */
static real    (*ObjFuncExt)(real *x, size_t nDims);
static real    *LowBndExt, *UpBndExt, *xExt;

/*******************************+++*******************************/
unsigned MinCont(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          real *LowBnd, real *UpBnd, size_t *Distrib,
          size_t nDims, int MinAlg, real *x, real *fx)
/*****************************************************************/
/*   Purpose:  Optimize continuous x's subject to simple lower   */
/*             and upper bounds by converting to an              */
/*             unconstrained minimization.                       */
/*                                                               */
/*   Args:     ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxFuncs  Maximum function evaluations.           */
/*             LowBnd    LowBnd[j] is the lower bound on x[j].   */
/*             UpBnd     UpBnd[j] is the upper bound on x[j].    */
/*             Distrib   Distrib[j] is the distribution of x[j]. */
/*                       (Only used if MinAlg is SIMPLEX.)       */
/*             nDims     Number of x dimensions.                 */
/*             MinAlg    The minimizer called for the            */
/*                       unconstrained minimization.             */
/*             x         Input:  Starting point.                 */
/*                       Output: "Optimal" point.                */
/*             fx        Input:  Function value at starting x.   */
/*                       Output: "Optimal" function value.       */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   96.03.08: Calls ObjFuncUncon to convert unconstrained       */
/*             ranges to constrained ranges.                     */
/*             Calls MinTryBounds.                               */
/*                                                               */
/*   Version:  1996.03.08                                        */
/*****************************************************************/
{
     real      *Obj, r;
     real      **Dir, **Simp;
     real      *LowBndExtCopy, *UpBndExtCopy, *xExtCopy;
     real      (*ObjFuncExtCopy)(real *x, size_t nDims);
     size_t    i, j;
     unsigned  nEvals;

     /* Save statics to local variables, */
     /* to enable recursive calling.     */
     ObjFuncExtCopy = ObjFuncExt;
     LowBndExtCopy  = LowBndExt;
     UpBndExtCopy   = UpBndExt;
     xExtCopy       = xExt;

     /* External equivalents. */
     ObjFuncExt = ObjFunc;
     LowBndExt  = LowBnd;
     UpBndExt   = UpBnd;

     xExt = AllocReal(nDims, NULL);

     /* Transform x's for unconstrained minimization. */
     for (j = 0; j < nDims; j++)
          x[j] = XToUncon(x[j], LowBnd[j], UpBnd[j]);

     if (MinAlg == SIMPLEXALG)
     {
          /* Allocate space for simplex and objectives. */
          Simp = AllocPtrReal(nDims + 1, NULL);
          for (i = 0; i < nDims + 1; i++)
               Simp[i] = AllocReal(nDims, NULL);
          Obj = AllocReal(nDims + 1, NULL);

          /* First row of simplex is starting point. */
          for (j = 0; j < nDims; j++)
               Simp[0][j] = x[j];
          Obj[0] = *fx;

          /* Remaining rows are random. */
          for (i = 1; i < nDims + 1; i++)
          {
               for (j = 0; j < nDims; j++)
               {
                    /* Random value in constrained range. */
                    r = RegTransformCont(RandUnif(), LowBnd[j],
                              UpBnd[j], Distrib[j]);
                    /* Transform to unconstrained. */
                    Simp[i][j] = XToUncon(r, LowBnd[j], UpBnd[j]);
               }

               Obj[i] = ObjFuncUncon(Simp[i], nDims);
          }
          nEvals = nDims;

          MaxFuncs = (nEvals <= MaxFuncs) ?
                    MaxFuncs - nEvals : 0;
          nEvals += Simplex(ObjFuncUncon, AbsTol, RelTol,
                    MaxFuncs, nDims, Simp, Obj);

          *fx = Obj[0];
          for (j = 0; j < nDims; j++)
               x[j] = Simp[0][j];

          for (i = 0; i < nDims + 1; i++)
               AllocFree(Simp[i]);
          AllocFree(Simp);
          AllocFree(Obj);
     }
     else if (MinAlg == POWELLALG)
     {
          /* Allocate space for directions and initialize them. */
          Dir = AllocPtrReal(nDims, NULL);
          for (i = 0; i < nDims; i++)
          {
               Dir[i] = AllocReal(nDims, NULL);
               for (j = 0; j < nDims; j++)
                    Dir[i][j] = 0.0;
               Dir[i][i] = 1.0;
          }

          nEvals = Powell(ObjFuncUncon, AbsTol, RelTol,
                    MaxFuncs, nDims, x, Dir, fx);

          for (i = 0; i < nDims; i++)
               AllocFree(Dir[i]);
          AllocFree(Dir);
     }

     /* Transform back to constrained x's. */
     for (j = 0; j < nDims; j++)
          x[j] = UnconToX(x[j], LowBnd[j], UpBnd[j]);

     AllocFree(xExt);

     /* Restore statics. */
     ObjFuncExt = ObjFuncExtCopy;
     LowBndExt  = LowBndExtCopy;
     UpBndExt   = UpBndExtCopy;
     xExt       = xExtCopy;

     nEvals += MinTryBounds(ObjFunc, nDims, LowBnd, UpBnd, x, fx);

     return nEvals;
}

/* Scale log transformation to suit optimizers. */
/* #define SHRINK_FACTOR    (500.0)             */

/*******************************+++*******************************/
real XToUncon(real x, real a, real b)
/*****************************************************************/
/* Purpose:    Return unconstrained u corresponding to x in      */
/*             [a, b].                                           */
/*                                                               */
/* 1995.04.03                                                    */
/* 1999.07.01: Restored u = sqrt(x - a) for x in [a, infinity].  */
/*****************************************************************/
{
     real u;

     if (a == -REAL_MAX && b == REAL_MAX)
          /* Already unconstrained. */
          u = x;
     else if (a > -REAL_MAX && b < REAL_MAX)
     {
          /* Closed interval [a, b]. */
          CodeCheck(b - a != 0.0);

          u = asin(2.0 * (x - a) / (b - a) - 1.0);

          /*
          p = (x - a) / (b - a);
          if (p == 1.0)
               u = LN_MAX / SHRINK_FACTOR;
          else
               u = SafeLog(p / (1.0 - p)) / SHRINK_FACTOR;
          */
     }
     else if (a == -REAL_MAX)
          /* (-infinity, b]. */
          u = sqrt(b - x);
          /* u = SafeLog(b - x) / SHRINK_FACTOR; */
     else
          /* [a, +infinity). */
          u = sqrt(x - a); 
          /* u = sqrt(log(x - a + 1.0)); */
          /* u = SafeLog(x - a) / SHRINK_FACTOR; */

     return u;
}

/*******************************+++*******************************/
real UnconToX(real u, real a, real b)
/*****************************************************************/
/*                                                               */
/* Purpose:    Return x in [a, b] corresponding to an            */
/*             unconstrained u.                                  */
/*                                                               */
/* 1994.06.20                                                    */
/* 1999.07.01: Restored x = a + u^2 for x in [a, infinity].      */
/*****************************************************************/
{
     real x;

     if (a == -REAL_MAX && b == REAL_MAX)
          /* Already unconstrained. */
          x = u;

     else if (a > -REAL_MAX && b < REAL_MAX)
          /* Closed interval [a, b]. */
          x = a + 0.5 * (sin(u) + 1.0) * (b - a);

          /* x = a + (b - a) / (1 + SafeExp(-SHRINK_FACTOR * u)); */

     else if (a == -REAL_MAX)
          /* (-infinity, b]. */
          x = b - u * u;
          /* x = b - SafeExp(SHRINK_FACTOR * u); */

     else
          /* [a, +infinity). */
          x = a + u * u;
          /* x = a + SafeExp(u * u) - 1.0; */
          /* x = a + SafeExp(SHRINK_FACTOR * u); */

     return x;
}

/*******************************+++*******************************/
real ObjFuncUncon(real *xUncon, size_t nDims)
/*****************************************************************/
/*   Purpose:  Convert unconstrained x's to their constrained    */
/*             ranges, and call objective.                       */
/*                                                               */
/*   Returns:  Objective.                                        */
/*                                                               */
/*   Version:  1996.03.08                                        */
/*****************************************************************/
{
     size_t    j;

     for (j = 0; j < nDims; j++)
          xExt[j] = UnconToX(xUncon[j], LowBndExt[j], UpBndExt[j]);

     return((*ObjFuncExt)(xExt, nDims));
}

