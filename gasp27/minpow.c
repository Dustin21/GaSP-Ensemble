/*****************************************************************/
/*   DIRECTION SET (POWELL'S) METHOD FOR MULTIDIMENSIONAL        */
/*   UNCONSTRAINED MINIMIZATION                                  */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "min.h"

static real sqrarg;
#define SQR(a) (sqrarg = (a), sqrarg*sqrarg)

/*******************************+++*******************************/
/*   unsigned  Powell(real (*ObjFunc)(real *x, size_t nDims),    */
/*                  real AbsTol, real RelTol, unsigned MaxFuncs, */
/*                  size_t nDims, real *x, real **d,             */
/*                  real *fx)                                   */
/*                                                               */
/*   Purpose:  Multidimensional minimization using Powell's      */
/*             algorithm.                                        */
/*                                                               */
/*   Args:     ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxFuncs  Maximum function evaluations.           */
/*             nDims     Number of dimensions.                   */
/*             x         Input:  Starting point.                 */
/*                       Output: "Optimal" point.                */
/*             d         Input:  nDims x nDims matrix with       */
/*                               columns containing the initial  */
/*                               search directions.              */
/*                       Output: Final direction set.            */
/*             fx        Input:  Function value at x.            */
/*                       Output: "Optimal" function value.       */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*****************************************************************/

unsigned Powell(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          size_t nDims, real *x, real **d, real *fx)
{
     size_t    i, iBig, j;
     real      ave, diff, OldObj, fx0, fx1, MaxDel;
     real      *di =NULL;
     real      *x0 =NULL;
     real      *x1 =NULL;
     unsigned  NumFuncs;

     di = AllocReal(nDims, NULL);
     x0 = AllocReal(nDims, NULL);
     x1 = AllocReal(nDims, NULL);

     /* Save the starting point. */
     for (j = 0; j < nDims; j++)
          x0[j] = x[j];

     NumFuncs = 0;
     while (NumFuncs < MaxFuncs)
     {
          fx0 = *fx;

          /* Loop over all directions in the set.          */
          /* MaxDel will be the biggest function decrease, */
          /* and iBig will index the associated direction. */
          MaxDel = 0.0;
          iBig = 0;
          for (i = 0; i < nDims && NumFuncs < MaxFuncs; i++)
          {
               OldObj = *fx;

               /* Copy direction i and minimize along it. */
               for (j = 0; j < nDims; j++)
                    di[j] = d[j][i];
               NumFuncs += MinLine(ObjFunc, AbsTol / nDims,
                    RelTol / nDims, MaxFuncs - NumFuncs, nDims,
                    x, di, fx);

               if (OldObj - *fx > MaxDel)
               {
                    /*  Largest decrease so far. */
                    MaxDel = OldObj - *fx;
                    iBig = i;
               }
          }

          diff =  fx0 - *fx;
          /* Have to be careful here to avoid overflow. */
          ave = 0.5 * fabs(fx0) + 0.5 * fabs(*fx);

          if (diff <= AbsTol || diff <= RelTol * ave ||
                    nDims == 1 || NumFuncs >= MaxFuncs)
               break;   /* Converged or too many evaluations. */

          /* Construct the extrapolated point and the average */
          /* direction moved, and save the current point.     */
          for (j = 0; j < nDims; j++)
          {
               x1[j] = 2.0 * x[j] - x0[j];
               di[j] = x[j] - x0[j];
               x0[j] = x[j];
          }

          /* Function value at the extrapolated point. */
          fx1 = (*ObjFunc)(x1, nDims);
          NumFuncs++;

          if (fx1 < fx0 && 2.0 * (fx0 - 2.0 * (*fx) + fx1) *
                    SQR(fx0 - *fx - MaxDel) < MaxDel * SQR(fx0 - fx1))
          {
               /* Minimize along new direction. */
               NumFuncs += MinLine(ObjFunc, AbsTol, RelTol,
                         MaxFuncs, nDims, x, di, fx);

               /* Include new direction. */
               for (j = 0; j < nDims; j++)
                    d[j][iBig] = di[j];
          }
     }

     AllocFree(di);
     AllocFree(x0);
     AllocFree(x1);

     return NumFuncs;
}
