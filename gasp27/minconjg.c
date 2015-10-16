/*****************************************************************/
/*   CONJUGATE GRADIENT (FLETCHER-REEVES-POLAK-RIBIERE) METHOD   */
/*   FOR MULTIDIMENSIONAL UNCONSTRAINED MINIMIZATION             */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "min.h"

static real (*ExtObjFunc)(size_t n, real *x, real *g);

/*******************************+++*******************************/
unsigned MinConjGrad(real (*ObjFunc)(size_t n, real *x, real *g),
     real AbsTol, real RelTol, unsigned MaxEvals,
     size_t n, real *x, real *fx)
/*****************************************************************/
/*   Purpose:  Multidimensional minimization using conjugate     */
/*             gradients.                                        */
/*                                                               */
/*   Args:     ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxEvals  Maximum function evaluations.           */
/*             n         Number of dimensions.                   */
/*             x         Input:  Starting point.                 */
/*                       Output: "Optimal" point.                */
/*             fx        Output: "Optimal" function value.       */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   Comments: fx need not be set before calling.  This          */
/*             should be changed for compatibility with other    */
/*             minimization routines.                            */
/*             Objective functions should be fixed up to deal    */
/*             with gradient/no-gradient methods.                */
/*                                                               */
/*   Version:  1996.03.04                                        */
/*****************************************************************/
{
     size_t   j;
     real     gg, gam, fp, dgg;
     real     *g, *h, *xi;
     unsigned nEvals;

     ExtObjFunc = ObjFunc;

     g  = AllocReal(n, NULL);
     h  = AllocReal(n, NULL);
     xi = AllocReal(n, NULL);
        
     fp = (*ObjFunc)(n, x, xi);

     for (j = 0; j < n; j++)
     {
          g[j] = -xi[j];
          xi[j] = h[j] = g[j];
     }

     nEvals = 0;
     while (nEvals < MaxEvals)
     {
          *fx = fp;
          nEvals += MinLine(ObjFuncNoGrad, AbsTol / n,
                    RelTol / n, MaxEvals - nEvals, n,
                    x, xi, fx);

          if (ApproxEq(fp, *fx, AbsTol, RelTol) ||
                    n == 1 || nEvals >= MaxEvals)
               break;   /* Converged or too many evaluations. */

          fp = (*ObjFunc)(n, x, xi);
          nEvals++;

          gg = VecSS(g, n);

          for (dgg = 0.0, j = 0; j < n; j++)
               /* Polak-Ribiere variation. */
               dgg += (xi[j] + g[j]) * xi[j];

          if (gg == 0.0)
               break;

          gam = dgg / gg;
          for (j = 0; j < n; j++)
          {
               g[j] = -xi[j];
               xi[j] = h[j] = g[j] + gam * h[j];
          }
     }

     AllocFree(g);
     AllocFree(h);
     AllocFree(xi);

     return nEvals;
}

real ObjFuncNoGrad(real *x, size_t n)
{
     return (*ExtObjFunc)(n, x, NULL);
}


