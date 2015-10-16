/*****************************************************************/
/*   ROUTINES FOR ONE-DIMENSIONAL MINIMIZATION                   */
/*****************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "min.h"

#define ZEPS 1.0e-10

#define GOLD   1.618034  /* Default ratio by which successive */
                         /* intervals are magnified.          */
#define CGOLD  0.3819660
#define GLIMIT 100.0     /* Maximum magnification allowed for */
                         /* parabolic-fit step.               */

#define TINY             1.0e-20
#define SIGN(a,b)        ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)    (a) = (b); (b) = (c); (c) = (d)

/*******************************+++*******************************/
unsigned MinBracket(real (*ObjFunc)(real x), real *a, real *x,
          real *b, real *fa, real *fx, real *fb)
/*****************************************************************/
/* Purpose:    Bracket a minimum.                                */
/*                                                               */
/* Args:       ObjFunc   The objective function.                 */
/*             a, x      Input:  Initial distinct points         */
/*                               (a < x).                        */
/*             a, x, b   Output: Points bracketing the minimum.  */
/*             fx        Input:  Function value at x.            */
/*             fa, fx,   Output: Function values at a, x, and b. */
/*             fb                                                */
/*                                                               */
/* Returns:    Number of function evaluations.                   */
/*                                                               */
/* 1990.12.27: Created.                                          */
/* 1999.06.16: [a, b] expanded to get away from an initial flat  */
/*             region.                                           */
/*****************************************************************/
{
     real      fu, q, r, temp, u, ulim;
     unsigned  NumFuncs;

     *fa = (*ObjFunc)(*a);
     NumFuncs = 1;

     while (*fa == *fx && NumFuncs < 10)
     {
          /* Function is flat here: Expand [a, x]. */
          *a = (*a) - GOLD * (*x - *a);
          *x = (*x) + GOLD * (*x - *a);
          *fa = (*ObjFunc)(*a);
          *fx = (*ObjFunc)(*x);
          NumFuncs += 2;
     } 

     if (*fx > *fa)
     {
          /*  Switch roles of a and x so that we can go downhill */
          /* in the direction from a to x.                       */
          temp =  *a;  *a = *x;   *x = temp;
          temp = *fa; *fa = *fx; *fx = temp;
     }

     /* First guess for b. */
     *b = (*x) + GOLD * (*x - *a);
     *fb = (*ObjFunc)(*b);
     NumFuncs++;

     /* Repeat until we bracket. */
     while (*fx > *fb)
     {
          /*  Compute u by parabolic extrapolation from a, x, b.  */
          /*  TINY is used to prevent any possible division by 0. */
          r = (*x - *a) * (*fx - *fb);
          q = (*x - *b) * (*fx - *fa);
          u = (*x) - ((*x - *b) * q - (*x - *a) * r) /
               (2.0 * SIGN(max(fabs(q - r), TINY), q - r));

          /* We won't go any farther than this. */
          ulim = (*x) + GLIMIT * (*b - *x);

          if ((*x - u) * (u - *b) > 0.0)
          {
               /* Parabolic u is between x and b: try it. */
               fu = (*ObjFunc)(u);
               NumFuncs++;
               if (fu < *fb)
               {
                    /* Minimum between x and b. */
                    *a = (*x);
                    *x = u;
                    *fa = (*fx);
                    *fx = fu;
                    break;
               }
               else if (fu > *fx)
               {
                    /* Minimum between a and u.  */
                    *b = u;
                    *fb = fu;
                    break;
               }

               /* Parabolic fit no use: default magnification. */
               u = (*b) + GOLD * (*b - *x);
               fu = (*ObjFunc)(u);
               NumFuncs++;
          }

          else if ((*b - u) * (u - ulim) > 0.0)
          {
               /* Parabolic fit is between b and its allowed limit. */
               fu = (*ObjFunc)(u);
               NumFuncs++;
               if (fu < *fb)
               {
                    SHFT(*x, *b, u, *b + GOLD * (*b - *x));
                    SHFT(*fx, *fb, fu, (*ObjFunc)(u));
                    NumFuncs++;
               }
          }
          else if ((u - ulim) * (ulim - *b) >= 0.0)
          {
               /* Limit parabolic u to maximum allowed value. */
               u = ulim;
               fu = (*ObjFunc)(u);
               NumFuncs++;
          }
          else
          {
               /* Reject parabolic u; use default magnification. */
               u = (*b) + GOLD * (*b - *x);
               fu  =  (*ObjFunc)(u);
               NumFuncs++;
          }

          /* Eliminate oldest point and continue. */
          SHFT(*a, *x, *b, u);
          SHFT(*fa, *fx, *fb, fu);
     }

     if (*a > *b)
     {
          /*  Restore order of a and b. */
          temp =  *a;  *a = *b;   *b = temp;
          temp = *fa; *fa = *fb; *fb = temp;
     }

     return NumFuncs;
}

/*******************************+++*******************************/
/*                                                               */
/*   unsigned  Brent(real (*ObjFunc)(real x), real AbsTol,       */
/*                  real RelTol, unsigned MaxFuncs, real *a,     */
/*                  real *x, real *b, real *fa, real *fx,        */
/*                  real *fb)                                    */
/*                                                               */
/*   Purpose:  One-dimensional minimization using Brent's        */
/*             algorithm.                                        */
/*                                                               */
/*   Args:     ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxFuncs  Maximum function evaluations.           */
/*             a, x, b   Input /   Bracketing triple of          */
/*             fa, fx, fb Output   abscissas and the             */
/*                                 corresponding function value. */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   Comment:  Function values for all of fa, fx, and fb must be */
/*             supplied; they are used for an initial parabolic  */
/*             interpolation.                                    */
/*                                                               */
/*   Version:  1991 February 3                                   */
/*                                                               */
/*****************************************************************/

unsigned Brent(real (*ObjFunc)(real x), real AbsTol, real RelTol,
          unsigned MaxFuncs, real *a, real *x, real *b, real *fa,
          real *fx, real *fb)
{
     real      ave, d, diff, e, etemp, fu, fv, fw;
     real      p, q, r, tol, tol1, tol2, u, v, w, xm;
     unsigned  NumFuncs;

     /* Distance moved on the step before last - force an  */
     /* initial attempt at parabolic interpolation.        */
     e = *b - *a;

     tol = sqrt(EPSILON);

     /* fw is the second smallest function value. */
     if (*fa < *fb)
     {
          w = *a; fw = *fa; v = *b; fv = *fb;
     }
     else
     {
          w = *b; fw = *fb; v = *a; fv = *fa;
     }

     for (NumFuncs = 0; ; NumFuncs++)
     {
          diff = fv - *fx;
          /* Have to be careful here to avoid overflow. */
          ave = 0.5 * fabs(fv) + 0.5 * fabs(*fx);

          if (diff <= AbsTol || diff <= RelTol * ave ||
                    NumFuncs >= MaxFuncs)
               break;   /* Converged or too many evaluations. */

          xm = *a * 0.5 + *b * 0.5;
          tol2 = 2.0 * (tol1 = tol * fabs(*x) + ZEPS);

          if (fabs(e) > tol1)
          {
               /* Construct a trial parabolic fit. */
               r = (*x - w) * (*fx - fv);
               q = (*x - v) * (*fx - fw);
               p = (*x - v) * q - (*x - w) * r;
               q = 2.0 * (q - r);
               if (q > 0.0)
                    p = -p;
               q = fabs(q);
               etemp = e;
               e = d;

               /* Is parabolic fit acceptable? */
               if (fabs(p) >= fabs(0.5 * q * etemp) ||
                         p <= q * (*a - *x) || p >= q * (*b - *x))
                    /* No: take the golden section step into */
                    /* the larger of the two segments.       */
                    d = CGOLD * (e = (*x >= xm ? *a - *x : *b - *x));
               else
               {
                    /* Take the parabolic step. */
                    d = p / q;
                    u = *x + d;
                    if (u - *a < tol2 || *b - u < tol2)
                         d = SIGN(tol1, xm - *x);
               }
          }
          else
               d = CGOLD * (e = (*x >= xm ? *a - *x : *b - *x));

          /* The one function evaluation per iteration. */
          u = (fabs(d) >= tol1 ? *x + d : *x + SIGN(tol1, d));
          fu = (*ObjFunc)(u);

          /* Update for next iteration. */
          if (fu <= *fx)
          {
               if (u >= *x)
                    {*a = *x; *fa = *fx;}
               else
                    {*b = *x, *fb = *fx;}
               SHFT(v, w, *x, u);
               SHFT(fv, fw, *fx, fu);
          }
          else
          {
               if (u < *x)
                    {*a = u; *fa = fu;}
               else
                    {*b = u; *fb = fu;}

               if (fu <= fw || w == *x)
               {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
               }
               else if (fu <= fv || v == *x || v == w)
               {
                    v = u;
                    fv = fu;
               }
          }
     }

     return NumFuncs;
}

/* These external variables communicate with f1dim. */
static real    *ExtD = NULL;
static real    *ExtX = NULL;
static real    *NewX = NULL;
static real    (*ExtObjFunc)(real *x, size_t nDims);
static size_t  ExtnDims;

/*******************************+++*******************************/
unsigned MinLine(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          size_t nDims, real *x, real *d, real *fx)
/*******************************+++*******************************/
/* Purpose:    Minimize a multi-dimensional function along the   */
/*             line from x in the direction d.                   */
/*                                                               */
/* Args:       ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxFuncs  Maximum function evaluations.           */
/*             nDims     Number of dimensions.                   */
/*             x         Input:  Starting point.                 */
/*                       Output: x giving minimum.               */
/*             d         Input:  Direction of line.              */
/*                       Output: Vector displacement that x was  */
/*                               moved.                          */
/*             fx        Input:  Function value at x.            */
/*                       Output: Minimum function value.         */
/*                                                               */
/* Returns:    Number of function evaluations.                   */
/*                                                               */
/* 1992.04.16: Created.                                          */
/* 1999.06.16: Initial check whether d is the zero vector.       */
/*****************************************************************/
{
     size_t    j;
     real      a, b, fa, fb, z;
     unsigned  NumFuncs;

     /* Check that d is not the zero vector. */
     for (j = 0; j < nDims; j++)
          if (d[j] != 0.0)
               break;
     if (j == nDims)
          return 0;          

     NewX = AllocReal(nDims, NULL);

     /* Copy to externals. */
     ExtnDims = nDims;
     ExtObjFunc = ObjFunc;
     ExtX       = x;
     ExtD       = d;

     /* Initial guess for brackets. */
     a = -0.1;
     z = 0.0;
     b = 0.1;

     /* Bracket minimum on line. */
     NumFuncs = MinBracket(f1dim, &a, &z, &b, &fa, fx, &fb);

     /* Minimize along line. */
     MaxFuncs = (NumFuncs <= MaxFuncs) ? MaxFuncs - NumFuncs : 0;
     NumFuncs += Brent(f1dim, AbsTol, RelTol, MaxFuncs, &a, &z, &b,
          &fa, fx, &fb);

     /* New x and d vectors for return. */
     for (j = 0; j < nDims; j++)
     {
          d[j] *= z;
          x[j] += d[j];
     }

     AllocFree(NewX);

     return NumFuncs;
}

real f1dim(real alpha)
{
     size_t    j;

     for (j = 0; j < ExtnDims; j++)
          NewX[j] = ExtX[j] + alpha * ExtD[j];

     return (*ExtObjFunc)(NewX, ExtnDims);
}

