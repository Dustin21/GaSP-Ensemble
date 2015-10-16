/*****************************************************************/
/*   ROUTINES FOR MISCELLANEOUS MATH                             */
/*                                                               */
/*   Copyright (c) William J. Welch 1992--95.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

extern int     ErrorSeverityLevel;

/*******************************+++*******************************/
boolean ApproxEq(real a, real b, real AbsTol, real RelTol)
/*****************************************************************/
/*   Purpose:  Are a and b approximately equal?                  */
/*                                                               */
/*   Returns:  0 or 1.                                           */
/*                                                               */
/*   Version:  1992 July 29                                      */
/*****************************************************************/
{
     real Ave, Diff;

     Diff =  fabs(a - b);

     /* Have to be careful here to avoid overflow. */
     Ave = 0.5 * fabs(a) + 0.5 * fabs(b);

     return (Diff <= AbsTol || Diff <= RelTol * Ave);
}

/*******************************+++*******************************/
real SafeExp(real x)
/*****************************************************************/
/*   Purpose:  Safe version of exp(x).                           */
/*                                                               */
/*   Returns:  exp(x) or NA_REAL.                                */
/*                                                               */
/*   Comment:  Change like SafeLog10.                            */
/*                                                               */
/*   Version:  1995 May 2                                        */
/*****************************************************************/
{
     if (x == NA_REAL)
     {
          ErrorSeverityLevel = SEV_WARNING;
          Error("exp(NA) gives NA.\n");
          ErrorSeverityLevel = SEV_ERROR;
          return NA_REAL;
     }
     else if (fabs(x) < EPSILON)
          /* Bug in exp(x) for certain small values of x! */
          return 1.0;
     else if (x <= LN_MIN)
          return 0.0;
     else if (x >= LN_MAX)
          return REAL_MAX;
     else
          return exp(x);
}

/*******************************+++*******************************/
void SafeLog10(size_t n, const real *x, real *log10x,
     size_t *nInputNA, size_t *nDomErr)
/*****************************************************************/
/*   Purpose:  Safe version of log10.                            */
/*             On return, log10x[i] = NA_REAL if x[i] = NA_REAL, */
/*                                    NA_REAL if x[i] <= 0.0,    */
/*                                    log10(x[i]) otherwise,     */
/*             *nInputNA counts the number of NA_REALs in x,     */
/*             and *nDomErr counts the number of x[i] <= 0       */
/*             occurrences.                                      */
/*                                                               */
/*   Version:  1994 October 12                                   */
/*****************************************************************/
{
     size_t    i, nDom, nNA;

     for (nDom = 0, nNA = 0, i = 0; i < n; i++)
     {
          if (x[i] == NA_REAL)
          {
               log10x[i] = NA_REAL;
               nNA++;
          }
          else if (x[i] <= 0.0)
          {
               log10x[i] = NA_REAL;
               nDom++;
          }
          else
               log10x[i] = log10(x[i]);
     }

     *nDomErr = nDom;
     *nInputNA = nNA;

     return;
}

/*******************************+++*******************************/
real RootMSE(size_t n, const real *y1, const real *y2,
          real *MaxErr, size_t *IndexMaxErr)
/*****************************************************************/
/*   Purpose:  Compute the root MSE and the maximum error        */
/*             between y1 and y2, ignoring any cases with NA's.  */
/*                                                               */
/*   Args:     n         Length of vectors.                      */
/*                       y1        First vector.                 */
/*             y2        Second vector.                          */
/*             MaxErr    Output: Maximum absolute error.         */
/*             IndexMaxErr                                       */
/*                       Output: Case with the maximum absolute  */
/*                       error.                                  */
/*                                                               */
/*   Returns:  The root MSE.                                     */
/*                                                               */
/*   Version:  1994 November 9                                   */
/*****************************************************************/
{
     real      AbsDiff, AbsMaxErr, Diff, rMSE;
     size_t    i, nNotMissing;

     nNotMissing = 0;
     rMSE = *MaxErr = AbsMaxErr = 0.0;
     *IndexMaxErr = 0;

     for (i = 0; i < n; i++)
     {
          if (y1[i] == NA_REAL || y2[i] == NA_REAL)
               continue;

          nNotMissing++;
          AbsDiff = fabs(Diff = y1[i] - y2[i]);
          rMSE += Diff * Diff;
          if (AbsDiff > AbsMaxErr)
          {
               AbsMaxErr   = AbsDiff;
               *MaxErr     = Diff;
               *IndexMaxErr = i;
          }
     }

     if (nNotMissing > 0)
          rMSE = sqrt(rMSE / nNotMissing);
     else
     {
          rMSE = *MaxErr = NA_REAL;
          *IndexMaxErr = INDEX_ERR;
     }

     return rMSE;
}

/*******************************+++*******************************/
int LevelLex(size_t n, const size_t *nLevels, size_t *Level)
/*****************************************************************/
/*   Purpose:  Generate the next level combination of the n      */
/*             levels in Level (lexicographic order), where      */
/*             variable j has levels 0,..., nLevels[j].          */
/*                                                               */
/*   Returns:  ALL_DONE if no more combinations remain;          */
/*             OK       otherwise.                               */
/*                                                               */
/*   Version:  1996.03.26                                        */
/*****************************************************************/
{
     size_t    j;

     if (n == 0)
          return ALL_DONE;

     j = n - 1;

     /* Increment last level. */
     Level[j]++;

     while (Level[j] >= nLevels[j] && j > 0)
     {
         Level[j] = 0;
         Level[j-1]++;
         j--;
     }

     return (Level[0] >= nLevels[0]) ? ALL_DONE : OK;
}

/*******************************+++*******************************/
size_t Combinations(size_t n, size_t m)
/*****************************************************************/
/*   Purpose:  Return n choose m.                                */
/*                                                               */
/*   Version:  1996.01.28                                        */
/*****************************************************************/
{
     size_t    j, nCombs;

     CodeCheck(m <= n);

     nCombs = 1;
     m = min(m, n - m);
     for (j = 0; j < m; j++)
          /* Can't use *= here, because of integer divide. */
          nCombs = nCombs * (n - j) / (j + 1);

     return nCombs;
}

/*******************************+++*******************************/
real Pythag(real a, real b)
/*****************************************************************/
/*   Purpose:  Return sqrt(a^2 + b^2).                           */
/*                                                               */
/*   Version:  1999.02.19                                        */
/*****************************************************************/
{
     real absa, absb, p;

     absa = fabs(a);
     absb = fabs(b);
     
     if (absa > absb)
          p = absa * sqrt(1.0 + sqr(absb / absa));
     else if (absb == 0.0)
          p = 0.0;
     else
          p = absb * sqrt(1.0 + sqr(absa / absb));

     return p;     
}

