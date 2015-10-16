/*****************************************************************/
/*   BASIC LINEAR ALGEBRA (SUB)ROUTINES                          */
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

/*******************************+++*******************************/
real VecDotProd(size_t n, const real *a, const real *b)
/*****************************************************************/
/* Purpose:    Return a[0] * b[0] + ... + a[n-1] * b[n-1].       */
/*                                                               */
/* Comment:    Use VecDotProdRange for dot product of elements i */
/*             through j.                                        */
/*                                                               */
/* 1991.05.17: Created (as DotProd).                             */
/* 1999.03.03: DotProd replaced by VecDotProd with reordering    */
/*             of parameters.                                    */
/*                                                               */
/* Version:    1999.03.03                                        */
/*****************************************************************/
{
     real      d;
     size_t    i;

     for (d = 0.0, i = 0; i < n; i++)
          d += a[i] * b[i];

     return d;
}

/*******************************+++*******************************/
real VecDotProdRange(size_t i, size_t j, const real *a,
     const real *b)
/*****************************************************************/
/* Purpose:    Return a[i] * b[i] + ... + a[j] * b[j].           */
/*                                                               */
/* Comment:    Use VecDotProd for dot product of first n         */
/*             elements.                                         */
/*                                                               */
/* 1999.03.03: Created.                                          */
/*                                                               */
/* Version:    1999.03.03                                        */
/*****************************************************************/
{
     real      d;
     size_t    k;

     for (d = 0.0, k = i; k <= j; k++)
          d += a[k] * b[k];

     return d;
}

/*******************************+++*******************************/
/*                                                               */
/*   real      VecSS(const real *a, size_t n)                    */
/*                                                               */
/*   Purpose:  Return a[0] * a[0] + ... + a[n-1] * a[n-1].       */
/*                                                               */
/*   Comment:  If sqrt(sum of squares) required, then Norm2      */
/*             may avoid under/overflow.                         */
/*                                                               */
/*   Version:  1991 May 17                                       */
/*                                                               */
/*****************************************************************/

real VecSS(const real *a, size_t n)
{
     real      SS;
     size_t    i;

     for (SS = 0.0, i = 0; i < n; i++)
          SS += a[i] * a[i];

     return SS;
}

/*******************************+++*******************************/
/*                                                               */
/*   real      VecSSCor(const real *a, size_t n)                 */
/*                                                               */
/*   Purpose:  Return sum of squares of a[0],..., a[n-1],        */
/*             corrected for the mean.                           */
/*                                                               */
/*   Version:  1991 July 10                                      */
/*                                                               */
/*****************************************************************/

real VecSSCor(const real *a, size_t n)
{
     real      Mean, SS;
     size_t    i;

     Mean = VecSum(a, n) / n;

     for (SS = 0.0, i = 0; i < n; i++)
          SS += (a[i] - Mean) * (a[i] - Mean);

     return SS;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      VecInit(real s, size_t n, real *a)                */
/*                                                               */
/*   Purpose:  Initialize a[0], ..., a[n-1] to the scalar s.     */
/*                                                               */
/*   Version:  1991 July 9                                       */
/*                                                               */
/*****************************************************************/

void VecInit(real s, size_t n, real *a)
{
     size_t i;

     for (i = 0; i < n; i++)
          a[i] = s;
}

/*******************************+++*******************************/
void VecStrInit(const string s, size_t n, string *a)
/*****************************************************************/
/*   Purpose:  Initialize a[0], ..., a[n-1] to the string s.     */
/*                                                               */
/*   Version:  1994 February 16                                  */
/*****************************************************************/
{
     size_t i;

     for (i = 0; i < n; i++)
     {
          if (a[i] != NULL)
               AllocFree(a[i]);
          a[i] = StrDup(s);
     }
}

/*******************************+++*******************************/
real VecSum(const real *a, size_t n)
/*****************************************************************/
/*   Purpose:  Return a[0] + ... + a[n-1].                       */
/*                                                               */
/*   Version:  1991 May 17                                       */
/*****************************************************************/
{
     real   sum;
     size_t i;

     for (sum = 0.0, i = 0; i < n; i++)
          sum += a[i];

     return sum;
}

/*******************************+++*******************************/
real VecSumAbs(size_t n, const real *a)
/*****************************************************************/
/*   Purpose:  Return |a[0]| + ... + |a[n-1]|.                   */
/*                                                               */
/*   Version:  1995 February 16                                  */
/*****************************************************************/
{
     real   sum;
     size_t i;

     for (sum = 0.0, i = 0; i < n; i++)
          sum += fabs(a[i]);

     return sum;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      VecAddScalar(real s, size_t n, real *a)           */
/*                                                               */
/*   Purpose:  Add scalar s to a[0], ..., a[n-1].                */
/*                                                               */
/*   Version:  1991 May 17                                       */
/*                                                               */
/*****************************************************************/

void VecAddScalar(real s, size_t n, real *a)
{
     size_t i;

     for (i = 0; i < n; i++)
          a[i] += s;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      VecMultScalar(real s, size_t n, real *a)          */
/*                                                               */
/*   Purpose:  Multiply a[0], ..., a[n-1] by scalar s.           */
/*                                                               */
/*   Version:  1991 May 17                                       */
/*                                                               */
/*****************************************************************/

void VecMultScalar(real s, size_t n, real *a)
{
     size_t i;

     for (i = 0; i < n; i++)
          a[i] *= s;
}

/*******************************+++*******************************/
void VecAddVec(real s, const real *a, size_t n, real *b)
/*****************************************************************/
/*   Purpose:  Add s * a[i] to b[i].                             */
/*                                                               */
/*   Version:  1995 April 22                                     */
/*****************************************************************/
{
     size_t    i;

     if (s == 1.0)
          for (i = 0; i < n; i++)
               b[i] += a[i];

     else if (s == -1.0)
          for (i = 0; i < n; i++)
               b[i] -= a[i];

     else
          for (i = 0; i < n; i++)
               b[i] += s * a[i];

     return;
}

/*******************************+++*******************************/
void VecAddVecNew(size_t n, const real *a, const real *b, real *c)
/*****************************************************************/
/*   Purpose:  Compute c[i] = a[i] + b[i] for i = 1,..., n.      */
/*                                                               */
/*   Comment:  Rename this as VecAddVec, with the old VecAddVec  */
/*             becoming VecIncVec.                               */
/*                                                               */
/*   Version:  1992 May 29                                       */
/*****************************************************************/
{
     size_t    i;

     for (i = 0; i < n; i++)
          c[i] = a[i] + b[i];

     return;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      VecMultVec(const real *a, size_t n, real *b)      */
/*                                                               */
/*   Purpose:  Replace b[i] by a[i] * b[i].                      */
/*                                                               */
/*   Version:  1991 July 15                                      */
/*                                                               */
/*****************************************************************/

void VecMultVec(real *a, size_t n, real *b)
{
     size_t    i;

     for (i = 0; i < n; i++)
          b[i] *= a[i];

     return;
}

/*******************************+++*******************************/
/*                                                               */
/*   real      VecMax(const real *a, size_t n)                   */
/*                                                               */
/*   Purpose:  Return maximum of a[0],..., a[n-1].               */
/*                                                               */
/*   Version:  1991 May 22                                       */
/*                                                               */
/*****************************************************************/

real VecMax(const real *a, size_t n)
{
     real      m;
     size_t    i;

     m = (n > 0) ? a[0] : 0.0;

     for (i = 1; i < n; i++)
          m = max(a[i], m);

     return m;
}

/*******************************+++*******************************/
real VecMin(const real *a, size_t n)
/*****************************************************************/
/*   Purpose:  Return minimum of a[0],..., a[n-1].               */
/*                                                               */
/*   Version:  1991.05.22                                        */
/*****************************************************************/
{
     real      m;
     size_t    i;

     m = (n > 0) ? a[0] : 0.0;

     for (i = 1; i < n; i++)
          m = min(a[i], m);

     return m;
}

/*******************************+++*******************************/
boolean VecHasNA(size_t n, real *a)
/*****************************************************************/
/*   Purpose:  Does a[0],..., a[n-1] contain at least one NA?    */
/*                                                               */
/*   Version:  1996.04.03                                        */
/*****************************************************************/
{
     size_t    i;

     for (i = 0; i < n; i++)
          if (a[i] == NA_REAL)
               break;

     return (i < n);
}

/*******************************+++*******************************/
/*                                                               */
/*   void      GivRot(real *a, real *b, real *c, real *s)        */
/*                                                               */
/*   Purpose:  Construct givens rotation and apply it to *a and  */
/*             *b.                                               */
/*                                                               */
/*   Version:  1990 February 22                                  */
/*                                                               */
/*****************************************************************/

void GivRot(real *a, real *b, real *c, real *s)
{
     real r, rho, scale, z;

     if (fabs(*a) > fabs(*b))
          rho = *a;
     else
          rho = *b;

     scale = fabs(*a) + fabs(*b);
     if (scale == 0.0)
     {
          *c = 1.0;
          *s = 0.0;
          r  = 0.0;
     }
     else
     {
          r = scale * sqrt( sqr(*a / scale) + sqr(*b / scale) );
          if (rho < 0.0)
               r = -r;
          *c = *a / r;
          *s = *b / r;
     }
     z = 1.0;
     if (fabs(*a) > fabs(*b))
          z = *s;
     if (fabs(*b) >= fabs(*a) && *c != 0.0)
          z = 1.0 / *c;
     *a = r;
     *b = z;
     return;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      VecIntCopy(const int *a, size_t n, int *b)        */
/*   void      VecCopy(const real *a, size_t n, real *b)         */
/*   void      VecSize_tCopy(const size_t *a, size_t n,          */
/*                  size_t *b)                                   */
/*   void      VecStrCopy(const string *a, size_t n, string *b)  */
/*                                                               */
/*   Purpose:  Copy a[i] to b[i] for i = 0,..., n - 1.           */
/*                                                               */
/*   Version:  1994 November 18                                  */
/*                                                               */
/*****************************************************************/

void VecIntCopy(const int *a, size_t n, int *b)
{
     size_t    i;

     for (i = 0; i < n; i++)
          b[i] = a[i];

     return;
}

void VecCopy(const real *a, size_t n, real *b)
{
     size_t    i;

     for (i = 0; i < n; i++)
          b[i] = a[i];

     return;
}

void VecSize_tCopy(const size_t *a, size_t n, size_t *b)
{
     size_t    i;

     for (i = 0; i < n; i++)
          b[i] = a[i];

     return;
}

void VecStrCopy(const string *a, size_t n, string *b)
{
     size_t    i;

     if (a == NULL || b == NULL)
          return;

     for (i = 0; i < n; i++)
     {
          if (b[i] != NULL)
               AllocFree(b[i]);

          b[i] = StrDup(a[i]);
     }

     return;
}

/*******************************+++*******************************/
void VecCopyStride(size_t n, size_t aStride, const real *a,
          size_t bStride, real *b)
/*****************************************************************/
/*   Purpose:  Copy a[i*aStride] to b[i*bStride] for             */
/*             i = 0,..., n - 1.                                 */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 18                                   */
/*****************************************************************/
{
     size_t    aIndex, bIndex, i;

     CodeCheck(aStride != 0 && bStride != 0);

     if (aStride == 1 && bStride == 1)
          for (i = 0; i < n; i++)
               b[i] = a[i];

     else
          for (aIndex = 0, bIndex = 0, i = 0; i < n;
                    i++, aIndex += aStride, bIndex += bStride)
               b[bIndex] = a[aIndex];

     return;
}

/*******************************+++*******************************/
void VecCopyIndex(size_t n, const size_t *aIndex, const real *a,
          const size_t *bIndex, real *b)
/*****************************************************************/
/*   Purpose:  Copy a[aIndex[i]] to b[b[Index[i]] for            */
/*             i = 0,..., n - 1.  If aIndex is NULL, then        */
/*             aIndex[i] = i is assumed; similarly bIndex.       */
/*                                                               */
/*   Version:  1992 May 27                                       */
/*****************************************************************/
{
     size_t    i, j, k;

     for (i = 0; i < n; i++)
     {
          j = (aIndex != NULL) ? aIndex[i] : i;
          k = (bIndex != NULL) ? bIndex[i] : i;
          b[k] = a[j];
     }

     return;
}

/*******************************+++*******************************/
/*                                                               */
/*   string    *VecIntStr(const int *a, size_t n, string *s)     */
/*   string    *VecRealStr(const real *a, size_t n, string *s)   */
/*   string    *VecSize_tStr(const size_t *a, size_t n,          */
/*                 string *s)                                    */
/*                                                               */
/*   Purpose:  Convert a[i] to string s[i] for i = 0,..., n - 1. */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 18                                   */
/*****************************************************************/

string *VecIntStr(const int *a, size_t n, string *s)
{
     size_t    i;

     if (n == 0)
          return NULL;

     CodeCheck(a != NULL && s != NULL);

     for (i = 0; i < n; i++)
         s[i] = StrDup(StrFromInt(a[i]));

     return s;
}

string *VecRealStr(const real *a, size_t n, string *s)
{
     size_t    i;

     if (n == 0)
          return NULL;

     CodeCheck(a != NULL && s != NULL);

     for (i = 0; i < n; i++)
         /* Default precision. */
         s[i] = StrDup(StrFromReal(a[i], "", PRECISION, 'g'));

     return s;
}

string *VecSize_tStr(const size_t *a, size_t n, string *s)
{
     size_t    i;

     if (n == 0)
          return NULL;

     CodeCheck(a != NULL && s != NULL);

     for (i = 0; i < n; i++)
         s[i] = StrDup(StrFromSize_t(a[i]));

     return s;
}

/*******************************+++*******************************/
/*                                                               */
/*   size_t    VecStrInt(const string *s, size_t n, int *a)      */
/*   size_t    VecStrReal(const string *s, size_t n, real *a)    */
/*   size_t    VecStrSize_t(const string *s, size_t n,           */
/*                  size_t *a)                                   */
/*                                                               */
/*   Purpose:  Convert string s[i] to a[i] for i = 0,..., n - 1. */
/*                                                               */
/*   Returns:  i        if s[i] cannot be converted;             */
/*             INDEX_OK otherwise.                               */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 18                                   */
/*****************************************************************/

size_t VecStrInt(const string *s, size_t n, int *a)
{
     size_t    i;

     if (n == 0)
          return INDEX_OK;

     CodeCheck(s != NULL && a != NULL);

     for (i = 0; i < n; i++)
         if (StrToInt(s[i], &a[i]) != OK)
               return i;

     return INDEX_OK;
}

size_t VecStrReal(const string *s, size_t n, real *a)
{
     size_t    i;

     if (n == 0)
          return INDEX_OK;

     CodeCheck(s != NULL && a != NULL);

     for (i = 0; i < n; i++)
         if (StrToReal(s[i], &a[i]) != OK)
               return i;

     return INDEX_OK;
}

size_t VecStrSize_t(const string *s, size_t n, size_t *a)
{
     size_t    i;

     if (n == 0)
          return INDEX_OK;

     CodeCheck(s != NULL && a != NULL);

     for (i = 0; i < n; i++)
         if (StrToSize_t(s[i], &a[i]) != OK)
               return i;

     return INDEX_OK;
}

/*******************************+++*******************************/
size_t VecSize_tIndex(size_t Target, size_t n, const size_t *a)
/*****************************************************************/
/*   Purpose:  Find first i such that Target matches a[i].       */
/*                                                               */
/*   Returns:  i          if a[i] matches Target;                */
/*             INDEX_ERR  otherwise.                             */
/*                                                               */
/*   Comments: Combine with StrIndex in libstr.c, etc.           */
/*                                                               */
/*   Version:  1994 February 2                                   */
/*****************************************************************/
{
     size_t    i;

     for (i = 0; i < n; i++)
          if (a[i] == Target)
               return i;

     return INDEX_ERR;
}
