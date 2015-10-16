/*****************************************************************/
/*   ROUTINES THAT MANUIPULATE (UPPER) TRIANGULAR MATRICES:      */
/*   (1) SOLVING TRIANGULAR SYSTEMS                              */
/*   (2) COMPUTING AND UPDATING THE CHOLESKY FACTOR R AS IN      */
/*   R'R = X'X OR X = QR.                                        */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--96.                    */
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
int TriForSolve(const Matrix *R, const real *b, size_t StartOff,
          real *x)
/*****************************************************************/
/*   Purpose:  Forward solve triangular system of equations      */
/*             R'x = b for x.                                    */
/*             Calling with TriForSolve(R, b, StartOff, b) will  */
/*             overwrite b with x.                               */
/*             Solving starts at element with offset StartOff    */
/*             (useful if new columns of R have been appended).  */
/*                                                               */
/*   Returns:  NONUNIQ_ERR if more than one solution is possible */
/*                         (zero is taken);                      */
/*             NUMERIC_ERR if equations cannot be solved because */
/*                         of a zero diagonal element of R;      */
/*             OK          otherwise.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for x         */
/*             (unless x = b).                                   */
/*                                                               */
/*   Version:  1991 May 20                                       */
/*****************************************************************/
{
     real      numer, *Rj, Rjj;
     int       ErrNum;
     size_t    j, n, nzero;

     n = MatNumCols(R);

     /* Leading zeros in b simplify calculations. */
     nzero = 0;
     while (nzero < n && b[nzero] == 0.0)
     {
          if (nzero >= StartOff)
                x[nzero] = 0.0;
          nzero++;
     }

     ErrNum = OK;
     for (j = max(nzero, StartOff); j < n && ErrNum != NUMERIC_ERR;
               j++)
     {
          Rj  = MatCol(R, j);
          Rjj = Rj[j];
          numer = b[j] - DotProd(Rj + nzero, x + nzero, j - nzero);
          if (Rjj != 0.0)
               x[j] = numer / Rjj;
          else if (numer == 0.0)
          {
               /* 0 * x = 0; take solution x = 0. */
               x[j] = 0.0;
               ErrNum = NONUNIQ_ERR;
          }
          else
               /* 0 * x = nonzero. */
               ErrNum = NUMERIC_ERR;
     }

     return ErrNum;
}

/*******************************+++*******************************/
/*                                                               */
/*   int       TriBackSolve(const Matrix *R, const real *b,      */
/*                  real *x)                                     */
/*                                                               */
/*   Purpose:  Backward solve the triangular system of equations */
/*             R x = b for x.                                    */
/*             Calling with TriBackSolve(R, b, b) will overwrite */
/*             b with x.                                         */
/*                                                               */
/*   Returns:  NONUNIQUE   if more than one solution is possible */
/*                         (zero is taken);                      */
/*             NUMERIC_ERR if equations cannot be solved because */
/*                         of a zero diagonal element of R;      */
/*             OK          otherwise.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for x         */
/*             (unless x = b).                                   */
/*                                                               */
/*   Version:  1991 May 20                                       */
/*                                                               */
/*****************************************************************/

int TriBackSolve(const Matrix *R, const real *b, real *x)
{
     int       ErrNum;
     real      Diag;
     size_t    j, jj, n;

     n = MatNumCols(R);

     if (x != b)
          for (j = 0; j < n; j++)
               x[j] = b[j];

     ErrNum = OK;
     for (j = 0; j < n && ErrNum != NUMERIC_ERR; j++)
     {
          jj = n - 1 - j;   /* Backwards solve. */
          Diag = MatElem(R, jj, jj);
          if (j > 0)
               VecAddVec(-x[jj+1], MatCol(R, jj + 1), jj + 1, x);
          if (Diag != 0.0)
               x[jj] /= Diag;
          else if (x[jj] == 0.0)
          {
               /* 0 * x = 0; take solution x = 0. */
               x[jj] = 0.0;
               ErrNum = NONUNIQ_ERR;
          }
          else
               /* 0 * x = nonzero */
               ErrNum = NUMERIC_ERR;
     }

     return ErrNum;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      TriDet(const Matrix *R, real *d1, int *d2)        */
/*                                                               */
/*   Purpose:  Compute determinant of triangular matrix R in the */
/*             form *d1 * 10 ** (*d2), where                     */
/*                  0.1 <= *d1 < 1 or *d1 = 0.                   */
/*                                                               */
/*   Version:  1991 May 20                                       */
/*                                                               */
/*****************************************************************/

void TriDet(const Matrix *R, real *d1, int *d2)
{
     int       p;
     real      d;
     size_t    j, n;

     n = MatNumCols(R);
     if (n == 0)
     {
         *d1 = 0.0;
         *d2 = 0;
         return;
     }

     d = 0.1;
     p = 1;
     for (j = 0; j < n && d != 0.0; j++)
     {
          d *= MatElem(R, j, j);
          while (fabs(d) >= 1.0)
          {
               d /= 10.0;
               p++;
          }
          if (fabs(d) > 0.0)
               while (fabs(d) < 0.1)
               {
                    d *= 10.0;
                    p--;
               }
          else
               p = 0;
     }
     *d1 = d;
     *d2 = p;
     return;
}

/*******************************+++*******************************/
real TriCond(const Matrix *R)
/*****************************************************************/
/* Purpose:    Return estimate of condition number of R.         */
/*                                                               */
/* Comment:    From LINPACK, STRCO, p. C.86.                     */
/*                                                               */
/* 1995.04.10: Created                                           */
/* 1999.06.17: Some braces added to avoid ambiguous else.        */
/*****************************************************************/
{
     int  j, k, kk, n;
     real e, Rcond, Rnorm, Rsave, ynorm, s, sm;
     real w, wk, wkm, ek ;
     real *z;

     if (MatEmpty(R))
          return 1.0;

     n = MatNumCols(R);

     z = AllocReal(n, NULL);       /* Used for workspace. */

     /* Compute the 1 norm of R */
     Rnorm = 0.0;
     for (j = 0; j < n; j++)
     {
          Rsave = VecSumAbs(j + 1, MatCol(R, j));
          if (Rnorm < Rsave)
               Rnorm = Rsave;
     }

     /* Rcond = (norm(R) * (estimate of norm(inverse(R)))) */
     ek = 1.0;
     VecInit(0.0, n, z);

     for (k = 0; k < n; k++)
     {
          if (z[k] != 0.0)
          {
               if (-z[k] < 0.0)
                    ek = -fabs(ek);
               else
                    ek = fabs(ek);
          }          

          if (fabs(ek - z[k]) >  fabs(MatElem(R, k, k)) )
          {
               /* Rescale z and adjust ek accordingly. */
               s = fabs(MatElem(R, k, k)) / fabs(ek-z[k]);
               VecMultScalar(s, n, z);
               ek = s * ek;
          }

          /* Set values for wk, wkm, s, and sm. */
          wk  =  ek - z[k];
          wkm = -ek - z[k];
          s  = fabs(wk);
          sm = fabs(wkm);
          if ( (e = MatElem(R, k, k)) != 0.0)
          {
              wk  = wk  / e;
              wkm = wkm / e;
          }
          else
          {
              wk  = 1.0;
              wkm = 1.0;
          }

          /* core steps of forward elimination */
          if (k < n - 1)
          {
               for (j = k + 1; j < n; j++)
               {
                    sm += fabs(z[j] + wkm * MatElem(R, k, j));
                    z[j] += wk * MatElem(R, k, j);
                    s += fabs(z[j]);
               }
               if (s < sm)
               {
                    w  = wkm - wk;
                    wk = wkm;
                    for (j = k + 1; j < n; j++)
                         z[j] += w * MatElem(R, k, j);
               }
          }
          z[k] = wk;
     } /* end of k loop */

     /* Scale z such that norm(z) = 1, i.e. ynorm = 1.  */
     /* Note that z contains the y-vector in R * z = y. */
     s = VecSumAbs(n, z);
     VecMultScalar(1.0 / s, n, z);
     ynorm = 1.0;

     /* Solve R * z = y. */
     for (kk = 0; kk < n; kk++)
     {
          k = n - 1 - kk;
          e = MatElem(R, k, k);
          if (fabs(z[k]) > fabs(e))
          {
               /* Scale z and adjust ynorm accordingly. */
               s = fabs(e) / fabs(z[k]);
               VecMultScalar(s, n, z);
               ynorm *= s;
          }
          /* Start backward elimination for the kth row. */
          if (e != 0.0)
               z[k] = z[k] / e;
          if (e == 0.0)
               z[k] = 1.0;
          if (kk < n - 1)
               VecAddVec(-z[k], MatCol(R, k), k, z);
     }

     /* Scale z such that znorm = 1 and adjust ynorm, */
     /* i.e. ynorm/znorm = ynorm.                     */
     s = VecSumAbs(n, z);
     VecMultScalar(1.0 / s, n, z);
     ynorm = ynorm / s;

     if (Rnorm != 0.0)
          Rcond = Rnorm / ynorm;
     else
          Rcond = REAL_MAX;

     AllocFree(z);

     return Rcond;
}

/*******************************+++*******************************/
size_t TriCholesky(const Matrix *S, size_t FirstOff, Matrix *R)
/*****************************************************************/
/*   Purpose:  Cholesky decomposition R'R of symmetric matrix S. */
/*             The decomposition starts at the column with       */
/*             offset FirstOff (useful if appending new columns).*/
/*                                                               */
/*   Returns:  Rank      if the rank (number of nonzero          */
/*                       diagonals of R) is less than the order  */
/*                       of R;                                   */
/*             OK        otherwise.                              */
/*                                                               */
/*   Comment:  See SPOFA, LINPACK, p. C.28.                      */
/*             The calling routine must allocate space for R.    */
/*             Adapted to complete decomposition if rank         */
/*             deficient.                                        */
/*             Calling with TriCholesky(S,..., S) will overwrite */
/*             S.                                                */
/*                                                               */
/*   96.02.20: Adapted to continue if zero diagonal encountered. */
/*                                                               */
/*   Version:  1996.02.20                                        */
/*****************************************************************/
{
     real      *Sj, *Ri, *Rj;
     real      r, t;
     size_t    i, j, n, Rank;

     n = MatNumCols(S);

     R->Shape = UP_TRIANG;

     for (j = FirstOff; j < n; j++)
     {
          /* Pointers to column j of S and R. */
          Sj = MatCol(S, j);
          Rj = MatCol(R, j);

          for (r = 0.0, i = 0; i < j; i++)
          {
               /* Pointer to column i of R. */
               Ri = MatCol(R, i);

               if (Ri[i] > 0.0)
                    t = (Sj[i] - DotProd(Rj, Ri, i)) / Ri[i];
               else
                    t = 0.0;

               Rj[i] = t;
               r += t * t;
          }

          r = Sj[j] - r;

          if (r > 0.0)
               Rj[j] = sqrt(r);
          else
               Rj[j] = 0.0;
     }

     for (Rank = 0, j = 0; j < n; j++)
          if (MatElem(R, j, j) > 0.0)
               Rank++;

     return (Rank == n) ? OK : Rank;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      TriRect(const Matrix *X, Matrix *R)               */
/*                                                               */
/*   Purpose:  Compute R in the decomposition X = QR.            */
/*                                                               */
/*   Comment:  The calling routine must allocate space for R.    */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*                                                               */
/*****************************************************************/

void TriRect(const Matrix *X, Matrix *R)
{
     real      *c, *Rj, *row, *s;
     size_t    i, j, NumCols, NumRows;

     NumCols = MatNumCols(X);
     NumRows = MatNumRows(X);

     /* Allocate space for the sines and cosines used in */
     /* TriUpdate, and to hold a row of X.              */
     c   = AllocReal(NumCols, NULL);
     s   = AllocReal(NumCols, NULL);
     row = AllocReal(NumCols, NULL);

     /* Initialize R. */
     for (j = 0; j < NumCols; j++)
          for (Rj = MatCol(R, j), i = 0; i <= j; i++)
               Rj[i] = 0.0;

     /* Update R for each row of X. */
     for (i = 0; i < NumRows; i++)
     {
          for (j = 0; j < NumCols; j++)
               row[j] = MatElem(X, i, j);
          TriUpdate(row, 1.0, R, c, s);
     }

     AllocFree(c);
     AllocFree(s);
     AllocFree(row);

     return;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      TriUpdate(const real *xrow, real wt, Matrix *R,   */
/*                  real *c, real *s)                            */
/*                                                               */
/*   Purpose:  Update Cholesky R of X'X when a new row xrow is   */
/*             appended with non-negative weight wt.             */
/*             On return c and s contain the cosines and sines   */
/*             of the transforming rotations.                    */
/*                                                               */
/*   Comment:  Calling routine must allocate space for c and s.  */
/*             See SCHUD, LINPACK, p. C.112.                     */
/*                                                               */
/*   Version:  1991 May 20                                       */
/*                                                               */
/*****************************************************************/

void TriUpdate(const real *xrow, real wt, Matrix *R, real *c,
          real *s)
{
     size_t    i, j, n;
     real      *Rj;
     real      sqrtwt, t, xj;

     sqrtwt = sqrt(wt);

     n  = MatNumCols(R);
     for (j = 0; j < n; j++)
     {
          xj = sqrtwt * xrow[j];

          /* Pointer to column j of R. */
          Rj = MatCol(R, j);

          /* Apply the previous rotations. */
          for (i = 0; i < j; i++)
          {
               t  = c[i] * Rj[i] + s[i] * xj;
               xj = c[i] * xj     - s[i] * Rj[i];
               Rj[i] = t;
          }

          /* Compute next rotation - also computes jth diagonal. */
          GivRot(&Rj[j], &xj, &c[j], &s[j]);
          /* Fix negative diagonal. */
          if (Rj[j] < 0.0)
          {
               Rj[j] = -Rj[j];
               c[j] = -c[j];
               s[j] = -s[j];
          }
     }

     return;
}
/*******************************+++*******************************/
/*                                                               */
/*   int       TriDownDate(const real *xrow, real wt, Matrix *R, */
/*                  real *c, real *s)                            */
/*                                                               */
/*   Purpose:  Downdate Cholesky R of X'X when row xrow is       */
/*             deleted with non-negative weight wt.              */
/*             On return c and s contain the cosines and sines   */
/*             of the transforming rotations.                    */
/*                                                               */
/*   Returns:  OK          if downdating is carried out;         */
/*             NUMERIC_ERR otherwise.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for c and s.  */
/*             See SCHDD, LINPACK, p. C.116.                     */
/*                                                               */
/*   Version:  1991 July 11                                      */
/*                                                               */
/*****************************************************************/

int TriDownDate(const real *xrow, real wt, Matrix *R, real *c, real *s)
{
     real      a, alpha, b, norm, scale, sqrtwt, t, xx;
     real      *Rj = NULL;
     size_t    i, ii, j, jj, n;

     n      = MatNumCols(R);
     sqrtwt = sqrt(wt);

     /* Solve R' s = sqrtwt * x, using s temporarily for soln. */
     for (j = 0; j < n; j++)
     {
          /* Pointer to column j of R (row j of L = R'). */
          Rj = MatCol(R, j);

          if (Rj[j] == 0.0)
               return NUMERIC_ERR;
          else
               s[j] = (sqrtwt * xrow[j] - DotProd(Rj, s, j))
                         / Rj[j];
     }

     /* Replace this with Norm2 to avoid under/overflow. */
     norm = sqrt(VecSS(s, n));

     if (norm >= 1.0)
          return NUMERIC_ERR;

     alpha = sqrt(1.0 - norm * norm);

     /* Determine the transformations. */
     for (jj = 0; jj < n; jj++)
     {
           j = n - jj - 1;   /* Necessary because of unsigned. */
           scale = alpha + fabs(s[j]);
           a = alpha / scale;
           b = s[j]  / scale;
           norm = sqrt(a * a + b * b);
           c[j] = a / norm;
           s[j] = b / norm;
           alpha = scale * norm;
     }

     /* Apply the transformations to R. */
     for (j = 0; j < n; j++)
     {
          xx = 0.0;
          /* Pointer to column j of R. */
          Rj = MatCol(R, j);
          for (ii = 0; ii <= j; ii++)
          {
              i = j - ii;
              t = c[i] * xx + s[i] * Rj[i];
              Rj[i] = c[i] * Rj[i] - s[i] * xx;
              xx = t;
          }
     }

     return OK;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      TriPerm(size_t FirstOff, size_t LastOff,          */
/*                  Matrix *R, real *c, real *s)                 */
/*                                                               */
/*   Purpose:  Permute the columns of R so that columns (with    */
/*             offsets) FirstOff,..., LastOff are cyclically     */
/*             permuted to columns                               */
/*                  FirstOff + 1,..., LastOff, FirstOff.         */
/*             Then re-triangularize R by applying an orthogonal */
/*             matrix Q (a product of Givens rotations).         */
/*                                                               */
/*   Args:     FirstOff, Offsets defining the cyclic permutation */
/*             LastOff   (see above).                            */
/*             R         Input:  UP_TRIANG matrix;               */
/*                       Output: Re-triangularized UP_TRIANG     */
/*                               matrix.                         */
/*             c, s      Output: c[i] and s[i] are the sines     */
/*                               and cosines of the rotations    */
/*                               in the ii, ii + 1 plane,        */
/*                               where ii = FirstOff + i,        */
/*                               and i = 0,...,                  */
/*                               LastOff - FirstOff - 1.         */
/*                                                               */
/*   Comment:  Calling routine must allocate space for c and s   */
/*             (at least LastOff - FirstOff elements).           */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*                                                               */
/*****************************************************************/

void TriPerm(size_t FirstOff, size_t LastOff, Matrix *R, real *c,
          real *s)
{
     real      cc, ElimElem, ss, t;
     real      *Col, *TempCol;
     size_t    j, k;

     /* Permute columns of FirstOff,..., LastOff of R */
     /* (only first FirstOff + 1 elements).           */
     TempCol = AllocReal(FirstOff + 1, NULL);
     VecCopy(MatCol(R, FirstOff), FirstOff + 1, TempCol);
     for (j = FirstOff; j < LastOff; j++)
          VecCopy(MatCol(R, j + 1), FirstOff + 1, MatCol(R, j));
     VecCopy(TempCol, FirstOff + 1, MatCol(R, LastOff));
     AllocFree(TempCol);

     for (j = FirstOff; j < LastOff; j++)
     {
          Col = MatCol(R, j);

          /* The element still in row j + 1, column j + 1         */
          /* (which is permuted to column j) must be eliminated.  */
          ElimElem = MatElem(R, j + 1, j + 1);
          GivRot(&Col[j], &ElimElem, &cc, &ss);
          c[j-FirstOff] = cc;
          s[j-FirstOff] = ss;

          /* We can now permute the elements in row j + 1. */
          for (k = j + 1; k < LastOff; k++)
          {
               t = MatElem(R, j + 1, k + 1);
               MatPutElem(R, j + 1, k, t);
          }
          MatPutElem(R, j + 1, LastOff, 0.0);

          /* Apply the rotation to the remaining columns of R. */
          for (k = j + 1; k < MatNumCols(R); k++)
          {
               Col = MatCol(R, k);
               t        =  cc * Col[j] + ss * Col[j+1];
               Col[j+1] = -ss * Col[j] + cc * Col[j+1];
               Col[j]   = t;
          }
     }

     return;
}

/*******************************+++*******************************/
void TriVec(const Matrix *R, const real *b, real *c)
/*****************************************************************/
/*   Purpose:  Compute c = R * b.                                */
/*                                                               */
/*   Comment:  Calling routine must allocate space for c.        */
/*             TriVec(R, b, b) overwrites b with R * b.          */
/*                                                               */
/*   Version:  1996.01.15                                        */
/*****************************************************************/
{
     real      *Rj;
     size_t    j, n;

     n = MatNumCols(R);

     for (j = 0; j < n; j++)
     {
          Rj = MatCol(R, j);
          VecAddVec(b[j], Rj, j, c);
          c[j] = b[j] * Rj[j];
     }
}
