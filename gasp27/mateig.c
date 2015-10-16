/*****************************************************************/
/*   ROUTINES FOR EIGEN ANALYSIS OF SYMMETRIC MATRICES           */
/*                                                               */
/*   Copyright (c) William J. Welch 1999.                        */
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
int MatEig(boolean SortValues, matrix *S, real *eVal, matrix *V)
/*****************************************************************/
/* Purpose:    Compute eigenvalue decomposition, V Diag(e) V',   */
/*             of symmetric matrix S, where the columns of V are */
/*             the eigenvectors and eVal contains the            */
/*             eigenvalues.                                      */
/*             If SortValues is TRUE, then the eigenvalues (and  */
/*             corresponding eigenvectors) are sorted largest to */
/*             smallest.                                         */
/*                                                               */
/* Returns:    NUMERIC_ERR if MatEigTriDiag fails to converge;   */
/*             OK          otherwise.                            */
/*                                                               */
/* Comments:   Calling routine must allocate space for eVal.     */
/*             V must be allocated RECT by the calling routine.  */
/*             If S and V are distinct, then S is unaltered.     */
/*             If V = S, then S is overwritten; in this case S   */
/*             must be allocated RECT.                           */
/*                                                               */
/* 1999.03.13: Created.                                          */
/*                                                               */
/* Version:    1999.03.13                                        */
/*****************************************************************/
{
     int       ErrNum;
     size_t    n;
     real      *w;
 
     CodeCheck(MatType(V) == REAL);
     CodeCheck(MatShape(V) == RECT);
     n = MatNumRows(V);
     CodeCheck(n == MatNumCols(V))

     if (S != V)
     {
          CodeCheck(MatType(S) == REAL);
          CodeCheck(MatShape(S) == SYM || MatShape(S) == RECT);
          CodeCheck(n == MatNumRows(S) && n == MatNumCols(S))
     }

     /* Workspace for subdiagonals of tridiagonal matrix. */
     w = AllocReal(n, NULL);

     MatTriDiag(S, eVal, w, V);
     ErrNum = MatEigTriDiag(SortValues, eVal, w, V);

     AllocFree(w);

     return ErrNum;
}

/*******************************+++*******************************/
void MatTriDiag(matrix *S, real *d, real *e, matrix *Z)
/*****************************************************************/
/* Purpose:    Reduce a symmetric matrix S to Z T Z', where T is */
/*             a symmetric tridiagonal matrix, and Z is the      */
/*             orthogonal transformation matrix.                 */
/*             On exit, d contains the diagonals of T, and e     */
/*             contains the subdiagonals (with e[0] = 0.0).      */ 
/*                                                               */
/* Comments:   Calling routine must allocate space for d and e.  */
/*             Z must be allocated RECT by the calling routine.  */
/*             If S and Z are distinct, then S is unaltered.     */
/*             If Z = S, then S is overwritten; in this case S   */
/*             must be allocated RECT.                           */
/*                                                               */
/*             Translation of tred2, Handbook for Automatic      */
/*             Computation, Volume II, pp. 212--226.             */
/*                                                               */
/* 1999.03.20: Created.                                          */
/*                                                               */
/* Version:    1999.03.20                                        */
/*****************************************************************/
{
     real     f, g, h, scale;
     real     *c;
     size_t   i, j, k, n;

     CodeCheck(MatType(Z) == REAL);
     CodeCheck(MatShape(Z) == RECT);
     n = MatNumRows(Z);
     CodeCheck(n == MatNumCols(Z))

     if (S != Z)
     {
          CodeCheck(MatType(S) == REAL);
          CodeCheck(MatShape(S) == SYM || MatShape(S) == RECT);
          CodeCheck(n == MatNumRows(S) && n == MatNumCols(S))
     }

     if (n == 0)
          return;

     /* Copy *upper* triangle of S to *lower* triangle of Z; */
     /* i.e., copy columns of S to rows of Z.                */
     for (j = 0; j < n; j++)
          for (c = MatCol(S, j), i = 0; i <= j; i++)
               MatPutElem(Z, j, i, c[i]);

     /* Copy last column (row) of S to d. */
     VecCopy(MatCol(S, n - 1), n, d);

     for (i = n - 1; i > 0; i--)
     {
          h = scale = 0.0;

          if (i > 1)
               /* Scale row using d[0],..., d[i-1]. */
               scale = VecSumAbs(i, d);

          if (i == 1 || scale == 0.0)
          {
               e[i] = d[i-1];

               for (j = 0; j < i; j++)
               {
                    d[j] = MatElem(Z, i - 1, j);
                    MatPutElem(Z, i, j, 0.0);
                    MatPutElem(Z, j, i, 0.0);
               }
          }
          else
          {
               VecMultScalar(1.0 / scale, i, d);
               h += VecSS(d, i);
 
               f = d[i-1];
               g = ((f > 0.0) ? -sqrt(h) : sqrt(h));

               e[i] = scale * g;
               h -= f * g;
               d[i-1] = f - g;

               /* Form S * u. */
               VecInit(0.0, i, e);
               for (j = 0; j < i; j++)
               {
                    f = d[j];
                    MatPutElem(Z, j, i, f);
                    g = e[j] + MatElem(Z, j, j) * f;
               
                    c = MatCol(Z, j);
                    g += VecDotProdRange(j + 1, i - 1, c, d);
                    VecAddVecRange(j + 1, i - 1, f, c, e);

                    e[j] = g;
               }

               /* Form p. */
               VecMultScalar(1.0 / h, i, e);
               f = VecDotProd(i, d, e);

               /* Form q. */
               VecAddVec(-f / (2.0 * h), d, i, e);

               /* Form reduced S. */
               for (j = 0; j < i; j++)
               {
                    f = d[j];
                    g = e[j];
                    c = MatCol(Z, j);
                    for (k = j; k < i; k++)
                         c[k] -= f * e[k] + g * d[k];
                    d[j] = c[i-1];
                    c[i] = 0.0;
               }
          }
          d[i] = h;
     }

     /* Accumulation of transformation matrices. */
     for (i = 1; i < n; i++)
     {
          c = MatCol(Z, i - 1);
          c[n-1] = c[i-1];
          c[i-1] = 1.0;
          h = d[i];

          c = MatCol(Z, i);
          
          if (h != 0.0)
          {
               for (k = 0; k < i; k++)
                    d[k] = c[k] / h;

               for (j = 0; j < i; j++)
               {
                    g = VecDotProd(i, c, MatCol(Z, j));
                    VecAddVec(-g, d, i, MatCol(Z, j));
               }
          }

          VecInit(0.0, i, c);
     }

     for (i = 0; i < n; i++)
     {
          d[i] = MatElem(Z, n - 1, i);
          MatPutElem(Z, n - 1, i, 0.0);
     }

     MatPutElem(Z, n - 1, n - 1, 1.0);
     e[0] = 0.0;
     
     return;
}

# define MAX_IT     30

/*******************************+++*******************************/
int MatEigTriDiag(boolean SortValues, real *d, real *e, matrix *Z)
/*****************************************************************/
/* Purpose:    Compute the eigenvalues and eigenvectors of the   */
/*             symmetric tridiagonal matrix with diagonal        */
/*             elements in d and subdiagonal elements in e (e[0] */
/*             is arbitrary).                                    */
/*             On exit, d contains the eigenvalues in ascending  */
/*             order, e has been destroyed, and the columns of Z */
/*             contain the eigenvectors.                         */
/*             If SortValues is TRUE, then the eigenvalues (and  */
/*             corresponding eigenvectors) are sorted largest to */
/*             smallest.                                         */
/*             Usually, on entry d, e, and Z are from            */
/*             MatTriDiag; this gives the eigen decomposition of */
/*             a symmetric matrix.  To obtain the eigenvectors   */
/*             of a tridiagonal matrix, Z should be the identity */
/*             on entry.                                         */
/*                                                               */
/* Returns:    NUMERIC_ERR if the algorithm failed to converge;  */
/*             OK          otherwise.                            */
/*                                                               */
/* Comments:   Z must be allocated RECT by the calling routine.  */
/*                                                               */
/*             Translation of imtql2, Handbook for Automatic     */
/*             Computation, Volume II, 241--248.                 */
/*                                                               */
/* 1999.03.30: Created.                                          */
/*                                                               */
/* Version:    1999.03.30                                        */
/*****************************************************************/
{
     boolean   Converged, Underflow;
     real      b, c, dd, f, g, p, r, s;
     real      *Zi, *Zi1;
     real      **ColPtr;
     size_t    i, ii, iter, j, k, m, n;
     size_t    *Index;
   
     CodeCheck(MatType(Z) == REAL);
     CodeCheck(MatShape(Z) == RECT);
     n = MatNumRows(Z);
     CodeCheck(n == MatNumCols(Z))

     if (n <= 1)
          return OK;     

     for (i = 1; i < n; i++)
          e[i-1] = e[i];
     e[n-1] = 0.0;

     Converged = TRUE;
     for (j = 0; j < n && Converged; j++)
     {
          iter = 0;
          Converged = FALSE;
          do 
          {
               /* Look for a small subdiagonal element */
               /* to split the matrix.                 */     
               for (m = j; m < n - 1; m++)
               {
                    dd = fabs(d[m]) + fabs(d[m+1]);
                    if (fabs(e[m]) + dd == dd)
                         break;
               }

               if (m == j)
                    Converged = TRUE;   /* Next j. */

               else if (iter++ >= MAX_IT)
                    break;

               else
               {
                    /* Form shift. */
                    g = (d[j+1] - d[j]) / (2.0 * e[j]);
                    r = Pythag(g, 1.0);
                    g = d[m] - d[j] + e[j] / (g + sign(r, g));
                    s = c = 1.0;
                    p = 0.0;

                    /* i from m - 1 to j (but i is size_t). */
                    for (Underflow = FALSE, ii = 1; ii <= m - j; ii++)
                    {
                         i = m - ii;
                         f = s * e[i];
                         b = c * e[i];
                         e[i+1] = r = Pythag(f, g);

                         if (r == 0.0)
                         {                         
                              /* Recover from underflow. */
                              d[i+1] -= p;
                              e[m] = 0.0;
                              Underflow = TRUE;
                              break;    /* Next iteration. */
                         }
                         else
                         {
                              s = f / r;
                              c = g / r;
                              g = d[i+1] - p;
                              r = (d[i] - g) * s + 2.0 * c * b;
                              p = s * r;
                              d[i+1] = g + p;
                              g = c * r - b;

                              /* Form eigenvector. */
                              Zi = MatCol(Z, i);
                              Zi1 = MatCol(Z, i + 1);
                              for (k = 0; k < n; k++)
                              {
                                   f = Zi1[k];
                                   Zi1[k] = s * Zi[k] + c * f;
                                   Zi[k]  = c * Zi[k] - s * f;
                              }
                         }
                    }

                    if (!Underflow)
                    {
                         d[j] -= p;
                         e[j] = g;
                         e[m] = 0.0;
                    }
               }
          } while (!Converged); 
     }

     if (SortValues && Converged)
     {
          /* Index for sorting eigenvalues and eigenvectors. */
          Index = AllocSize_t(n, NULL);
          QuickIndex(d, n, Index);

          /* Reverse the sort; i.e., largest to smallest. */
          for (j = 0; j < n / 2; j++)
          {
               i = Index[j];
               Index[j] = Index[n-1-j];
               Index[n-1-j] = i;
          }

          /* Copy eigenvalues to e, then back to d in ascending order. */
          VecCopy(d, n, e);
          VecCopyIndex(n, Index, e, NULL, d);

          /* Copy eigenvector pointers, then back in ascending order. */
          ColPtr = AllocPtrReal(n, NULL);
          for (j = 0; j < n; j++)
               ColPtr[j] = MatCol(Z, j);
          for (j = 0; j < n; j++)
               MatPutCol(Z, j, ColPtr[Index[j]]);
               
          AllocFree(Index);
          AllocFree(ColPtr);
     }
     
     return Converged ? OK : NUMERIC_ERR;
}     

     
