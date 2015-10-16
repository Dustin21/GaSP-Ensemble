/*****************************************************************/
/*                                                               */
/*   ROUTINES FOR COMPUTING A QR DECOMPOSITION                   */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--91.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

/*******************************+++*******************************/
/*                                                               */
/*   size_t    QRLS(Matrix *F, real *y, Matrix *Q, Matrix *R,    */
/*                  real *c, real *res)                          */
/*                                                               */
/*   Purpose:  Gram-Schmidt QR decomposition for least squares.  */
/*                                                               */
/*   Args:     F         Input:  Expanded design matrix.         */
/*             y         Input:  Response vector.                */
/*             Q, R      Output: QR decompostion of F.           */
/*             c         Output: R betahat = c can be solved     */
/*                                 for betahat.                  */
/*             res       Output: Least squares residuals.        */
/*                                                               */
/*   Comments: Calling with F = Q will overwrite F with Q.       */
/*             Calling with y = res will overwrite y with res.   */
/*             The calling routine must allocate space for Q, R, */
/*             c, and res.                                       */
/*                                                               */
/*   Returns:  j + 1 if R[j, j] becomes zero;                    */
/*             OK    otherwise.                                  */
/*                                                               */
/*   Version:  1991 July 11                                      */
/*                                                               */
/*****************************************************************/

size_t QRLS(Matrix *F, real *y, Matrix *Q, Matrix *R, real *c,
          real *res)
{
     real      r_jj, temp;
     real      *Q_j;
     size_t    i, j, k, n, p;

     n = Q->NumRows;
     p = Q->NumCols;

     if (F != Q)
          /* Copy F to Q. */
          for (j = 0; j < p; j++)
               VecCopy(MatCol(F, j), n, MatCol(Q, j));

     /* The Gram-Schmidt method iterates on Q a column at a    */
     /* time, producing at step j the jth column of Q, the jth */
     /* row of R, and the jth element of c.                    */

     for (j = 0; j < p; j++)
     {
          Q_j = MatCol(Q, j);

          /* Replace this with SNRM2. */
          r_jj = sqrt(VecSS(Q_j, n));
          MatPutElem(R, j, j, r_jj);
          if (r_jj <= 0.0)
               return j + 1;

          VecMultScalar(1.0 / r_jj, n, Q_j);

          for (k = j + 1; k < R->NumCols; k++)
               MatPutElem(R, j, k, DotProd(Q_j, MatCol(Q, k), n));

          c[j] = DotProd(y, Q_j, n);

          /* Update the remaining columns of the Q matrix. */
          for (k = j + 1; k < p; k++)
               VecAddVec(-MatElem(R, j, k), Q_j, n, MatCol(Q, k));
     }

     /* Compute the residuals from the LS fit. */
     /* (Is this the right method?)            */
     for (i = 0; i < n; i++)
     {
          for (temp = 0.0, j = 0; j < Q->NumCols; j++)
               temp += MatElem(Q, i, j) * c[j];
          res[i] = y[i] - temp;
     }

     return OK;
}

