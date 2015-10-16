/*****************************************************************/
/*   ROUTINES FOR REAL, SYMMETRIC MATRICES                       */
/*                                                               */
/*   Copyright (c) William J. Welch 1995.                        */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

/*******************************+++*******************************/
void MatSymCol(const Matrix *S, size_t ColIndex, real *col)
/*****************************************************************/
/*   Purpose:  Load column ColIndex of S into col.               */
/*                                                               */
/*   Comment:  Calling routine must allocate space for col.      */
/*                                                               */
/*   Version:  1995 April 22                                     */
/*****************************************************************/
{
     size_t    i;

     CodeCheck(MatType(S) == REAL);
     CodeCheck(MatShape(S) == SYM);

     for (i = 0; i <= ColIndex; i++)
          col[i] = MatElem(S, i, ColIndex);
     for (i = ColIndex + 1; i < MatNumRows(S); i++)
          col[i] = MatElem(S, ColIndex, i);

     return;
}

/*******************************+++*******************************/
void MatSymUpdate(real w, const real *v, Matrix *S)
/*****************************************************************/
/*   Purpose:  Add w * v * v' to S.                              */
/*                                                               */
/*   Version:  1995 April 22                                     */
/*****************************************************************/
{
     size_t    j;

     CodeCheck(MatType(S) == REAL);
     CodeCheck(MatShape(S) == SYM);

     for (j = 0; j < MatNumCols(S); j++)
          VecAddVec(w * v[j], v, j + 1, MatCol(S, j));

     return;
}

/*******************************+++*******************************/
real MatSymQuadForm(const real *v, const Matrix *S)
/*****************************************************************/
/*   Purpose:  Return a' S a.                                    */
/*                                                               */
/*   Version:  1995 April 22                                     */
/*****************************************************************/
{
     real      q;
     real      *c;
     size_t    j;

     CodeCheck(MatType(S) == REAL);
     CodeCheck(MatShape(S) == SYM);

     for (q = 0.0, j = 0; j < MatNumCols(S); j++)
     {
          c = MatCol(S, j);
          q += 2.0 * v[j] * DotProd(c, v, j) + c[j] * v[j] * v[j];
     }

     return q;
}
