/*****************************************************************/
/*   MATRIX UTILITY ROUTINES                                     */
/*                                                               */
/*   Copyright (c) William J. Welch 1991--96.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

extern int     ErrorSeverityLevel;

/*******************************+++*******************************/
/*   int       *MatIntCol(const Matrix *M, size_t j)             */
/*   real      *MatCol(const Matrix *M, size_t j)                */
/*   size_t    *MatSize_tCol(const Matrix *M, size_t j)          */
/*   string    *MatStrCol(const Matrix *M, size_t j)             */
/*                                                               */
/*   Purpose:  Return a pointer of the appropriate type to       */
/*             column j.                                         */
/*                                                               */
/*   Comment:  Terminates program if the type is illegal or j is */
/*             out of range.                                     */
/*             94.03.25: rearranged to avoid compiler warning    */
/*                       messages.                               */
/*                                                               */
/*   Version:  1996 January 10                                   */
/*****************************************************************/

int *MatIntCol(const Matrix *M, size_t j)
{
     CodeCheck(j < MatNumCols(M) && MatColType(M, j) == INTEGER)

     return M->IntElem[j];
}

real *MatCol(const Matrix *M, size_t j)
{
     CodeCheck(j < MatNumCols(M) && MatColType(M, j) == REAL)

     return M->Elem[j];
}

size_t *MatSize_tCol(const Matrix *M, size_t j)
{
     CodeCheck(j < MatNumCols(M) && MatColType(M, j) == SIZE_T)

     return M->Size_tElem[j];
}

string *MatStrCol(const Matrix *M, size_t j)
{
     CodeCheck(j < MatNumCols(M) && MatColType(M, j) == STRING)

     return M->StrElem[j];
}

/*******************************+++*******************************/
/*   int       *MatIntColFind(const Matrix *M, const string Name,*/
/*                  boolean HardFail)                            */
/*   real      *MatColFind(const Matrix *M, const string Name,   */
/*                  boolean HardFail)                            */
/*   size_t    *MatSize_tColFind(const Matrix *M,                */
/*                  const string Name, boolean HardFail)         */
/*   string    *MatStrColFind(const Matrix *M, const string Name,*/
/*                  boolean HardFail)                            */
/*                                                               */
/*   Purpose:  Return a pointer of the appropriate type to       */
/*             a named column.                                   */
/*                                                               */
/*   Comment:  If the named column does not exist or the type is */
/*             illegal, then either:                             */
/*             (i)  the program terminates (HardFail == TRUE) or */
/*             (ii) NULL is returned (HardFail == FALSE).        */
/*             94.03.25: rearranged to avoid compiler warning    */
/*                       messages.                               */
/*                                                               */
/*   Version:  1996 January 10                                    */
/*****************************************************************/

int *MatIntColFind(const Matrix *M, const string Name,
          boolean HardFail)
{
     size_t j;

     j = MatColIndex(M, Name);

     if (j < MatNumCols(M) && MatColType(M, j) == INTEGER)
          return M->IntElem[j];
     else
          CodeCheck(!HardFail);

     return NULL;
}

real *MatColFind(const Matrix *M, const string Name,
          boolean HardFail)
{
     size_t j;

     j = MatColIndex(M, Name);

     if (j < MatNumCols(M) && MatColType(M, j) == REAL)
          return M->Elem[j];
     else
          CodeCheck(!HardFail);

     return NULL;
}

size_t *MatSize_tColFind(const Matrix *M, const string Name,
          boolean HardFail)
{
     size_t j;

     j = MatColIndex(M, Name);

     if (j < MatNumCols(M) && MatColType(M, j) == SIZE_T)
          return M->Size_tElem[j];
     else
          CodeCheck(!HardFail);

     return NULL;
}

string *MatStrColFind(const Matrix *M, const string Name,
          boolean HardFail)
{
     size_t j;

     j = MatColIndex(M, Name);

     if (j < MatNumCols(M) && MatColType(M, j) == STRING)
          return M->StrElem[j];
     else
          CodeCheck(!HardFail);

     return NULL;
}

/*******************************+++*******************************/
void *MatVoidCol(const Matrix *M, size_t j)
/*****************************************************************/
/*   Purpose:  Return a pointer to void for column j; e.g., if   */
/*             column j is INTEGER, then returns M->IntElem[j].  */
/*                                                               */
/*   Comment:  Terminates program if MatColType(M, j) is illegal */
/*             or if j is out of range.                          */
/*                                                               */
/*   94.03.25: Rearranged to avoid compiler warning messages.    */
/*   96.01.22: Repeated ifs replaced by switch and CodeBug       */
/*             parameters changed.                               */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     CodeCheck(j < MatNumCols(M))

     switch (MatColType(M, j))
     {
          case INTEGER:
               return M->IntElem[j];

          case REAL:
               return M->Elem[j];

          case SIZE_T:
               return M->Size_tElem[j];

          case STRING:
               return M->StrElem[j];

          default:
               CodeBug("Illegal type");
     }

     return NULL;   /* Never gets to here. */
}

/*******************************+++*******************************/
void MatInitValue(real s, Matrix *M)
/*****************************************************************/
/*   Purpose:  Initialize all elements of the REAL matrix M to   */
/*             the value s.                                      */
/*                                                               */
/*   Version:  1994 April 22                                     */
/*****************************************************************/
{
     size_t j;

     CodeCheck(MatType(M) == REAL);

     for (j = 0; j < MatNumCols(M); j++)
          VecInit(s, MatColLen(M,j), MatCol(M, j));

     return;
}

/*******************************+++*******************************/
void MatRow(const Matrix *M, size_t RowIndex, real *row)
/*****************************************************************/
/*   Purpose:  Load row RowIndex of M into row.                  */
/*                                                               */
/*   Comment:  Calling routine must allocate space for row.      */
/*                                                               */
/*   96.01.10  TriRow subsumed.                                  */
/*   96.01.22  CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     size_t    j, jFirst;

     CodeCheck(MatType(M) == REAL);

     switch(MatShape(M))
     {
          case RECT:
               jFirst = 0;
               break;

          case SYM:
          case UP_TRIANG:
               VecInit(0.0, RowIndex, row);
               jFirst = RowIndex;
               break;

          default:
                CodeBug("Illegal shape");
     }

     for (j = jFirst; j < MatNumCols(M); j++)
          row[j] = MatElem(M, RowIndex, j);

     return;
}

/*******************************+++*******************************/
void MatRowPut(const real *row, size_t RowIndex, Matrix *M)
/*****************************************************************/
/*   Purpose:  Load row into row RowIndex of M.                  */
/*                                                               */
/*   Comment:  Replaced TriRowPut.                               */
/*                                                               */
/*   96.01.22  CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     size_t    jFirst, j;

     CodeCheck(MatType(M) == REAL);

     switch(MatShape(M))
     {
          case RECT:
               jFirst = 0;
               break;

          case SYM:
          case UP_TRIANG:
               jFirst = RowIndex;
               break;

          default:
                CodeBug("Illegal shape");
     }

     for (j = jFirst; j < MatNumCols(M); j++)
          MatPutElem(M, RowIndex, j, row[j]);

     return;
}

/*******************************+++*******************************/
void MatVec(const Matrix *M, const real *x, real *y)
/*****************************************************************/
/*   Purpose:  Compute y = M'x.                                  */
/*                                                               */
/*   Comment:  Calling routine must allocate space for y.        */
/*                                                               */
/*   Version:  1995 February 14                                  */
/*****************************************************************/
{
     size_t    j, n;

     n = MatNumRows(M);
     for (j = 0; j < MatNumCols(M); j++)
          y[j] = DotProd(MatCol(M, j), x, n);
}

/*******************************+++*******************************/
void MatMultElemWise(const Matrix *A, Matrix *B)
/*****************************************************************/
/*   Purpose:  Replace B_ij by A_ij * B_ij.                      */
/*                                                               */
/*   Version:  1995 April 22                                     */
/*****************************************************************/
{
     size_t    j;

     CodeCheck(MatType(A) == REAL);
     CodeCheck(MatType(B) == REAL);
     CodeCheck(MatNumRows(A) == MatNumRows(B));
     CodeCheck(MatNumCols(A) == MatNumCols(B));

     for (j = 0; j < MatNumCols(A); j++)
          VecMultVec(MatCol(A, j), MatColLen(A, j), MatCol(B, j));
}

/*******************************+++*******************************/
/*                                                               */
/*   void      MatStack(const Matrix *M, boolean ByCols,         */
/*                  real *v)                                     */
/*                                                               */
/*   Purpose:  Stack the columns of M into the vector v.         */
/*                                                               */
/*   Version:  1992 February 9                                   */
/*                                                               */
/*****************************************************************/

void MatStack(const Matrix *M, boolean ByCols, real *v)
{
     size_t    j, NumCols, NumRows;

     NumRows = MatNumRows(M);
     NumCols = MatNumCols(M);

     for (j = 0; j < NumCols; j++)
          if (ByCols)
               VecCopy(MatCol(M, j), NumRows, v + j * NumRows);
          else
               VecCopyStride(NumRows, 1, MatCol(M, j), NumCols,
                         v + j);

     return;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      MatUnStack(const real *v, boolean ByCols,         */
/*                  Matrix *M)                                   */
/*                                                               */
/*   Purpose:  Unstack the vector v into the columns of M.       */
/*                                                               */
/*   Version:  1992 February 9                                   */
/*                                                               */
/*****************************************************************/

void MatUnStack(const real *v, boolean ByCols, Matrix *M)
{
     size_t    j, NumCols, NumRows;

     NumRows = MatNumRows(M);
     NumCols = MatNumCols(M);

     for (j = 0; j < NumCols; j++)
          if (ByCols)
               VecCopy(v + j * NumRows, NumRows, MatCol(M, j));
          else
               VecCopyStride(NumRows, NumCols, v + j, 1,
                         MatCol(M, j));

     return;
}

/*******************************+++*******************************/
int MatMerge(Matrix *M1, Matrix *M2)
/*****************************************************************/
/*   Purpose:  Merge matrices M1 and M2 so that, on return, M1   */
/*             contains the columns of both matrices, and M2 is  */
/*             0 x 0.                                            */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if M1 and M2 are incompatible        */
/*                          (the matrices are not merged);       */
/*             OK           otherwise.                           */
/*                                                               */
/*   Comments: Only labelled, rectangular matrices may be merged.*/
/*             A row-label mismatch causes a warning error, but  */
/*             merging continues and OK is returned.             */
/*             Change to allow unlabelled matrices.              */
/*                                                               */
/*   Version:  1995 May 2                                        */
/*****************************************************************/
{
     size_t    i, j, NumRows;

     if (MatNumCols(M2) == 0)
          return OK;

     if (MatShape(M1) != RECT || MatShape(M2) != RECT)
     {
          Error("Only rectangular matrices can be merged.\n");
          return INCOMPAT_ERR;
     }

     if ( (NumRows = MatNumRows(M1)) != MatNumRows(M2))
     {
          Error("Cannot merge matrices with different numbers "
                    "of rows.\n");
          return INCOMPAT_ERR;
     }

     for (j = 0; j < MatNumCols(M2); j++)
          if (StrIndex(MatColName(M2, j), MatColNames(M1),
                    MatNumCols(M1)) != INDEX_ERR)
          {
               Error("Cannot merge matrices with repeated column "
                         "names.\n");
               return INCOMPAT_ERR;
          }


     if ( (i = StrVecCmp(MatRowNames(M1), MatRowNames(M2),
               NumRows)) < NumRows)
     {
          ErrorSeverityLevel = SEV_WARNING;
          Error("Merging matrices with different row labels: "
                    "%s versus %s.\n",
                    MatRowName(M1, i), MatRowName(M2, i));
          ErrorSeverityLevel = SEV_ERROR;
     }

     for (j = 0; j < MatNumCols(M2); j++)
     {
          MatColumnAdd(MatColName(M2, j), MatColType(M2, j), M1);

          /* Copy column j of M2 to the last column of M1. */
          MatCopyCol(j, M2, MatNumCols(M1) - 1, M1);
     }

     MatFree(M2);

     return OK;
}

/*******************************+++*******************************/
size_t MatColConvert(size_t j, int NewColType, Matrix *M)
/*****************************************************************/
/*   Purpose:  Convert column j of M to a new column type.       */
/*                                                               */
/*   Returns:  i        if the element i cannot be converted     */
/*                      (M is unchanged);                        */
/*             INDEX_OK otherwise.                               */
/*                                                               */
/*   96.01.22: IllegalType replaced by CodeCheck or CodeBug.     */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     int       *i;
     real      *r;
     size_t    BadIndex, n, *z;
     string    *s;

     CodeCheck(MatColType(M, j) == STRING);

     s = MatStrCol(M, j);
     n = MatColLen(M, j);

     switch (NewColType)
     {
          case INTEGER:
               i = AllocInt(n, NULL);
               if ( (BadIndex = VecStrInt(s, n, i)) == INDEX_OK)
                    MatPutIntCol(M, j, i);
               else
                    AllocFree(i);
               break;

          case REAL:
               r = AllocReal(n, NULL);
               if ( (BadIndex = VecStrReal(s, n, r)) == INDEX_OK)
                    MatPutCol(M, j, r);
               else
                    AllocFree(r);
               break;

          case SIZE_T:
               z = AllocSize_t(n, NULL);
               if ( (BadIndex = VecStrSize_t(s, n, z)) == INDEX_OK)
                    MatPutSize_tCol(M, j, z);
               else
                    AllocFree(z);
               break;

          default:
               CodeBug("Illegal type");
     }

     if (BadIndex == INDEX_OK)
     {
          AllocStrFree(n, 0, s);
          MatPutStrCol(M, j, NULL);
          MatPutColType(M, j, NewColType);
          if (NewColType != MatType(M))
               MatPutType(M, MIXED);
     }

     return BadIndex;
}
