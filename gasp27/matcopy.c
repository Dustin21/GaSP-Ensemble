/*****************************************************************/
/*   ROUTINES FOR COPYING MATRICES                               */
/*                                                               */
/*   Copyright (c) William J. Welch 1992--94.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

extern size_t  nPointers;

/*******************************+++*******************************/
void MatCopySub(size_t m, size_t n, size_t SrcRowOffset,
          size_t SrcColOffset, const Matrix *Src,
          size_t DestRowOffset, size_t DestColOffset, Matrix *Dest)
/*****************************************************************/
/*   Purpose:  Copy m x n submatrix from Src, starting at        */
/*             (SrcRowOffset, SrcColOffset), to Dest, starting   */
/*             at (DestRowOffset, DestColOffset).                */
/*                                                               */
/*   Comments: Calling routine must allocate space for Dest.     */
/*             Labels are not copied.                            */
/*             MatCopy() is a macro with a simplified            */
/*             parameter list for copying an entire matrix.      */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 18                                   */
/*****************************************************************/
{
     size_t    j, Len;

     CodeCheck(MatType(Src)  == MatType(Dest));
     CodeCheck(MatShape(Src) == MatShape(Dest));
     CodeCheck(SrcRowOffset  + m <= MatNumRows(Src));
     CodeCheck(DestRowOffset + m <= MatNumRows(Dest));
     CodeCheck(SrcColOffset  + n <= MatNumCols(Src));
     CodeCheck(DestColOffset + n <= MatNumCols(Dest));

     /* Copy elements. */
     for (j = 0; j < n; j++)
     {
          Len = MatColLen(Src, SrcColOffset + j);
          if (Len > SrcRowOffset)
          {
               Len = min(m, Len - SrcRowOffset);

               MatCopyColSub(Len, SrcColOffset + j,
                         SrcRowOffset, Src, DestColOffset + j,
                         DestRowOffset, Dest);
          }
     }
}

/*******************************+++*******************************/
void MatDup(const Matrix *Src, Matrix *Dest)
/*****************************************************************/
/*   Purpose:  Duplicates matrix Src to matrix Dest.             */
/*                                                               */
/*   Comments: Space is allocated here for Dest.                 */
/*                                                               */
/*   Version:  1994 November 24                                  */
/*****************************************************************/
{
     MatAllocate(MatNumRows(Src), MatNumCols(Src), MatShape(Src),
               MatType(Src), MatColTypes(Src), MatLabelled(Src),
               Dest);

     /* Copy labels. */
     if (MatLabelled(Src))
     {
          MatPutText(Dest, MatText(Src));
          VecStrCopy(MatRowNames(Src), MatNumRows(Src),
                    MatRowNames(Dest));
          VecStrCopy(MatColNames(Src), MatNumCols(Src),
                    MatColNames(Dest));
     }

     MatCopy(Src, Dest);
}

/*******************************+++*******************************/
void MatDupIndex(size_t nIn, const size_t *In, const Matrix *M,
          Matrix *NewM)
/*****************************************************************/
/*   Purpose:  Select rows indexed by In from M and put them in  */
/*             NewM.                                             */
/*                                                               */
/*   Comment:  NewM is allocated here.                           */
/*                                                               */
/*   Version:  1994 November 24                                  */
/*****************************************************************/
{
     size_t    i;

     MatAllocate(nIn, MatNumCols(M), RECT, MatType(M),
               MatColTypes(M), MatLabelled(M), NewM);

     if (MatLabelled(M))
     {
           MatPutText(NewM, MatText(M));
           VecStrCopy(MatColNames(M), MatNumCols(M),
                     MatColNames(NewM));
     }

     for (i = 0; i < nIn; i++)
          MatCopyRow(In[i], M, i, NewM);

     return;
}

/*******************************+++*******************************/
void MatCopyColSub(size_t m, size_t j, size_t SrcOffset,
          const Matrix *Src, size_t k, size_t DestOffset,
          Matrix *Dest)
/*****************************************************************/
/*   Purpose:  Copy m elements of column j of Src, starting at   */
/*             SrcOffset, to column k of Dest, starting at       */
/*             DestOffset.                                       */
/*                                                               */
/*   Comment:  MatCopyCol() is a macro with a simplified         */
/*             parameter list for copying an entire column.      */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*   96.01.22: IllegalType replaced by CodeBug.                  */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     CodeCheck(MatColType(Src, j) == MatColType(Dest, k));
     CodeCheck(SrcOffset  + m <= MatColLen(Src, j));
     CodeCheck(DestOffset + m <= MatColLen(Dest, k));

     switch (MatColType(Src, j))
     {
          case INTEGER:
               VecIntCopy(MatIntCol(Src, j) + SrcOffset, m,
                         MatIntCol(Dest, k) + DestOffset);
               break;

          case REAL:
               VecCopy(MatCol(Src, j) + SrcOffset, m,
                         MatCol(Dest, k) + DestOffset);
               break;

          case SIZE_T:
               VecSize_tCopy(MatSize_tCol(Src, j) + SrcOffset, m,
                         MatSize_tCol(Dest, k) + DestOffset);
               break;

          case STRING:
               VecStrCopy(MatStrCol(Src, j) + SrcOffset, m,
                         MatStrCol(Dest, k) + DestOffset);
               break;

          default:
               CodeBug("Illegal type");
     }

     return;
}

/*******************************+++*******************************/
void MatCopyRow(size_t i, const Matrix *Src, size_t k,
          Matrix *Dest)
/*****************************************************************/
/*   Purpose:  Copy row i of Src to row k of Dest.               */
/*             If both matrices are labelled, then the row label */
/*             is also copied.                                   */
/*                                                               */
/*   Comments: Calling routine must allocate space for Dest.     */
/*             Both Src and Dest must be RECT matrices.          */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*   96.01.22: IllegalType replaced by CodeBug.                  */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     size_t    j;

     CodeCheck(MatNumCols(Src) == MatNumCols(Dest));
     CodeCheck(MatShape(Src) == RECT && MatShape(Dest) == RECT);
     CodeCheck(i < MatNumRows(Src) && k < MatNumRows(Dest));

     if (MatLabelled(Src) && MatLabelled(Dest))
          MatPutRowName(Dest, k, MatRowName(Src, i));

     for (j = 0; j < MatNumCols(Src); j++)
     {
          CodeCheck(MatColType(Src, j) == MatColType(Dest, j));

          switch (MatColType(Src, j))
          {
               case INTEGER:
                    MatPutIntElem(Dest, k, j,
                              MatIntElem(Src, i, j));
                    break;

               case REAL:
                    MatPutElem(Dest, k, j, MatElem(Src, i, j));
                    break;

               case SIZE_T:
                    MatPutSize_tElem(Dest, k, j,
                              MatSize_tElem(Src, i, j));
                    break;

               case STRING:
                    MatPutStrElem(Dest, k, j,
                              MatStrElem(Src, i, j));
                    break;

               default:
                    CodeBug("Illegal type");
          }
     }

     return;
}
