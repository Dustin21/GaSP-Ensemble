/*****************************************************************/
/*                                                               */
/*   MATRIX (RE)ALLOCATION ROUTINES                              */
/*                                                               */
/*   Copyright (c) William J. Welch 1991--94.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

/*******************************+++*******************************/
void MatInit(int Shape, int Type, boolean Labelled, Matrix *M)
/*****************************************************************/
/*   Purpose:  Initialize a matrix as 0 x 0.                     */
/*                                                               */
/*   Version:  1995 January 13                                   */
/*****************************************************************/
{
     M->NumRows = M->NumCols = 0;

     M->Shape     = Shape;

     M->Type = Type;
     M->ColType = NULL;

     M->IntElem    = NULL;
     M->Elem       = NULL;
     M->Size_tElem = NULL;
     M->StrElem    = NULL;

     M->Labelled = Labelled;

     M->Text = NULL;
     M->RowName = M->ColName = NULL;

     M->Initialized = YES;

     M->Next = NULL;
}

/*******************************+++*******************************/
void MatReAllocate(size_t NewNumRows, size_t NewNumCols,
          const int *NewColType, Matrix *M)
/*****************************************************************/
/*   Purpose:  Re-allocate a matrix.                             */
/*                                                               */
/*   Args:     NewNumRows     Number of rows.                    */
/*             NewNumCols     Number of columns.                 */
/*             NewColType     Types for each column - only used  */
/*                            if the matrix Type is MIXED.       */
/*                            NewColType will not be accessed    */
/*                            for "old" columns.                 */
/*                            NewColType[j] always refers to     */
/*                            column j, not "new" column j.      */
/*             M              The matrix to be re-allocated.     */
/*                                                               */
/*   Comments: Use MatAllocate() for first allocation.           */
/*             Shape and Labelled cannot be changed.             */
/*             Types of "old" columns cannot be changed.         */
/*                                                               */
/*   Version:  1995 March 9                                      */
/*****************************************************************/
{
     size_t    j, NewLen, OldNumCols, OldNumRows;

     CodeCheck(MatInitialized(M));

     OldNumRows = MatNumRows(M);
     OldNumCols = MatNumCols(M);

     /* Free excess columns. */
     for (j = NewNumCols; j < OldNumCols; j++)
          MatColReAlloc(0, j, M);

     if (NewNumCols != OldNumCols)
     {
          /* Reallocate column pointers and column types. */
          M->Elem       = AllocPtrReal(NewNumCols, M->Elem);
          M->IntElem    = AllocPtrInt(NewNumCols, M->IntElem);
          M->Size_tElem = AllocPtrSize_t(NewNumCols, M->Size_tElem);
          M->StrElem    = AllocPtrStr(NewNumCols, M->StrElem);
          M->ColType    = AllocInt(NewNumCols, M->ColType);
     }

     /* Initialize pointers and types for new columns. */
     for (j = OldNumCols; j < NewNumCols; j++)
     {
          M->IntElem[j]    = NULL;
          M->Elem[j]       = NULL;
          M->Size_tElem[j] = NULL;
          M->StrElem[j]    = NULL;

          if (MatType(M) != MIXED)
               M->ColType[j] = MatType(M);
          else if (NewColType != NULL)
               M->ColType[j] = NewColType[j];
          else
               Fatal("Code bug: NewColType not assigned in "
                          "MatReAllocate.\n");

          if (M->ColType[j] != MatType(M))
               MatPutType(M, MIXED);
     }

     /* Reallocate all surviving columns. */
     for (j = 0; j < NewNumCols; j++)
     {
          NewLen = (M->Shape == RECT) ? NewNumRows : j + 1;
          MatColReAlloc(NewLen, j, M);
     }

     if (MatLabelled(M))
     {
          /* Reallocate labels. */
          M->RowName = AllocStrFree(OldNumRows, NewNumRows,
                    M->RowName);
          M->ColName = AllocStrFree(OldNumCols, NewNumCols,
                    M->ColName);
     }

     MatPutNumRows(M, NewNumRows);
     MatPutNumCols(M, NewNumCols);

     return;
}

/*******************************+++*******************************/
void MatFree(Matrix *M)
/*****************************************************************/
/*   Purpose:  Free a matrix and any further matrices in the     */
/*             list.                                             */
/*                                                               */
/*   Version:  1995 March 9                                      */
/*****************************************************************/
{
     do
     {
          /* Make M 0 x 0. */
          MatReAllocate(0, 0, NULL, M);

          if (MatLabelled(M))
          {
               AllocFree(M->Text);
               M->Text = NULL;
          }

     } while ( (M = M->Next) != NULL);
}

/*******************************+++*******************************/
/*                                                               */
/*   void      MatColReAlloc(size_t NewLen, size_t j, Matrix *M) */
/*                                                               */
/*   Purpose:  Re-allocate a column of a matrix.                 */
/*                                                               */
/*   Args:     NewLen    New length.                             */
/*             j         Index of the column.                    */
/*             M         The matrix to be re-allocated.          */
/*                                                               */
/*   Version:  1992 January 27                                   */
/*                                                               */
/*****************************************************************/

void MatColReAlloc(size_t NewLen, size_t j, Matrix *M)
{
     int       *IntElem;
     real      *Elem;
     size_t    i, OldLen;
     size_t    *Size_tElem;

     OldLen = (M->IntElem[j] || M->Elem[j] || M->Size_tElem[j]
               || M->StrElem[j]) ? MatColLen(M, j) : 0;

     if (NewLen == OldLen)
          return;

     switch (MatColType(M, j))
     {
          case REAL:
               Elem = M->Elem[j] = AllocReal(NewLen, M->Elem[j]);
               for (i = OldLen; i < NewLen; i++)
                    Elem[i] = 0.0;
               break;
          case INTEGER:
               IntElem = M->IntElem[j]
                         = AllocInt(NewLen, M->IntElem[j]);
               for (i = OldLen; i < NewLen; i++)
                    IntElem[i] = 0;
               break;
          case SIZE_T:
               Size_tElem = M->Size_tElem[j]
                         = AllocSize_t(NewLen, M->Size_tElem[j]);
               for (i = OldLen; i < NewLen; i++)
                    Size_tElem[i] = 0;
               break;
          case STRING:
               M->StrElem[j] = AllocStrFree(OldLen, NewLen,
                         M->StrElem[j]);
               break;

          default:
               Fatal("Code bug: Illegal type in MatColReAlloc.\n");
     }

     return;
}

/*******************************+++*******************************/
size_t MatColumnAdd(const string ColName, int NewColType,
          Matrix *M)
/*****************************************************************/
/*   Purpose:  Find, or add, a named column.                     */
/*                                                               */
/*   Args:     ColName        Cannot be NULL.                    */
/*             NewColType     The type of the column.            */
/*             M              A matrix.                          */
/*                                                               */
/*   Returns:  Index of the new column.                          */
/*                                                               */
/*   Comments: See macros MatIntColAdd(), etc.                   */
/*             M must be labelled.                               */
/*                                                               */
/*   96.01.18: CodeBug parameters changed and CodeBug replaced   */
/*             by CodeCheck.                                     */
/*   96.01.22: CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     int       *ColType;
     size_t    ColIndex, i, nRows, NewNumCols;

     CodeCheck(MatLabelled(M));

     if ( (ColIndex = MatColIndex(M, ColName)) != INDEX_ERR)
     {
          /* Column already exists. */
          CodeCheck(MatColType(M, ColIndex) == NewColType);
     }
     else
     {
          /* New column. */
          ColIndex = MatNumCols(M);

          nRows = MatNumRows(M);
          NewNumCols = MatNumCols(M) + 1;

          ColType = AllocInt(NewNumCols, NULL);
          ColType[NewNumCols-1] = NewColType;

          MatReAllocate(nRows, NewNumCols, ColType, M);

          AllocFree(ColType);

          MatPutColName(M, NewNumCols - 1, ColName);

          /* Initialize as Not Available. */
          switch (NewColType)
          {
               case REAL:
                    VecInit(NA_REAL, nRows,
                              MatCol(M, NewNumCols - 1));
                    break;

               case INTEGER:
                    for (i = 0; i < nRows; i++)
                         MatPutIntElem(M, i, NewNumCols - 1, NA_INT);
                    break;

               case SIZE_T:
                    for (i = 0; i < nRows; i++)
                         MatPutSize_tElem(M, i, NewNumCols - 1,
                                   NA_SIZE_T);
                    break;

               case STRING:
                    VecStrInit(NOT_AVAIL, nRows,
                              MatStrCol(M, NewNumCols - 1));
                    break;

               default:
                    CodeBug("Illegal type");
          }
     }

     return ColIndex;
}

/*******************************+++*******************************/
void MatRowAdd(size_t nVars, const string *VarName,
          const real *row, matrix *M)
/*****************************************************************/
/*   Purpose:  Put the nVars variables in VarNames into a new    */
/*             row of M.                                         */
/*                                                               */
/*   Version:  1996.04.07                                        */
/*****************************************************************/
{
     real      *Col;
     size_t    j, n;

     n = MatNumRows(M);
     MatReAlloc(n + 1, MatNumCols(M), M);

     for (j = 0; j < nVars; j++)
     {
          Col = MatColFind(M, VarName[j], YES);
          Col[n] = row[j];
     }

     return;
}
