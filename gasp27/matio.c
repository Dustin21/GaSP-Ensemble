/*****************************************************************/
/*   MATRIX INPUT-OUTPUT ROUTINES                                */
/*                                                               */
/*   Copyright (c) William J. Welch 1991--95.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

#define BETWEEN_SPACES   2

/*******************************+++*******************************/
int MatRead(FILE *InpFile, int Type, Matrix *M)
/*******************************+++*******************************/
/*   Purpose:  Read a matrix from a file.                        */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  Only RECT, Labelled, REAL or STRING matrices can  */
/*             be read.                                          */
/*                                                               */
/*   96.01.22: IllegalType replaced by CodeCheck.                */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     boolean   Finished;
     int       ErrNum;
     Matrix    Block;
     string    Buf;

     CodeCheck(Type == REAL || Type == STRING);

     /* Initialize M and Block as 0 x 0 matrices. */
     MatInit(RECT, Type, YES, M);
     MatInit(RECT, Type, YES, &Block);

     /* Get matrix name. */
     while ( (Buf = BufRead(InpFile)) != NULL &&
               strstr(Buf, "---") == NULL &&
               strstr(Buf, "___") == NULL)
          MatPutText(M, StrCatAlloc(MatText(M), Buf));

     if (Buf == NULL)
     {
          Error("Found nothing following the description.\n");
          ErrNum = INPUT_ERR;
     }
     else
     {
           /* Read first block into M. */
           ErrNum = MatReadABlock(InpFile, Type, M, &Finished);
     }

     while (ErrNum == OK && !Finished)
     {
          ErrNum = MatReadABlock(InpFile, Type, &Block, &Finished);

          if (ErrNum == OK)
               ErrNum = MatMerge(M, &Block);

          if (ErrNum != OK)
               MatFree(&Block);
     }

     if (ErrNum != OK)
          MatFree(M);

     fclose(InpFile);

     return ErrNum;
}

#define ROW_ALLOC   100  /* Was 32. */

/*******************************+++*******************************/
/*                                                               */
/*   int MatReadABlock(FILE *InpFile, int Type, Matrix *Block,   */
/*             boolean *Finished)                                */
/*                                                               */
/*   Purpose:  Read a block of a matrix from a file.             */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Version:  1993 October 11                                   */
/*                                                               */
/*****************************************************************/

int MatReadABlock(FILE *InpFile, int Type, Matrix *Block,
          boolean *Finished)
{
     int       ErrNum;
     real      r;
     size_t    CaseLabels, i, j, NumCols;
     string    Buf, Token;

     /* Start a new line. */
     Buf = BufRead(InpFile);

     /* Are case labels supplied?. */
     if (stricmp(Token = BufForceTok(InpFile, &Buf), "Case") == 0)
     {
          CaseLabels = YES;
          Token = BufForceTok(InpFile, &Buf);
     }
     else
          CaseLabels = NO;

     /* Read column labels. */
     NumCols = 0;
     while (Token != NULL && strstr(Token, "---") == NULL
               && strstr(Token, "___") == NULL)
     {
          MatReAllocate(0, NumCols + 1, NULL, Block);
          MatPutColName(Block, NumCols, Token);
          NumCols++;
          Token = BufForceTok(InpFile, &Buf);
     }

     if (Token == NULL)
     {
          Error("Column labels should terminate with a line "
                    "containing \"---\" or \"___\".\n");
          ErrNum = INPUT_ERR;
     }

     else if (NumCols == 0)
          ErrNum = OK;

     else
     {
          ErrNum = OK;

          /* Start a new line for the data. */
          Buf = BufRead(InpFile);

          /* Read tokens until EOF or next block. */
          i = 0;
          j = 0;
          while (ErrNum == OK &&
                    (Token = BufForceTok(InpFile, &Buf)) != NULL &&
                    strstr(Token, "---") == NULL &&
                    strstr(Token, "___") == NULL)
          {
               if (j == 0 && i % ROW_ALLOC == 0)
                    /* Reallocate ROW_ALLOC new rows at once */
                    /* to avoid fragmentation of memory.     */
                    MatReAllocate(i + ROW_ALLOC, NumCols, NULL, Block);

               if (j == 0 && CaseLabels)
                    MatPutRowName(Block, i, Token);
               else if (Type == STRING)
                    MatPutStrElem(Block, i, j - CaseLabels, Token);
               else if (StrToReal(Token, &r) != OK)
               {
                    Error("\"%s\" at row %d, column \"%s\" should "
                              "be a (real) number.\n", Token, i + 1,
                              MatColName(Block, j - CaseLabels));
                    ErrNum = INPUT_ERR;
                    break;
               }
               else
                    MatPutElem(Block, i, j - CaseLabels, r);

               if (++j == NumCols + CaseLabels)
               {
                    i++;
                    j = 0;
               }
          }

          if (ErrNum == OK && j != 0)
          {
               Error("Row %d is incomplete.\n", i + 1);
               ErrNum = INPUT_ERR;
          }
     }

     *Finished = (Token == NULL) ? YES : NO;

     if (ErrNum != OK)
          MatFree(Block);
     else
          /* Reallocate for actual number of rows. */
          MatReAllocate(i, NumCols, NULL, Block);

     return ErrNum;
}

/*******************************+++*******************************/
void MatWrite(Matrix *M, FILE *OutFile)
/*****************************************************************/
/*   Purpose:  Write a matrix to a file.                         */
/*                                                               */
/*   Comment:  Only RECT, Labelled matrices can be written.      */
/*             A matrix row may be written on several output     */
/*             lines.  Use MatWriteBlock() for "nice" output.    */
/*                                                               */
/*   Version:  1995 October 25                                   */
/*****************************************************************/
{
     boolean   RightAdj;
     char      *Conversion;
     int       *Precision;
     size_t    CaseWidth, ColsPerLine, ColWidth, i, j;
     size_t    LineWidth, MaxColWidth, NumCols, NumRows;
     string    *ColName, s;

     NumRows = MatNumRows(M);
     NumCols = MatNumCols(M);
     ColName = MatColNames(M);

     /* Output the text. */
     FileOutput(OutFile, "%s", (MatText(M) != NULL) ? MatText(M) :
               "Unnamed matrix.\n\n");

     /* Column width for case labels. */
     CaseWidth = MatCaseWidth(M, &RightAdj);

     /* Get maximum column width and precisions for real columns. */
     Conversion = AllocChar(NumCols, NULL);
     Precision  = AllocInt(NumCols, NULL);
     for (MaxColWidth = 0, j = 0; j < NumCols; j++)
     {
          ColWidth = MatColWidth(M, j, &Precision[j],
                    &Conversion[j]);
          MaxColWidth = max(ColWidth + BETWEEN_SPACES, MaxColWidth);
     }

     /* Always output at least one column per line. */
     ColsPerLine = max(1, (OUTPUT_COLS - CaseWidth) / MaxColWidth);
     LineWidth = CaseWidth + min(NumCols, ColsPerLine) * MaxColWidth;

     /* Output "Case" and column names. */

     for (j = 0; j < LineWidth; j++)
          FileOutput(OutFile, "-");
     FileOutput(OutFile, "\n");

     FileOutput(OutFile, "%*s", (int) CaseWidth, "Case");
     for (j = 0; j < NumCols; j++)
     {
          if (j % ColsPerLine == 0 && j != 0)
               FileOutput(OutFile, "\n%*s", (int) CaseWidth, "");
          FileOutput(OutFile, "%*s", (int) MaxColWidth, ColName[j]);
     }
     FileOutput(OutFile, "\n");

     for (j = 0; j < LineWidth; j++)
          FileOutput(OutFile, "-");
     FileOutput(OutFile, "\n\n");

     /* Output the data. */
     for (i = 0; i < NumRows; i++)
     {
          FileOutput(OutFile, "%*s", (int) CaseWidth,
                    MatRowName(M, i));
          for (j = 0; j < NumCols; j++)
          {
               if (j % ColsPerLine == 0 && j != 0)
                    FileOutput(OutFile, "\n%*s",
                              (int) CaseWidth, "");

               s = MatElemToStr(M, i, j, Precision[j],
                         Conversion[j]);
               FileOutput(OutFile, "%*s", (int) MaxColWidth, s);
          }
          FileOutput(OutFile, "\n");
     }

     AllocFree(Conversion);
     AllocFree(Precision);

     return;
}

/*******************************+++*******************************/
void MatWriteBlock(Matrix *M, boolean CaseLabels, FILE *OutFile)
/*****************************************************************/
/*   Purpose:  Write a matrix to a file in blocks of columns     */
/*             that fit on a single output line.                 */
/*                                                               */
/*   Comment:  Only RECT, Labelled matrices can be written.      */
/*                                                               */
/*   Version:  1995 October 25                                   */
/*****************************************************************/
{
     boolean   RightAdj;
     char      *Conversion;
     size_t    CaseWidth, FirstCol, i, j, LastCol;
     size_t    LineWidth, NumCols, NumRows;
     size_t    *ColWidth;
     int       *Precision;
     string    *ColName, s;

     NumRows = MatNumRows(M);
     NumCols = MatNumCols(M);
     ColName = MatColNames(M);

     /* Output the matrix name. */
     FileOutput(OutFile, "%s", (MatText(M) != NULL) ? MatText(M) :
               "Unnamed matrix.\n\n");

     /* Get column widths and precisions. */
     ColWidth   = AllocSize_t(NumCols, NULL);
     Conversion = AllocChar(NumCols, NULL);
     Precision = AllocInt(NumCols, NULL);
     for (j = 0; j < NumCols; j++)
          ColWidth[j] = MatColWidth(M, j, &Precision[j],
                    &Conversion[j]);

     /* Column width for case labels. */
     if (CaseLabels)
          CaseWidth = MatCaseWidth(M, &RightAdj);
     else
          CaseWidth = 0;

     /* Output in blocks of columns that will fit on a line. */
     FirstCol = 0;
     while (FirstCol < NumCols)
     {
          /* Always output at least one column, whatever. */
          LastCol = FirstCol;
          LineWidth = ((CaseLabels) ? CaseWidth + BETWEEN_SPACES : 0)
                    + ColWidth[FirstCol];

          while (LastCol + 1 < NumCols &&
                    LineWidth + BETWEEN_SPACES
                    + ColWidth[LastCol+1] <= OUTPUT_COLS)
               LineWidth += BETWEEN_SPACES + ColWidth[++LastCol];

          for (j = 0; j < LineWidth; j++)
               FileOutput(OutFile, "-");
          FileOutput(OutFile, "\n");

          if (CaseLabels)
               FileOutput(OutFile, (RightAdj) ? "%*s" : "%-*s",
                         (int) CaseWidth, "Case");
          for (j = FirstCol; j <= LastCol; j++)
          {
               if (CaseLabels || j > FirstCol)
                    FileOutput(OutFile, "%*s", BETWEEN_SPACES, "");
               FileOutput(OutFile, (MatColType(M, j) == STRING)
                         ? "%-*s" : "%*s", (int) ColWidth[j],
                         ColName[j]);
          }
          FileOutput(OutFile, "\n");

          for (j = 0; j < LineWidth; j++)
               FileOutput(OutFile, "-");
          FileOutput(OutFile, "\n\n");

          for (i = 0; i < NumRows; i++)
          {
               if (CaseLabels)
                    FileOutput(OutFile, (RightAdj) ? "%*s" : "%-*s",
                              (int) CaseWidth, MatRowName(M, i));
               for (j = FirstCol; j <= LastCol; j++)
               {
                    if (CaseLabels || j > FirstCol)
                         FileOutput(OutFile, "%*s", BETWEEN_SPACES, "");
                    s = MatElemToStr(M, i, j, Precision[j],
                              Conversion[j]);
                    FileOutput(OutFile, (MatColType(M, j) == STRING)
                              ? "%-*s" : "%*s", (int) ColWidth[j], s);
               }
               FileOutput(OutFile, "\n");
          }

          if ( (FirstCol = LastCol + 1) < NumCols)
               FileOutput(OutFile, "\n");
     }

     AllocFree(ColWidth);
     AllocFree(Conversion);
     AllocFree(Precision);

     return;
}

/*******************************+++*******************************/
size_t MatColWidth(const Matrix *M, size_t j, int *Precision,
               char *Conversion)
/*****************************************************************/
/*   Purpose:  Return the column width necessary for outputting  */
/*             column j of M.                                    */
/*             If the column is of type real, then, on return,   */
/*             *Precision will be the minimum precision that     */
/*             does not lose accuracy, and *Conversion will be   */
/*             'e' or 'f'.                                       */
/*                                                               */
/*   Version:  1995 October 25                                   */
/*****************************************************************/
{
     size_t    ExponLen, i, Width;
     string    DecPoint, Expon, s, StrEnd;

     Width = 0;

     if (MatColType(M, j) == REAL)
     {
          *Conversion = 'g';
          *Precision = 0;

          /* Are there any e conversions? */
          for (i = 0; i < MatNumRows(M); i++)
               if (strchr(StrFromReal(MatElem(M, i, j), "",
                         PRECISION, 'g'), 'e') != NULL)
               {
                    *Conversion = 'e';
                    break;
               }

          for (i = 0; i < MatNumRows(M); i++)
          {
               /* Conversion is either e or g (not f) to maintain */
               /* the minimum number of significant digits.       */
               s = StrFromReal(MatElem(M, i, j), "#",
                         PRECISION, *Conversion);

               if (stricmp(s, NOT_AVAIL) == 0)
                    Width = max(strlen(s), Width);
               else
               {
                    ExponLen = 0;
                    if ( (Expon = strchr(s, 'e')) != NULL)
                    {
                         ExponLen = strlen(Expon);

                         /* Delete exponent part of string. */
                         *Expon = NULL;
                    }

                    /* Find decimal point.  This assumes that the */
                    /* real has been converted to a string with a */
                    /* # flag, so that '.' is  always present.    */
                    DecPoint = strchr(s, '.');
                    CodeCheck(DecPoint != NULL);

                    /* Delete trailing zeros to get precision. */
                    StrEnd = s + strlen(s) - 1;
                    while (StrEnd > DecPoint && *StrEnd == '0')
                         *StrEnd-- = NULL;
                    *Precision = max(StrEnd - DecPoint, *Precision);

                    /* Delete digits after decimal point;        */
                    /* if none, delete the decimal point itself. */
                    if (StrEnd == DecPoint)
                         *DecPoint = NULL;
                    else
                         *(DecPoint + 1) = NULL;

                    /* Everything but precision. */
                    Width = max(strlen(s) + ExponLen, Width);
               }
          }

          Width += (size_t) *Precision;
          if (*Conversion == 'g')
               *Conversion = 'f';
     }
     else
          for (i = 0; i < MatNumRows(M); i++)
               Width = max(strlen(MatElemToStr(M, i, j, 'x', -1)),
                         Width);

     Width = max(Width, strlen(MatColName(M, j)));

     return Width;
}

/*******************************+++*******************************/
size_t MatCaseWidth(const Matrix *M, boolean *RightAdj)
/*****************************************************************/
/*   Purpose:  Return the column width necessary for outputting  */
/*             the case column of M.                             */
/*                                                               */
/*   Version:  1995 October 25                                   */
/*****************************************************************/
{
     size_t    CaseNumber, i, Width;
     string    s;

     *RightAdj = TRUE;
     for (Width = 0, i = 0; i < MatNumRows(M); i++)
     {
          s = MatRowName(M, i);

          Width = max(strlen(s), Width);

          if (StrToSize_t(s, &CaseNumber) != OK)
               /* Non-numeric case labels. */
               *RightAdj = FALSE;
     }

     Width = max(Width, strlen("Case"));

     return Width;
}

/*******************************+++*******************************/
string MatElemToStr(const Matrix *M, size_t i, size_t j,
          int Precision, char Conversion)
/*****************************************************************/
/*   Purpose:  Return a string representing element i, j of M.   */
/*             Precision and Conversion are only used if the     */
/*             element is real.                                  */
/*                                                               */
/*   96.01.22: IllegalType replaced by CodeBug.                  */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     string    s;

     switch (MatColType(M, j))
     {
          case INTEGER:
               s = StrFromInt(MatIntElem(M, i, j));
               break;

          case REAL:
               s = StrFromReal(MatElem(M, i, j), "", Precision,
                         Conversion);
               break;

          case SIZE_T:
               s = StrFromSize_t(MatSize_tElem(M, i, j));
               break;

          case STRING:
               s = MatStrElem(M, i, j);
               break;

          default:
               CodeBug("Illegal type");
     }

     return s;
}

/*******************************+++*******************************/
void MatPrint(const Matrix *M)
/*****************************************************************/
/*   Purpose:  Print a (small) real matrix of any shape to       */
/*             stdin.  Intended for debugging.                   */
/*                                                               */
/*   Version:  1993 December 26                                  */
/*****************************************************************/
{
     size_t    i, j;

     for (i = 0; i < MatNumRows(M); i++)
     {
          for (j = 0; j < MatNumCols(M); j++)
               if (j < MatColLen(M, j))
                    printf("%15e", MatElem(M, i, j));
               else
                    printf("%15s", "");
          printf("\n");
     }
}
