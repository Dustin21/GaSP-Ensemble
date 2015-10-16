/*****************************************************************/
/*   ROUTINES TO OVERSEE MATRIX AND SCALAR DATABASE OBJECTS      */
/*                                                               */
/*   Copyright (c) William J. Welch 1994.                        */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"
#include "kriging.h"
#include "optdes.h"
#include "alex.h"

extern int     Seed;

extern Matrix  Cand;
extern Matrix  ExpReg;
extern Matrix  PredReg;
extern Matrix  SPModMat;
extern Matrix  XDescrip;
extern Matrix  Y;
extern size_t  CritNum;
extern string  FuncName;

Matrix         DbStatus;

static string  DbStatusColName[] = DB_STATUS_COL_NAMES;
static int     DbStatusColType[] = DB_STATUS_COL_TYPES;

/******************************+++********************************/
int DbCheck(const string Name)
/*****************************************************************/
/*   Purpose:  Return status of database object and add a row to */
/*             DbStatus.                                         */
/*                                                               */
/*   Returns:  INPUT_ERR, INCOMPAT_ERR, or OK.                   */
/*                                                               */
/*   Version:  1996.02.11                                        */
/*****************************************************************/
{
     DbMatrix  *D;
     int       ErrNum;
     size_t    j, Row, Index;

     if (stricmp(Name, EXP_REG) == 0)
     {
          if ( (ErrNum = RegExtract(&XDescrip, X_DESCRIP, "." FIT,
                    &ExpReg)) == OK)
               ErrNum = RegCandCompat(&Cand, &ExpReg);

     }

     else if (stricmp(Name, PRED_REG) == 0)
     {
          if ( (ErrNum = RegExtract(&XDescrip, X_DESCRIP, "." PRED,
                    &PredReg)) == OK)
               ErrNum = RegCandCompat(&Cand, &PredReg);
     }

     if (stricmp(Name, EXP_REG) == 0 ||
               stricmp(Name, PRED_REG) == 0)
     {
          if (ErrNum != OK)
          {
               Row = StrIndex(X_DESCRIP,
                         MatStrCol(&DbStatus, DB_STATUS_OBJ_COL),
                         MatNumRows(&DbStatus));
               MatPutStrElem(&DbStatus, Row, DB_STATUS_OK_COL, NO_STR);
          }
          return ErrNum;
     }

     if (MatNumRows(&DbStatus) == 0)
     {
          MatAllocate(1, NumStr(DbStatusColName), RECT, MIXED,
                    DbStatusColType, YES, &DbStatus);
          for (j = 0; j < MatNumCols(&DbStatus); j++)
               MatPutColName(&DbStatus, j, DbStatusColName[j]);
          MatPutText(&DbStatus, "Database status:\n");
     }
     else
          /* Add a new row to DbStatus. */
          MatReAlloc(MatNumRows(&DbStatus) + 1,
                    MatNumCols(&DbStatus), &DbStatus);

     /* Status information goes in this row. */
     Row = MatNumRows(&DbStatus) - 1;

     MatPutStrElem(&DbStatus, Row, DB_STATUS_OBJ_COL, Name);

     if ( (Index = DbScalIndex(Name, NO)) != INDEX_ERR)
     {
          MatPutStrElem(&DbStatus, Row, DB_STATUS_CLASS_COL,
                    "Scalar");
          ErrNum = DbScalCheck(Index);
          MatPutStrElem(&DbStatus, Row, DB_STATUS_VAL_COL,
                    DbScalValue(Index));
          MatPutSize_tElem(&DbStatus, Row, DB_STATUS_ROWS_COL,
                    NA_SIZE_T);
          MatPutSize_tElem(&DbStatus, Row, DB_STATUS_COLS_COL,
                    NA_SIZE_T);

          if (ErrNum == OK && stricmp(Name, DESIGN_CRIT) == 0)
               /* Errors in subsequent input ignored. */
               CritInput(CritNum);
     }
     else
     {
          MatPutStrElem(&DbStatus, Row, DB_STATUS_CLASS_COL,
                    "Matrix");
          D = DbMatFind(Name, YES);

          ErrNum = DbMatCheck(D,
                    &MatStrElem(&DbStatus, Row, DB_STATUS_VAL_COL));
          MatPutSize_tElem(&DbStatus, Row, DB_STATUS_ROWS_COL,
                    D->M->NumRows);
          MatPutSize_tElem(&DbStatus, Row, DB_STATUS_COLS_COL,
                    D->M->NumCols);

          if (D->IsOutput == YES ||
                    ((stricmp(Name, REG_MOD) == 0 ||
                    stricmp(Name, SP_MOD) == 0) &&
                    stricmp(FuncName, "Fit") == 0))
               /* Store the matrix title as the row name:  */
               /* Will be used later in DbOutputMatStatus. */
               MatPutRowName(&DbStatus, Row, D->Title);
     }

     MatPutStrElem(&DbStatus, Row, DB_STATUS_OK_COL, (ErrNum == OK) ?
               YES_STR : NO_STR);

     /* Always override defaults. */
     if (stricmp(Name, RAN_NUM_SEED) == 0)
          RandInit(Seed, Seed, Seed);

     return ErrNum;
}

/*******************************+++*******************************/
void DbOutputMatStatus(void)
/*****************************************************************/
/*   Purpose:  Output status of output matrices.                 */
/*                                                               */
/*   Version:  1995 May 12                                       */
/*****************************************************************/
{
     Matrix    OutMatStatus;
     size_t    i, ii;

     MatAllocate(0, 2, RECT, STRING, NULL, YES, &OutMatStatus);
     MatPutColName(&OutMatStatus, 0, "Matrix");
     MatPutColName(&OutMatStatus, 1, "Contains");
     MatPutText(&OutMatStatus, "Output-matrix status:\n");

     for (ii = 0, i = 0; i < MatNumRows(&DbStatus); i++)
     {
          if (DbStatus.RowName[i] == NULL)
               /* Not an output matrix. */
               continue;

          MatReAlloc(ii + 1, 2, &OutMatStatus);
          MatPutStrElem(&OutMatStatus, ii, 0,
                    MatStrElem(&DbStatus, i, DB_STATUS_OBJ_COL));
          MatPutStrElem(&OutMatStatus, ii, 1,
                    MatRowName(&DbStatus, i));
          ii++;
     }

     /* Do not write case labels. */
     MatWriteBlock(&OutMatStatus, NO, stdout);
     Output("\n");
     MatFree(&OutMatStatus);
}
