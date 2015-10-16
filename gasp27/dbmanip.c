/*****************************************************************/
/*   ROUTINES TO MANIPULATE DATA                                 */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--6.                     */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"
#include "alex.h"

extern boolean DesignJob;

extern Matrix  X;
extern Matrix  YMat;
extern Matrix  YDescrip;
extern Matrix  YTrueMat;

extern size_t  nPointers;
extern size_t  nProtected;
extern string  ErrorVar;

size_t    nCasesX  = 0;
size_t    *IndexX  = NULL;

/* Set for for each y variable analyzed. */
real      *y       = NULL;
real      *yTrue   = NULL;
size_t    nCasesXY = 0;
size_t    *IndexXY = NULL;
string    yName    = NULL;

string TransName[] = TRANSFORM_NAMES;

/*******************************+++*******************************/
int DbProcessDescrip(const string DescripName,
     const Matrix *Descrip, DbMatrix *D)
/*****************************************************************/
/*   Purpose:  Process a data matrix by applying Descrip, i.e.,  */
/*             including only variables that are to be analyzed, */
/*             replacing cases out of range by NA's, and         */
/*             applying transformations.                         */
/*                                                               */
/*   96.03.07: IsFit/IsPred used for MinCol and MaxCol.          */
/*                                                               */
/*   Version:  1996.03.07                                        */
/*****************************************************************/
{
     boolean   Fixed;
     int       ErrNum;
     Matrix    *Data;
     Matrix    NewData;
     real      VarMax, VarMin;
     real      *MaxCol, *MinCol, *NewCol, *OldCol;
     size_t    i, j, n, nDomErr, nInputNA, nVars, TransNum;
     string    NewColName;
     string    *Analyze, *Support, *Trans, *Variable;

     ErrNum = OK;
     Data = D->M;

     nVars    = MatNumRows(Descrip);
     Variable = MatStrColFind(Descrip, VARIABLE,       YES);
     Support  = MatStrColFind(Descrip, SUPPORT,        NO);
     Analyze  = MatStrColFind(Descrip, ANALYZE,        NO);
     Trans    = MatStrColFind(Descrip, TRANSFORMATION, NO);

     if (D->IsFit)
     {
          MinCol = MatColFind(Descrip, MIN "." FIT,  NO);
          MaxCol = MatColFind(Descrip, MAX "." FIT,  NO);
     }
     else if (D->IsPred)
     {
          MinCol = MatColFind(Descrip, MIN "." PRED,  NO);
          MaxCol = MatColFind(Descrip, MAX "." PRED,  NO);
     }
     else
          MinCol = MaxCol = NULL;

     if (MinCol == NULL)
          MinCol = MatColFind(Descrip, MIN, NO);
     if (MaxCol == NULL)
          MaxCol = MatColFind(Descrip, MAX, NO);

     if (nVars == MatNumCols(Data) &&
               StrVecCmp(Variable, MatColNames(Data), nVars)
               == INDEX_OK &&
               Analyze == NULL && Trans == NULL &&
               MinCol == NULL && MaxCol == NULL)
          return OK;

     n = MatNumRows(Data);
     MatAllocate(n, 0, RECT, REAL, NULL, YES, &NewData);
     VecStrCopy(MatRowNames(Data), n, MatRowNames(&NewData));

     /*
     Output("Processing %s:\n", D->Name);
     Output("\n%15s%11s%20s%6s\n", VARIABLE, "InputNAs",
               "TransformationNAs", "Cases");
     */

     for (j = 0; j < nVars && ErrNum == OK; j++)
     {
          if (Analyze != NULL && stricmp(Analyze[j], NO) == 0)
               continue;

          Fixed = (Support != NULL &&
                    stricmp(Support[j], FIXED) == 0);

          if ( (OldCol = MatColFind(Data, Variable[j], NO)) == NULL)
          {
               ErrNum = INCOMPAT_ERR;

               if (!DesignJob && stricmp(D->Name, CAND) != 0)
                    Incompatibility(DB_COL_VAR, D->Name,
                              Variable[j], DescripName);

               else if (DesignJob && stricmp(D->Name, X_MAT) == 0 &&
                         nProtected > 0)
                    Incompatibility(DB_X_COL_PROT, Variable[j]);

               else if (DesignJob && stricmp(D->Name, X_MAT) == 0 &&
                         Fixed)
                    Incompatibility(DB_X_COL_FIXED, Variable[j]);

               else
                    ErrNum = OK;

               continue;
          }

          TransNum = (Trans == NULL) ? NONE :
                    StrIndex(Trans[j], TransName, NumStr(TransName));

          NewColName = (TransNum == NONE) ? StrDup(Variable[j]) :
                    StrPaste(2, TransName[TransNum], Variable[j]);

          NewCol = MatCol(&NewData,
                    MatColumnAdd(NewColName, REAL, &NewData));

          AllocFree(NewColName);

          VarMin = (MinCol != NULL) ? MinCol[j] : -REAL_MAX;
          VarMax = (MaxCol != NULL) ? MaxCol[j] :  REAL_MAX;

          /* nInputNA counts cases with NA's already in OldCol, */
          /* or where NA's are produced by out of range.        */
          /* nDomErr counts NA's produced by transformation     */
          /* domain errors.                                     */
          nInputNA = nDomErr = 0;

          /* Count NA's already in input or from out of range. */
          /* Range check not applied in design.                */
          for (i = 0; i < n; i++)
               if (OldCol[i] == NA_REAL || (!DesignJob &&
                         (OldCol[i] < VarMin || OldCol[i] > VarMax)))
               {
                    NewCol[i] = NA_REAL;
                    nInputNA++;
               }
               else
                     NewCol[i] = OldCol[i];

          if (TransNum == LOG)
              /* Apply log10 transformation. */
              SafeLog10(n, NewCol, NewCol, &nInputNA, &nDomErr);

          /*
          Output("%15s%11d%20d%6d\n", Variable[j], nInputNA,
                    nDomErr, n - nInputNA - nDomErr);
          */
     }

     if (ErrNum == OK)
     {
          MatFree(Data);
          MatDup(&NewData, Data);
          D->IsProcessed = YES;
     }
     else
          MatFree(&NewData);

     return ErrNum;
}

/*******************************+++*******************************/
size_t DbnActiveY(void)
/*****************************************************************/
/*   Purpose:  Return the number of y variables to be analyzed.  */
/*                                                               */
/*   Version:  1994 March 9                                      */
/*****************************************************************/
{
     size_t    j, nAnalyze;
     string    *Analyze;

     if ( (Analyze = MatStrColFind(&YDescrip, ANALYZE, NO)) == NULL)
          nAnalyze = MatNumRows(&YDescrip);

     else
     {
          nAnalyze = 0;
          for (j = 0; j < MatNumRows(&YDescrip); j++)
               if (stricmp(Analyze[j], YES_STR) == 0)
                   nAnalyze++;
     }

     return nAnalyze;
}

/*******************************+++*******************************/
size_t DbIndexXY(size_t j)
/*****************************************************************/
/*   Purpose:  Put all cases with no NA's in X or in y variable  */
/*             j in IndexXY.                                     */
/*             On exit, y points to the YMat column for variable */
/*             j, yTrue points to the YTrueMat column, and yName */
/*             is the name.                                      */
/*                                                               */
/*   Returns:  The number of cases with no NA's.                 */
/*                                                               */
/*   96.05.28: yName allocated to fix memory leakage.            */
/*   96.06.09: Y (global) renamed YMat to avoid conflict with y. */
/*             YTrue (global) renamed YTrueMat to avoid conflict */
/*             with yTrue.                                       */
/*                                                               */
/*   Version:  1996.06.09                                        */
/*****************************************************************/
{
     size_t    i, ii, n, TransNum;
     size_t    *Cases;
     string    *Analyze, *Trans, *Variable;

     Variable = MatStrColFind(&YDescrip, VARIABLE,       YES);
     Analyze  = MatStrColFind(&YDescrip, ANALYZE,        NO);
     Trans    = MatStrColFind(&YDescrip, TRANSFORMATION, NO);

     if (Analyze == NULL || stricmp(Analyze[j], YES_STR) == 0)
     {
          /* Variable j to be analyzed. */
          TransNum = (Trans == NULL) ? NONE :
                    StrIndex(Trans[j], TransName, NumStr(TransName));

          if (TransNum == NONE)
               yName = StrReplace(Variable[j], yName);
          else
          {
               AllocFree(yName);
               yName = StrPaste(2, TransName[TransNum],
                         Variable[j]);
          }

          y = MatColFind(&YMat, yName, YES);

          n = DbIndexX();

          IndexXY = AllocSize_t(n, IndexXY);

          for (nCasesXY = 0, ii = 0; ii < n; ii++)
          {
               i = IndexX[ii];
               if (y[i] != NA_REAL)
                    IndexXY[nCasesXY++] = i;
          }

          if (nCasesXY == 0)
          {
               ErrorVar = yName;
               Error("No data cases");
          }

          Cases = MatSize_tColAdd(CASES, &YDescrip);
          Cases[j] = nCasesXY;

          /* True y's for this y variable (if there are any). */
          yTrue = MatColFind(&YTrueMat, yName, NO);
     }
     else
     {
          nCasesXY = 0;
          yName = NULL;
          y     = NULL;
          yTrue = NULL;
     }


     return nCasesXY;
}

/*******************************+++*******************************/
size_t DbIndexX(void)
/*****************************************************************/
/*   Purpose:  Put all cases with no NA's in X in IndexX.        */
/*                                                               */
/*   Returns:  The number of cases with no NA's.                 */
/*                                                               */
/*   Version:  1996.04.14                                        */
/*****************************************************************/
{
     real      *xRow;
     size_t    i, n, nXVars;

     n      = MatNumRows(&X);
     nXVars = MatNumCols(&X);

     IndexX = AllocSize_t(n, IndexX);

     xRow = AllocReal(nXVars, NULL);

     for (nCasesX = 0, i = 0; i < n; i++)
     {
          MatRow(&X, i, xRow);

          if (!VecHasNA(nXVars, xRow))
               IndexX[nCasesX++] = i;
     }

     AllocFree(xRow);

     return nCasesX;
}
