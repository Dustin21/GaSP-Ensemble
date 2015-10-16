/*****************************************************************/
/*   ROUTINES TO EXECUTE VISUALIZATION COMPUTATIONS              */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--6.                     */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"
#include "kriging.h"
#include "alex.h"

extern boolean      ErrorSave;
extern int          ErrorSeverityLevel;
extern string       ErrorVar;

extern boolean      RanErr;

extern LinModel     RegMod;
extern LinModel     SPMod;

extern Matrix       ANOVAPerc;
extern Matrix       MainEff;
extern Matrix       JointEff;
extern Matrix       PredReg;
extern Matrix       SPModMat;
extern Matrix       T;
extern Matrix       X;
extern Matrix       YDescrip;

extern real         InterPerc;
extern real         MainPerc;
extern real         *y;

extern size_t       CorFamNum;
extern size_t       nCasesXY;
extern size_t       *IndexXY;

extern string       yName;


static string       SummaryStats[] = {VARIABLE, TRANSFORMATION,
                         CASES, ANOVA_TOTAL_PERC};

/*******************************+++*******************************/
int Visualize(void)
/*****************************************************************/
/*   Purpose:  Compute effects.                                  */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*                                                               */
/*   1996.03.25: Average and its standard error computed.        */
/*               Averaging w.r.t. groups of variables.           */
/*   1996.04.12: Completed removed.                              */
/*   1996.04.04: X and y include NA's.                           */
/*   1996.04.14: KrigModAlloc/KrigModData.                       */
/*   1996.04.14: KrigModAlloc/KrigModData; DbIndexXY.            */
/*   2009.05.07: Multiple correlation families                   */
/*****************************************************************/
{
     int            ErrNum, ErrReturn;
     KrigingModel   KrigMod;
     Matrix         GroupVarIndex;
     real           *ANOVATot, *Average, *ErrVar, *Perc;
     real           *SEAve, *SPVar;
     size_t         i, j, nNotNA;
     size_t         *GroupSize;

     ErrVar = MatColFind(&YDescrip, ERR_VAR, NO);
     SPVar  = MatColFind(&YDescrip, SP_VAR, YES);

     /* Determine group structure (GroupSize and  */
     /* GroupVarIndex allocated in RegGroupings). */
     RegGroupings(&PredReg, &GroupSize, &GroupVarIndex);

     /* Add columns to the YDescrip matrix. */
     ANOVATot = MatColAdd(ANOVA_TOTAL_PERC,    &YDescrip);
     Average  = MatColAdd(AVERAGE,             &YDescrip);
     SEAve    = MatColAdd(STD_ERR "." AVERAGE, &YDescrip);

     /* Compute effects for each response. */
     ErrReturn = OK;
     ErrorSave = YES;
     for (j = 0; j < MatNumRows(&YDescrip); j++)
     {
          if (DbIndexXY(j) == 0)
               continue;

          /* Delete all previous main and joint effects */
          /* for this y variable.                       */
          RemEffectsRows(yName, &MainEff);
          RemEffectsRows(yName, &JointEff);

          ErrorVar = yName;

          /* Set up kriging model. */
          KrigModAlloc(nCasesXY, MatNumCols(&X), yName, &T, &RegMod,
                    &SPMod, CorFamNum, RanErr, &KrigMod);
          KrigModData(nCasesXY, IndexXY, &X, y, &KrigMod);

          /* SPModMat contains the correlation parameters. */
          ErrNum = KrigModSetUp(&SPModMat, yName, SPVar[j],
                    (ErrVar != NULL) ? ErrVar[j] : 0.0, &KrigMod);

          if (ErrNum == OK)
          {
               Perc = MatColAdd(yName, &ANOVAPerc);

               ErrNum = CompEffects(&KrigMod, yName, &PredReg,
                         GroupSize, &GroupVarIndex, MainPerc, InterPerc,
                         Perc, &Average[j], &SEAve[j]);
          }

          if (ErrNum != OK)
               ErrReturn = ErrNum;

          ANOVATot[j] = 0.0;
          for (nNotNA = 0, i = 0; i < MatNumRows(&ANOVAPerc); i++)
               if (Perc[i] != NA_REAL)
               {
                    nNotNA++;
                    ANOVATot[j] += Perc[i];
               }
          if (nNotNA == 0)
               ANOVATot[j] = NA_REAL;

          KrigModFree(&KrigMod);
     }

     OutputSummary(&YDescrip, NumStr(SummaryStats), SummaryStats);

     AllocFree(GroupSize);
     MatFree(&GroupVarIndex);

     return ErrReturn;
}

/*******************************+++*******************************/
void RemEffectsRows(const string yName, Matrix *Eff)
/*****************************************************************/
/*   Purpose:  Remove any rows in Eff corresponding to           */
/*             y-variable yName.                                 */
/*                                                               */
/*   96.04.04: Call to DbSelectCases replaced with code here.    */
/*                                                               */
/*   Version:  1996.04.04                                        */
/*****************************************************************/
{
     Matrix    NewEff;
     size_t    i, n, nIn;
     size_t    *In;
     string    *yNameCol;

     if ( (n = MatNumRows(Eff)) == 0)
          return;

     In = AllocSize_t(n, NULL);

     yNameCol = MatStrColFind(Eff, VARIABLE ".y", YES);

     for (nIn = 0, i = 0; i < n; i++)
          if (stricmp(yNameCol[i], yName) != 0)
               In[nIn++] = i;

     if (nIn < n)
     {
          MatDupIndex(nIn, In, Eff, &NewEff);
          MatFree(Eff);
          *Eff = NewEff;
     }

     AllocFree(In);
}

/*******************************+++*******************************/
int CompEffects(KrigingModel *KrigMod, const string yName,
     const Matrix *PredReg, const size_t *GroupSize,
     const Matrix *GroupVarIndex, real MainPerc, real InterPerc,
     real *Perc, real *Average, real *SEAve)
/*****************************************************************/
/*   Purpose:    Compute ANOVA percentage contributions and put  */
/*               important effects in MAIN_EFF and JOINT_EFF     */
/*               matrices.                                       */
/*                                                               */
/*   Returns:    OK or an error number.                          */
/*                                                               */
/*   Comment:    Calling routine must allocate space for Perc,   */
/*               Average, and SEAve.                             */
/*                                                               */
/*   1996.02.12: Main effects generated for fixed x variables;   */
/*               checks for zero SSTot changed.                  */
/*   1996.02.19: AllocFree(SE) added.                            */
/*   1996.03.20: Average and its standard error computed.        */
/*   1996.03.25: Averaging w.r.t. groups of variables.           */
/*   1996.04.12: Temporary output showing progress.              */
/*   1996.08.31: Plotting coordinates for grouped variables      */
/*               appended.                                       */
/*   1996.02.17: Weights for integration in variance of an       */
/*               effect from RegLevelWt.                         */
/*                                                               */
/*   Version:    1996.02.17                                      */
/*****************************************************************/
{
     boolean   *ActiveGroup;
     int       ErrNum, SevSave;
     Matrix    IndexSP;
     real      RAve, SSTot, VarEff, VarEffRow;
     real      *Eff = NULL, *fAve, *rAve;
     real      *SE = NULL;
     size_t    c, kSP, i, i1, i2, j, jj, j1, j2, m, m1, m2;
     size_t    n, nGroups, x1Index, x2Index;
     size_t    IndexGroup[2];
     size_t    *nSPTerms, *xIndex;

     kSP     = ModDF(KrigSPMod(KrigMod));
     n       = MatNumRows(KrigChol(KrigMod));
     nGroups = MatNumCols(GroupVarIndex);

     /* Workspace in KrigMod. */
     fAve = KrigMod->fRow;
     rAve = KrigMod->r;

     /* Allocations. */
     ActiveGroup = (boolean *) AllocGeneric(nGroups,
               sizeof(boolean), NULL);
     nSPTerms = AllocSize_t(nGroups, NULL);
     MatAllocate(kSP, nGroups, RECT, SIZE_T, NULL, NO, &IndexSP);

     for (j = 0; j < nGroups; j++)
          nSPTerms[j] = KrigSPActiveTerms(KrigMod, GroupSize[j],
                         MatSize_tCol(GroupVarIndex, j),
                         MatSize_tCol(&IndexSP, j));

     /* Average predictor w.r.t. all x variables. */
     AvePred(KrigMod, PredReg, 0, NULL, GroupSize, GroupVarIndex,
               nSPTerms, &IndexSP, fAve, rAve, &RAve);

     if ( (ErrNum = KrigYHatSE(KrigMod, RAve, fAve, rAve, Average,
               SEAve)) == OK)
     {
          /* Subtract average prediction from data,   */
          /* so that predictions will have average 0. */
          VecAddScalar(-(*Average), n, KrigY(KrigMod));

          /* Get new decompositions and total sum of squares */
          /* for all variables.                              */
          KrigCorMat(0, NULL, KrigMod);
          if ( (ErrNum = KrigDecompose(KrigMod)) == OK)
                ErrNum = CompSSTot(KrigMod, PredReg, GroupSize,
                          GroupVarIndex, nSPTerms, &IndexSP, &SSTot);
     }

     if (ErrNum == OK)
     {
          /* Determine active groups. */
          for (j = 0; j < nGroups; j++)
          {
               xIndex = MatSize_tCol(GroupVarIndex, j);
               for (jj = 0; jj < GroupSize[j]; jj++)
                    if (KrigIsXActive(KrigMod, xIndex[jj]))
                         /* Active variable found. */
                         break;
               ActiveGroup[j] = (jj < GroupSize[j]);
          }
     }

     if (SSTot < sqrt(EPSILON))
     {
          SSTot = 0.0;
          SevSave = ErrorSeverityLevel;
          ErrorSeverityLevel = SEV_WARNING;
          Error("No variation in predictor");
          ErrorSeverityLevel = SevSave;
     }

     /* Main effects and contributions. */
     for (j = 0; j < nGroups && ErrNum == OK; j++)
     {
          OutputTemp("Variable: %s  Main effect: %d", yName, j + 1);

          if (!ActiveGroup[j] && MainPerc > 0.0)
          {
               Perc[j] = 0.0;
               continue;
          }

          x1Index = MatSize_tElem(GroupVarIndex, 0, j);
          m = RegNumLevels(PredReg, x1Index);

          Eff = AllocReal(m, Eff);
          SE  = AllocReal(m, SE);

          /* Compute main effect of group j. */
          AnyEffect(KrigMod, PredReg, 1, &j, GroupSize,
                    GroupVarIndex, nSPTerms, &IndexSP, Eff, SE);

          if (SSTot == 0.0)
               Perc[j] = NA_REAL;
          else
          {
               for (VarEff = 0.0, i = 0; i < m; i++)
                    VarEff += RegLevelWt(PredReg, x1Index, i)
                              * Eff[i] * Eff[i];
                    
               Perc[j] = VarEff / SSTot * 100.0;
               if (Perc[j] < sqrt(EPSILON))
                    Perc[j] = 0.0;
          }

          /* Previously required GroupSize[j] == 1. */
          if (Perc[j] != NA_REAL && Perc[j] >= MainPerc)
          {
               /* Generate plotting coordinates. */

               /* Add average prediction back in. */
               VecAddScalar(*Average, m, Eff);

               /* Append effect to MAIN_EFF. */
               AppendEffect(yName, 1, &j, PredReg, GroupSize,
                         GroupVarIndex, Eff, SE, &MainEff);
          }
     }

     /* Joint effects and interaction contributions. */
     for (c = nGroups, j1 = 0; j1 < nGroups - 1; j1++)
     {
          IndexGroup[0] = j1;
          x1Index = MatSize_tElem(GroupVarIndex, 0, j1);
          m1 = RegNumLevels(PredReg, x1Index);
          
          for (j2 = j1 + 1; j2 < nGroups; j2++, c++)
          {
               OutputTemp("Variable: %s  Joint effect: %d", yName,
                         c - nGroups + 1);

               IndexGroup[1] = j2;

               if ( (!ActiveGroup[j1] || !ActiveGroup[j2]) &&
                         InterPerc > 0.0)
               {
                    Perc[c] = 0.0;
                    continue;
               }

               x2Index = MatSize_tElem(GroupVarIndex, 0, j2);
               m2 = RegNumLevels(PredReg, x2Index);
                  
               Eff = AllocReal(m1 * m2, Eff);
               SE  = AllocReal(m1 * m2, SE);

               /* Compute joint effect of variables j1 and j2. */
               AnyEffect(KrigMod, PredReg, 2, IndexGroup, GroupSize,
                         GroupVarIndex, nSPTerms, &IndexSP, Eff, SE);

               /* Contribution of the *interaction* effect. */
               if (SSTot == 0.0)
                    Perc[c] = NA_REAL;
               else
               {
                    for (VarEff = 0.0, i1 = 0; i1 < m1; i1++)
                    {
                         for (VarEffRow = 0.0, i2 = 0; i2 < m2; i2++)
                              VarEffRow += RegLevelWt(PredReg,
                                        x2Index, i2)
                                        * Eff[i1 * m2 + i2]
                                        * Eff[i1 * m2 + i2];
                         VarEff += RegLevelWt(PredReg, x1Index, i1)
                                   * VarEffRow;
                    }
                                                  
                    Perc[c] = VarEff / SSTot * 100.0
                              - ((Perc[j1] != NA_REAL) ? Perc[j1]
                              : 0.0)
                              - ((Perc[j2] != NA_REAL) ? Perc[j2]
                              : 0.0);
                    if (Perc[c] < sqrt(EPSILON))
                         Perc[c] = 0.0;
               }

               /* Previously required:                         */
               /* GroupSize[j1] == 1 && and GroupSize[j2] == 1 */
               if (Perc[c] != NA_REAL && Perc[c] >= InterPerc)
               {
                    /* Generate plotting coordinates. */

                    /* Add average back in. */
                    VecAddScalar(*Average, m1 * m2, Eff);

                    /* Append joint effect to JOINT_EFF. */
                    AppendEffect(yName, 2, IndexGroup, PredReg,
                              GroupSize, GroupVarIndex, Eff, SE,
                              &JointEff);
               }
          }
     }

     OutputTemp("");

     AllocFree(ActiveGroup);
     AllocFree(nSPTerms);
     MatFree(&IndexSP);
     AllocFree(Eff);
     AllocFree(SE);

     return ErrNum;
}

/*******************************+++*******************************/
void AvePred(KrigingModel *KrigMod, const Matrix *PredReg,
     size_t nGroups, const size_t *IndexGroup,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     real *fAve, real *rAve, real *RAve)
/*****************************************************************/
/* Purpose:    Average f and r w.r.t. the groups *not* in        */
/*             IndexGroup.                                       */
/*                                                               */
/* 1996.03.25: Averaging w.r.t. groups of variables.             */
/* 1996.04.10: RAve calculation changed to deal correctly        */
/*             with random error.                                */
/* 1996.02.17: Wt from RegLevelWt instead of 1/m.                */
/* 1996.02.18: All elements in R multiplied by                   */
/*             KrigMod->SPVarProp, including diagonals, to       */
/*             predict average of f(x) beta + Z *without*        */
/*             epsilon.                                          */
/*                                                               */
/* 1996.02.18: Created?                                          */
/* 2009.05.13: KrigCorVec arguments changed                      */
/*****************************************************************/
{
     LinModel  *RegMod, *SPMod;
     Matrix    GAve;
     real      wRw, wRwj, SPVarPropSave;
     real      *f, *fj, *g, *r, *rj, *Rj = NULL, *Wt = NULL, *xRow;
     size_t    i, j, kReg, kSP, m, n;
     size_t    *IndexSPCol, *xIndex;

     n      = MatNumRows(KrigChol(KrigMod));
     RegMod = KrigRegMod(KrigMod);
     SPMod  = KrigSPMod(KrigMod);
     kReg   = ModDF(RegMod);
     kSP    = ModDF(SPMod);

     /* Allocations. */
     f  = AllocReal(kReg, NULL);
     fj = AllocReal(kReg, NULL);
     g  = AllocReal(kSP,  NULL);
     r  = AllocReal(n,    NULL);
     rj = AllocReal(n,    NULL);
     MatInit(RECT, REAL, NO, &GAve);

     /* Workspace in KrigMod. */
     xRow = KrigMod->xRow;

     VecInit(1.0, kReg, fAve);
     VecInit(KrigMod->SPVarProp, n, rAve);

     wRw = KrigMod->SPVarProp;

     /* SPVarProp is applied only once. */
     SPVarPropSave = KrigMod->SPVarProp;
     KrigMod->SPVarProp = 1.0;

     for (j = 0; j < MatNumCols(GroupVarIndex); j++)
     {
          if (VecSize_tIndex(j, nGroups, IndexGroup) != INDEX_ERR)
               continue;

          xIndex = MatSize_tCol(GroupVarIndex, j);
          IndexSPCol = MatSize_tCol(IndexSP, j);

          m = RegNumLevels(PredReg, xIndex[0]);
          CodeCheck(m > 0);
          Rj = AllocReal(m, Rj);
          Wt = AllocReal(m, Wt);
          MatReAlloc(m, kSP, &GAve);

          VecInit(0.0, kReg, fj);
          VecInit(0.0, n,    rj);

          for (wRwj = 0.0, i = 0; i < m; i++)
          {
               fgrGroup(KrigMod, PredReg, GroupSize[j], xIndex, i,
                        nSPTerms[j], IndexSPCol, xRow, f, g, r);

               Wt[i] = RegLevelWt(PredReg, xIndex[0], i);
               
               VecAddVec(Wt[i], f, kReg, fj);
               VecAddVec(Wt[i], r, n,    rj);

               MatRowPut(g, i, &GAve);
               KrigCorVec(g, &GAve, i, nSPTerms[j], IndexSPCol, YES,
                         KrigMod, Rj);
               wRwj += Wt[i] * (Wt[i] + 2.0 * DotProd(Wt, Rj, i));
          }

          VecMultVec(fj, kReg, fAve);
          VecMultVec(rj, n,    rAve);
          
          wRw *= wRwj;
     }

     *RAve = wRw;

     KrigMod->SPVarProp = SPVarPropSave;

     AllocFree(f);
     AllocFree(fj);
     AllocFree(g);
     AllocFree(r);
     AllocFree(Rj);
     AllocFree(rj);
     AllocFree(Wt);

     MatFree(&GAve);
}

/*******************************+++*******************************/
int CompSSTot(KrigingModel *KrigMod, const Matrix *PredReg,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     real *SSTot)
/*****************************************************************/
/* Purpose:    Compute SS(Total) for predictor.                  */
/*                                                               */
/* Returns:    OK or an error number.                            */
/*                                                               */
/* 1996.01.17: MatRow and MatRowPut replaced TriRow and          */
/*             TriRowPut.                                        */
/* 1996.02.20: TriCholesky returns rank.                         */
/* 1996.03.25: Averaging w.r.t. groups of variables.             */
/* 1999.03.29: Cholesky decomposition of frfr replaced by eigen  */
/*             decomposition.                                    */
/*                                                               */
/* Version:    1999.03.29                                        */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    frfr, frfrj;
     real      a;
     real      *eVal, *v;
     size_t    j, k, n;

     ErrNum = OK;

     n = MatNumRows(KrigChol(KrigMod));
     k = ModDF(KrigRegMod(KrigMod));

     /* Allocations:                                            */
     /* frfr is RECT because it is overwritten by eigenvectors. */
     MatAlloc(k + n, k + n, RECT, &frfr);
     MatAlloc(k + n, k + n, SYM,  &frfrj);

     /* Workspace in KrigMod. */
     eVal = KrigMod->fr;

     /* frfr must be symmetric for frfrAve. */
     MatPutShape(&frfr, SYM);
     frfrAve(KrigMod, PredReg, GroupSize, GroupVarIndex, nSPTerms,
               IndexSP, &frfrj, &frfr);
     MatPutShape(&frfr, RECT);          

     /* Eigen decomposition, overwriting frfr with eigenvectors. */
     if ( (ErrNum = MatEig(YES, &frfr, eVal, &frfr)) != OK)
          Error("Eigen decomposition of averaging moment matrix failed.");

     for (*SSTot = 0.0, j = 0; j < k + n && ErrNum == OK; j++)
     {

          if (eVal[j] < EPSILON * eVal[0])
               break;
               
          /* Eigenvector j is row j of V'. */
          v = MatCol(&frfr, j);

          /* Overwrite first k elements of v with solution. */
          if ( (ErrNum = TriForSolve(KrigR(KrigMod), v, 0, v))
                    != OK)
               Error("Ill-conditioned expanded-design matrix.\n");

          /* Overwrite next n elements of v with solution. */
          else if ( (ErrNum = TriForSolve(KrigChol(KrigMod),
                    v + k, 0, v + k)) != OK)
               Error("Ill-conditioned correlation matrix.\n");

          else
          {
               a = VecDotProd(k, v, KrigMod->RBeta)
                         + VecDotProd(n, v + k, KrigMod->ResTilde);
               *SSTot += eVal[j] * a * a;
          }     
     }                    

     MatFree(&frfr);
     MatFree(&frfrj);

     return ErrNum;
}

/*******************************+++*******************************/
void AnyEffect(KrigingModel *KrigMod, const Matrix *PredReg,
     size_t nGroups, const size_t *IndexGroup,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     real *Eff, real *SE)
/*****************************************************************/
/*   Purpose:  Compute any (arbitrary-degree) effect.            */
/*                                                               */
/*   Comments: Calling routine must allocate space for Eff and   */
/*             SE.                                               */
/*                                                               */
/*   96.03.25: Averaging w.r.t. groups of variables.             */
/*                                                               */
/*   Version:  1996.04.02                                        */
/*****************************************************************/
{
     real      RAve, SPVarPropSave;
     real      *fAve, *f, *fj, *g, *r, *rj, *rAve, *xRow;
     size_t    i, j, jj, k, n;
     size_t    *Level, *nLevels, *xIndex;

     n = MatNumRows(KrigChol(KrigMod));
     k = ModDF(KrigRegMod(KrigMod));

     /* Workspace in KrigMod. */
     fAve = KrigMod->fRow;
     g    = KrigMod->gRow;
     rAve = KrigMod->r;
     xRow = KrigMod->xRow;

     /* Allocations. */
     f  = AllocReal(k, NULL);
     fj = AllocReal(k, NULL);
     r  = AllocReal(n, NULL);
     rj = AllocReal(n, NULL);
     Level   = AllocSize_t(nGroups, NULL);
     nLevels = AllocSize_t(nGroups, NULL);

     AvePred(KrigMod, PredReg, nGroups, IndexGroup, GroupSize,
               GroupVarIndex, nSPTerms, IndexSP, fAve, rAve, &RAve);

     /* SPVarProp was applied in AvePred. */
     SPVarPropSave = KrigMod->SPVarProp;
     KrigMod->SPVarProp = 1.0;

     /* Adjust averages for variables in the effect. */

     /* Set up initial group levels. */
     for (j = 0; j < nGroups; j++)
     {
          Level[j] = 0;
          nLevels[j] = RegNumLevels(PredReg,
                    MatSize_tElem(GroupVarIndex, 0, IndexGroup[j]));
     }

     /* For each level combination. */
     i = 0;
     do
     {
          VecCopy(fAve, k, f);
          VecCopy(rAve, n, r);

          /* Adjust f and r for combination i */
          /* of variables in effect.          */
          for (jj = 0; jj < nGroups; jj++)
          {
               j = IndexGroup[jj];

               xIndex = MatSize_tCol(GroupVarIndex, j);

               fgrGroup(KrigMod, PredReg, GroupSize[j], xIndex,
                         Level[jj], nSPTerms[j],
                         MatSize_tCol(IndexSP, j), xRow, fj, g, rj);

               VecMultVec(fj, k, f);
               VecMultVec(rj, n, r);
          }

          /* Need to transform r with T? */

          KrigYHatSE(KrigMod, RAve, f, r, &Eff[i], &SE[i]);

          i++;
     } while (LevelLex(nGroups, nLevels, Level) != ALL_DONE);

     KrigMod->SPVarProp = SPVarPropSave;

     AllocFree(f);
     AllocFree(fj);
     AllocFree(r);
     AllocFree(rj);
     AllocFree(Level);
     AllocFree(nLevels);

     return;
}

/*******************************+++*******************************/
void AppendEffect(const string yName, size_t DegreeEff,
          const size_t *IndexGroup, const Matrix *PredReg,
          const size_t *GroupSize, const Matrix *GroupVarIndex,
          real *Eff, real *SE, Matrix *EffMat)
/*****************************************************************/
/*   Purpose:  Append an effect to the bottom of EffMat.         */
/*                                                               */
/*   96.03.25: Generates combinations instead of XEffect.        */
/*   96.08.31: Effect can include grouped variables.             */
/*                                                               */
/*   Version:  1996.08.31                                        */
/*****************************************************************/
{
     real      l;
     size_t    i, j, jj, nRowsNew, nRowsOld, xIndex;
     size_t    *Level, *nLevels;
     string    s;

     nRowsOld = MatNumRows(EffMat);

     /* Allocations. */
     Level   = AllocSize_t(DegreeEff, NULL);
     nLevels = AllocSize_t(DegreeEff, NULL);

     /* Set up initial variable levels. */
     for (nRowsNew = 1, j = 0; j < DegreeEff; j++)
     {
          Level[j] = 0;
          nLevels[j] = RegNumLevels(PredReg,
                    MatSize_tElem(GroupVarIndex, 0, IndexGroup[j]));
          nRowsNew *= nLevels[j];
     }

     MatReAlloc(nRowsOld + nRowsNew, MatNumCols(EffMat), EffMat);

     /* Put x-variable names in first new row of EffMat. */
     /* Names will be copied for further rows.           */
     for (jj = 0; jj < DegreeEff; jj++)
     {
          j = IndexGroup[jj];
          xIndex = MatSize_tElem(GroupVarIndex, 0, j);

          if (GroupSize[j] == 1)
               MatPutStrElem(EffMat, nRowsOld, jj,
                         RegVar(PredReg, xIndex));
          else
          {
               s = StrPaste(2, GROUP, StrFromSize_t(
                         RegCandGroup(PredReg, xIndex)));
               MatPutStrElem(EffMat, nRowsOld, jj, s);
               AllocFree(s);
          }
     }

     i = 0;
     do
     {
          for (jj = 0; jj < DegreeEff; jj++)
          {
               if (i > 0)
                    MatPutStrElem(EffMat, nRowsOld + i, jj,
                              MatStrElem(EffMat, nRowsOld + i - 1, jj));

               j = IndexGroup[jj];
               xIndex = MatSize_tElem(GroupVarIndex, 0, j);

               l = (GroupSize[j] == 1) ?
                         RegLevel(PredReg, xIndex, Level[jj]) :
                         (real) (Level[jj] + 1);
               MatPutElem(EffMat, nRowsOld + i, DegreeEff + jj + 1,
                         l);
          }

          MatPutStrElem(EffMat, nRowsOld + i, DegreeEff, yName);
          MatPutElem(EffMat, nRowsOld + i, 2 * DegreeEff + 1, Eff[i]);
          MatPutElem(EffMat, nRowsOld + i, 2 * DegreeEff + 2, SE[i]);
          i++;

     } while (LevelLex(DegreeEff, nLevels, Level) != ALL_DONE);

     AllocFree(Level);
     AllocFree(nLevels);

     return;
}

