/*****************************************************************/
/*   BASIC ROUTINES FOR THE MODEL                                */
/*             Y = REGRESSION + STOCHASTIC PROCESS               */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--2009.                  */
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
#include "kriging.h"
#include "alex.h"

/*******************************+++*******************************/
void KrigModAlloc(size_t nCases, size_t nXVars, const string yName,
     const Matrix *T, const LinModel *RegMod,
     const LinModel *SPMod, size_t CorFam, boolean RanErr,
     KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Initialize KrigMod, and allocate F, G, etc.       */
/*                                                               */
/*   1996.04.14: Some code moved to KrigModData.                 */
/*   2009.05.14: Multiple correlation families                   */
/*****************************************************************/
{
     size_t    kReg, kSP;

     KrigMod->Y = AllocReal(nCases, NULL);

     KrigMod->yName = yName;

     KrigMod->T = T;

     KrigMod->RegMod = RegMod;
     KrigMod->SPMod  = SPMod;
     KrigMod->CorFam = CorFam;
     KrigMod->RanErr = RanErr;

     kReg = ModDF(RegMod);
     kSP  = ModDF(SPMod);

     MatAlloc(nCases, kReg, RECT, KrigF(KrigMod));
     MatAlloc(nCases, kSP,  RECT, KrigG(KrigMod));

     MatAllocate(nCases, kSP, RECT, SIZE_T, NULL, NO,
               KrigSteps(KrigMod));
     KrigMod->MaxSteps = AllocSize_t(kSP, NULL);
     MatAlloc(nCases, kSP, RECT, KrigDist(KrigMod));

     CorParAlloc(CorFam, kSP, ModTermNames(SPMod), KrigCorPar(KrigMod));

     MatAlloc(nCases, nCases, UP_TRIANG, KrigChol(KrigMod));
     MatAlloc(nCases, kReg,   RECT,      KrigQ(KrigMod));
     MatAlloc(kReg,   kReg,   UP_TRIANG, KrigR(KrigMod));

     KrigMod->RBeta    = AllocReal(kReg, NULL);
     KrigMod->Beta     = AllocReal(kReg, NULL);
     KrigMod->ResTilde = AllocReal(nCases, NULL);

     KrigMod->xRow = AllocReal(nXVars, NULL);
     KrigMod->fRow = AllocReal(kReg, NULL);
     KrigMod->fr   = AllocReal(kReg + nCases, NULL);
     KrigMod->gRow = AllocReal(kSP, NULL);
     KrigMod->r    = AllocReal(nCases, NULL);
     KrigMod->w1   = AllocReal(nCases, NULL);
     KrigMod->w2   = AllocReal(nCases, NULL);

     /* Further initializations, etc. for T. */
     KrigModAllocT(KrigMod);

     return;
}

/*******************************+++*******************************/
void KrigModAllocT(KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Further initializations and allocations from      */
/*             transformations T.                                */
/*                                                               */
/*   96.04.04: KrigMod->Y reallocated instead of allocated.      */
/*                                                               */
/*   Version:  1996.04.04                                        */
/*****************************************************************/
{
     Matrix    *C, *Chol, *F, *T;
     real      *y;
     size_t    j, k, t;

     T = KrigT(KrigMod);

     if (T == NULL || (t = MatNumCols(T)) == 0)
          /* No transformations. */
          return;

     /* y becomes T'y. */
     y = KrigY(KrigMod);
     KrigMod->Y = AllocReal(t, KrigMod->Y);
     MatVec(T, y, KrigMod->Y);

     /* F becomes T'F. */
     F = KrigF(KrigMod);
     k = MatNumCols(F);
     for (j = 0; j < k; j++)
     {
          MatVec(T, MatCol(F, j), KrigMod->w1);
          VecCopy(KrigMod->w1, t, MatCol(F, j));
     }
     MatReAlloc(t, k, F);
     MatReAlloc(t, k, KrigQ(KrigMod));

     /* Chol was allocated as n x n: use the space for C. */
     Chol = KrigChol(KrigMod);
     C    = KrigC(KrigMod);
     *C = *Chol;

     /* Chol is t x t. */
     MatAlloc(t, t, UP_TRIANG, Chol);
}

/*******************************+++*******************************/
void KrigModFree(KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Free kriging model.                               */
/*                                                               */
/*   96.04.04: KrigMod->Y freed.                                 */
/*                                                               */
/*   Version:  1996.04.04                                        */
/*****************************************************************/
{
     AllocFree(KrigY(KrigMod));

     MatFree(KrigF(KrigMod));
     MatFree(KrigG(KrigMod));

     MatFree(KrigSteps(KrigMod));
     AllocFree(KrigMod->MaxSteps);
     MatFree(KrigDist(KrigMod));

     MatFree(KrigCorPar(KrigMod));

     MatFree(KrigChol(KrigMod));
     MatFree(KrigQ(KrigMod));
     MatFree(KrigR(KrigMod));

     AllocFree(KrigMod->RBeta);
     AllocFree(KrigMod->Beta);
     AllocFree(KrigMod->ResTilde);

     AllocFree(KrigMod->xRow);
     AllocFree(KrigMod->fRow);
     AllocFree(KrigMod->fr);
     AllocFree(KrigMod->gRow);
     AllocFree(KrigMod->r);
     AllocFree(KrigMod->w1);
     AllocFree(KrigMod->w2);

     KrigModFreeT(KrigMod);
}

/*******************************+++*******************************/
void KrigModFreeT(KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Free space allocated because of transformations T.*/
/*                                                               */
/*   96.04.04: KrigMod->Y not freed here.                        */
/*                                                               */
/*   Version:  1996.04.04                                        */
/*****************************************************************/
{
     Matrix    *T;

     T = KrigT(KrigMod);

     if (T != NULL && MatNumCols(T) > 0)
          MatFree(KrigC(KrigMod));

     return;
}

/*******************************+++*******************************/
void KrigModData(size_t nCases, const size_t *RowIndex,
     const Matrix *X, const real *y, KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Set up y, F, G, and call KrigGSpacing.            */
/*                                                               */
/*   Version:  1996.04.14                                        */
/*****************************************************************/
{
     /* CritAMSE etc. do not have y at design stage. */
     if (y != NULL)
          VecCopyIndex(nCases, RowIndex, y, NULL, KrigY(KrigMod));

     ModFMatRowIndex(KrigRegMod(KrigMod), nCases, RowIndex, X,
               KrigF(KrigMod));
     ModFMatRowIndex(KrigSPMod(KrigMod),  nCases, RowIndex, X,
               KrigG(KrigMod));

     KrigGSpacing(KrigMod);

     return;
}

/*******************************+++*******************************/
int KrigModSetUp(const Matrix *CorPar, const string yName,
     real SPVar, real ErrVar, KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Set up parameters and decompositions.             */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  This could be eliminated: CorParSetUp does most   */
/*             most of the work.                                 */
/*                                                               */
/*   Version:  1996.01.20                                        */
/*****************************************************************/
{
     int  ErrNum;

     ErrNum = CorParSetUp(CorPar, yName, SPVar, ErrVar, KrigMod);

     if (ErrNum == OK)
     {
          KrigCorMat(0, NULL, KrigMod);
          ErrNum = KrigDecompose(KrigMod);
     }

     return ErrNum;
}

/*******************************+++*******************************/
void KrigGSpacing(KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Set up Steps and MaxSteps.                        */
/*                                                               */
/*   Version:  1995 January 16                                   */
/*****************************************************************/
{
     Matrix    *G;
     real      Eps, Gap, s, ss, StepLen;
     real      *DistCol, *GCol, *r;
     size_t    i, j, n;
     size_t    *StepsCol, *MaxSteps;

     Eps = sqrt(EPSILON);

     G = KrigG(KrigMod);
     n = MatNumRows(G);

     MaxSteps = KrigMod->MaxSteps;

     /* Workspace. */
     r = KrigMod->r;

     for (j = 0; j < MatNumCols(G); j++)
     {
          GCol = MatCol(G, j);
          VecCopy(GCol, n, r);
          QuickReal(n, r);

          /* Find StepLen, the minimum nonzero gap. */
          for (StepLen = r[n-1] - r[0], i = 1; i < n; i++)
               if ( (Gap = r[i] - r[i-1]) > 0.0)
                    StepLen = min(StepLen, Gap);

          MaxSteps[j] = 0;
          if (StepLen == 0.0 ||
                    StepLen < (r[n-1] - r[0]) / (n - 1) - Eps)
               continue;

          StepsCol = MatSize_tCol(KrigSteps(KrigMod), j);
          for (i = 0; i < n; i++)
          {
               s  = (GCol[i] - r[0]) / StepLen;
               ss = floor(s + Eps);
               if (ApproxEq(s, ss, Eps, 0.0))
               {
                     StepsCol[i] = (size_t) ss;
                     MaxSteps[j] = max(StepsCol[i], MaxSteps[j]);
               }
               else
               {
                     MaxSteps[j] = 0;
                     break;
               }
          }

          if (MaxSteps[j] > 0)
          {
               DistCol = MatCol(KrigDist(KrigMod), j);
               for (i = 0; i < MaxSteps[j]; i++)
                    DistCol[i] = (i + 1) * StepLen;
          }
     }
}

/*******************************+++*******************************/
void KrigCorMat
(
     size_t       nActive,   /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     KrigingModel *KrigMod
)
/*****************************************************************/
/*   Purpose:  Put the correlation matrix into Chol.             */
/*                                                               */
/*   Version:  1995 February 14                                  */
/*****************************************************************/
{
     Matrix    *C, *Chol, *T;
     real      *w1, *w2;
     size_t    i, j, k, n, t;

     T = KrigT(KrigMod);

     if (T == NULL || (t = MatNumCols(T)) == 0)
          /* Correlation matrix can go directly into Chol. */
          KrigCorC(nActive, Active, KrigMod, KrigChol(KrigMod));

     else
     {
          /* Put correlation matrix for original responses in C. */
          KrigCorC(nActive, Active, KrigMod, KrigC(KrigMod));

          /* Chol = T'CT is symmetric t x t. */
          C    = KrigC(KrigMod);
          Chol = KrigChol(KrigMod);
          n    = MatNumRows(C);
          w1   = KrigMod->w1;
          w2   = KrigMod->w2;
          for (j = 0; j < t; j++)
          {
               for (k = 0; k < n; k++)
               {
                    /* Row k of C equals column k: load into w1. */
                    MatSymCol(C, k, w1);

                    /* Row k of C * column j of T. */
                    w2[k] = DotProd(w1, MatCol(T, j), n);
               }

               for (i = 0; i <= j; i++)
                    /* Row i of T' * w2. */
                    MatPutElem(Chol, i, j, DotProd(MatCol(T, i),
                              w2, n));
          }
     }

     return;
}

/*******************************+++*******************************/
void KrigCorC
(
     size_t       nActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     KrigingModel *KrigMod,
     Matrix       *C
)
/*****************************************************************/
/*   Purpose:  Put the correlation matrix for the original       */
/*             responses into C (which could be a member of      */
/*             KrigMod).                                         */
/*                                                               */
/* 1995.07.27: Created?                                          */
/* 2009.05.14: KrigCorVec replaces PECor, and CorParIsActive     */
/*             replaces PEIsActive (multiple correlation         */
/*             families)                                         */                       
/*****************************************************************/
{
     Matrix    *CorPar;
     real      *Cor, *CCol, *gRow;
     size_t    i, j, k, kk, n, NumActiveIrreg, StepsDiff;
     size_t    *ActiveIrreg, *MaxSteps, *StepsCol;

     CorPar = KrigCorPar(KrigMod);
     gRow   = KrigMod->gRow;
     n      = MatNumRows(KrigG(KrigMod));

     MaxSteps = KrigMod->MaxSteps;

     if (Active == NULL)
          nActive = MatNumCols(KrigG(KrigMod));

     /* Allocations. */
     ActiveIrreg = AllocSize_t(nActive, NULL);

     /* Workspace of length n. */
     Cor = KrigMod->w1;

     /* ActiveIrreg indexes the G columns with irregular spacing. */
     for (NumActiveIrreg = 0, kk = 0; kk < nActive; kk++)
     {
          k = (Active == NULL) ? kk : Active[kk];
          if (MaxSteps[k] == 0)
               /* Irregular spacing. */
               ActiveIrreg[NumActiveIrreg++] = k;
     }

     MatPutElem(C, 0, 0, 1.0);

     /* Irregularly spaced G columns. */
     for (j = 1; j < n; j++)
     {
          MatRow(KrigG(KrigMod), j, gRow);

          /* Only have j correlations in column j. */
          KrigCorVec(gRow, KrigG(KrigMod), j, NumActiveIrreg,
                    ActiveIrreg, YES, KrigMod, MatCol(C, j));
          MatPutElem(C, j, j, 1.0);
     }

     /* Even if there are no irregularly-spaced columns, */
     /* C has now been initialized to 1's, and SPVarProp */
     /* has been applied.                                */

     /* For the regularly spaced columns of G, distances  */
     /* are from zero (min distance has been subtracted). */
     VecInit(0.0, MatNumCols(KrigG(KrigMod)), gRow);

     for (kk = 0; kk < nActive; kk++)
     {
          k = (Active == NULL) ? kk : Active[kk];

          if (MaxSteps[k] == 0 ||
                    !CorParIsActive(KrigCorFam(KrigMod), CorPar, k))
               /* Irregularly spaced or inactive column. */
               continue;

          /* Compute the set of possible correlations. */
          Cor[0] = 1.0;
          KrigCorVec(gRow, KrigDist(KrigMod), MaxSteps[k], 1,
                    &k, NO, KrigMod, Cor + 1);
 
          /* Put appropriate correlations in correlation matrix. */
          StepsCol = MatSize_tCol(KrigSteps(KrigMod), k);
          for (j = 1; j < n; j++)
          {
               CCol = MatCol(C, j);
               /* Only have j correlations in column j. */
               for (i = 0; i < j; i++)
               {
                    StepsDiff = max(StepsCol[i], StepsCol[j])
                              - min(StepsCol[i], StepsCol[j]);
                    CCol[i] *= Cor[StepsDiff];
               }
          }
     }

     AllocFree(ActiveIrreg);

     return;
}

/*******************************+++*******************************/
int KrigDecompose(KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  (Re-)compute decompositions for kriging model.    */
/*                                                               */
/*   Return:   NUMERIC_ERR if C or Inverse(Chol') * F are not    */
/*                         full rank;                            */
/*             OK          otherwise.                            */
/*                                                               */
/*   Comment:  Chol must already hold the correlation matrix and */
/*             is overwritten.                                   */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 18                                   */
/*****************************************************************/
{
     Matrix    *Chol, *F, *Q, *R;
     real      *Beta, *RBeta, *ResTilde, *Y;

     Y = KrigMod->Y;
     F = KrigF(KrigMod);

     Chol = KrigChol(KrigMod);

     Q = KrigQ(KrigMod);
     R = KrigR(KrigMod);

     RBeta    = KrigMod->RBeta;
     ResTilde = KrigMod->ResTilde;
     Beta     = KrigMod->Beta;

     /* Overwrite correlation matrix with Cholesky decomposition. */
     if (TriCholesky(Chol, 0, Chol) != OK)
     {
          Error("Ill-conditioned Cholesky factor.\n");
          return NUMERIC_ERR;
     }

     /* Put FTilde in Q, and YTilde in ResTilde. */
     /* As TriCholesky did not return an error, KrigSolve cannot. */
     CodeCheck(KrigSolve(Chol, F, Y, Q, ResTilde) == OK);

     /* Gram-Schmidt QR orthogonalization of FTilde. */
     if (QRLS(Q, ResTilde, Q, R, RBeta, ResTilde) != OK)
     {
          Error("Cannot perform QR decomposition.\n");
          return NUMERIC_ERR;
     }

     /* Compute regression-model beta's. */
     if (TriBackSolve(R, RBeta, Beta) != OK)
     {
          Error("Cannot compute regression beta's.\n");
          return NUMERIC_ERR;
     }

     return OK;
}


int KrigSolve(const Matrix *Chol, const Matrix *F, const real *Y,
          Matrix *FTilde, real *YTilde)
{
     int       ErrNum;
     size_t    j;

     /* Solve Chol' FTilde = F for FTilde */
     ErrNum = OK;
     for (j = 0; j < MatNumCols(F) && ErrNum == OK; j++)
          ErrNum = TriForSolve(Chol, MatCol(F, j), 0,
                    MatCol(FTilde, j));

     /* Solve Chol' YTilde = Y for YTilde; */
     if (ErrNum == OK)
           ErrNum = TriForSolve(Chol, Y, 0, YTilde);

     return ErrNum;
}


/*******************************+++*******************************/
void KrigTranformVec(KrigingModel *KrigMod, real *r)
/*****************************************************************/
/*   Purpose:  Overwrite r with T'r.                             */
/*                                                               */
/*   Version:  1995 July 27                                      */
/*****************************************************************/
{
     Matrix    *T;

     T = KrigT(KrigMod);
     if (T != NULL && MatNumCols(T) > 0)
     {
          /* Overwrite r with T'r. */
          MatVec(T, r, KrigMod->w1);
          VecCopy(KrigMod->w1, MatNumCols(T), r);
     }
}

/*******************************+++*******************************/
boolean KrigIsXActive(const KrigingModel *KrigMod, size_t j)
/*****************************************************************/
/*   Purpose:  Determine whether x variable j is active,         */
/*             in either the regression model or the stochastic- */
/*             process model.                                    */
/*                                                               */
/*   Version:  1996.03.22                                        */
/*****************************************************************/
{
     LinModel  *RegMod, *SPMod;
     real      *Beta, *Theta;

     RegMod = KrigRegMod(KrigMod);
     SPMod  = KrigSPMod(KrigMod);

     Beta = KrigMod->Beta;
     Theta = MatCol(KrigCorPar(KrigMod), 0);

     return ModIsXActive(RegMod, Beta, j)
               || ModIsXActive(SPMod, Theta, j);
}

/*******************************+++*******************************/
size_t KrigSPActiveTerms(const KrigingModel *KrigMod,
          size_t nActiveX, const size_t *xIndex, size_t *IndexTerm)
/*****************************************************************/
/*   Purpose:  Find the indices for terms in the stochastic      */
/*             process model that involve one or more x          */
/*             variables in xIndex and are active (theta > 0).   */
/*             On return, IndexTerm contains the indices.        */
/*                                                               */
/*   Returns:  The number of terms found.                        */
/*                                                               */
/*   Comment:  Calling routine must allocate space for IndexTerm */
/*             (of length ModDF(SPMod)).                         */
/*                                                               */
/*   Version:  1996.03.28                                        */
/*****************************************************************/
{
     real      *Theta;

     Theta = MatCol(KrigCorPar(KrigMod), 0);

     return ModActiveTerms(KrigSPMod(KrigMod), Theta,
               nActiveX, xIndex, IndexTerm);
}

/*******************************+++*******************************/
void frfrAve(KrigingModel *KrigMod, const Matrix *PredReg,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     matrix *frfrj, matrix *frfr)
/*****************************************************************/
/*   Purpose:    Compute average fr(fr)^T over a region.         */
/*                                                               */
/*   Comment:    Calling routine must allocate space for         */
/*               (k + n) * (k + n) matrices frfrj (used for      */
/*               workspace) and frfr.                            */
/*                                                               */
/*   1996.02.16: Wt from RegLevelWt instead of 1/m.              */
/*                                                               */
/*   Version:    1996.02.16                                      */
/*****************************************************************/
{
     real      SPVarPropSave, Wt;
     real      *fr, *g, *xRow;
     size_t    i, j, k, m, n;
     size_t    *xIndex;

     n = MatNumRows(KrigChol(KrigMod));
     k = ModDF(KrigRegMod(KrigMod));

     /* Workspace in KrigMod. */
     fr   = KrigMod->fr;
     g    = KrigMod->gRow;
     xRow = KrigMod->xRow;

     /* Initialize ff part of integral matrix to 1.0,     */
     /* fr part to SPVarProp, and rr part to SPVarProp^2. */
     MatInitValue(0.0, frfr);
     VecInit(1.0, k, fr);
     VecInit(KrigMod->SPVarProp, n, fr + k);
     MatSymUpdate(1.0, fr, frfr);

     /* SPVarProp is applied only once. */
     SPVarPropSave = KrigMod->SPVarProp;
     KrigMod->SPVarProp = 1.0;

     for (j = 0; j < MatNumCols(GroupVarIndex); j++)
     {
          xIndex = MatSize_tCol(GroupVarIndex, j);

          MatInitValue(0.0, frfrj);

          m = RegNumLevels(PredReg, xIndex[0]);
          for (i = 0; i < m; i++)
          {
               fgrGroup(KrigMod, PredReg, GroupSize[j], xIndex, i,
                         nSPTerms[j], MatSize_tCol(IndexSP, j),
                         xRow, fr, g, fr + k);

               Wt = RegLevelWt(PredReg, xIndex[0], i);
               MatSymUpdate(Wt, fr, frfrj);
          }

          MatMultElemWise(frfrj, frfr);
     }

     /* Need to transform with T? */

     KrigMod->SPVarProp = SPVarPropSave;

     return;
}

/*******************************+++*******************************/
void fgrGroup(const KrigingModel *KrigMod, const Matrix *PredReg,
     size_t nXVars, const size_t *xIndex, size_t Level,
     size_t nSPTerms, const size_t *IndexSP, real *xRow, real *f,
     real *g, real *r)
/*****************************************************************/
/*   Purpose:  Compute f, g, and r for a group.                  */
/*                                                               */
/*   Comment:  f, g, r, and xRow (used for workspace) must be    */
/*             allocated by the calling routine.                 */
/*                                                               */
/* 1996.03.28: Created                                           */
/* 2009.05.13: KrigCorVec arguments changed                      */
/*****************************************************************/
{
     size_t n;

     RegLevelsGroup(PredReg, nXVars, xIndex, Level, xRow);

     XToFActive(KrigRegMod(KrigMod), nXVars, xIndex, xRow, f);
     XToFActive(KrigSPMod(KrigMod),  nXVars, xIndex, xRow, g);

     n = MatNumRows(KrigChol(KrigMod));

     KrigCorVec(g, KrigG(KrigMod), n, nSPTerms, IndexSP, YES,
          KrigMod, r);
}




