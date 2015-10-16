/*****************************************************************/
/*   MINIMIZE OVER CONTINUOUS (SIMPLE BOUNDS) AND/OR DISCRETE    */
/*   VARIABLES                                                   */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--96.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "min.h"

/* These external variables communicate with ObjCont. */
static real    (*ObjFuncExt)(real *x, size_t nDims);
static real    *xExt;
static size_t  *IndexCont, nDimsExt;

/*******************************+++*******************************/
unsigned MinAnyX(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          const Matrix *XReg, size_t nDims, int MinAlg, real *x,
          real *Obj)
/*****************************************************************/
/*   Purpose:  Optimize any objective function of x, where x may */
/*             include fixed, discrete, or continuous            */
/*             variables.                                        */
/*                                                               */
/*   Args:     ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxFuncs  Maximum function evaluations.           */
/*             XReg      X-variable feasibility region.          */
/*             nDims     Number of dimensions.                   */
/*             MinAlg    The minimizer called for unconstrained  */
/*                       minimization.                           */
/*             x         Input:  Starting point (including any   */
/*                               fixed variables);               */
/*                       Output: "Optimal" point.                */
/*             Obj       Input:  Objective at x;                 */
/*                       Output: "Optimal" objective.            */
/*                                                               */
/*   Returns:  Total number of function evaluations.             */
/*                                                               */
/*   96.01.18: CodeBug parameters changed and CodeBug replaced   */
/*             by CodeCheck.                                     */
/*   96.01.18: CodeBug parameters changed.                       */
/*   96.03.07: Extrapolation at end of each iteration.           */
/*   96.03.08: MinConverged replaced by ApproxEq.                */
/*             Extrapolation removed.                            */
/*   96.03.09: Extrapolation at end of each iteration.           */
/*                                                               */
/*   Version:  1996.03.09                                        */
/*****************************************************************/
{
     real      ObjOld;
     real      *ContMax, *ContMin, *xCont, *xOld;
     real      *xExtCopy;
     real      (*ObjFuncExtCopy)(real *x, size_t nDims);
     size_t    i, ThisGroup;
     size_t    *ContDistrib, **GroupIndex, *Group;
     size_t    *GroupSize, *UngroupedIndex;
     size_t    *IndexContCopy, nDimsExtCopy;
     size_t    j, nContVars, nGroups, nUngroupedVars;
     unsigned  nEvals, NumOpts;

     /* Save statics to local variables, */
     /* to enable recursive calling.     */
     ObjFuncExtCopy = ObjFuncExt;
     xExtCopy = xExt;
     IndexContCopy = IndexCont;
     nDimsExtCopy = nDimsExt;

     CodeCheck(RegNumVars(XReg) == nDims);

     ContDistrib = AllocSize_t(nDims, NULL);
     IndexCont   = AllocSize_t(nDims, NULL);
     ContMax     = AllocReal(nDims, NULL);
     ContMin     = AllocReal(nDims, NULL);
     xCont       = AllocReal(nDims, NULL);
     xOld        = AllocReal(nDims, NULL);

     UngroupedIndex = AllocSize_t(nDims, NULL);

     GroupIndex  = AllocPtrSize_t(nDims, NULL);
     Group       = AllocSize_t(nDims, NULL);
     GroupSize   = AllocSize_t(nDims, NULL);

     /* Number of continuous variables. */
     nContVars = 0;

     /* Number of GRID or ungrouped DISCRETE variables. */
     nUngroupedVars = 0;

     /* Number of grouped DISCRETE variables. */
     nGroups = 0;

     /* Count the continuous variables, etc. */
     for (j = 0; j < nDims; j++)
          switch (RegSupport(XReg, j))
          {
               case FIXED:
                    break;

               case CONTINUOUS:
                    IndexCont[nContVars] = j;
                    ContMin[nContVars]   = RegMin(XReg, j);
                    ContMax[nContVars]   = RegMax(XReg, j);
                    ContDistrib[nContVars++]
                              = RegDistrib(XReg, j);
                    break;

               case GRID:
                    UngroupedIndex[nUngroupedVars++] = j;
                    break;

               case DISCRETE:
                    ThisGroup = RegCandGroup(XReg, j);

                    if (ThisGroup == 0)
                         /* Ungrouped. */
                         UngroupedIndex[nUngroupedVars++] = j;

                    else
                    {
                         for (i = 0; i < nGroups; i++)
                              if (ThisGroup == Group[i])
                                   break;

                         if (i == nGroups)
                         {
                              /* New candidate group. */
                              Group[i] = ThisGroup;
                              GroupSize[i] = 1;
                              GroupIndex[i] = AllocSize_t(nDims, NULL);
                              GroupIndex[i][0] = j;
                              nGroups++;
                         }
                         else
                              /* Existing candidate group. */
                              GroupIndex[i][GroupSize[i]++] = j;
                    }

                    break;

               default:
                    CodeBug("Illegal support");
          }

     /* External equivalents. */
     ObjFuncExt = ObjFunc;
     nDimsExt   = nDims;
     xExt       = x;

     /* Iterate until converged. */
     nEvals = 0;
     do
     {
          /* Used in test for convergence. */
          ObjOld = *Obj;
          NumOpts = 0;

          VecCopy(x, nDims, xOld);

          if (nContVars > 0)
          {
               /* Load continuous x's into xCont. */
               VecCopyIndex(nContVars, IndexCont, x, NULL, xCont);

               nEvals += MinCont(ObjCont, AbsTol, RelTol, MaxFuncs,
                         ContMin, ContMax, ContDistrib, nContVars,
                         MinAlg, xCont, Obj);
               NumOpts++;

               /* Put best continuous levels back in x. */
               VecCopyIndex(nContVars, NULL, xCont, IndexCont, x);
          }

          /* Ungrouped variables. */
          for (j = 0; j < nUngroupedVars; j++)
          {
               /* Optimize ungrouped-variable j. */
               nEvals += MinDisc(1, &UngroupedIndex[j], XReg, x,
                         Obj);
               NumOpts++;
          }

          /* Grouped variables. */
          for (j = 0; j < nGroups; j++)
          {
               /* Optimize group j. */
               nEvals += MinDisc(GroupSize[j], GroupIndex[j],
                         XReg, x, Obj);
               NumOpts++;
          }

          /* Try extrapolating. */
          nEvals += MinExtrap(ObjFunc, XReg, nDims, xOld, x, Obj);

     } while (NumOpts > 1 &&
               !ApproxEq(ObjOld, *Obj, AbsTol, RelTol));

     AllocFree(ContDistrib);
     AllocFree(IndexCont);
     AllocFree(ContMax);
     AllocFree(ContMin);
     AllocFree(xCont);
     AllocFree(xOld);

     AllocFree(UngroupedIndex);

     for (j = 0; j < nGroups; j++)
          AllocFree(GroupIndex[j]);
     AllocFree(GroupIndex);
     AllocFree(Group);
     AllocFree(GroupSize);

     /* Restore statics. */
     ObjFuncExt = ObjFuncExtCopy;
     xExt = xExtCopy;
     IndexCont = IndexContCopy;
     nDimsExt = nDimsExtCopy;

     return nEvals;
}

/*******************************+++*******************************/
unsigned MinDisc(size_t NumVars, const size_t *VarIndex,
     const Matrix *XReg, real *x, real *Obj)
/*****************************************************************/
/*   Purpose:  Minimize over the NumVars DISCRETE or GRID        */
/*             variables with indices in VarIndex.               */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   Version:  1992 May 28                                       */
/*****************************************************************/
{
     real      TrialObj;
     real      *xBest;
     size_t    i, j, Index, NumLevels;

     xBest = AllocReal(NumVars, NULL);

     /* Save the current levels. */
     VecCopyIndex(NumVars, VarIndex, x, NULL, xBest);

     /* The number of levels is the same for all variables. */
     NumLevels = RegNumLevels(XReg, VarIndex[0]);

     for (i = 0; i < NumLevels; i++)
     {
          /* Load level i for each variable into x. */
          for (j = 0; j < NumVars; j++)
          {
               Index = VarIndex[j];
               x[Index] = RegLevel(XReg, Index, i);
          }

          TrialObj = (*ObjFuncExt)(x, nDimsExt);

          if (TrialObj < *Obj)
          {
               /* Best point found so far. */
               *Obj = TrialObj;
               VecCopyIndex(NumVars, VarIndex, x, NULL, xBest);
          }
     }

     /* Put best levels back into x. */
     VecCopyIndex(NumVars, NULL, xBest, VarIndex, x);

     AllocFree(xBest);

     return NumLevels;   /* Number of function evaluations. */
}

/*******************************+++*******************************/
real ObjCont(real *xCont, size_t nContVars)
/*****************************************************************/
/*   Purpose:  Put continuous variables in xExt (x), and call    */
/*             objective function.                               */
/*                                                               */
/*   Returns:  Objective.                                        */
/*                                                               */
/*   Version:  1996.03.08                                        */
/*****************************************************************/
{
     size_t c, j;

     /* Put continuous x's into xExt (x). */
     for (c = 0; c < nContVars; c++)
     {
          j = IndexCont[c];
          xExt[j] = xCont[c];
     }

     return((*ObjFuncExt)(xExt, nDimsExt));
}

/*******************************+++*******************************/
unsigned MinMultiStart(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          const Matrix *XReg, size_t nDims, int MinAlg,
          const Matrix *StartPt, real *x, real *Obj)
/*****************************************************************/
/*   Purpose:  Multiple local searches (via MinAnyX) from each   */
/*             row in StartPt.                                   */
/*                                                               */
/*   Returns:  Total number of function evaluations.             */
/*                                                               */
/*   Version:  1996.04.14                                        */
/*****************************************************************/
{
     real      ObjTry;
     real      *xTry;
     size_t    i;
     unsigned  nEvals;

     CodeCheck(MatNumCols(StartPt) == nDims);

     xTry = AllocReal(nDims, NULL);

     nEvals = 0;
     for (*Obj = REAL_MAX, i = 0; i < MatNumRows(StartPt); i++)
     {
          /* Local search starting at row i of StartPt. */
          MatRow(StartPt, i, xTry);

          ObjTry = ObjFunc(xTry, nDims);
          nEvals++;

          nEvals += MinAnyX(ObjFunc, AbsTol, RelTol, MaxFuncs,
                    XReg, nDims, MinAlg, xTry, &ObjTry);

          if (ObjTry < *Obj)
          {
               *Obj = ObjTry;
               VecCopy(xTry, nDims, x);
          }
     }

     AllocFree(xTry);

     return nEvals;
}
