/*****************************************************************/
/*   ROUTINES TO MANAGE THE LINEAR MODEL                         */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--94.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"

/*******************************+++*******************************/
void XToFActive(const LinModel *Mod, size_t nActive,
          const size_t *xActive, const real *x, real *f)
/*****************************************************************/
/*   Purpose:  Compute the coefficients of the linear model      */
/*             parameters.                                       */
/*                                                               */
/*   Args:     Mod       Linear model.                           */
/*             nActive   Number of active x variables (only used */
/*                       if xActive != NULL).                    */
/*             xActive   If xActive != NULL, then contains the   */
/*                       indices of the x variables that are     */
/*                       active.                                 */
/*             x         Vector of explanatory variables.        */
/*             f         Linear-model coefficients.              */
/*                                                               */
/*   Version:  1994 February 2                                   */
/*****************************************************************/
{
     Matrix    *Term;
     real      fi, xx;
     size_t    i, j, xIndex;

     /* For each linear-model term. */
     for (i = 0; i < Mod->nTerms; i++)
     {
          Term = Mod->Term + i;

          /* For each component of this term. */
          for (fi = 1.0, j = 0; j < MatNumRows(Term); j++)
          {
               xIndex = ModxIndex(Term, j);

               if (xActive != NULL &&
                         VecSize_tIndex(xIndex, nActive, xActive)
                         == INDEX_ERR)
                    continue;

               xx = x[xIndex];

               if (ModCatLevel(Term, j) == 0)
                    /* Quantitative variable. */
                    fi *= ModFn(xx, ModFunc(Term, j));
               else if ((size_t) xx != ModCatLevel(Term, j))
               {
                    /* Categorical dummy takes value 0. */
                    fi = 0.0;
                    break;
               }

               /* Remaining case, dummy takes value 1, */
               /* needs no action.                     */
          }

          f[i] = fi;
     }

     return;
}

/*******************************+++*******************************/
void ModFMatRowIndex(const LinModel *Mod, size_t nRows,
    const size_t *RowIndex, const Matrix *X, Matrix *F)
/*****************************************************************/
/*   Purpose:  Expand, according to model Mod, the nRows rows of */
/*             X with indices RowIndex to give F.  If RowIndex   */
/*             is NULL then all rows of X are used.              */
/*                                                               */
/*   Comment:  The calling routine must allocate space for F     */
/*             (which has nRows rows).                           */
/*                                                               */
/*   Version:  1996.03.31                                        */
/*****************************************************************/
{
     real      *fRow, *xRow;
     size_t    i, ii, nXVars;

     if (RowIndex == NULL)
          nRows = MatNumRows(X);

     nXVars = MatNumCols(X);

     fRow = AllocReal(Mod->nTerms, NULL);
     xRow = AllocReal(nXVars, NULL);

     for (i = 0; i < nRows; i++)
     {
          ii = (RowIndex == NULL) ? i : RowIndex[i];

          MatRow(X, ii, xRow);
          XToF(Mod, xRow, fRow);
          MatRowPut(fRow, i, F);
     }

     AllocFree(fRow);
     AllocFree(xRow);

     return;
}

/*******************************+++*******************************/
boolean ModIsXActive(const LinModel *Mod, const real *Beta,
          size_t xIndex)
/*****************************************************************/
/*   Purpose:  Is a specified x variable active in a linear      */
/*             model?                                            */
/*                                                               */
/*   Args:     Mod       Linear model.                           */
/*             Beta      Linear-model coefficients.              */
/*             xIndex    Index of the x variable of interest.    */
/*                                                               */
/*   Version:  1992 February 28                                  */
/*****************************************************************/
{
     size_t    i;

     /* Look for the x variable in each linear-model term. */
     for (i = 0; i < Mod->nTerms; i++)
          if (ModIsXActiveInTerm(Mod, Beta, xIndex, i))
               /* x variable need be active in only one term. */
               return YES;

     return NO;
}

/*******************************+++*******************************/
size_t ModActiveTerms(const LinModel *Mod, const real *Beta,
          size_t nActiveX, const size_t *xIndex, size_t *IndexTerm)
/*****************************************************************/
/*   Purpose:  Find the indices for terms that involve one or    */
/*             more x variables in xIndex and are active         */
/*             (beta > 0).  On return, IndexTerm contains the    */
/*             indices.                                          */
/*                                                               */
/*   Returns:  The number of term indices found.                 */
/*                                                               */
/*   Comment:  Calling routine must allocate space for IndexTerm */
/*             (of length ModDF(Mod)).                           */
/*                                                               */
/*   Version:  1995 July 27                                      */
/*****************************************************************/
{
     size_t    i, j, nActiveTerms;

     for (nActiveTerms = 0, i = 0; i < Mod->nTerms; i++)
          for (j = 0; j < nActiveX; j++)
               if (ModIsXActiveInTerm(Mod, Beta, xIndex[j], i))
               {
                    /* Only one x variable need */
                    /* be active in term i.     */
                    IndexTerm[nActiveTerms++] = i;
                    break;
               }

     return nActiveTerms;
}

/*******************************+++*******************************/
boolean ModIsXActiveInTerm(const LinModel *Mod, const real *Beta,
          size_t xIndex, size_t TermIndex)
/*****************************************************************/
/*   Purpose:  Is a specified x variable active in a specified   */
/*             linear model term?                                */
/*                                                               */
/*   Args:     Mod       Linear model.                           */
/*             Beta      Linear-model coefficients.              */
/*             xIndex    Index of the x variable of interest.    */
/*             TermIndex Index of the term interest.             */
/*                                                               */
/*   Version:  1994 February 28                                  */
/*****************************************************************/
{
     Matrix    *Term;
     size_t    j;

     if (Beta[TermIndex] == 0.0)
          /* Whole term is inactive. */
          return NO;

     Term = Mod->Term + TermIndex;

     /* For each component of this term. */
     for (j = 0; j < MatNumRows(Term); j++)
          if (ModxIndex(Term, j) == xIndex)
               /* Found the x variable. */
               return YES;

     return NO;
}

