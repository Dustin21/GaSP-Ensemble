/*****************************************************************/
/*   MINIMIZE BY EXTRAPOLATING A PROMISING DIRECTION             */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--92.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "min.h"

extern size_t  nPointers;

/* Extrapolation schedule. */
real Gamma[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0};

#define NUMEXTRAP   (sizeof(Gamma) / sizeof(real))

/*******************************+++*******************************/
unsigned MinExtrap(real (*ObjFunc)(real *x, size_t nDims),
          const Matrix *Reg, size_t nDims, const real *xOld,
          real *xNew, real *Obj)
/*****************************************************************/
/*   Purpose:  Find optimal extrapolation of xOld through       */
/*             xNew, a better point.                            */
/*                                                               */
/*   Args:     ObjFunc   Calulates objective function.           */
/*             Reg       Feasiblity region.                      */
/*             nDims     Dimension of points.                    */
/*             xOld      The old point.                          */
/*             xNew      Input:  The new point.                  */
/*                       Output: The "optimal" point.            */
/*             Obj       Input:  The objective for xNew.        */
/*                       Output: The "optimal" objective.        */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   Version:  1996.03.07                                        */
/*****************************************************************/
{
     real      BestGamma, NewObj;
     real      *xExtrap;
     size_t    i;
     unsigned  nEvals;

     xExtrap = AllocReal(nDims, NULL);

     nEvals = 0;
     for (BestGamma = 0.0, i = 0; i < NUMEXTRAP; i++)
     {
          if (Extrap(nDims, xOld, xNew, Gamma[i], Reg, xExtrap)
                    == NO)
               break;

          NewObj = (*ObjFunc)(xExtrap, nDims);
          nEvals++;

          if (NewObj < *Obj)
          {
               BestGamma = Gamma[i];
               *Obj = NewObj;
          }
          else
               break;
     }

     Extrap(nDims, xOld, xNew, BestGamma, Reg, xNew);

     AllocFree(xExtrap);

     return nEvals;
}

/*******************************+++*******************************/
boolean Extrap(size_t nDims, const real *xOld, const real *xNew,
          real Gamma, const Matrix *Reg, real *xExtrap)
/*****************************************************************/
/*   Purpose:  Extrapolate xOld through xNew by factor Gamma,  */
/*             keeping the extrapolated point in the feasibility */
/*             region, Reg.                                      */
/*                                                               */
/*   Args:     nDims     Dimension of points.                    */
/*             xOld     Old point.                              */
/*             xNew     New point.                              */
/*             Gamma     Extrapolation factor.                   */
/*             Reg       Feasiblity region.                      */
/*             xExtrap  Output: Extrapolated point.             */
/*                                                               */
/*   Returns:  YES if there is at least one continuous variable  */
/*                 that can be extrapolated;                     */
/*             NO  otherwise.                                    */
/*                                                               */
/*   Version:  1996.03.18                                        */
/*****************************************************************/
{
     boolean   CanExtrap;
     size_t    i;

     for (CanExtrap = NO, i = 0; i < nDims; i++)
     {
          if (RegSupport(Reg, i) != CONTINUOUS)
               xExtrap[i] = xNew[i];
          else
          {
               CanExtrap = YES;

               xExtrap[i] = xNew[i]
                         + Gamma * (xNew[i] - xOld[i]);

               if (xExtrap[i] < RegMin(Reg, i))
                    xExtrap[i] = RegMin(Reg, i);
               else if (xExtrap[i] > RegMax(Reg, i))
                    xExtrap[i] = RegMax(Reg, i);
          }
     }

     return CanExtrap;
}

/* Points within TOL_RANGE * Range of a bound are "close" */
/* to the bound.                                          */
#define TOL_RANGE   0.01

/*******************************+++*******************************/
unsigned MinTryBounds(real (*ObjFunc)(real *x, size_t nDims),
          size_t nDims, const real *LowBnd, const real *UpBnd,
          real *x, real *Obj)
/*****************************************************************/
/*   Purpose:  For each (continuous) variable close to one of    */
/*             its bounds, try putting it at the bound.          */
/*                                                               */
/*   Args:     ObjFunc   Calulates objective function.           */
/*             nDims     Dimension of x.                         */
/*             LowBnd    Lower bounds on x.                      */
/*             UpBnd     Upper bounds on x.                      */
/*             x         Input:  Starting point.                 */
/*                       Output: A possibly better point.        */
/*             Obj       Input:  The objective for x.            */
/*                       Output: A possibly better objective.    */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   Version:  1996.03.08                                        */
/*****************************************************************/
{
     real      NewObj;
     real      Tol, xSave;
     size_t    j;
     unsigned  nEvals;

     nEvals = 0;

     for (j = 0; j < nDims; j++)
     {
          if (x[j] == LowBnd[j] || x[j] == UpBnd[j])
               continue;

          Tol = TOL_RANGE * (UpBnd[j] - LowBnd[j]);

          if (x[j] - LowBnd[j] < Tol)
          {
               /* Close to lower bound. */
               xSave = x[j];
               x[j] = LowBnd[j];
               NewObj = (*ObjFunc)(x, nDims);
               nEvals++;
               if (NewObj < *Obj)
                    *Obj = NewObj;
               else
                    x[j] = xSave;
          }
          else if (UpBnd[j] - x[j] < Tol)
          {
               /* Close to upper bound. */
               xSave = x[j];
               x[j] = UpBnd[j];
               NewObj = (*ObjFunc)(x, nDims);
               nEvals++;
               if (NewObj < *Obj)
                    *Obj = NewObj;
               else
                    x[j] = xSave;
          }
     }

     return nEvals;
}


