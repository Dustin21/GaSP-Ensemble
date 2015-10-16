/*****************************************************************/
/*   ROUTINES FOR MANAGING GENERIC CORRELATION-PARAMETER         */
/*   MATRICES AND CORRELATIONS                                   */
/*                                                               */   
/*   IMPLEMENTED CORRELATION FUNCTIONS SHOULD BE ACCESSED VIA    */
/*   ONLY THESE FUNCTIONS.                                       */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--2009.                  */
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
#include "alex.h"

/*******************************+++*******************************/
void CorParAlloc
(
     size_t         CorFam,        /* Correlation family         */
     size_t         NumTerms,      /* Number of terms.           */
     const string   *TermName,     /* Term names ("x1", etc.)    */
     Matrix         *CorPar        /* Output: allocated and      */
                                   /* labelled correlation-      */
                                   /* parameter matrix.          */
)
/*****************************************************************/
/* Purpose: Allocate and label a correlation-parameter matrix.   */
/*                                                               */
/* 2009.05.14: Created                                           */
/*****************************************************************/
{
     size_t    i;

     if (CorFam == COR_FAM_POW_EXP)
          PEAlloc(NumTerms, CorPar);
     else if (CorFam == COR_FAM_MATERN)
          MaternAlloc(NumTerms, CorPar);

     /* Label rows (columns labelled for specific family) */     
     for (i = 0; i < NumTerms; i++)
          MatPutRowName(CorPar, i, TermName[i]);

     return;
}


/*******************************+++*******************************/
void CorParStart
(
     size_t CorFam,      /* Correlation family                   */
     const Matrix *G,    /* Expanded-design matrix for the       */
                         /* stochastic-process model.            */
     Matrix *CorPar,     /* Output: Starting values of the       */
                         /* correlation parameters.              */
     Matrix *CorReg      /* Output: Feasibility region for the   */
                         /* correlation parameters.              */
)
/*****************************************************************/
/* Purpose: Return starting values for the correlation           */
/*          parameters and their optimization region.            */
/*                                                               */
/* 2009.05.14: Created.                                          */
/*****************************************************************/
{
     if (CorFam == COR_FAM_POW_EXP)
          PEStart(G, CorPar, CorReg);
     else if (CorFam == COR_FAM_MATERN)
          MaternStart(G, CorPar, CorReg);

     return;
}


/*******************************+++*******************************/
unsigned CorParTest
(
     size_t CorFam,      /* Correlation family                   */
     Matrix *RegCorPar,  /* Feasibility region for the           */
                         /* correlation parameters.              */
     size_t TermIndex,   /* Index of the tested term.            */
     real   AbsTol,
     real   CritLogLikeDiff,
     Matrix *CorPar,     /* Input:  Correlation parameters;      */
                         /* Output: Row TermIndex may change.    */
     real   *NegLogLike  /* Input: Negative log likelihood;      */
                         /* Output: New value.                   */
)
/*****************************************************************/
/* Purpose: Test parameters against "null" values for a single   */
/*          term.                                                */
/*                                                               */
/* Returns: Number of function evaluations.                      */
/*                                                               */
/* 2009.05.14: Created.                                          */
/*****************************************************************/
{
     unsigned  NumFuncs;
     
     if (CorFam == COR_FAM_POW_EXP)
          NumFuncs = PETest(RegCorPar, TermIndex, AbsTol, CritLogLikeDiff,
                    CorPar, NegLogLike);
     else if (CorFam == COR_FAM_MATERN)
          NumFuncs = MaternTest(RegCorPar, TermIndex, AbsTol, CritLogLikeDiff,
                    CorPar, NegLogLike);

     return NumFuncs;
}

/*******************************+++*******************************/
int CorParSetUp(const Matrix *CorPar, const string yName,
     real SPVar, real ErrVar, KrigingModel *KrigMod)
/*****************************************************************/
/*   Purpose:  Set up correlation parameters.                    */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Version:  1996.01.20                                        */
/*****************************************************************/
{
     int  ErrNum;

     if ( (ErrNum = CorParExtract(CorPar, yName, YES,
               KrigCorPar(KrigMod))) == OK)
     {
          if (KrigMod->RanErr)
          {
               KrigMod->SigmaSq   = SPVar + ErrVar;
               KrigMod->SPVarProp = SPVar / KrigMod->SigmaSq;
          }
          else
          {
               KrigMod->SigmaSq   = SPVar;
               KrigMod->SPVarProp = 1.0;
          }
     }

     return ErrNum;
}

/*******************************+++*******************************/
int CorParExtract(const Matrix *CorPar, const string YName,
          boolean ErrorMessage, Matrix *CorParYName)
/*****************************************************************/
/*   Purpose:  Extract the correlation parameters for YName from */
/*             CorPar and put them into CorParYName.             */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if CorPar does not include the       */
/*                          the columns for YName's parameters;  */
/*             OK           otherwise.                           */
/*                                                               */
/*   Comment:  Calling routine must allocate space for           */
/*             CorParYName.                                      */
/*             Legality of the values in CorPar and              */
/*             compatibility of CorPar with the stochastic-      */
/*             process model are assumed.                        */
/*                                                               */
/*   96.03.12: ErrorMessage added.                               */
/*                                                               */
/*   Version:  1996.03.12                                        */
/*****************************************************************/
{
     int       ErrNum;
     real      *TargetCol;
     size_t    j;
     string    TargetName;

     ErrNum = OK;
     for (j = 0; j < MatNumCols(CorParYName) && ErrNum == OK; j++)
     {
          TargetName = StrPaste(3, MatColName(CorParYName, j), ".",
                    YName);

          if ( (TargetCol = MatColFind(CorPar, TargetName, NO))
                    == NULL)
          {
               if (ErrorMessage)
                    Error(DB_COMPULSORY, SP_MOD, TargetName);
               ErrNum = INCOMPAT_ERR;
          }
          else
          {
               MatPutColName(CorParYName, j, TargetName);
               VecCopy(TargetCol, MatNumRows(CorParYName),
                         MatCol(CorParYName, j));
          }

          AllocFree(TargetName);
     }

     return ErrNum;
}

/*******************************+++*******************************/
boolean CorParIsActive
(
     size_t       CorFam,    /* Correlation family               */
     const Matrix *CorPar,  /* Correlation parameters           */
     size_t       TermIndex  /* Index of the term of interest    */
)
/*****************************************************************/
/* Purpose: Is term TermIndex active in the correlation          */
/*             function?                                         */
/*                                                               */
/* 2009.05.14: Created                                           */
/*****************************************************************/
{
     if (CorFam == COR_FAM_POW_EXP)
          return PEIsActive(CorPar,TermIndex);
     else if (CorFam == COR_FAM_MATERN)
          return MaternIsActive(CorPar,TermIndex);
     else
     {
          CodeBug("Illegal correlation family\n");
          return 0;  /* So compiler always has a return type. */
     }
}

/*******************************+++*******************************/
void KrigCorVec
(
     const real   *g,        /* An arbitrary point in g space.   */
     const Matrix *G,        /* Matrix of arbitrary points.      */
     size_t       n,         /* Only the correlations for the    */
                             /* first n rows of G are computed.  */
     size_t       nActive,   /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     boolean      applySPVarProp,  /* Should correlations be     */
                                   /* multiplied by SPVarProp?   */
     const KrigingModel *KrigMod,
     real         *r
)
/*****************************************************************/
/* Purpose: Compute correlations between the point g and the     */
/*             points in the first n rows of G.                  */
/*                                                               */
/* 2009.05.14: Created                                           */
/*****************************************************************/
{
     if (KrigCorFam(KrigMod) == COR_FAM_POW_EXP)
          PECor(g, G, n, nActive, Active, KrigCorPar(KrigMod), r);
     else if (KrigCorFam(KrigMod) == COR_FAM_MATERN)
          MaternCor(g, G, n, nActive, Active, KrigCorPar(KrigMod), r);

     if (applySPVarProp && KrigMod->SPVarProp < 1.0)
          VecMultScalar(KrigMod->SPVarProp, n, r);
}


