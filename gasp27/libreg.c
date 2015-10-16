/*****************************************************************/
/*   ROUTINES FOR MANAGING REGIONS                               */
/*                                                               */
/*   Copyright (c) William J. Welch 1992--96.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

static string RegColName[] = REG_COL_NAMES;
static int    RegColType[] = REG_COL_TYPES;

static string DistribName[] = DISTRIB_NAMES;
static string SupportName[] = SUPPORT_NAMES;

/* Number of levels for approximating a continuous variable. */
/* 1996.02.17: changed from NUM_CONT_LEVELS 10.              */
/* 2002.04.12: NUM_CONT_LEVELS_EXC increased from 10 to 100. */
/* 2002.04.12: NUM_CONT_LEVELS_INC increased from 11 to 101. */
/* 2002.05.14: Back again to 10 and 11.                      */
#define NUM_CONT_LEVELS_EXC  10
#define NUM_CONT_LEVELS_INC  11


/*******************************+++*******************************/
void RegAlloc(size_t nVars, Matrix *Reg)
/*****************************************************************/
/*   Purpose:    Allocate space for a region matrix with nVars   */
/*               variables.                                      */
/*                                                               */
/*   1996.02.17: Default changed to inclusive range.             */
/*                                                               */
/*   Version:    1996.02.17                                      */
/*****************************************************************/
{
     size_t    i, j;
     string    VarName;

     MatAllocate(nVars, NumStr(RegColName), RECT, MIXED,
               RegColType, YES, Reg);

     for (j = 0; j < MatNumCols(Reg); j++)
          MatPutColName(Reg, j, RegColName[j]);

     /* Defaults: */
     for (i = 0; i < nVars; i++)
     {
          VarName = StrPaste(2, "x", StrFromSize_t(i + 1));
          RegPutVar(Reg, i, VarName);
          AllocFree(VarName);

          RegPutSupport(       Reg, i, CONTINUOUS);
          RegPutMin(           Reg, i, 0.0);
          RegPutMax(           Reg, i, 1.0);
          RegPutNumLevels(     Reg, i, 0);
          RegPutNumCats(       Reg, i, 0);
          RegPutDistrib(       Reg, i, UNIFORM);
          RegPutInclusive(     Reg, i, YES);
          RegPutCandGroup(     Reg, i, 0);
          RegPutStep(          Reg, i, 0.0);
          RegPutCandMatIndex(  Reg, i, 0);
          RegPutWt(            Reg, i, 1.0);
          RegPutAnalyze(       Reg, i, YES);
          RegPutTransformation(Reg, i, NONE);
     }

     return;
}

/*******************************+++*******************************/
int RegExtract(const Matrix *XDescrip, const string XDescripName,
     const string ColExt, Matrix *Reg)
/*****************************************************************/
/*   Purpose:  Extract region Reg from XDescrip.                 */
/*                                                               */
/*   Returns:  INPUT_ERR    if user's matrix is illegal;         */
/*             OK           otherwise.                           */
/*                                                               */
/*   96.01.18: CodeBug parameters changed.                       */
/*   96.01.22: CodeBug parameters changed.                       */
/*   96.10.21: Empty (no VARIABLE column) XDescrip trapped.      */
/*             (Previously generated hard fail.)                 */
/*   96.02.18: */
/*                                                               */
/*   Version:  1996.02.03/1996.10.21                             */
/*****************************************************************/
{
     boolean   IsCat;
     int       ErrNum;
     real      Step, TempMax, TempMin;
     real      *Max, *Min, *Wt;
     size_t    i, SupportNum;
     size_t    *CandGroup, *nCats, *nLevels;
     string    *Distrib, *Inclusive, *Support, *Variable;

     Variable  = MatStrColFind(XDescrip, VARIABLE, NO);

     if (Variable == NULL)
          return INPUT_ERR;

     nCats     = MatSize_tColFind(XDescrip, NUM_CATS,  NO);
     Inclusive = MatStrColFind(   XDescrip, INCLUSIVE, NO);
     Wt        = MatColFind(      XDescrip, WEIGHT,    NO);

     if ( (Support = MatStrColFind(XDescrip,
               StrPaste(2, SUPPORT, ColExt),  NO)) == NULL)
          Support = MatStrColFind(XDescrip, SUPPORT, NO);

     if ( (Min = MatColFind(XDescrip,
               StrPaste(2, MIN, ColExt),  NO)) == NULL)
          Min = MatColFind(XDescrip, MIN, NO);

     if ( (Max = MatColFind(XDescrip,
               StrPaste(2, MAX, ColExt),  NO)) == NULL)
          Max = MatColFind(XDescrip, MAX, NO);

     if ( (nLevels = MatSize_tColFind(XDescrip,
               StrPaste(2, NUM_LEVELS, ColExt),  NO)) == NULL)
          nLevels = MatSize_tColFind(XDescrip, NUM_LEVELS, NO);

     if ( (Distrib = MatStrColFind(XDescrip,
               StrPaste(2, DISTRIBUTION, ColExt), NO)) == NULL)
          Distrib = MatStrColFind(XDescrip, DISTRIBUTION, NO);

     if ( (CandGroup = MatSize_tColFind(XDescrip,
               StrPaste(2, CAND_GROUP, ColExt), NO)) == NULL)
          CandGroup = MatSize_tColFind(XDescrip, CAND_GROUP, NO);

     RegAlloc(MatNumRows(XDescrip), Reg);

     ErrNum = OK;
     for (i = 0; i < MatNumRows(XDescrip) && ErrNum == OK; i++)
     {
          RegPutVar(Reg, i, Variable[i]);

          if (nCats != NULL)
               RegPutNumCats(Reg, i, nCats[i]);

          if (Wt != NULL)
               RegPutWt(Reg, i, Wt[i]);

          IsCat = (RegNumCats(Reg, i) > 0);

          /* Default in RegAlloc is CONTINUOUS. */
          if (Support != NULL)
          {
               SupportNum = StrIndex(Support[i], SupportName,
                         NumStr(SupportName));
               RegPutSupport(Reg, i, SupportNum);
          }

          if (Min != NULL || Max != NULL)
          {
               TempMin = (Min != NULL) ? Min[i] : RegMin(Reg, i);
               TempMax = (Max != NULL) ? Max[i] : RegMax(Reg, i);

               if (TempMin > TempMax)
               {
                    Error(REG_MINMAX, XDescripName, Variable[i]);
                    ErrNum = INPUT_ERR;
                    break;
               }
               else
               {
                    RegPutMin(Reg, i, TempMin);
                    RegPutMax(Reg, i, TempMax);
               }
          }

          switch (RegSupport(Reg, i))
          {
               case CONTINUOUS:
                    if (IsCat)
                    {
                         Error(REG_CAT_CONT, XDescripName, Variable[i]);
                         ErrNum = INPUT_ERR;
                         break;
                    }

                    if (Distrib != NULL)
                         RegPutDistrib(Reg, i,
                                   StrIndex(Distrib[i], DistribName,
                                   NumStr(DistribName)));

                    if (RegDistrib(Reg, i) == NORMAL)
                         /* Change default range  */
                         /* to avoid bad scaling. */
                         RegPutInclusive(Reg, i, NO);
                         
                    if (Inclusive != NULL)
                         RegPutInclusive(Reg, i,
                                   (stricmp(Inclusive[i], YES_STR) == 0) ?
                                   YES : NO);

                    if (RegInclusive(Reg, i) &&
                             RegDistrib(Reg, i) == NORMAL)
                    {
                         Error(REG_INC_NORM, XDescripName, Variable[i]);
                         ErrNum = INPUT_ERR;
                    }

                    /* Approximation by discrete levels. */
                    if (RegMin(Reg, i) == RegMax(Reg, i))
                         RegPutNumLevels(Reg, i, 1);
                    else if (RegInclusive(Reg, i))
                         RegPutNumLevels(Reg, i, NUM_CONT_LEVELS_INC);
                    else
                         RegPutNumLevels(Reg, i, NUM_CONT_LEVELS_INC);

                    break;

               case DISCRETE:
                    if (CandGroup != NULL)
                         RegPutCandGroup(Reg, i, CandGroup[i]);
                    break;

               case FIXED:
                    RegPutNumLevels(Reg, i, 1);
                    break;

               case GRID:
                    if (nLevels == NULL || nLevels[i] < 1)
                    {
                         Error(REG_NUM_LEVELS, XDescripName,
                                   Variable[i]);
                         ErrNum = INPUT_ERR;
                    }

                    else if (nLevels[i] == 1 &&
                              RegMin(Reg, i) != RegMax(Reg, i))
                    {
                         Error(REG_NUM_LEVELS_MINMAX, XDescripName,
                                   Variable[i]);
                         ErrNum = INPUT_ERR;
                    }

                    else if (IsCat && (RegMin(Reg, i) < 1.0 ||
                              fmod(RegMin(Reg, i), 1.0) != 0.0))
                    {
                         Error(REG_CAT_MIN, XDescripName, Variable[i]);
                         ErrNum = INPUT_ERR;
                    }

                    else if (IsCat &&
                              (RegMax(Reg, i) > (double) nCats[i] ||
                              fmod(RegMax(Reg, i), 1.0) != 0.0))
                    {
                         Error(REG_CAT_MAX, XDescripName, Variable[i]);
                         ErrNum = INPUT_ERR;
                    }

                    else
                    {
                         RegPutNumLevels(Reg, i, nLevels[i]);
                         Step = (nLevels[i] == 1) ? 0.0 :
                                   (RegMax(Reg, i) - RegMin(Reg, i)) /
                                   (nLevels[i] - 1.0);
                         RegPutStep(Reg, i, Step);
                    }

                    if (ErrNum == OK && IsCat && nLevels[i] > 1 &&
                             (RegStep(Reg, i) < 0.99 ||
                             RegStep(Reg, i) > 1.01))
                    {
                         Error(REG_CAT_STEP, XDescripName, Variable[i]);
                         ErrNum = INPUT_ERR;
                    }

                    break;

               default:
                    CodeBug("Illegal support");
          }

          if (CandGroup != NULL && CandGroup[i] > 0 &&
                    RegSupport(Reg, i) != DISCRETE)
          {
               Error(REG_CAND_GROUP, Variable[i], CandGroup[i]);
               ErrNum = INPUT_ERR;
          }
     }

     if (ErrNum != OK)
          MatFree(Reg);

     return ErrNum;
}

/*******************************+++*******************************/
int RegCandCompat(const Matrix *Cand, Matrix *Reg)
/*****************************************************************/
/*   Purpose:  Update Reg for DISCRETE variables.                */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if there is an incompatibility;      */
/*             OK           otherwise.                           */
/*                                                               */
/*   Version:  1995 February 20                                  */
/*****************************************************************/
{
     int       ErrNum;
     real      *Col;
     size_t    CandIndex, i, nCands;

     ErrNum = OK;
     for (i = 0; i < MatNumRows(Reg) && ErrNum == OK; i++)
     {
          if (RegSupport(Reg, i) == DISCRETE)
          {
               CandIndex = MatColIndex(Cand, RegVar(Reg, i));

               if (CandIndex != INDEX_ERR)
               {
                    RegPutCandMatIndex(Reg, i, CandIndex);
                    nCands = MatNumRows(Cand);
                    RegPutNumLevels(Reg, i, nCands);
                    Col = MatCol(Cand, CandIndex);
                    RegPutMin(Reg, i, VecMin(Col, nCands));
                    RegPutMax(Reg, i, VecMax(Col, nCands));
               }
               else
               {
                    Incompatibility(REG_CAND, RegVar(Reg, i));
                    ErrNum = INCOMPAT_ERR;
               }
          }
     }

     Reg->Next = Cand;

     return ErrNum;
}

/*******************************+++*******************************/
void RegRandPt(const Matrix *Reg, real *x)
/*****************************************************************/
/*   Purpose:  Generate a random point x in region Reg.          */
/*                                                               */
/*   Comment:  Calling routine must allocate space for x.        */
/*                                                               */
/*   Version:  1995 February 22                                  */
/*****************************************************************/
{
     real      u;
     size_t    j, jj, Group;

     for (j = 0; j < MatNumRows(Reg); j++)
     {
          if (RegSupport(Reg, j) == FIXED)
               continue;

          u = RandUnif();
          x[j] = RegTransform(u, Reg, j);

          if (RegSupport(Reg, j) == DISCRETE &&
                    (Group = RegCandGroup(Reg, j)) > 0)
               /* Apply same u to all previous variables */
               /* in the group.                          */
               for (jj = 0; jj < j; jj++)
               if (RegSupport(Reg, jj) == DISCRETE &&
                         RegCandGroup(Reg, jj) == Group)
                    x[jj] = RegTransform(u, Reg, jj);
     }

     return;
}

/*******************************+++*******************************/
real RegTransform(real u, const Matrix *Reg, size_t j)
/*****************************************************************/
/*   Purpose:  Transform Uniform[0, 1] u into valid value for    */
/*             variable j of Reg.                                */
/*                                                               */
/*   Returns:  Transformed value.                                */
/*                                                               */
/*   96.01.18: CodeBug parameters changed.                       */
/*   96.01.22: CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     real r;

     switch (RegSupport(Reg, j))
     {
          case CONTINUOUS:
               r = RegTransformCont(u, RegMin(Reg, j),
                         RegMax(Reg, j), RegDistrib(Reg, j));
               break;

          case GRID:
          case DISCRETE:
               r = RegLevel(Reg, j,
                         (size_t) floor(u * RegNumLevels(Reg, j)));
               break;

          default:
               CodeBug("Illegal support");
     }

     return r;
}

/*******************************+++*******************************/
real RegTransformCont(real u, real a, real b, size_t DistribNum)
/*****************************************************************/
/* Purpose:    Transform Uniform[0, 1] u into the continuous     */
/*             region [a, b].                                    */
/*                                                               */
/* Returns:    Transformed value.                                */
/*                                                               */
/* 1996.01.18: CodeBug parameters changed.                       */
/*             InvNormal replaced by CDFInvNorm.                 */
/* 1996.01.22: CodeBug parameters changed.                       */
/* 1999.06.25: case ARCTAN added.                                */
/*****************************************************************/
{
     real r;

     switch (DistribNum)
     {
          case UNIFORM:
               r = a + u * (b - a);
               break;

          case NORMAL:
               r = a + (0.5 + CDFInvNorm(u) / 6.0) * (b - a);
               break;

          case LOGUNIFORM:
               r = a * pow(b / a, u);
               break;

          case ARCSIN:
               r = a + (sin((2.0 * u - 1.0) * HALFPI) + 1.0)
                         * 0.5 * (b - a);
               break;

          case ARCTAN:
               r = a + tan(u * atan(b - a));
               break;

          /* EXPONENTIAL was previously also allowed: */
          /*   er = exp(-range);                      */
          /*   r = a - log(er + u * (1 - er));     */
          /*   break;                              */

          default:
               CodeBug("Illegal distribution");
     }

     return r;
}

/*******************************+++*******************************/
real RegLevel(const Matrix *Reg, size_t j, size_t LevelIndex)
/*****************************************************************/
/*   Purpose:  Return level LevelIndex of variable j.            */
/*                                                               */
/*   96.01.18: CodeBug parameters changed.                       */
/*   96.01.22: CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     real r, u;

     switch (RegSupport(Reg, j))
     {
          case CONTINUOUS:
               if (RegInclusive(Reg, j))
                    u = (real) LevelIndex /
                              (real) (RegNumLevels(Reg, j) - 1);
               else
                    u = (LevelIndex + 0.5) / RegNumLevels(Reg, j);

               r = RegTransformCont(u, RegMin(Reg, j),
                         RegMax(Reg, j), RegDistrib(Reg, j));
               break;

          case DISCRETE:
               r = MatElem(Reg->Next, LevelIndex,
                         RegCandMatIndex(Reg, j));
               break;

          case FIXED:
               CodeCheck(LevelIndex == 0);
               r = RegMin(Reg, j);
               break;

          case GRID:
               r = RegMin(Reg, j) + LevelIndex * RegStep(Reg, j);
               break;

          default:
               CodeBug("Illegal support");
     }

     return r;
}

/*******************************+++*******************************/
real RegLevelWt(const Matrix *Reg, size_t j, size_t LevelIndex)
/*****************************************************************/
/*   Purpose:  Return (integration) weight for level LevelIndex  */
/*             of variable j.                                    */
/*                                                               */
/*   Version:  1996.02.17                                        */
/*****************************************************************/
{
     real      Wt;
     size_t    m;

     m = RegNumLevels(Reg, j);
     CodeCheck(m > 0);

     switch (RegSupport(Reg, j))
     {
          case CONTINUOUS:
               if (m == 1)
                    Wt = 1.0;
               else if (RegInclusive(Reg, j))
               {
                    CodeCheck(is_odd(m));

                    if (LevelIndex == 0 || LevelIndex == m - 1)
                         Wt = 1.0;
                    else if (is_odd(LevelIndex))
                         Wt = 4.0;
                    else
                         Wt = 2.0;

                    /* Number of intervals is 1 + (m - 3) / 2. */
                    Wt = Wt / 6.0 / (1 + (m - 3) / 2);
               }     
               else
                    Wt = 1.0 / m;
               break;

          case DISCRETE:
               Wt = 1.0 / m;
               break;

          case FIXED:
               CodeCheck(m == 1 && LevelIndex == 0);
               Wt = 1.0;
               break;

          case GRID:
               Wt = 1.0 / m;
               break;

          default:
               CodeBug("Illegal support");
     }

     return Wt;
}

/*******************************+++*******************************/
size_t RegGroupIndices(const Matrix *Reg, size_t j,
     size_t *Index)
/*****************************************************************/
/*   Purpose:  On return, Index contains the indices of the      */
/*             variables in the same group as variable j         */
/*             (including j).                                    */
/*                                                               */
/*   Returns:  Number of variables in the group                  */
/*             (1 for ungrouped).                                */
/*                                                               */
/*   Comments: Calling routine must allocate space for Index.    */
/*                                                               */
/*   Version:  1996.03.28                                        */
/*****************************************************************/
{
     size_t    Group, jj, GroupSize;

     if ( (Group = RegCandGroup(Reg, j)) == 0)
     {
          /* Variable j is not grouped. */
          GroupSize = 1;
          Index[0] = j;
     }
     else
     {
          /* Variable j is grouped. */
          for (GroupSize = 0, jj = 0; jj < MatNumRows(Reg); jj++)
               if (RegCandGroup(Reg, jj) == Group)
                    Index[GroupSize++] = jj;
     }

     return GroupSize;
}

/*******************************+++*******************************/
size_t RegGroupings(const Matrix *Reg, size_t **GroupSize,
          Matrix *Index)
/*****************************************************************/
/*   Purpose:  Returns the number of groups.  On return,         */
/*             *GroupSize[i] is the number of variables in group */
/*             i, and column i of Index contains the indices of  */
/*             the variables in group i.                         */
/*                                                               */
/*   Returns:  Number of groups.                                 */
/*                                                               */
/*   Comments: GroupSize and Index allocated here.               */
/*                                                               */
/*   Version:  1996.03.28                                        */
/*****************************************************************/
{
     size_t    j, nGroups, nXVars;
     size_t    *IndexCol;

     nXVars = MatNumRows(Reg);

     *GroupSize = AllocSize_t(nXVars, NULL);

     MatAllocate(nXVars, nXVars, RECT, SIZE_T, NULL, NO, Index);

     for (nGroups = 0, j = 0; j < nXVars; j++)
     {
          IndexCol = MatSize_tCol(Index, nGroups);
          (*GroupSize)[nGroups] = RegGroupIndices(Reg, j,
                    IndexCol);
          if (IndexCol[0] == j)
               /* A new group or an ungrouped variable. */
               nGroups++;
     }

     MatReAlloc(nXVars, nGroups, Index);

     return nGroups;
}

/*******************************+++*******************************/
void RegLevelsGroup(const Matrix *Reg, size_t GroupSize,
          const size_t *Index, size_t LevelIndex, real *x)
/*****************************************************************/
/*   Purpose:  Load level LevelIndex of the GroupSize variables  */
/*             with indices Index into x.                        */
/*                                                               */
/*   Comments: Calling routine must allocate space for x.        */
/*                                                               */
/*   Version:  1996.03.28                                        */
/*****************************************************************/
{
     size_t    j, jj;

     for (jj = 0; jj < GroupSize; jj++)
     {
          j = Index[jj];
          x[j] = RegLevel(Reg, j, LevelIndex);
     }

     return;
}

/*******************************+++*******************************/
boolean RegIsCand(const Matrix *Reg)
/*****************************************************************/
/*   Purpose:  Is Reg a candidate set (i.e., all discrete        */
/*             variables in the same group)?                     */
/*                                                               */
/*   Version:  1995 February 4                                   */
/*****************************************************************/
{
     size_t    Group, j;

     if (RegSupport(Reg, 0) != DISCRETE ||
               (Group = RegCandGroup(Reg, 0)) == 0)
          return NO;

     for (j = 1; j < MatNumRows(Reg); j++)
          if (RegSupport(Reg, j) != DISCRETE ||
                    RegCandGroup(Reg, j) != Group)
               return NO;

     return YES;
}
