/*****************************************************************/
/*   ROUTINES TO PARSE THE LINEAR MODEL                          */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--95.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"

extern int     ErrorSeverityLevel;

int            TermColType[] = TERM_COL_TYPES;

/*******************************+++*******************************/
int ModParse1(size_t nTerms, const string *TermStr,
          const string ModName, LinModel *Mod)
/*****************************************************************/
/*   Purpose:  First parse of linear-model.                      */
/*             TermStr strings are parsed and LinModel is set up.*/
/*             There is no checking for compatibility with the   */
/*             x-variable names or categorical information,      */
/*             so each term's xIndex column is not assigned.     */
/*             On return, Mod points to the partially-parsed     */
/*             linear model.                                     */
/*                                                               */
/*   Returns:  INPUT_ERR (Mod is ModFree'd) or OK.               */
/*                                                               */
/*   Version:  1994 February 9                                   */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    *Term;
     size_t    i;
     string    Comp, TermStrCopy;

     Mod->nTerms    = nTerms;
     Mod->TermNames = TermStr;

     /* Allocate a matrix for each term. */
     Mod->Term = (Matrix *) AllocGeneric(nTerms, sizeof(Matrix),
               NULL);

     for (ErrNum = OK, i = 0; i < nTerms && ErrNum == OK; i++)
     {
          Term = Mod->Term + i;

          /* Initialize term to have no rows (components). */
          MatAllocate(0, MOD_NUM_COLS, RECT, MIXED, TermColType,
                    NO, Term);

          TermStrCopy = StrDup(TermStr[i]);

          Comp = strtok(TermStrCopy, "*");
          while (Comp != NULL && ErrNum == OK)
          {
               if (strcmp(Comp, "1") != 0 &&
                         (ErrNum = ModAddComp(Comp, Term)) != OK)
                    Error(MOD_TERM, i + 1, ModName);

               Comp = strtok(NULL, "*");
          }

          AllocFree(TermStrCopy);
     }

     if (ErrNum != OK)
          ModFree(Mod);

     return ErrNum;
}

/*******************************+++*******************************/
int ModAddComp(string Comp, Matrix *Term)
/*****************************************************************/
/*   Purpose:  Add a component to a linear model term.           */
/*                                                               */
/*   Returns:  INPUT_ERR or OK.                                  */
/*                                                               */
/*   Version:  1992 February 15                                  */
/*****************************************************************/
{
     int       Fn;
     size_t    CatLevel, Last;

     /* Parse the component. */
     if (ModParseComp(Comp, &CatLevel, &Fn) != OK)
          return INPUT_ERR;

     /* Check x variable has not already appeared in this term. */
     if (StrIndex(Comp, ModxNames(Term), MatNumRows(Term))
               != INDEX_ERR)
     {
          Error("%s should not appear more than once in a term.\n",
                    Comp);
          return INPUT_ERR;
     }

     /* Allocate the new component. */
     MatReAllocate(MatNumRows(Term) + 1, MatNumCols(Term), NULL,
               Term);
     Last = MatNumRows(Term) - 1;
     ModPutxName(Term, Last, Comp);
     ModPutFunc(Term, Last, Fn);
     ModPutCatLevel(Term, Last ,CatLevel);

     return OK;
}

/*******************************+++*******************************/
int ModParseComp(string Comp, size_t *CatLevel, int *Fn)
/*****************************************************************/
/*   Purpose:  Parse a linear-model component of the form        */
/*             Name, Name[Level], or Name^Fn.                    */
/*                                                               */
/*   Args:     Comp      Input:  the component to be parsed;     */
/*                       Output: the variable name terminated    */
/*                               to exclude "[.]" or "^.".       */
/*             CatLevel  Output: the level.                      */
/*             Fn        Output: function (i.e., exponent).      */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*                                                               */
/*   Version:  1993 June 4                                       */
/*****************************************************************/
{
     int       ErrNum;
     string    LevelStr, NextToken;

     *CatLevel = 0;
     *Fn       = 0;

     if (StrBrackets(Comp, &LevelStr, &NextToken) != OK)
     {
          Error("Mismatching brackets.\n");
          ErrNum = INPUT_ERR;
     }

     else if (NextToken != NULL && *NextToken != NULL)
     {
          Error("Characters after \"[]\".\n");
          ErrNum = INPUT_ERR;
     }

     else if (LevelStr != NULL)
     {
          /* Categorical variable. */
          if (StrToSize_t(LevelStr, CatLevel) != OK ||
                    *CatLevel == 0)
          {
               Error("Level must be an integer > 0.\n");
               ErrNum = INPUT_ERR;
          }
          else
               ErrNum = OK;
     }
     else
     {
          *CatLevel = 0;
          ErrNum = ModFnParse(Comp, Fn);
     }

     return ErrNum;
}

/*******************************+++*******************************/
int ModParse2(size_t nXVars, const string *xName,
          const size_t *nCats, const string ModName, LinModel *Mod)
/*****************************************************************/
/*   Purpose:  Second parse of linear-model.                     */
/*             Set up x-variable indices and check the           */
/*             categorical level for each component.             */
/*                                                               */
/*   Returns:  INCOMPAT_ERR or OK.                               */
/*                                                               */
/*   Version:  1995 May 2                                        */
/*****************************************************************/
{
     int       ErrNum;
     Matrix    *Term;
     size_t    i, j, *CatLevel, *xIndex;
     string    *xNameTerm;

     ErrNum = OK;

     /* For each term. */
     for (i = 0; i < Mod->nTerms && ErrNum == OK; i++)
     {
          Term = Mod->Term + i;

          if (MatNumRows(Term) > 0)
          {
               xNameTerm = ModxNames(Term);
               xIndex   = ModxIndexes(Term);
               CatLevel = ModCatLevels(Term);
          }

          /* For each component. */
          for (j = 0; j < MatNumRows(Term) && ErrNum == OK; j++)
          {
               if ( (xIndex[j] = StrIndex(xNameTerm[j], xName,
                         nXVars)) == INDEX_ERR)
               {
                    Error("%s must appear as an x variable.\n",
                              xNameTerm[j]);
                    ErrNum = INCOMPAT_ERR;
               }

               else if (CatLevel[j] > 0)
               {
                    if (nCats == NULL || nCats[xIndex[j]] == 0)
                    {
                         Error("%s has a categorical level so "
                                   "must have " NUM_CATS " > 0.\n",
                                   xNameTerm[j]);
                         ErrNum = INCOMPAT_ERR;
                    }

                    else if (CatLevel[j]
                              > nCats[xIndex[j]])
                    {
                         Error("The level of %s cannot exceed "
                                   NUM_CATS ".\n", xNameTerm[j]);
                         ErrNum = INCOMPAT_ERR;
                    }
               }

               else if (nCats != NULL && nCats[xIndex[j]] > 0)
               {
                    ErrorSeverityLevel = SEV_WARNING;
                    Error("%s has " NUM_CATS " > 0, but is "
                              "appearing linearly.\n",
                              xNameTerm[j]);
                    ErrorSeverityLevel = SEV_ERROR;
                    Output(MOD_TERM, i + 1, ModName);
               }
          }

          if (ErrNum != OK)
               Error(MOD_TERM, i + 1, ModName);
     }

     return ErrNum;
}

/*******************************+++*******************************/
void ModFree(LinModel *Mod)
/*****************************************************************/
/*   Purpose:  Free space allocated for the linear model.        */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*****************************************************************/
{
     size_t    i;

     for (i = 0; i < Mod->nTerms; i++)
          MatFree(Mod->Term + i);

     AllocFree(Mod->Term);

     Mod->nTerms = 0;
     Mod->Term   = NULL;

     return;
}

