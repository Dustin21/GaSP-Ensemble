/*****************************************************************/
/*   Copyright (c) William J. Welch 1990--96.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*   Version: 1996.03.31                                         */
/*****************************************************************/

#define TERM_COL_TYPES   {STRING, SIZE_T, INTEGER, SIZE_T}
                        /* "xName", "xIndex", "Func", "CatLevel" */
#define MOD_NUM_COLS     4

#define ModxName(   Term, j)  MatStrElem(   Term, j, 0)
#define ModxIndex(  Term, j)  MatSize_tElem(Term, j, 1)
#define ModFunc(    Term, j)  MatIntElem(   Term, j, 2)
#define ModCatLevel(Term, j)  MatSize_tElem(Term, j, 3)

#define ModxNames(    Term)   MatStrCol(   Term, 0)
#define ModxIndexes(  Term)   MatSize_tCol(Term, 1)
#define ModFuncs(     Term)   MatIntCol(   Term, 2)
#define ModCatLevels( Term)   MatSize_tCol(Term, 3)

#define ModPutxName(   Term, j, s) MatPutStrElem(   Term, j, 0, s)
#define ModPutxIndex(  Term, j, z) MatPutSize_tElem(Term, j, 1, z)
#define ModPutFunc(    Term, j, i) MatPutIntElem(   Term, j, 2, i)
#define ModPutCatLevel(Term, j, z) MatPutSize_tElem(Term, j, 3, z)

/* New type for a linear model. */
typedef struct
{
     size_t    nTerms;        /* Number of terms.        */
     string    *TermNames;    /* Names of the terms.     */
     Matrix    *Term;         /* A matrix for each term. */
} LinModel;

#define ModDF(Mod)            ((Mod)->nTerms)
#define ModTermNames(Mod)     ((Mod)->TermNames)

#define MOD_TERM    "At term %d of %s.\n"

/* model.c: */

/*****************************************************************/
void XToFActive(const LinModel *Mod, size_t nActive,
          const size_t *xActive, const real *x, real *f);
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
/*****************************************************************/

/*****************************************************************/
void ModFMatRowIndex(const LinModel *Mod, size_t nRows,
    const size_t *RowIndex, const Matrix *X, Matrix *F);
/*****************************************************************/
/*   Purpose:  Expand, according to model Mod, the nRows rows of */
/*             X with indices RowIndex to give F.  If RowIndex   */
/*             is NULL then all rows of X are used.              */
/*                                                               */
/*   Comment:  The calling routine must allocate space for F     */
/*             (which has nRows rows).                           */
/*****************************************************************/

/*****************************************************************/
boolean ModIsXActive(const LinModel *Mod, const real *Beta,
          size_t xIndex);
/*****************************************************************/
/*   Purpose:  Is a specified x variable active in a linear      */
/*             model?                                            */
/*                                                               */
/*   Args:     Mod       Linear model.                           */
/*             Beta      Linear-model coefficients.              */
/*             xIndex    Index of the x variable of interest.    */
/*****************************************************************/

/*****************************************************************/
size_t ModActiveTerms(const LinModel *Mod, const real *Beta,
          size_t nActiveX, const size_t *xIndex,
          size_t *IndexTerm);
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
/*****************************************************************/

/*****************************************************************/
boolean ModIsXActiveInTerm(const LinModel *Mod, const real *Beta,
          size_t xIndex, size_t TermIndex);
/*****************************************************************/
/*   Purpose:  Is a specified x variable active in a specified   */
/*             linear model term?                                */
/*                                                               */
/*   Args:     Mod       Linear model.                           */
/*             Beta      Linear-model coefficients.              */
/*             xIndex    Index of the x variable of interest.    */
/*             TermIndex Index of the term interest.             */
/*****************************************************************/

#define XToF(M,x,f)      XToFActive(M, 0, NULL, x, f)
#define ModFMat(M,X,F)   ModFMatRowIndex(M, 0, NULL, X, F)


/* modfn.c: */

/*****************************************************************/
int ModFnParse(string Comp, int *Fn);
/*****************************************************************/
/*   Purpose:  Parse Comp for a linear-model function.           */
/*             On return, *Fn is the function found, and Comp    */
/*             has all function-related characters removed.      */
/*                                                               */
/*   Returns:  INPUT_ERR or OK.                                  */
/*****************************************************************/

/*****************************************************************/
real ModFn(real x, int fn);
/*****************************************************************/
/*   Purpose:  Return linear-model component when function fn is */
/*             applied to x.                                     */
/*****************************************************************/


/* modparse.c: */

/*****************************************************************/
int ModParse1(size_t nTerms, const string *TermStr,
          const string ModName, LinModel *Mod);
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
/*****************************************************************/

/*****************************************************************/
int ModAddComp(string Comp, Matrix *Term);
/*****************************************************************/
/*   Purpose:  Add a component to a linear model term.           */
/*                                                               */
/*   Returns:  INPUT_ERR or OK.                                  */
/*****************************************************************/

/*****************************************************************/
int ModParseComp(string Comp, size_t *CatLevel, int *Fn);
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
/*   Version:  1993 June 4                                       */
/*****************************************************************/

/*****************************************************************/
int ModParse2(size_t nXVars, const string *xName,
          const size_t *nCats, const string ModName, LinModel *Mod);
/*****************************************************************/
/*   Purpose:  Second parse of linear-model.                     */
/*             Set up x-variable indices and check the           */
/*             categorical level for each component.             */
/*                                                               */
/*   Returns:  INPUT_ERR or OK.                                  */
/*****************************************************************/

/*****************************************************************/
void ModFree(LinModel *Mod);
/*****************************************************************/
/*   Purpose:  Free space allocated for the linear model.        */
/*****************************************************************/
