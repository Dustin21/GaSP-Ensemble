/*****************************************************************/
/*   Copyright (c) William J. Welch 1994--6.                     */
/*   All rights reserved.                                        */
/*                                                               */
/*   Comment:  Can function prototype for DCritVal, etc. be      */
/*             avoided?  They are all the same.                  */
/*                                                               */
/*   Version:  1996.04.20                                        */
/*****************************************************************/

#define DES_MAX_CHANGES  2    /* Latin hypercube algorithm. */

#define ALL_COORD        "AllCoordinates"
#define FEDOROV          "Fedorov"
#define LATIN_HYPERCUBE  "LatinHypercube"
#define SEQUENTIAL       "Sequential"
#define DES_ALG_NAMES    {ALL_COORD, FEDOROV, LATIN_HYPERCUBE, \
                          SEQUENTIAL}
#define ALL_COORD_ALG         0
#define FEDOROV_ALG           1
#define LATIN_HYPERCUBE_ALG   2
#define SEQUENTIAL_ALG        3

#define EVALUATE         "Evaluate"

#define LINK_FN_NAMES \
     {"Identity", "Logit", "Probit", "CLogLog", "LogLog", \
      "Log", "Inverse", "InverseSq", "No1", "No2", "No3"}
#define LINK_IDENTITY    0
#define LINK_LOGIT       1
#define LINK_PROBIT      2
#define LINK_C_LOG_LOG   3
#define LINK_LOG_LOG     4
#define LINK_LOG         5
#define LINK_INVERSE     6
#define LINK_INVERSE_SQ  7
#define LINK_NO_1        8
#define LINK_NO_2        9
#define LINK_NO_3       10

#define VAR_FN_NAMES \
     {"Identity", "Binary", "Linear", "Quadratic", "Cubic", \
      "No1", "No2", "No3"}
#define VAR_FN_IDENTITY  0
#define VAR_FN_BINARY    1
#define VAR_FN_LINEAR    2
#define VAR_FN_QUADRATIC 3
#define VAR_FN_CUBIC     4
#define VAR_FN_NO_1      5
#define VAR_FN_NO_2      6
#define VAR_FN_NO_3      7

#define LIKE_NAMES \
     {"Normal", "InverseGaussian", "Binomial", "Poisson", \
      "Gamma", "No1", "No2", "No3"}
#define LIKE_NORMAL                0
#define LIKE_INVERSE_GAUSSIAN      1
#define LIKE_BINOMIAL              2
#define LIKE_POISSON               3
#define LIKE_GAMMA                 4
#define LIKE_NO_1                  5
#define LIKE_NO_2                  6
#define LIKE_NO_3                  7


/* Used by distance-based criteria. */
#define MAX_METRIC  100.0


/* crit.c: */

/*****************************************************************/
string CritName(size_t CritNum);
/*****************************************************************/
/*   Purpose:  Return the name of criterion CritNum.             */
/*                                                               */
/*   Returns:  Name or NULL.                                     */
/*****************************************************************/

/*****************************************************************/
int CritInput(size_t CritNum);
/*****************************************************************/
/*   Purpose:  Make sure all input for criterion CritNum is in   */
/*             the database.                                     */
/*                                                               */
/*   Returns:  INPUT_ERR, INCOMPAT_ERR, or OK.                   */
/*****************************************************************/

/*****************************************************************/
int CritSetUp(size_t CritNum, size_t nCases, size_t nXVars);
/*****************************************************************/
/*   Purpose:  Set up for for criterion CritNum.                 */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
int CritCleanUp(size_t CritNum, size_t nCases,
     const size_t *RowIndex, const Matrix *X);
/*****************************************************************/
/*   Purpose:  For criterion CritNum, clean up allocations,      */
/*             compute futher properties of the best design X,   */
/*             etc.                                              */
/*****************************************************************/

/*****************************************************************/
int CritVal(size_t CritNum, size_t nCases, const size_t *RowIndex,
          const matrix *X, real *Crit);
/*****************************************************************/
/*   Purpose:  For criterion CritNum, compute the criterion for  */
/*             the nCases X rows in RowIndex.                    */
/*             On return, *Crit is the criterion value.          */
/*****************************************************************/

/*****************************************************************/
int CritUpdate(size_t CritNum, boolean WhatIf, size_t nRowChanges,
     const size_t *RowIndexChange, const real *WtChange,
     size_t nColChanges, const size_t *ColIndexChange,
     const Matrix *X, real *Crit);
/*****************************************************************/
/*   Purpose:  For criterion CritNum, update the criterion when  */
/*             the nRowChanges X rows in RowIndexChange are      */
/*             changed by WtChange (may have negative elements). */
/*             Only the nColChanges columns in ColIndexChange    */
/*             are changed during a downdate-whatif-update cycle.*/
/*             On return, *Crit is the criterion value.          */
/*             If WhatIf is true, then do *not* update working   */
/*             arrays.                                           */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
int CritFit(size_t CritNum, boolean FirstCall);
/*****************************************************************/
/*   Purpose:  For criterion CritNum, fit model parameters.      */
/*****************************************************************/

/*****************************************************************/
real CritCond(size_t CritNum);
/*****************************************************************/
/*   Purpose:  For criterion CritNum, return the condition       */
/*             number.                                           */
/*****************************************************************/


/* critamse.c: */

int  AMSESetUp(const vector *IndexRow, const matrix *X);
int  AMSECleanUp(const vector *IndexRow, const matrix *X);
int  AMSECritVal(const vector *IndexRow, const real *WtAll,
          const matrix *X, real *Crit);
int  AMSEFit(boolean FirstCall);
real AMSECond(void);


/* critcens.c: */

int DCensSetUp(const matrix *X);
int DCensCritVal(size_t n, const Matrix *X, const real *Wt,
          real *Crit);
int DCensUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nChanges, const size_t *RowIndex,
     const real *WtChange, size_t nColChanges,
     const size_t *ColIndex, real *Crit);
real DCensCond(void);
void XToFTMat(const Matrix *X, size_t nChanges,
          const size_t *RowIndex);
real ScaleY(real e, size_t k, real *fRow, real *ParVec);
real LogLikeExp(size_t k, real *Beta, real *dLike);
real LogLikeWei(size_t nPars, real *par, real *dLike);
void AdjustWtExp(size_t n, const real *Wt, real *WtAdj);
void AdjustWtWei(size_t n, const real *Wt, real *WtAdj);
real DetWeiAdj(size_t n);


/* critd.c: */

int DSetUp(const matrix *X);
int DCritVal(size_t n, const Matrix *X, const real *Wt, real *Crit);
int DUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, size_t nColChanges,
     const size_t *ColIndex, real *Crit);

/*****************************************************************/
int DSafeWhatIf(size_t n, const Matrix *X, const real *NewWt,
          size_t nRowChanges, const size_t *RowIndex,
          const real *WtChange, const LinModel *LocalM,
          const Matrix *LocalR, real LocalDet, real *Crit);
/*****************************************************************/
/*   Purpose:  A safe WhatIf calculation for D optimality.       */
/*             It acts on the local variables LocalM, LocalR,    */
/*             and LocalDet (rather than the global RegMod, R,   */
/*             and Det) so that it may be used by other criteria.*/
/*             On return *Crit is the criterion value.           */
/*                                                               */
/*   Comment:  Implemented only for nRowChanges = 1.             */
/*                                                               */
/*   Returns:  OK (if the quick whatif calculation cannot be     */
/*             performed, slower updating is used, so            */
/*             NUMERIC_ERR is trapped and overcome).             */
/*****************************************************************/

/*****************************************************************/
real DCond(void);
/*****************************************************************/
/*   Purpose:  Return condition number of R.                     */
/*****************************************************************/

/*****************************************************************/
real DetFF(Matrix *R, size_t n);
/*****************************************************************/
/*   Purpose:  Return det(first q rows/cols of R) ** (-2/q),     */
/*             where q = min(n, k) and k is the size of R.       */
/*                                                               */
/*   Comment:  The first min(n, k) rows/cols are considered to   */
/*             help sequential algorithms get started.           */
/*             R's numbers of rows and columns may be changed    */
/*             temporarily.                                      */
/*                                                               */
/*   Returns:  det ** (-2/q) if first q diagonal elements > 0.0; */
/*             REAL_MAX      otherwise.                          */
/*****************************************************************/


/* critdglm.c: */

int DGLMSetUp(const matrix *X);
int DGLMCritVal(size_t n, const Matrix *X, const real *Wt,
          real *Crit);
int DGLMUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, size_t nColChanges,
     const size_t *ColIndex, real *Crit);
real DGLMCond(void);
void DGLMAdjustWt(const real *Wt, size_t nWts, const Matrix *X,
          const size_t *RowIndex, real *AdjWt);


/* critds.c: */

int DsSetUp(size_t nMax, size_t NumXVars);
int DsCritVal(size_t n, const Matrix *X, const real *Wt,
     real *Crit);


/* crite.c: */

int ESetUp(const matrix *X);
int ECritVal(size_t n, const Matrix *X, const real *Wt, real *Crit);


/* critl2.c: */

int L2SetUp(size_t nMax, size_t NumXVars);
int L2CritVal(size_t n, const Matrix *X, const real *Wt,
     real *Crit);
int L2Update(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, real *Crit);


/* critg.c */

int GSetUp(const matrix *X);
int GCritVal(size_t n, const Matrix *X, const real *Wt, real *Crit);
int GUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, size_t nColChanges,
     const size_t *ColIndex, real *Crit);
real GCond(void);

/*****************************************************************/
real GMaxfFFf(const Matrix *PredReg);
/*****************************************************************/
/*   Purpose:  calls MinAnyX to maximize the unscaled prediction */
/*             variance, f'Inv(F'F)f over the prediction region. */
/*                                                               */
/*   Returns:  maximum value of f'Inv(F'F)f if F is full rank;   */
/*             MAX_VAR otherwise.                                */
/*****************************************************************/

/*****************************************************************/
real GPredVar(real *x, size_t nXVars);
/*****************************************************************/
/*   Purpose:  Called by MinAnyX to compute -f(x)'Inv(F'F)f(x)   */
/*             for any x.                                        */
/*                                                               */
/*   Returns:  f'Inv(F'F)f if F is full rank;                    */
/*             MAX_VAR otherwise.                                */
/*****************************************************************/


/* critgcd.c: */

int GCDSetUp(size_t nMax, size_t NumXVars);
int GCDCritVal(size_t n, const Matrix *X, const real *Wt,
     real *Crit);
int GCDUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, real *Crit);


/* critcov.c: */

int CoverSetUp(const matrix *X);
int CoverCritVal(size_t n, const Matrix *X, const real *Wt,
     real *Crit);
int CoverUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, size_t nColChanges,
     const size_t *ColIndex, real *Crit);

/*****************************************************************/
void CoverChange(const Matrix *X, size_t nRowChanges,
          const size_t *RowIndex, const real *WtChange,
          boolean WhatIf, Matrix *CoverCount, real *Crit);
/*****************************************************************/
/*   Purpose:  Called by Update to update the criterion.         */
/*                                                               */
/*   Returns:  OK (no error conditions).                         */
/*****************************************************************/

/*****************************************************************/
void CoverProjUpdate(real WtChange, size_t k, size_t p,
          boolean UpdateCoverCount, Matrix *CoverCount);
/*****************************************************************/
/*   Purpose:  Update Noncover for projection p when candidate k */
/*             is added (WtChange = 1) or deleted                */
/*             (WtChange = -1).  If UpdateCoverCount is true,    */
/*             then CoverCount is also updated for projection p. */
/*****************************************************************/

/*****************************************************************/
size_t CoverIdentifyCand(size_t nXVars, real *xRow);
/*****************************************************************/
/*   Purpose:  Return the candidate number of xRow.              */
/*****************************************************************/

/*****************************************************************/
size_t CoverFindCands(size_t k, size_t p, boolean QuickBreak,
          size_t *CoverCountCol);
/*****************************************************************/
/*   Purpose:  Find the candidates covered by candidate k in     */
/*             projection p.                                     */
/*             On return, CoveredCand contains the indices of    */
/*             the covered candidates.                           */
/*                                                               */
/*   Returns:  Number of covered candidates.                     */
/*****************************************************************/


/* critmaxd.c: */

int MaxDistSetUp(const matrix *X);
int MaxDistCritVal(size_t n, const Matrix *X, const real *Wt,
     real *Crit);
int MaxDistUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, size_t nColChanges,
     const size_t *ColIndex, real *Crit);

/*****************************************************************/
void MaxDistInit(void);
/*****************************************************************/
/*   Purpose:  Initialize SumRecip and NumZeroDists.             */
/*****************************************************************/

/*****************************************************************/
real MaxDistChange(size_t n, const Matrix *X, const real *NewWt,
     size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, Matrix *LocalSumRecip,
     Matrix *LocalNumZeroDists);
/*****************************************************************/
/*   Purpose:  Called by MaxDistUpdate to carry out the update.  */
/*                                                               */
/*   Returns:  OK (no error conditions).                         */
/*****************************************************************/


/* critmind.c: */

int MinDistSetUp(const matrix *X);
int MinDistCritVal(size_t n, const Matrix *X, const real *Wt,
     real *Crit);
int MinDistUpdate(boolean WhatIf, size_t n, const Matrix *X,
     const real *NewWt, size_t nRowChanges, const size_t *RowIndex,
     const real *WtChange, size_t nColChanges,
     const size_t *ColIndex, real *Crit);

/*****************************************************************/
void MinDistChange(size_t n, size_t i, size_t ProjDim,
          const size_t *xIndex, real WtChange);
/*****************************************************************/
/*   Purpose:  Update/downdate for point i and a given           */
/*             projection.                                       */
/*****************************************************************/


/* critrff.c: */

/*****************************************************************/
int RFFInit(size_t n, const Matrix *X, const real *Wt,
          const LinModel *M, Matrix *R);
/*****************************************************************/
/*   Purpose:  Initialize Cholesky R of F'F.                     */
/*                                                               */
/*   Returns:  NUMERIC_ERR if downdate cannot be performed;      */
/*             OK          otherwise.                            */
/*                                                               */
/*   Comment:  NUMERIC_ERR can only result from negative         */
/*             weights, which should not be used for             */
/*             initialization.                                   */
/*             Calling routine must allocate space for R.        */
/*****************************************************************/

/*****************************************************************/
int RFFUpdate(size_t n, const Matrix *X, const real *NewWt,
          size_t nRowChanges, const size_t *RowIndex,
          const real *WtChange, const LinModel *M, Matrix *R);
/*****************************************************************/
/*   Purpose:  Update Cholesky R of F'F when the nRowChanges X   */
/*             rows in RowIndex change by WtChange (may have     */
/*             negative elements).                               */
/*                                                               */
/*   Returns:  OK (if a downdate cannot be performed, R is       */
/*             initialized from scratch, so NUMERIC_ERR is       */
/*             trapped and overcome).                            */
/*                                                               */
/*   Comment:  An infinite recursion could arise if NewWt        */
/*             contains negative elements - this is              */
/*             theoretically impossible.                         */
/*****************************************************************/


/* critutil.c: */

/*****************************************************************/
real *NormalizeWts(real Metric, Matrix *Reg);
/*****************************************************************/
/*   Purpose:  Normalize weights in Reg so that each x variable  */
/*             has a distance of 1 between min and max.          */
/*****************************************************************/

/*****************************************************************/
size_t ProjGenerate(size_t nXVars, size_t ProjDimMin,
          size_t ProjDimMax, size_t **ProjDim, Matrix *Proj);
/*****************************************************************/
/*   Purpose:  Generate projections of dimensions                */
/*             ProjDimMin,..., ProjDimMax.                       */
/*             On return, *ProjDim[p] contains the dimension of  */
/*             projection p, and column p of Proj contains the   */
/*             indices of the x variables in projection p.       */
/*                                                               */
/*   Returns:  Number of projections.                            */
/*                                                               */
/*   Comment:  ProjDim and Proj are allocated here.              */
/*****************************************************************/

/*****************************************************************/
size_t ProjGenSingleDim(size_t nXVars, size_t d, Matrix *Proj);
/*****************************************************************/
/*   Purpose:  Add projections of dimension d to Proj.           */
/*                                                               */
/*   Returns:  Number of new projections.                        */
/*****************************************************************/

/*****************************************************************/
void OneDimDist(size_t nCats, real x, const real *xVec, size_t n,
     real Wt, real Metric, real *Dist);
/*****************************************************************/
/*   Purpose:  Compute the one-dimensional distances             */
/*                  Dist[i] = [Wt * (x - xVec[i])] ** Metric     */
/*             for i = 0,...,n.  Note that Wt is applied after   */
/*             exponentiating.                                   */
/*****************************************************************/

/*****************************************************************/
void MultiDist(size_t nCases, size_t nCols, const size_t *ColIndex,
     const Matrix *MargDist, boolean MetricIsMax, real *Dist);
/*****************************************************************/
/*   Purpose:  Put multidimensional distances in Dist.           */
/*                                                               */
/*   Comment:  Dist must be allocated by the calling routine.    */
/*****************************************************************/

/*****************************************************************/
int CritCorpar(boolean FirstCall, KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Set up correlation parameters in KrigMod, by      */
/*             reading them or executing a fit command.          */
/*****************************************************************/


/* desall.c: */

/*****************************************************************/
int AllCoord(size_t CritNum, size_t n, size_t nProtected,
          const Matrix *ExpReg, real TolAbs, real TolRel,
          Matrix *X, real *BestCrit, real *CondNum);
/*****************************************************************/
/*   Purpose:  Generate an n-point design by direct              */
/*             minimization over all coordinates.                */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
real AllObj(real *StackedX, size_t NumCoords);
/*****************************************************************/
/*   Purpose:  Called by optimizer to evaluate a design.         */
/*                                                               */
/*   Returns:  New criterion.                                    */
/*****************************************************************/


/* desfed.c: */

/*****************************************************************/
int Fedorov(size_t CritNum, size_t n, size_t nProtected,
          const Matrix *ExpReg, real TolAbs, real TolRel,
          Matrix *X, real *BestCrit, real *CondNum);
/*****************************************************************/
/*   Purpose:  Generate an n-point design using modified         */
/*             "Fedorov" exchanges.                              */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/* deslhs.c: */

/*****************************************************************/
int LHSOpt(size_t CritNum, size_t n, size_t nProtected,
          const Matrix *ExpReg, real TolAbs, real TolRel,
          Matrix *X, real *BestCrit, real *CondNum);
/*****************************************************************/
/*   Purpose:  Generate an "optimal" n-point design within the   */
/*             class of Latin hypercubes.                        */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/* desseq.c: */

/*****************************************************************/
int Seq(size_t CritNum, size_t n, size_t nProtected, size_t Tries,
          const Matrix *ExpReg, real TolAbs, real TolRel,
          Matrix *X, real *Crit, real *CondNum);
/*****************************************************************/
/*   Purpose:  Generate an n-point design sequentially.          */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/* desutil.c: */

/*****************************************************************/
void DesInit(size_t n, const Matrix *ExpReg, Matrix *X);
/*****************************************************************/
/*   Purpose:  Allocate/reallocate the X matrix.                 */
/*****************************************************************/

/*******************************+++*******************************/
int DesStart(size_t CritNum, size_t nProtected, size_t n,
     const Matrix *ExpReg, Matrix *X, real *Crit, real *CondNum);
/*****************************************************************/
/*   Purpose:  Put random points in rows nProtected,..., n - 1   */
/*             of X such that the criterion is well conditioned. */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*******************************+++*******************************/
void DesRand(size_t nProtected, size_t n, const Matrix *ExpReg,
          Matrix *X);
/*****************************************************************/
/*   Purpose:  Put random points in rows nProtected,..., n - 1   */
/*             of X.                                             */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
int DesBestAdd(size_t CritNum, size_t i, size_t n, Matrix *X,
          const real *NewWt, const Matrix *ExpReg, real TolAbs,
          real TolRel, real *BestCrit);
/*****************************************************************/
/*   Purpose:  Find the best design point for row i of X, and    */
/*             update working arrays for the change.             */
/*             On entry, row i of X is the starting point, and   */
/*             *BestCrit is the corresponding criterion value.   */
/*             On return, the best new point for row i of X has  */
/*             been put in X, and *BestCrit is the best          */
/*             criterion value found.                            */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
void DesFindBestAdd(size_t CritNum, size_t i, size_t n, Matrix *X,
          const real *NewWt, const Matrix *ExpReg, real TolAbs,
          real TolRel, real *BestCrit);
/*****************************************************************/
/*   Purpose:  Find the best design point for row i of X,        */
/*             without updating working arrays for the change.   */
/*             On entry, row i of X is the starting point, and   */
/*             *BestCrit is the corresponding criterion value.   */
/*             On return, the best new point for row i of X has  */
/*             been put in X, and *BestCrit is the best          */
/*             criterion value found.                            */
/*****************************************************************/

/*****************************************************************/
real DesAddObj(real *x, size_t nXVars);
/*****************************************************************/
/*   Purpose:  Called by optimizer to evaluate an added point.   */
/*                                                               */
/*   Returns:  New criterion.                                    */
/*****************************************************************/

