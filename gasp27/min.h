/*****************************************************************/
/*   Copyright (c) William J. Welch 1990--96.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*   Version: 1996.06.16                                        */
/*****************************************************************/

#define SIMPLEXALG 0
#define POWELLALG  1

/* min.c: */

/*****************************************************************/
unsigned MinAnyX(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          const Matrix *XReg, size_t nDims, int MinAlg, real *x,
          real *Obj);
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
/*             nDims   Number of dimensions.                   */
/*             MinAlg    The minimizer called for unconstrained  */
/*                       minimization.                           */
/*             x         Input:  Starting point (including any   */
/*                               fixed variables);               */
/*                       Output: "Optimal" point.                */
/*             Obj       Input:  Objective at x;                 */
/*                       Output: "Optimal" objective.            */
/*                                                               */
/*   Returns:  Total number of function evaluations.             */
/*****************************************************************/

/*****************************************************************/
unsigned MinDisc(size_t nDims, const size_t *VarIndex,
     const Matrix *XReg, real *x, real *Obj);
/*****************************************************************/
/*   Purpose:  Minimize over the nDims DISCRETE or GRID        */
/*             variables with indices in VarIndex.               */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*****************************************************************/

/*****************************************************************/
real ObjCont(real *xCont, size_t nContVars);
/*****************************************************************/
/*   Purpose:  Put continuous variables in ExtX (x), and call    */
/*             objective function.                               */
/*                                                               */
/*   Returns:  Objective.                                        */
/*****************************************************************/

/*****************************************************************/
unsigned MinMultiStart(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          const Matrix *XReg, size_t nDims, int MinAlg,
          const Matrix *StartPt, real *x, real *Obj);
/*****************************************************************/
/*   Purpose:  Multiple local searches (via MinAnyX) from each   */
/*             row in StartPt.                                   */
/*                                                               */
/*   Returns:  Total number of function evaluations.             */
/*****************************************************************/


/* minconjg.c: */

/*****************************************************************/
unsigned MinConjGrad(real (*ObjFunc)(size_t n, real *x, real *g),
     real AbsTol, real RelTol, unsigned MaxEvals,
     size_t n, real *x, real *fx);
/*****************************************************************/
/*   Purpose:  Multidimensional minimization using conjugate     */
/*             gradients.                                        */
/*                                                               */
/*   Args:     ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxEvals  Maximum function evaluations.           */
/*             n         Number of dimensions.                   */
/*             x         Output: "Optimal" point.                */
/*             fx        Output: "Optimal" function value.       */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*                                                               */
/*   Comments: x and fx need not be set before calling.  This    */
/*             should be changed for compatibility with other    */
/*             minimization routines.                            */
/*             Objective functions should be fixed up to deal    */
/*             with gradient/no-gradient methods.                */
/*****************************************************************/

real ObjFuncNoGrad(real *x, size_t n);


/* mincont.c: */

unsigned  MinCont(real (*ObjFunc)(real *x, size_t nDims),
               real AbsTol, real RelTol, unsigned MaxFuncs,
               real *LowBnd, real *UpBnd, size_t *Distrib,
               size_t nDims, int MinAlg, real *x, real *fx);
real      XToUncon(real x, real a, real b);
real      UnconToX(real u, real a, real b);

/*****************************************************************/
real ObjFuncUncon(real *xUncon, size_t nDims);
/*****************************************************************/
/*   Purpose:  Convert unconstrained x's to their constrained    */
/*             ranges, and call objective.                       */
/*                                                               */
/*   Returns:  Objective.                                        */
/*****************************************************************/

/* minextra.c: */

unsigned  MinExtrap(real (*ObjFunc)(real *x, size_t nDims),
               const Matrix *Reg, size_t nDims, const real *xOld,
               real *xNew, real *Obj);

/*****************************************************************/
boolean Extrap(size_t nDims, const real *xOld, const real *xNew,
          real Gamma, const Matrix *Reg, real *xExtrap);
/*****************************************************************/
/*   Purpose:  Extrapolate xOld through xNew by factor Gamma,    */
/*             keeping the extrapolated point in the feasibility */
/*             region, Reg.                                      */
/*                                                               */
/*   Args:     nDims     Dimension of points.                    */
/*             xOld      Old point.                              */
/*             xNew      New point.                              */
/*             Gamma     Extrapolation factor.                   */
/*             Reg       Feasiblity region.                      */
/*             xExtrap   Output: Extrapolated point.             */
/*                                                               */
/*   Returns:  YES if there is at least one continuous variable  */
/*                 that can be extrapolated;                     */
/*             NO  otherwise.                                    */
/*****************************************************************/

/*****************************************************************/
unsigned MinTryBounds(real (*ObjFunc)(real *x, size_t nDims),
          size_t nDims, const real *LowBnd, const real *UpBnd,
          real *x, real *Obj);
/*****************************************************************/
/*   Purpose:  For each (continuous) variable close to one of    */
/*             its bounds, try putting it at the bound.          */
/*                                                               */
/*   Args:     ObjFunc   Calulates objective function.           */
/*             LowBnd    Lower bounds on x.                      */
/*             UpBnd     Upper bounds on x.                      */
/*             nDims     Dimension of points.                    */
/*             x         Input:  Starting point.                 */
/*                       Output: A possibly better point.        */
/*             Obj       Input:  The objective for x.            */
/*                       Output: A possibly better objective.    */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*****************************************************************/


/* minone.c: */

/*****************************************************************/
unsigned MinBracket(real (*ObjFunc)(real x), real *a, real *x,
          real *b, real *fa, real *fx, real *fb);
/*****************************************************************/
/* Purpose:    Bracket a minimum.                                */
/*                                                               */
/* Args:       ObjFunc   The objective function.                 */
/*             a, x      Input:  Initial distinct points         */
/*                               (a < x).                        */
/*             a, x, b   Output: Points bracketing the minimum.  */
/*             fx        Input:  Function value at x.            */
/*             fa, fx,   Output: Function values at a, x, and b. */
/*             fb                                                */
/*                                                               */
/* Returns:    Number of function evaluations.                   */
/*****************************************************************/

unsigned  Brent(real (*ObjFunc)(real x), real AbsTol, real RelTol,
               unsigned MaxFuncs, real *a, real *x, real *b,
               real *fa, real *fx, real *fb);

/*****************************************************************/
unsigned MinLine(real (*ObjFunc)(real *x, size_t nDims),
          real AbsTol, real RelTol, unsigned MaxFuncs,
          size_t nDims, real *x, real *d, real *fx);
/*******************************+++*******************************/
/* Purpose:    Minimize a multi-dimensional function along the   */
/*             line from x in the direction d.                   */
/*                                                               */
/* Args:       ObjFunc   The objective function.                 */
/*             AbsTol    Absolute tolerance on function value    */
/*                       for convergence.                        */
/*             RelTol    Relative tolerance on function value    */
/*                       for convergence.                        */
/*             MaxFuncs  Maximum function evaluations.           */
/*             nDims     Number of dimensions.                   */
/*             x         Input:  Starting point.                 */
/*                       Output: x giving minimum.               */
/*             d         Input:  Direction of line.              */
/*                       Output: Vector displacement that x was  */
/*                               moved.                          */
/*             fx        Input:  Function value at x.            */
/*                       Output: Minimum function value.         */
/*                                                               */
/* Returns:    Number of function evaluations.                   */
/*****************************************************************/

real      f1dim(real alpha);


/* minpow.c: */

unsigned  Powell(real (*ObjFunc)(real *x, size_t nDims),
               real AbsTol, real RelTol, unsigned MaxFuncs,
               size_t nDims, real *x, real **d, real *fx);


/* minsimp.c: */

unsigned  Simplex(real (*ObjFunc)(real *x, size_t nDims),
               real AbsTol, real RelTol, unsigned MaxFuncs,
               size_t nDims, real **p, real *y);
real      SimpTry(real **p, real *y, real *psum, size_t nDims,
               real (*ObjFunc)(real *x, size_t nDims), size_t ihi,
               unsigned *NumFuncs, real fac);
