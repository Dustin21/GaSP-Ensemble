/*****************************************************************/
/*   Copyright (c) William J. Welch 1994--2009.                  */
/*   All rights reserved.                                        */
/*                                                               */
/*   2009.05.14: Multiple correlation families                   */
/*****************************************************************/

#define KRIG_MOD_DEFINED

typedef struct
{
     /* X not saved as available via F and G. */
     Matrix    *T;       /* Each of t columns is a transformation. */
     real      *Y;       /* Responses. */
     string    yName;

     LinModel  *RegMod;       /* Regression model. */
     LinModel  *SPMod;        /* Model for the stochastic process. */
     size_t    CorFam;        /* Correlation family */
     boolean   RanErr;        /* Is random error present? */

     Matrix    F;             /* Expanded design matrix for the */
                              /* regression model.              */

     Matrix    G;             /* Expanded design matrix for the */
                              /* stochastic-process model.      */

     /* Used to exploit regular spacing in G. */

     Matrix    Steps;         /* Steps[i,j] is the number of    */
                              /* steps from the min in column j */
                              /* to G[i,j].                     */
     size_t    *MaxSteps;     /* MaxSteps[j] is the maximum */
                              /* of column j of Steps.      */
     Matrix    Dist;          /* Dist[i,j] = (i + 1) * (step length) */
                              /* for i = 0,..., MaxSteps[j] - 1.     */

     /* If the spacing in column j of G is irregular, */
     /* then MaxSteps[j] = 0.                         */

     /* Fitted values of kriging-model parameters. */

     Matrix    CorPar;        /* Correlation parameters. */
     real      SigmaSq;       /* Total variance. */
     real      SPVarProp;     /* Proportion of SigmaSq due to */
                              /* the stochastic process.      */

     /* Correlation matrix and decompositions change */
     /* when fitted parameters change.               */

     Matrix    C;             /* Original correlation matrix. */
     Matrix    Chol;          /* Upper-triangular t x t Cholesky */
                              /* factor (Chol'Chol = T'CT).   */

     Matrix    Q;             /* QR = Inverse(Chol') F. */
     Matrix    R;

     real      *RBeta;        /* Solve R * Beta = RBeta */
                              /* to get Beta.           */
     real      *Beta;         /* Regression model beta's */
     real      *ResTilde;     /* Inverse(Chol') * generalized LS */
                              /* residuals.                      */

     /* Workspace. */
     real      *xRow;
     real      *fRow;
     real      *fr;
     real      *gRow;
     real      *r;
     real      *w1;
     real      *w2;
} KrigingModel;


#define KrigT(M)         ((M)->T)
#define KrigY(M)         ((M)->Y)
#define KrigYName(M)     ((M)->yName)
#define KrigRegMod(M)    ((M)->RegMod)
#define KrigSPMod(M)     ((M)->SPMod)
#define KrigCorFam(M)    ((M)->CorFam)
#define KrigRanErr(M)    ((M)->RanErr)
#define KrigF(M)         (&(M)->F)
#define KrigG(M)         (&(M)->G)
#define KrigSteps(M)     (&(M)->Steps)
#define KrigDist(M)      (&(M)->Dist)

#define KrigCorPar(M)    (&(M)->CorPar)

#define KrigC(M)         (&(M)->C)
#define KrigChol(M)      (&(M)->Chol)
#define KrigQ(M)         (&(M)->Q)
#define KrigR(M)         (&(M)->R)

#define COR_FAM_NAMES    {POW_EXP, MATERN}
#define COR_FAM_POW_EXP  0
#define COR_FAM_MATERN   1

/* krcorpar.c: */

/*******************************+++*******************************/
void CorParAlloc
(
     size_t         CorFam,        /* Correlation family         */
     size_t         NumTerms,      /* Number of terms.           */
     const string   *TermName,     /* Term names ("x1", etc.)    */
     Matrix         *CorPar        /* Output: allocated and      */
                                   /* labelled correlation-      */
                                   /* parameter matrix.          */
);
/*****************************************************************/
/*   Purpose:  Allocate and label a correlation-parameter matrix.*/
/*****************************************************************/

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
);
/*****************************************************************/
/* Purpose: Return starting values for the correlation           */
/*          parameters and their optimization region.            */
/*****************************************************************/

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
);
/*****************************************************************/
/* Purpose: Test parameters against "null" values for a single   */
/*          term.                                                */
/*                                                               */
/* Returns: Number of function evaluations.                      */
/*****************************************************************/

/*****************************************************************/
int CorParSetUp(const Matrix *CorPar, const string yName,
     real SPVar, real ErrVar, KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Set up correlation parameters.                    */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
int CorParExtract(const Matrix *CorPar, const string YName,
          boolean ErrorMessage, Matrix *CorParYName);
/*****************************************************************/
/*   Purpose:  Extract the correlation parameters for YName from */
/*             CorPar and put them into CorParYName.             */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if CorPar does not include the       */
/*                          the columns for YName's parameters   */
/*                          (CorParYName is MatFree'd);          */
/*             OK           otherwise.                           */
/*                                                               */
/*   Comment:  Calling routine must allocate space for           */
/*             CorParYName.                                      */
/*             Legality of the values in CorPar and              */
/*             compatibility of CorPar with the stochastic-      */
/*             process model are assumed.                        */
/*****************************************************************/

/*******************************+++*******************************/
boolean CorParIsActive
(
     size_t       CorFam,    /* Correlation family               */
     const Matrix *Corpar,  /* Correlation parameters           */
     size_t       TermIndex  /* Index of the term of interest    */
);
/*****************************************************************/
/* Purpose: Is term TermIndex active in the correlation          */
/*             function?                                         */
/*****************************************************************/

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
);
/*****************************************************************/
/* Purpose: Compute correlations between the point g and the     */
/*             points in the first n rows of G.                  */
/*****************************************************************/


/* kriging.c: */

/*******************************+++*******************************/
void KrigModAlloc(size_t nCases, size_t nXVars, const string yName,
     const Matrix *T, const LinModel *RegMod,
     const LinModel *SPMod, size_t CorFam, boolean RanErr,
     KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Initialize KrigMod, and allocate F, G, etc.       */
/*****************************************************************/

/*****************************************************************/
void KrigModAllocT(KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Further initializations and allocations from      */
/*             transformations T.                                */
/*****************************************************************/

/*****************************************************************/
void KrigModFree(KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Free kriging model.                               */
/*****************************************************************/

/*****************************************************************/
void KrigModFreeT(KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Free space allocated because of transformations T.*/
/*****************************************************************/

/*****************************************************************/
void KrigModData(size_t nCases, const size_t *RowIndex,
     const Matrix *X, const real *y, KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Set up y, F, G, and call KrigGSpacing.            */
/*****************************************************************/

/*****************************************************************/
int KrigModSetUp(const Matrix *CorPar, const string yName,
     real SPVar, real ErrVar, KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Set up parameters and decompositions.             */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
void KrigGSpacing(KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  Set up Steps and MaxSteps.                        */
/*****************************************************************/

/*****************************************************************/
void KrigCorMat
(
     size_t       NumActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     KrigingModel *KrigMod
);
/*****************************************************************/
/*   Purpose:  Put the correlation matrix into Chol.             */
/*****************************************************************/

/*****************************************************************/
void KrigCorC
(
     size_t       NumActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     KrigingModel *KrigMod,
     Matrix       *C
);
/*****************************************************************/
/*   Purpose:  Put the correlation matrix for the original       */
/*             responses into C (which could be a member of      */
/*             KrigMod).                                         */
/*****************************************************************/

/*****************************************************************/
int KrigDecompose(KrigingModel *KrigMod);
/*****************************************************************/
/*   Purpose:  (Re-)compute decompositions for kriging model.    */
/*                                                               */
/*   Return:   NUMERIC_ERR if C or Inverse(Chol') * F are not    */
/*                         full rank;                            */
/*             OK          otherwise.                            */
/*                                                               */
/*   Comment:  Chol must already hold the correlation matrix and */
/*             is overwritten.                                   */
/*****************************************************************/

int KrigSolve(const Matrix *Chol, const Matrix *F, const real *Y,
          Matrix *FTilde, real *YTilde);

/*****************************************************************/
void KrigTranformVec(KrigingModel *KrigMod, real *r);
/*****************************************************************/
/*   Purpose:  Overwrite r with T'r.                             */
/*****************************************************************/

/*****************************************************************/
boolean KrigIsXActive(const KrigingModel *KrigMod, size_t j);
/*****************************************************************/
/*   Purpose:  Determine whether x variable j is active,         */
/*             in either the regression model or the stochastic- */
/*             process model.                                    */
/*****************************************************************/

/*****************************************************************/
size_t KrigSPActiveTerms(const KrigingModel *KrigMod,
          size_t nActiveX, const size_t *xIndex, size_t *IndexTerm);
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
/*****************************************************************/

/*****************************************************************/
void frfrAve(KrigingModel *KrigMod, const Matrix *PredReg,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     matrix *frfrj, matrix *frfr);
/*****************************************************************/
/*   Purpose:  Compute average fr(fr)^T over a region.           */
/*                                                               */
/*   Comment:  Calling routine must allocate space for           */
/*             (k + n) * (k + n) matrices frfrj (used for        */
/*             workspace) and frfr.                              */
/*****************************************************************/

/*****************************************************************/
void fgrGroup(const KrigingModel *KrigMod, const Matrix *PredReg,
     size_t nXVars, const size_t *xIndex, size_t Level,
     size_t nSPTerms, const size_t *IndexSP, real *xRow, real *f,
     real *g, real *r);
/*****************************************************************/
/*   Purpose:  Compute f, g, and r for a group.                  */
/*                                                               */
/*   Comment:  f, g, r, and xRow (used for workspace) must be    */
/*             allocated by the calling routine.                 */
/*****************************************************************/


/* krmle.c: */

/*****************************************************************/
void MLEStart(KrigingModel *KrigMod, Matrix *RegCorPar);
/*****************************************************************/
/*   Purpose:  Put random values of correlation parameters and   */
/*             SPVarProp in KrigMod, and set up optimization     */
/*             region, RegCorPar (which includes SPVarProp).     */
/*                                                               */
/*   Returns:  OK or error condition.                            */
/*****************************************************************/

/*****************************************************************/
int MLEFit
(
     Matrix       *RegCorPar, /* Unchanged.                      */
     KrigingModel *KrigMod,
     real         LogLikeTol,
     real         CritLogLikeDiff,
     size_t       Try,
     real         *NegLogLike,/* Output: -log(likelihood).       */
     real         *CondNum,   /* Output: condition number.       */
     unsigned     *TotFuncs   /* Output: likelihood evaluations. */
);
/*****************************************************************/
/*   Purpose:  Fit regression and correlation parameters.        */
/*                                                               */
/*   Returns:  OK or error condition.                            */
/*****************************************************************/

/*****************************************************************/
unsigned MLEOneTerm
(
     real AbsTol,             /* Absolute convergence tolerance  */
                              /* on -log(likelihood).            */
     const Matrix *RegCorPar, /* Optimization region for all     */
                              /* parameters.                     */
     real  *NegLogLike        /* Input: Current -log(likelihood);*/
                              /* Output: Optimized value.        */
);
/*****************************************************************/
/*   Purpose:  Optimize correlation parameters for term          */
/*             TermIndex (external), other parameters remaining  */
/*             fixed.                                            */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*****************************************************************/

/*****************************************************************/
real MLELike(void);
/*****************************************************************/
/*   Purpose:  Computes the negative log likelihood given by     */
/*             1/2 [log det (C) + n log sigma hat squared].      */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*                                                               */
/*   Comment   Correct correlation matrix must be already loaded */
/*             into ExtKrigMod->Chol.                            */
/*****************************************************************/

/*****************************************************************/
real MLELikeUpdate
(
     real   *CorParRow,  /* Input: Correlation parameters for a  */
                         /* single term (one row of CorPar).     */
     size_t NumPars      /* Number of parameters for a single    */
                         /* term.                                */
);
/*****************************************************************/
/*   Purpose:  Called by optimizer to update -log(likelihood)    */
/*             when the correlation parameters change only for a */
/*             single stochastic-process term.                   */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*****************************************************************/

/*****************************************************************/
real MLELikeObj
(
     real   *CorParVec,  /* Input: Correlation parameters        */
                         /* arranged as a vector for all terms,  */
                         /* possibly followed by SPVarProp.      */
     size_t NumPars      /* Number of parameters.                */
);
/*****************************************************************/
/*   Purpose:  Called by optimizer to compute -log likelihood.   */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*****************************************************************/

/*****************************************************************/
real MLELikeScale
(
     real   *SPVarProp,  /* Scaling proportion for the   */
                         /* stochastic-process variance. */
     size_t OneInput     /* Not used.                    */
);
/*****************************************************************/
/*   Purpose:  Called by optimizer to compute -log likelihood    */
/*             when the off-diagonal elements of CPartial are    */
/*             re-scaled for *SPVarProp.                         */
/*                                                               */
/*   Returns:  -log(likelihood)                                  */
/*****************************************************************/


/* krmatern.c: */

/*******************************+++*******************************/
void MaternAlloc
(
     size_t         NumTerms,      /* Number of terms.           */
     Matrix         *CorPar        /* Output: allocated and      */
                                   /* labelled correlation-      */
                                   /* parameter matrix.          */
);
/*****************************************************************/
/* Purpose: Allocate correlation matrix and label columns.       */
/*****************************************************************/

/*******************************+++*******************************/
void MaternStart
(
     const Matrix *G,    /* Expanded-design matrix for the       */
                         /* stochastic-process model.            */
     Matrix *CorPar,     /* Output: Starting values of the       */
                         /* correlation parameters.              */
     Matrix *CorReg      /* Output: Feasibility region for the   */
                         /* correlation parameters.              */
);
/*****************************************************************/
/* Purpose:    Return starting values for the correlation        */
/*             parameters and their optimization region.         */
/*****************************************************************/

/*******************************+++*******************************/
void MaternCor(
     const real   *g,        /* A point.                         */
     const Matrix *G,        /* Matrix of points.                */
     size_t       n,         /* The correlations for only the    */
                             /* first n rows of G are computed.  */
     size_t       NumActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     const Matrix *CorPar,   /* Correlation parameters.          */
     real         *Cor       /* Output: correlations.            */
);
/*****************************************************************/
/* Purpose:  Compute correlations between the point g and the    */
/*           points in the first n rows of G.                    */
/*****************************************************************/

/*******************************+++*******************************/
void MaternCorOneDim(real h, const real *g, size_t n, real theta,
          real deriv, real *Cor);
/*****************************************************************/
/* Purpose:  Multiply the correlations Cor[0],...,Cor[n-1] by    */
/*           the 1-d Matern correlations from the distances      */
/*            between h and g[0],...,g[n-1].                     */
/*****************************************************************/

/*******************************+++*******************************/
unsigned MaternTest
(
     Matrix *CorReg,     /* Feasibility region for the           */
                         /* correlation parameters.              */
     size_t TermIndex,   /* Index of the tested term.            */
     real   AbsTol,
     real   CritLogLikeDiff,
     Matrix *CorPar,     /* Input:  Correlation parameters;      */
                         /* Output: Row TermIndex may change.    */
     real   *NegLogLike  /* Input: Negative log likelihood;      */
                         /* Output: New value.                   */
);
/*****************************************************************/
/* Purpose:  Test whether deriv = derivMax and/or theta = 0 for  */
/*           a single term.                                      */
/*                                                               */
/* Returns:  Number of function evaluations.                     */
/*****************************************************************/

/*******************************+++*******************************/
boolean MaternIsActive
(
     const Matrix *CorPar,  /* Correlation parameters           */
     size_t       TermIndex  /* Index of the term of interest    */

);
/*****************************************************************/
/* Purpose: Is term TermIndex active in the correlation          */
/*             function?                                         */
/*****************************************************************/


/* krpowexp.c: */

/*****************************************************************/
void PEStart
(
     const Matrix *G,    /* Expanded-design matrix for the       */
                         /* stochastic-process model.            */
     Matrix *CorPar,     /* Output: Starting values of the       */
                         /* correlation parameters.              */
     Matrix *CorReg      /* Output: Feasibility region for the   */
                         /* correlation parameters.              */
);
/*****************************************************************/
/* Purpose:    Return starting values for the correlation        */
/*             parameters and their optimization region.         */
/*                                                               */
/* Comment:    CorReg is allocated here.                         */
/*****************************************************************/

/*****************************************************************/
void PEAlloc
(
     size_t         NumTerms,      /* Number of terms.           */
     Matrix         *CorPar        /* Output: allocated and      */
                                   /* labelled correlation-      */
                                   /* parameter matrix.          */
);
/*****************************************************************/
/*   Purpose:  Allocate and label a correlation-parameter matrix.*/
/*****************************************************************/

/*****************************************************************/
void PECor(
     const real   *g,        /* A point.                         */
     const Matrix *G,        /* Matrix of points.                */
     size_t       n,         /* The correlations for only the    */
                             /* first n rows of G are computed.  */
     size_t       NumActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     const Matrix *CorPar,   /* Correlation parameters.          */
     real         *Cor       /* Output: correlations.            */
);
/*****************************************************************/
/*   Purpose:  Compute correlations between the point g and the  */
/*             points in the first n rows of G.                  */
/*****************************************************************/

/*****************************************************************/
void PEDist
(
     const real   *g,        /* A point.                         */
     const Matrix *G,        /* Matrix of points.                */
     size_t       n,         /* The distances for only the first */
                             /* n rows of G are computed.        */
     size_t       NumActive, /* Number of active terms           */
                             /* (only used if Active != NULL).   */
     const size_t *Active,   /* If != NULL, then contains the    */
                             /* indices of the active terms.     */
     const Matrix *CorPar,   /* Correlation parameters.          */
     real         *Dist      /* Output: distances.               */
);
/*****************************************************************/
/*   Purpose:  Compute distances from the point g to the points  */
/*             in the first n rows of G.                         */
/*****************************************************************/

/*****************************************************************/
void PEDistInc(real h, const real *g, size_t n, real Theta,
          real Alpha, real *Dist);
/*****************************************************************/
/*   Purpose:  Increment the distances Dist[0],...,Dist[n-1] for */
/*             the distances between h and g[0],...,g[n-1].      */
/*****************************************************************/

/*****************************************************************/
unsigned PETest
(
     Matrix *CorReg,     /* Feasibility region for the           */
                         /* correlation parameters.              */
     size_t TermIndex,   /* Index of the tested term.            */
     real   AbsTol,
     real   CritLogLikeDiff,
     Matrix *CorPar,     /* Input:  Correlation parameters;      */
                         /* Output: Row TermIndex may change.    */
     real   *NegLogLike  /* Input: Negative log likelihood;      */
                         /* Output: New value.                   */
);
/*****************************************************************/
/*   Purpose:  Test whether Alpha = 0 and/or Theta = 0 for a     */
/*             single term.                                      */
/*                                                               */
/*   Returns:  Number of function evaluations.                   */
/*****************************************************************/

/*******************************+++*******************************/
boolean PEIsActive
(
     const Matrix *Corpar,  /* Correlation parameters           */
     size_t       TermIndex  /* Index of the term of interest    */

);
/*****************************************************************/
/* Purpose: Is term TermIndex active in the correlation          */
/*             function?                                         */
/*****************************************************************/


/* krpred.c: */

/*****************************************************************/
int KrigPredSetUp
(
     const KrigingModel  *KrigMod, /* Kriging model and          */
                                   /* decompositions.            */
     real           *ResTildeTilde /* Output: Inv(C) * GLS       */
                                   /* residuals.                 */
);
/*****************************************************************/
/*   Purpose:  Set up for computing predictions *without*        */
/*             standard errors.                                  */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for the       */
/*             vector ResTildeTilde (length n).                  */
/*             For generating prediction coefficients, this      */
/*             seems to be unstable.                             */
/*****************************************************************/

/*****************************************************************/
void KrigPred(KrigingModel *KrigMod, const Matrix *XPred,
      const real *ResTildeTilde, real *YHat);
/*****************************************************************/
/*   Purpose:  Compute kriging predictions *without* standard    */
/*             errors for a matrix of points.                    */
/*                                                               */
/*   Comment:  Calling routine must allocate space for YHat.     */
/*****************************************************************/

/*****************************************************************/
int KrigPredSE(KrigingModel *KrigMod, const Matrix *XPred,
          real *YHat, real *SE);
/*****************************************************************/
/*   Purpose:  Compute kriging predictions *with* standard       */
/*             errors for a matrix of points.                    */
/*                                                               */
/*   Args:     KrigMod   Kriging model and decompositions.       */
/*             XPred     Prediction points.                      */
/*             YHat      Output: vector of predictions.          */
/*             SE        Output: vector of standard errors.      */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for YHat and  */
/*             SE.                                               */
/*****************************************************************/

/*****************************************************************/
int KrigYHatSE(KrigingModel *KrigMod, real RAve, real *f, real *r,
          real *YHat, real *SE);
/*****************************************************************/
/*   Purpose:  Compute a kriging prediction and a standard error */
/*             at a single point.                                */
/*                                                               */
/*   Args:     KrigMod   Kriging model and decompositions.       */
/*             RAve      Average R (1.0 for prediction at a      */
/*                       point; used for s.e. of average).       */
/*             f         Input: vector of linear model terms for */
/*                       the new point.                          */
/*                       Output: used for workspace.             */
/*             r         Input: vector of correlations between   */
/*                       the new point and the design points.    */
/*                       Output: used for workspace.             */
/*             YHat      Output: the prediction.                 */
/*             SE        Output: the standard error (computed    */
/*                       only if SE != NULL).                    */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
int KrigTilde(const KrigingModel *KrigMod, real *f, real *r);
/*****************************************************************/
/*   Purpose:  Overwrite f with fTilde and r with rTilde.        */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/
