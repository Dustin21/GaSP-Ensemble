/*****************************************************************/
/*   Copyright (c) William J. Welch 1994--2009.                  */
/*   All rights reserved.                                        */
/*                                                               */
/*   Version: 2011.07.06                                         */
/*****************************************************************/

#define COPYRIGHT   "(c) Copyright William J. Welch 1996-2011.  "\
"All rights reserved.\n\n"

#ifndef KRIG_MOD_DEFINED
     #include "kriging.h"
#endif

/* Type for implemented functions: */
typedef struct
{
     string         FuncName;

     /* Pointer to function carrying out the operation. */
     int            (*Func)(void);

     /* Inputs required. */
     const string   *Check;
} Function;

/* Type for a matrix in the database. */
typedef struct DbMatStruct
{
     string    Name;
     Matrix    *M;
     int       Type;
     string    Title;
     string    FileName;

     /* Compulsory and optional columns. */
     size_t    nCompCols;
     string    *CompCol;
     size_t    nOptCols;
     string    *OptCol;

     /* For compatibility, the row names must be identical */
     /* to those in RowComp.                               */
     struct DbMatStruct  *RowLabelsComp1;
     struct DbMatStruct  *RowLabelsComp2;

     /* Is it an "output" matrix.  If so, the default will */
     /* set up row labels as in RowLabelsComp1/2.          */
     boolean   IsOutput;

     /* Does the matrix hold x variables?  If so, the column    */
     /* names must include all active variables in XDescrip.    */
     boolean   IsX;
     /* Does the matrix hold y variables?  If so, ... YDescrip. */
     boolean   IsY;

     boolean   IsFit;
     boolean   IsPred;

     /* Has a data matrix (IsX or IsY) been processed for */
     /* transformations, etc?                             */
     boolean   IsProcessed;
} DbMatrix;

/* Column names and types for database-status matrix. */
#define DB_STATUS_COL_NAMES   {"Object", "OK?", "Class", "Value/File", \
                                   "Rows", "Columns"}
#define DB_STATUS_COL_TYPES   {STRING, STRING, STRING, STRING, \
                                   SIZE_T, SIZE_T}
#define DB_STATUS_OBJ_COL     0
#define DB_STATUS_OK_COL      1
#define DB_STATUS_CLASS_COL   2
#define DB_STATUS_VAL_COL     3
#define DB_STATUS_ROWS_COL    4
#define DB_STATUS_COLS_COL    5

/* Criteria for comparing models. */
#define MOD_COMP_CRIT_NAMES   {LIKELIHOOD, CROSS_VALIDATION}
#define MOD_COMP_CRIT_LIKE    0
#define MOD_COMP_CRIT_CV      1     


/* run.c: */

/*****************************************************************/
void Run(int argc, char *argv[], size_t nFns, Function *ImpFn,
          const string Banner, const string Prompt);
/*****************************************************************/
/* Purpose:    Run an application.                               */
/*****************************************************************/

/*****************************************************************/
void EventLoop(string FileName, FILE *InpFile, size_t nFns,
          Function *ImpFn, const string Prompt);
/*****************************************************************/
/*   Purpose:  Execute events.                                   */
/*****************************************************************/

/*****************************************************************/
int ExecuteFunc(int (*Func)(void));
/*****************************************************************/
/*   Purpose:  Execute a function (verb).                        */
/*****************************************************************/

/* gaspcv.c: */

/*****************************************************************/
int CrossValidate(void);
/*****************************************************************/
/*   Purpose:  Compute cross validations, CV.                    */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
int CalcCV(KrigingModel *KrigMod, real *YHatCV, real *SE);
/*****************************************************************/
/*   Purpose:  Compute kriging cross-validation predictions and, */
/*             optionally, their standard errors.                */
/*                                                               */
/*   Args:     KrigMod   Input: Kriging model without            */
/*                       decompositions.                         */
/*                       Output: Decompositions are garbage.     */
/*             YHatCV    Output: Cross-validation predictions.   */
/*             SE        Output: Standard errors (computed only  */
/*                       if SE != NULL).                         */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for YHatCV    */
/*             and SE.                                           */
/*             Better matrix updating for doing this?            */
/*             KrigMod decompositions are changed.               */
/*****************************************************************/


/* gaspseq.c: */

/*****************************************************************/
int DataAdaptSeqDes(void);
/*****************************************************************/
/*   Purpose:  Data-adaptive sequential design.                  */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
int GetUserObj(size_t nXVars, real *x, real *y);
/*****************************************************************/
/*   Purpose:  Get y(x) by calling the user's objective          */
/*             function (which may also change x).               */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
boolean MaximizeExpectedImprovement(size_t nXVars, real *xBest);
/*****************************************************************/
/*   Purpose:  Find xBest that maximizes the negative expected   */
/*             improvement in the minimum.                       */
/*                                                               */
/*   Returns:  Yes or NO, to indicate whether the expected       */
/*             improvement has converged to zero.                */
/*****************************************************************/

/*****************************************************************/
real ExpectedImprovement(real *x, size_t nXVars);
/*****************************************************************/
/*   Purpose:  Return the expected improvement in the minimum    */
/*             if the function is evaluated at x.                */
/*****************************************************************/

/*****************************************************************/
boolean MaximizeUncertainty(size_t nXVars, real *xBest);
/*****************************************************************/
/*   Purpose:  Find xBest, point with maximum uncertainty.       */
/*                                                               */
/*   Returns:  Yes or NO, to indicate whether the maximum        */
/*             uncertainty has converged to zero.                */
/*****************************************************************/

/*****************************************************************/
real Uncertainty(real *x, size_t nXVars);
/*****************************************************************/
/*   Purpose:  Return p(1 - p), where p is Pr(y > yCritical), a  */
/*             measure of uncertainty whether y < or > yCritical.*/
/*****************************************************************/

/* gaspvis.c: */

/*****************************************************************/
int Visualize(void);
/*****************************************************************/
/*   Purpose:  Compute effects.                                  */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
void RemEffectsRows(const string yName, Matrix *Eff);
/*****************************************************************/
/*   Purpose:  Remove any rows in Eff corresponding to           */
/*             y-variable yName.                                 */
/*****************************************************************/

/*****************************************************************/
int CompEffects(KrigingModel *KrigMod, const string yName,
     const Matrix *PredReg, const size_t *GroupSize,
     const Matrix *GroupVarIndex, real MainPerc, real InterPerc,
     real *Perc, real *Average, real *SEAve);
/*****************************************************************/
/*   Purpose:  Compute ANOVA percentage contributions and put    */
/*             important effects in MAIN_EFF and JOINT_EFF       */
/*             matrices.                                         */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for Perc,     */
/*             Average, and SEAve.                               */
/*****************************************************************/

/*****************************************************************/
void AvePred(KrigingModel *KrigMod, const Matrix *PredReg,
     size_t nGroups, const size_t *IndexEffectGroup,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     real *fAve, real *rAve, real *RAve);
/*****************************************************************/
/*   Purpose:  Average f and r w.r.t. the groups *not* in        */
/*             IndexEffectGroup.                                       */
/*****************************************************************/

/*****************************************************************/
int CompSSTot(KrigingModel *KrigMod, const Matrix *PredReg,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     real *SSTot);
/*****************************************************************/
/*   Purpose:  Compute SS(Total) for predictor.                  */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*****************************************************************/

/*****************************************************************/
void AnyEffect(KrigingModel *KrigMod, const Matrix *PredReg,
     size_t nGroups, const size_t *IndexEffectGroup,
     const size_t *GroupSize, const Matrix *GroupVarIndex,
     const size_t *nSPTerms, const Matrix *IndexSP,
     real *Eff, real *SE);
/*****************************************************************/
/*   Purpose:  Compute any (arbitrary-degree) effect.            */
/*                                                               */
/*   Comments: Calling routine must allocate space for Eff and   */
/*             SE.                                               */
/*****************************************************************/

/*****************************************************************/
void AppendEffect(const string yName, size_t DegreeEff,
          const size_t *IndexGroup, const Matrix *PredReg,
          const size_t *GroupSize, const Matrix *GroupVarIndex,
          real *Eff, real *SE, Matrix *EffMat);
/*****************************************************************/
/*   Purpose:  Append an effect to the bottom of EffMat.         */
/*****************************************************************/


/* acedeval.c: */

/*****************************************************************/
int EvalDes(void);
/*****************************************************************/
/*   Purpose:  Evaluate a design, X.                             */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/


/* gaspfit.c: */

/*****************************************************************/
int Fit(void);
/*****************************************************************/
/*   Purpose:  Fit the model parameters.                         */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
int FitBest(KrigingModel *KrigMod, size_t Tries, real *Beta,
     Matrix *CorPar, real *SPVar, real *ErrVar, real *NegLogLike,
     real *CVRootMSE, unsigned *nEvals, real *CondNum);
/*****************************************************************/
/*   Purpose:  Choose best of several MLE tries.                 */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/


/* acedoptd.c: */

int OptDes(void);


/* gasppred.c: */

/*****************************************************************/
int Predict(void);
/*****************************************************************/
/*   Purpose:  Compute predictions, their standard errors, and   */
/*             prediction coefficients.                          */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
void OutputSummary(Matrix *Summ, size_t nCols,
     const string *ColName);
/*****************************************************************/
/*   Purpose:  Output summary information.                       */
/*****************************************************************/


/* acedlhs.c: */

/*****************************************************************/
int RunLHS(void);
/*****************************************************************/
/*   Purpose:  Compute a Latin-hypercube design with specified   */
/*             correlations (default identity).                  */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
int ExpandCor(const Matrix *X, const Matrix *ExpReg,
          const Matrix *NonzeroCor, Matrix *TargetCor);
/*****************************************************************/
/*   Purpose:  Expand the correlation matrix NonzeroCor of       */
/*             nonzero correlations to a full correlation matrix */
/*             for all variables, TargetCor.                     */
/*                                                               */
/*   Returns:  INPUT_ERR if one of the row/column names in       */
/*                       NonzeroCor is not an X column;          */
/*             OK        otherwise.                              */
/*****************************************************************/


/* db.c: */

/*****************************************************************/
int DbCheck(const string Name);
/*****************************************************************/
/*   Purpose:  Return status of database object and add a row to */
/*             DbStatus.                                         */
/*                                                               */
/*   Returns:  INPUT_ERR, INCOMPAT_ERR, or OK.                   */
/*****************************************************************/

/*****************************************************************/
void DbOutputMatStatus(void);
/*****************************************************************/
/*   Purpose:  Output status of output matrices.                 */
/*****************************************************************/


/* dbmanip.c: */

/*****************************************************************/
int DbProcessDescrip(const string DescripName,
     const Matrix *Descrip, DbMatrix *D);
/*****************************************************************/
/*   Purpose:  Process a data matrix by applying Descrip, i.e.,  */
/*             including only variables that are to be analyzed, */
/*             replacing cases out of range by NA's, and         */
/*             applying transformations.                         */
/*****************************************************************/

/*****************************************************************/
size_t DbnActiveY(void);
/*****************************************************************/
/*   Purpose:  Return the number of y variables to be analyzed.  */
/*****************************************************************/

/*****************************************************************/
size_t DbIndexXY(size_t j);
/*****************************************************************/
/*   Purpose:  Put all cases with no NA's in X or in y variable  */
/*             j in IndexXY.                                     */
/*             On exit, y points to the Y column for variable j, */
/*             yTrue points to the YTtrue column, and yName is   */
/*             the name.                                         */
/*                                                               */
/*   Returns:  The number of cases with no NA's.                 */
/*****************************************************************/

/*****************************************************************/
size_t DbIndexX(void);
/*****************************************************************/
/*   Purpose:  Put all cases with no NA's in X in IndexX.        */
/*                                                               */
/*   Returns:  The number of cases with no NA's.                 */
/*****************************************************************/


/* dbmat.c: */

/*****************************************************************/
void DbMatInit(void);
/*****************************************************************/
/*   Purpose:  Initialize the matrix database.                   */
/*****************************************************************/

/*****************************************************************/
DbMatrix *DbMatFind(const string Name, boolean HardFail);
/*****************************************************************/
/*   Purpose:  Return the DbMat entry for Name.                  */
/*                                                               */
/*   Returns:  Pointer to DbMat[?] or NULL.                      */
/*****************************************************************/

/*****************************************************************/
int DbMatCheck(DbMatrix *D, string *FileName);
/*****************************************************************/
/*   Purpose:  Check whether matrix D in DbMat is legal and      */
/*             compatible.                                       */
/*             On return, *FileName contains the input file name.*/
/*                                                               */
/*   Returns:  INPUT_ERR, INCOMPAT_ERR, or OK.                   */
/*****************************************************************/

/*****************************************************************/
boolean DbMatDefault(DbMatrix *D);
/*****************************************************************/
/*   Purpose:  Assign a default to a matrix.                     */
/*                                                               */
/*   Returns:  YES if a default is assigned;                     */
/*             NO  otherwise.                                    */
/*****************************************************************/

/*****************************************************************/
void DbMatDefDescrip(size_t nXVars, const string *xName,
     Matrix *Descrip);
/*****************************************************************/
/*   Purpose:  Assign a default XDescrip or YDescrip.            */
/*****************************************************************/

/*****************************************************************/
int DbMatRead(DbMatrix *D, const string FileName);
/*****************************************************************/
/*   Purpose:  Read a matrix from a file and, if legal, put it   */
/*             in the database.                                  */
/*                                                               */
/*   Returns:  INPUT_ERR    if the matrix is illegal;            */
/*             INCOMPAT_ERR if the matrix is incompatible        */
/*                          (with XDescrip or YDescrip);         */
/*             OK           otherwise.                           */
/*****************************************************************/

/*****************************************************************/
void DbNewDescrip(const string DescripName);
/*****************************************************************/
/*   Purpose:  Free data matrices for a new Descrip matrix.      */
/*****************************************************************/

/*****************************************************************/
int DbMatWrite(DbMatrix *D, const string FileName,
          const string BlockingOption);
/*****************************************************************/
/*   Purpose:  Output a matrix to a file.                        */
/*             BlockingOption should be BLOCKED or UNBLOCKED;    */
/*             a NULL value is the same as BLOCKED.              */
/*                                                               */
/*   Returns:  OK        if successful;                          */
/*             INPUT_ERR otherwise.                              */
/*****************************************************************/


/* dbmatcom.c: */

/*****************************************************************/
int DbMatCompat(DbMatrix *D);
/*****************************************************************/
/*   Purpose:  Check compatibility of a matrix with other        */
/*             information in the database.                      */
/*                                                               */
/*   Returns:  INCOMPAT_ERR or OK.                               */
/*****************************************************************/

/*****************************************************************/
int DbRowLabels(const DbMatrix *D);
/*****************************************************************/
/*   Purpose:  Check consistency of row labels with those of up  */
/*             to two other database matrices.                   */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if an inconsistency is found;        */
/*             OK           otherwise.                           */
/*****************************************************************/

/*****************************************************************/
int DbRowLabelsComp(const DbMatrix *D1, const DbMatrix *D2);
/*****************************************************************/
/*   Purpose:  Compare the row labels of D1 and D2.              */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if an inconsistency is found;        */
/*             OK           otherwise.                           */
/*****************************************************************/

/*****************************************************************/
int DbMatX(const Matrix *X, Matrix *XReg);
/*****************************************************************/
/*   Purpose:  Check that X contains any protected runs and any  */
/*             fixed variables in XReg.                          */
/*                                                               */
/*   Returns:  OK or INCOMPAT_ERR.                               */
/*****************************************************************/

/*****************************************************************/
int DbMatCatCols(size_t nXVars, const string *xName,
          const size_t *nCats, const Matrix *M);
/*****************************************************************/
/*   Purpose:  Check that any categorical columns have legal     */
/*             values.                                           */
/*                                                               */
/*   Returns:  INCOMPAT_ERR or OK.                               */
/*****************************************************************/


/* dbmatleg.c: */

/*****************************************************************/
void DbColTemplInit(void);
/*****************************************************************/
/*   Purpose:  Initialize column templates.                      */
/*****************************************************************/

/*****************************************************************/
int DbMatLegal(DbMatrix *D);
/*****************************************************************/
/*   Purpose:  Check that D is legal.                            */
/*                                                               */
/*   Returns:  INPUT_ERR    if the matrix is illegal;            */
/*             OK           otherwise.                           */
/*                                                               */
/*   Comment:  May change the matrix, e.g., by setting up a      */
/*             default.                                          */
/*****************************************************************/

/*****************************************************************/
int DbMatColChk(DbMatrix *D);
/*****************************************************************/
/*   Purpose:  Check each column independently.                  */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*****************************************************************/

/*****************************************************************/
int DbMatCor(const Matrix *M);
/*****************************************************************/
/*   Purpose:  Check legality of a correlation matrix.           */
/*                                                               */
/*   Returns:  INPUT_ERR    if illegal;                          */
/*             OK           otherwise.                           */
/*****************************************************************/


/* dbscalar.c: */

/*****************************************************************/
void DbScalarInit(void);
/*****************************************************************/
/*   Purpose:  Initialize the scalar database.                   */
/*****************************************************************/

/*****************************************************************/
size_t DbScalIndex(const string Name, boolean ErrorMessage);
/*****************************************************************/
/*   Purpose:  Return the Scal index for Name.                   */
/*                                                               */
/*   Returns:  INDEX_ERR if not found;                           */
/*             OK        otherwise.                              */
/*****************************************************************/

/*****************************************************************/
int DbScalCheck(size_t Index);
/*****************************************************************/
/*   Purpose:  Check whether scalar Index in Scal is legal and   */
/*             compatible.                                       */
/*                                                               */
/*   Returns:  INPUT_ERR, INCOMPAT_ERR, or OK.                   */
/*****************************************************************/

/*****************************************************************/
string DbScalValue(size_t Index);
/*****************************************************************/
/*   Purpose:  Return the value of the scalar Index in Scal as a */
/*             string.                                           */
/*****************************************************************/

/*****************************************************************/
int DbScalParse(size_t Index, const string Token);
/*****************************************************************/
/*   Purpose:  Parse Token and, if legal, put it in the scalar   */
/*             database.                                         */
/*                                                               */
/*   Returns:  INPUT_ERR    if Token is illegal;                 */
/*             OK           otherwise.                           */
/*****************************************************************/

/*****************************************************************/
int DbScalCompat(size_t Index);
/*****************************************************************/
/*   Purpose:  Check Scal[Index] for compatibility with other    */
/*             database objects.                                 */
/*                                                               */
/*   Returns:  INCOMPAT_ERR or OK.                               */
/*****************************************************************/


/* lhs.c: */

/*****************************************************************/
void LHSRand(size_t nNew, size_t nProtected, const Matrix *XReg,
          Matrix *X);
/*****************************************************************/
/*   Purpose:  Add a completely random Latin hypercube of nNew   */
/*             points to the nProtected points already in X.     */
/*****************************************************************/

/*****************************************************************/
void LHSMargin(const size_t *Rank, size_t n, const Matrix *XReg,
          size_t j, real *Margin);
/*****************************************************************/
/*   Purpose:  Generate a margin (i.e., column) for a Latin-     */
/*             hypercube design.                                 */
/*****************************************************************/

/*****************************************************************/
int LHS(const Matrix *XReg, const Matrix *TargetCor, Matrix *X,
          real *Objective, real *CondNum);
/*****************************************************************/
/*   Purpose:  Generate a Latin hypercube X.                     */
/*             If there are no categorical or grouped-candidate  */
/*             variables, the correlations will be               */
/*             (approximately) matched to TargetCor; otherwise,  */
/*             a completely random Latin hypercube will be       */
/*             generated.                                        */
/*                                                               */
/*   Comment:  Assumes that any fixed variables are first in     */
/*             XReg, X, etc.                                     */
/*****************************************************************/

/*****************************************************************/
int LHSTargetCor(const Matrix *XReg, const Matrix *TargetCor,
          Matrix *X, real *Objective, real *CondNum);
/*****************************************************************/
/*   Purpose:  Generate a Latin hypercube X with correlations    */
/*             that (approximately) match TargetCor.             */
/*                                                               */
/*   Comment:  Assumes that any fixed variables are first in     */
/*             XReg, X, etc.                                     */
/*****************************************************************/

/*****************************************************************/
void LHSCorStats(const Matrix *R, const Matrix *XReg,
          const Matrix *TargetCor, real *MaxCor, real *RMSCor,
          size_t *iMax, size_t *jMax);
/*****************************************************************/
/*   Purpose:  Extract maximum and RMS deviation from the target */
/*             correlations.                                     */
/*             Correlations between fixed variables are ignored. */
/*****************************************************************/

/*****************************************************************/
int LHSImprove(const Matrix *R, const Matrix *TargetR,
          const Matrix *XReg, Matrix *X);
/*****************************************************************/
/*   Purpose:  Transform X to target correlation structure.      */
/*                                                               */
/*   Returns:  NUMERIC_ERR if R equations cannot be solved or do */
/*                         not have a unique solution;           */
/*             OK          otherwise.                            */
/*                                                               */
/*   Comment:  X columns must have means of zero.                */
/*****************************************************************/

/*****************************************************************/
void LHSValidX(const Matrix *XReg, const real *Mean, Matrix *X);
/*****************************************************************/
/*   Purpose:  Transform columns of X to valid LHS margins.      */
/*****************************************************************/



/* The following names, text, etc. are gathered here to */
/* facilitate consistency and possible translation into */
/* foreign languages (see also define.h).               */

/* Names of int scalars: */

#define RAN_NUM_SEED          "RandomNumberSeed"

/* Names of real scalars: */

#define CENSORING_LIMIT       "CensoringLimit"
#define COVER_DIST            "CoverageDistance"
#define DIST_METRIC           "DistanceMetric"
#define CRIT_COR              "CriticalCorrelation"
#define CRIT_LOG_LIKE_DIFF    "CriticalLogLikelihoodDifference"
#define INTER_EFF_PERC        "InteractionEffectPercentage"
#define LAMBDA                "Lambda"
#define LIFE_DIST             "LifetimeDistribution"
#define LOG_LIKE_TOL          "LogLikelihoodTolerance"
#define MAIN_EFF_PERC         "MainEffectPercentage"
#define METRIC                "Metric"
#define TOL                   "Tolerance"

/* Names of size_t scalars: */

#define PROJ_DIM         "ProjectionDimension"
#define PROTECTED_RUNS   "ProtectedRuns"
#define REFIT_RUNS       "RefitRuns"
#define RUNS             "Runs"
#define TRIES            "Tries"
#define N_X_VARS         "xVariables"

/* Names of string scalars: */

#define COR_FAM               "CorrelationFamily"
#define DESIGN_ALG            "DesignAlgorithm"
#define DESIGN_CRIT           "DesignCriterion"
#define GEN_PRED_COEF         "GeneratePredictionCoefficients"
#define IN_DIR                "InputDirectory"
#define MOD_COMP_CRIT         "ModelComparisonCriterion"
#define NORMALIZED_RANGES     "NormalizedRanges"
#define RAN_ERR               "RandomError"
#define SEQ_CRIT              "SequentialCriterion"
#define RESP_FUNC             "ResponseFunction"
#define OUT_DIR               "OutputDirectory"

/* Values taken by string scalars: */

#define LIKELIHOOD       "Likelihood"
#define CROSS_VALIDATION "CrossValidation"
#define MATERN           "Matern"
#define POW_EXP          "PowerExponential"

/* Names of matrices: */

#define ANOVA_PERC       "ANOVAPercentages"
#define CAND             "Candidates"
#define CV_MAT           "CrossValidations"
#define EXP_REG          "ExperimentalRegion"
#define JOINT_EFF        "JointEffects"
#define K_PHI            "KPhi"    /* Kiefer's Phi. */
#define MAIN_EFF         "MainEffects"
#define PRED_COEF        "PredictionCoefficients"
#define PRED_REG         "PredictionRegion"
#define PRIOR_SAMP       "PriorSample"
#define REG_MOD          "RegressionModel"
#define SP_MOD           "StochasticProcessModel"
#define T_MAT            "T"
#define X_AVERAGE        "XAverage"
#define X_COR            "XCorrelations"
#define X_DESCRIP        "XDescription"
#define X_MAT            "X"
#define X_PRED           "XPrediction"
#define Y_DESCRIP        "YDescription"
#define Y_MAT            "Y"
#define Y_PRED           "YPrediction"
#define Y_TRUE           "YTrue"

/* Matrix titles, output before column names: */

#define ANOVA_PERC_TITLE "ANOVA percentage contributions."
#define CV_TITLE         "CV predictions and standard errors."
#define LHS_TITLE        "Latin-hypercube experimental design."
#define MAIN_EFF_TITLE   "Important main effects."
#define JOINT_EFF_TITLE  "Important joint effects."
#define PRED_COEF_TITLE  "Coefficients for prediction."
#define REG_MOD_TITLE    "Estimated regression parameters."
#define SP_MOD_TITLE     "Estimated correlation parameters."
#define X_TITLE          "Experimental design."
#define X_DESCRIP_TITLE  "X variable descriptions."
#define Y_DESCRIP_TITLE  "Y variable descriptions and summary statistics."
#define Y_PRED_TITLE     "Predictions and standard errors."

/* Column names or extensions (e.g., Y.Pred): */

#define ALPHA            "Alpha"
#define ANOVA_TOTAL_PERC "ANOVATotal%"
#define BETA             "Beta"
#define CASE_MAX_ERR     "CaseMaxErr"
#define CASE_CV_MAX_ERR  "CaseCVMaxErr"
#define CASES            "Cases"
#define COEF             "Coef"
#define COND_NUM         "ConditionNumber"
#define CV_MAX_ERR       "CVMaxErr"
#define CV_ROOT_MSE      "CVRootMSE"
#define ERR_VAR          "ErrorVariance"
#define FIT              "Fit"
#define LOG_LIKE         "LogLikelihood"
/* MIN and MAX in define.h. */
#define MAX_ERR          "MaxErr"
#define PRED             "Pred"
#define ROOT_MSE         "RootMSE"
#define STANDARDIZED     "Standardized"
#define STD_ERR          "SE"
#define SP_VAR           "StochasticProcessVariance"
#define SP_VAR_PROP      "StochasticProcessVarianceProportion"
#define THETA            "Theta"


/* Information messages: */

/* #define DB_MAT_EMPTY     "%s matrix is empty.\n" */
#define DB_ROWS_COLS     "%d row%s and %d column%s read.\n"

#define LHS_CAT          "No correlation matching because of " \
                              "categorical variables.\n"
#define LHS_GROUP        "No correlation matching because of " \
                              "grouped-candidate variables.\n"
#define LHS_FIXED        "Put fixed variables first in " X_DESCRIP \
                              " for correlation matching.\n"
#define LHS_NA           "No correlation matching because of " \
                              "NA's in fixed variables.\n"
#define LHS_POS_DEF      "Target correlation matrix is not " \
                              "positive definite.\n"

/* Information put in database-status matrix: */
#define EMPTY_STR             "Empty"
#define DEFAULT_STR           "Default"
#define DEFAULT_X_AVERAGE     "Same as " X_MAT
#define OUTPUT_STR            "Output"

/* Old messages, no longer used: */
/*
#define NO_SP_MOD        "Stochastic-process model has a (linear) " \
                              "term for each x variable.\n"
#define NO_X_COR         "All x-variable correlations are zero.\n"
#define NO_X_PRED        "There are no new x values for prediction.\n"
#define NO_X_DESCRIP     "All columns in " X_MAT " are x variables.\n"
#define NO_Y_DESCRIP     "All columns in " Y_MAT " are y variables.\n"
*/

/* Error messages: */

#define DB_ANOVA         "The row names of " ANOVA_PERC  \
                               " must agree with " PRED_REG ".\n"
#define DB_BLOCK         "The last token should be \"Blocked\" " \
                              "or \"Unblocked\".\n"
#define DB_CAT           "%s is a categorical variable: It must "\
                              "take integer values between 1 and %d.\n"
#define DB_COL           "%s is not a valid column name.\n"
#define DB_COL_VAR       "%s must have a column %s\n" \
                              "(a Variable in %s).\n"
#define DB_COMPULSORY    "%s must have a column %s.\n"
#define DB_COR_DIAG      "Diagonal element %d must be 1.\n"
#define DB_COR_LABEL     "Row and column labels must be " \
                              "identical: reconcile %s and %s.\n"
#define DB_COR_ONE       "Element (%d, %d) must be between " \
                               "-1 and 1.\n"
#define DB_COR_SYM       "Elements (%d, %d) and (%d, %d) must " \
                              "be equal.\n"
#define DB_COR_VAR       "%s in " X_COR " must appear as a Variable " \
                              "in the\n" X_DESCRIP " matrix.\n"
#define DB_DUP_COL       "Column name %s should not be repeated.\n"
#define DB_DUP_TERM      TERM " %s should not be repeated.\n"
#define DB_DUP_VAR       VARIABLE " %s should not be repeated " \
                               "(with the same " TRANSFORMATION ").\n"
#define DB_MATRIX        "%s is not a database matrix.\n"
#define DB_MAT_ERR_AT    "At column %s, row (label) %s.\n"
#define DB_N_ROWS        "Matrices %s and %s must have the "\
                              "same number of rows.\n"
#define DB_POS_DEF       "Must be positive definite.\n"
#define DB_MIN_MAX       "%s must be less than or equal to %s.\n"
#define DB_PROJ_XVARS    "%s must be less than or equal to\n" \
                              "the number of x variables.\n"
#define DB_PROT_RUNS     PROTECTED_RUNS " must be less than " RUNS ".\n"
#define DB_SCALAR        "%s is not a database scalar.\n"
#define DB_SQ            "Must be square.\n"
#define DB_ROW_NAME      "Matrices %s and %s must have the same row " \
                              "labels: Reconcile %s and %s.\n"
#define DB_TEMPLATE      "Code bug: No template for %s.\n"

#define DB_X_COL_PROT    X_MAT " must have a column %s containing " \
                              "the protected runs.\n"
#define DB_X_COL_FIXED   X_MAT " must have a column " \
                              "for Fixed x variable %s.\n"
#define DB_X_N           X_MAT " must have %d rows (= " RUNS ") " \
                              "for the fixed x variables.\n"
#define DB_X_PROT        X_MAT " must have %d rows for the " \
                              "protected runs.\n"

/* Miscellaneous: */

#define BLOCKED     "Blocked"
#define SCREEN      "Screen"
#define UNBLOCKED   "Unblocked"
