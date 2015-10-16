/*****************************************************************/
/*   GASP MAIN PROGRAM                                           */
/*                                                               */
/*   Copyright (c) William J. Welch 1995--1999.                  */
/*   All rights reserved.                                        */
/*                                                               */
/* 2009.05.08: Derivatives.Min and Derivatives.Max added to      */
/*             FitCheck                                          */
/* 2011.07.06: StochasticProcessVarianceProportion.Min and       */
/*             StochasticProcessVarianceProportion.Max added to  */
/*             FitCheck                                          */
/*****************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "implem.h"
#include "lib.h"
#include "model.h"
#include "alex.h"

#define BANNER      "GaSP (Gaussian Stochastic Process) "\
"Version 2.7\n" COPYRIGHT

#define PROMPT      "GaSP> "

string Prompt = PROMPT;  /* Communicates with dbmat.c. */

boolean DesignJob = NO;

/* Inputs: */

const string CVCheck[]   = {IN_DIR, OUT_DIR,
                              X_DESCRIP, X_MAT, Y_DESCRIP, Y_MAT,
                              REG_MOD, SP_MOD, COR_FAM, RAN_ERR,
                              CV_MAT, NULL};

const string FitCheck[]  = {IN_DIR, OUT_DIR,
                              X_DESCRIP, X_MAT, Y_DESCRIP, Y_MAT,
                              REG_MOD, SP_MOD, COR_FAM, RAN_ERR,
                              THETA "." STANDARDIZED "." MIN,
                              THETA "." STANDARDIZED "." MAX,
                              ALPHA "." MIN, ALPHA "." MAX,
                              "Derivatives" "." MIN, "Derivatives" "." MAX,
                              SP_VAR_PROP "." MIN, SP_VAR_PROP "." MAX,
                              CRIT_LOG_LIKE_DIFF, LOG_LIKE_TOL,
                              TRIES, RAN_NUM_SEED, NULL};

/*
const string SeqDesCheck[] = {IN_DIR, OUT_DIR,
                              X_DESCRIP, CAND, EXP_REG, X_MAT,
                              Y_DESCRIP, Y_MAT,
                              REG_MOD, SP_MOD, COR_FAM, RAN_ERR,
                              THETA "." STANDARDIZED "." MIN,
                              THETA "." STANDARDIZED "." MAX,
                              ALPHA "." MIN, ALPHA "." MAX,
                              CRIT_LOG_LIKE_DIFF, LOG_LIKE_TOL,
                              REFIT_RUNS,
                              RESP_FUNC, SEQ_CRIT,
                              TOL "." ABSOLUTE, TOL "." RELATIVE,
                              RAN_NUM_SEED, NULL};
*/

const string PredCheck[] = {IN_DIR, OUT_DIR,
                              X_DESCRIP, X_MAT, Y_DESCRIP, Y_MAT,
                              REG_MOD, SP_MOD, COR_FAM, RAN_ERR,
                              X_PRED, Y_PRED, Y_TRUE,
                              GEN_PRED_COEF, PRED_COEF, NULL};

const string VisCheck[]  = {IN_DIR, OUT_DIR,
                              X_DESCRIP, CAND, PRED_REG, X_MAT,
                              Y_DESCRIP, Y_MAT,
                              REG_MOD, SP_MOD, COR_FAM, RAN_ERR,
                              MAIN_EFF_PERC, INTER_EFF_PERC,
                              ANOVA_PERC, MAIN_EFF, JOINT_EFF, NULL};

/* Implemented functions: */
static Function ImpFn[] =
{
     {"CrossValidate",          CrossValidate,    CVCheck  },
     {"Fit",                    Fit,              FitCheck },
     /*
     {"SequentialDesign",       DataAdaptSeqDes,  SeqDesCheck },
     */
     {"Predict",                Predict,          PredCheck},
     {"Visualize",              Visualize,        VisCheck }
};

#define NUM_FUNCS   (sizeof(ImpFn) / sizeof(Function))

int main(int argc, char *argv[])
{
     Run(argc, argv, NUM_FUNCS, ImpFn, BANNER, PROMPT);

     exit(0);
}
