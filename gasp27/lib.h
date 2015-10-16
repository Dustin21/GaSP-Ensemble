/*****************************************************************/
/*   Copyright (c) William J. Welch 1991--99.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*   Version: 1999.06.24                                         */
/*****************************************************************/

#include <stdarg.h>

#ifndef MAT_DEFINED
#include "matrix.h"
#endif

/* The following names, text, etc. are gathered here to */
/* facilitate consistency and possible translation into */
/* foreign languages.                                   */

/* libin.c error messages: */

#define CANNOT_READ_SIZE_T    "Cannot read integer number >= 0.\n"
#define CANNOT_READ_REAL      "Cannot read (real) number.\n"


/* libreg.c error messages: */

#define REG_CAND         DISCRETE_STR  " variable %s must appear " \
                             "in the Candidates matrix.\n"
#define REG_CAND_GROUP   "Variable %s has CandidateGroup = %d: " \
                         "Support must be Discrete.\n"
#define REG_CAT          "%s variable %s is categorical: "
#define REG_CAT_CONT     REG_CAT "It cannot also be continuous.\n"
#define REG_CAT_MAX      REG_CAT "Max must be an integer " \
                              "< or = the number of categorical " \
                              "levels.\n"
#define REG_CAT_MIN      REG_CAT "Min must be an integer > 0.\n"
#define REG_CAT_STEP     REG_CAT "Min, Max, and " NUM_LEVELS "\n" \
                              "must give a step length of 1 for " \
                              "variable %s.\n"
#define REG_MINMAX       "%s variable %s must have\n" \
                              "Min (default 0) < or = Max (default 1).\n"
#define REG_INC_NORM     "%s variable %s cannot have inclusive ranges\n" \
                              "because it has a normal distribution.\n"
#define REG_NUM_LEVELS   "%s variable %s must have " NUM_LEVELS " > 0.\n"
#define REG_NUM_LEVELS_MINMAX \
                         "%s variable %s (Grid with NumberLevels = 1)\n" \
                               "must have Min (default 0.5) = " \
                               "Max (default 0.5).\n"
#define REG_SUPPORT      "%s must have a column labelled " \
                              SUPPORT ".\n"

/* Type for a scalar. */
typedef struct
{
     string    Name;
     int       Type;     /* INTEGER, REAL, SIZE_T, or STRING. */
     int       IntVal;
     real      RealVal;
     size_t    Size_tVal;
     string    Str;
     boolean   Legal;
} scalar;

#define ScalName(s)                ((s)->Name)
#define ScalPutName(s, NewName)    ((s)->Name = StrDup(NewName))
#define ScalType(s)                ((s)->Type)
#define ScalPutType(s, NewType)    ((s)->Type = NewType)
#define ScalInt(s)                 ((s)->IntVal)
#define ScalPutInt(s, i)           ((s)->IntVal = i)
#define ScalReal(s)                ((s)->RealVal)
#define ScalPutReal(s, r)          ((s)->RealVal = r)
#define ScalSize_t(s)              ((s)->Size_tVal)
#define ScalPutSize_t(s, z)        ((s)->Size_tVal = z)
#define ScalStr(s)                 ((s)->Str)
#define ScalPutStr(s, NewStr)      ((s)->Str = StrDup(NewStr))
#define ScalLegal(s)               ((s)->Legal)
#define ScalPutLegal(s, NewLegal)  ((s)->Legal = NewLegal)


/* Type for a vector. */
typedef struct
{
     string    Name;
     int       Type;          /* INTEGER, REAL, SIZE_T, or STRING. */
     size_t    Length;        /* Actual length. */
     size_t    AllocLength;   /* Allocated (maximum) length. */
     int       *IntVal;
     real      *RealVal;
     size_t    *Size_tVal;
     string    *Str;
     boolean   Legal;
} vector;

#define VecName(v)                 ((v)->Name)
#define VecPutName(v, NewName)     ((v)->Name = StrDup(NewName))
#define VecType(v)                 ((v)->Type)
#define VecPutType(v, NewType)     ((v)->Type = NewType)
#define VecLength(v)               ((v)->Length)
#define VecAllocLength(v)          ((v)->AllocLength)
#define VecPutLength(v, NewLength) ((v)->Length = NewLength)
#define VecInts(v)                 ((v)->IntVal)
#define VecPutInts(v, i)           ((v)->IntVal = i)
#define VecInt(v, j)               ((v)->IntVal[j])
#define VecPutInt(v, j, i)         ((v)->IntVal[j] = i)
#define VecReals(v)                ((v)->RealVal)
#define VecPutReals(v, r)          ((v)->RealVal = r)
#define VecReal(v, j)              ((v)->RealVal[j])
#define VecPutReal(v, j, r)        ((v)->RealVal[j] = r)
#define VecSize_ts(v)              ((v)->Size_tVal)
#define VecPutSize_ts(v, z)        ((v)->Size_tVal = z)
#define VecSize_t(v, j)            ((v)->Size_tVal[j])
#define VecPutSize_t(v, j, z)      ((v)->Size_tVal[j] = z)
#define VecStrs(v)                 ((v)->Str)
#define VecPutStrs(v, NewStrs)     ((v)->Str = NewStrs)
#define VecStr(v, j)               ((v)->Str[j])
#define VecPutStr(v, j, NewStr)    ((v)->Str[j] = StrDup(NewStr))
#define VecLegal(v)                ((v)->Legal)
#define VecPutLegal(v, NewLegal)   ((v)->Legal = NewLegal)


/* Type for a generic singly-linked list. */
typedef struct ListStruct
{
     string              Name;
     int                 Type;
     void                *Data;
     struct ListStruct   *Next;
} List;


/* Type for a data template. */
typedef struct
{
     string    Name;
     int       Type;
     int       MinInt;
     int       MaxInt;
     real      MinReal;
     real      MaxReal;
     size_t    MinSize_t;
     size_t    MaxSize_t;
     string    *LegalStr;
     size_t    NumLegalStr;
} template;

#define TemplType(T)          ((T)->Type)
#define TemplStrIndex(s, T)   (StrIndex(s, T->LegalStr, T->NumLegalStr))

/* liballoc.c: */

char      *AllocChar(size_t n, char *p);
string    *AllocStr(size_t n, string *p);
string    **AllocPtrStr(size_t n, string **p);
real      *AllocReal(size_t n, real *p);
real      **AllocPtrReal(size_t n, real **p);
int       *AllocInt(size_t n, int *p);
int       **AllocPtrInt(size_t n, int **p);
size_t    *AllocSize_t(size_t n, size_t *p);
size_t    **AllocPtrSize_t(size_t n, size_t **p);

/*****************************************************************/
void *AllocGeneric(size_t n, size_t Size, void *p);
/*****************************************************************/
/*   Purpose:  Allocate or reallocate items of arbitrary size.   */
/*                                                               */
/*   Args:     n         Input     Number of items.              */
/*             Size      Input     Size (from sizeof) of an      */
/*                                 item.                         */
/*             p         Input     Pointer to items (NULL for    */
/*                                 first allocation).            */
/*                       Output    Points to the n items.        */
/*                                                               */
/*   Return:   Pointer to the n items.                           */
/*                                                               */
/*   Comment:  exit(1) is called if there is insufficient memory.*/
/*             This is the only difference from calloc() and     */
/*             realloc().                                        */
/*****************************************************************/

/*****************************************************************/
size_t AllocFindPtr(void *p);
/*****************************************************************/
/*   Purpose:  Return i such that Pointer[i] = p.                */
/*****************************************************************/

/*****************************************************************/
void AllocFree(void *p);
/*****************************************************************/
/*   Purpose:  Free p.                                           */
/*****************************************************************/

string    *AllocStrFree(size_t OldLen, size_t NewLen, string *s);
size_t    AllocMax(size_t Size);


/* libbufin.c: */

#define BufRead(InpFile) \
          BufReadPrompt(InpFile, NULL, NO)
#define BufTok(BufPtr)   BufToken(NULL, BufPtr, NULL)
#define BufForceTok(InpFile,BufPtr) \
          BufForceToken(InpFile, NULL, NO, NULL, BufPtr, NULL)

string    BufReadPrompt(FILE *InpFile, const string Prompt,
               boolean Echo);
string    BufToken(const string Delim, string *BufPtr,
               char *TermChar);
string    BufPromptToken(const string AppPrompt,
               const string TokPrompt, const string Assign);
string    BufForceToken(FILE *InpFile, const string Prompt,
               boolean Echo, const string Delim, string *BufPtr,
               char *TermChar);


/* libfile.c: */

FILE *FileOpen(const string FileName, const string Mode);
void FileBackup(const string FileName, FILE *Source);


/* libin.c: */

/*****************************************************************/
int InputSize_t(FILE *FilePtr, size_t *z);
/*****************************************************************/
/*   Purpose:  Read size_t z from file FilePtr.                  */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/

/*****************************************************************/
int InputReal(FILE *FilePtr, real *r);
/*****************************************************************/
/*   Purpose:  Read real r from file FilePtr.                    */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*****************************************************************/


/* liblist.c: */

/*****************************************************************/
List *ListAdd(const string Name, int Type, void *Data, List *L);
/*****************************************************************/
/*   Purpose:  Add a new member to the beginning of a list.      */
/*                                                               */
/*   Returns:  Pointer to the new list head.                     */
/*                                                               */
/*   Comment:  The calling routine must allocate void *Data for  */
/*             the new member.                                   */
/*****************************************************************/

/*****************************************************************/
void *ListData(const string Name, int Type, const List *L);
/*****************************************************************/
/*   Purpose:  Return a pointer to the data for a named member   */
/*             of a list.                                        */
/*                                                               */
/*   Returns:  Pointer if Name is found and Type is correct;     */
/*             NULL    otherwise.                                */
/*****************************************************************/


/* libmath.c: */

/*****************************************************************/
boolean ApproxEq(real a, real b, real AbsTol, real RelTol);
/*****************************************************************/
/*   Purpose:  Are a and b approximately equal?                  */
/*                                                               */
/*   Returns:  0 or 1.                                           */
/*****************************************************************/

/*****************************************************************/
real SafeExp(real x);
/*****************************************************************/
/*   Purpose:  Safe version of exp(x).                           */
/*                                                               */
/*   Returns:  exp(x) or NA_REAL.                                */
/*                                                               */
/*   Comment:  Change like SafeLog10.                            */
/*****************************************************************/

/*****************************************************************/
void SafeLog10(size_t n, const real *x, real *log10x,
     size_t *nInputNA, size_t *nDomErr);
/*****************************************************************/
/*   Purpose:  Safe version of log10.                            */
/*             On return, log10x[i] = NA_REAL if x[i] = NA_REAL, */
/*                                    NA_REAL if x[i] <= 0.0,    */
/*                                    log10(x[i]) otherwise,     */
/*             *nInputNA counts the number of NA_REALs in x,     */
/*             and *nDomErr counts the number of x[i] <= 0       */
/*             occurrences.                                      */
/*****************************************************************/

/*****************************************************************/
real RootMSE(size_t n, const real *y1, const real *y2,
          real *MaxErr, size_t *IndexMaxErr);
/*****************************************************************/
/*   Purpose:  Compute the root MSE and the maximum error        */
/*             between y1 and y2, ignoring any cases with NA's.  */
/*                                                               */
/*   Args:     n         Length of vectors.                      */
/*                       y1        First vector.                 */
/*             y2        Second vector.                          */
/*             MaxErr    Output: Maximum absolute error.         */
/*             IndexMaxErr                                       */
/*                       Output: Case with the maximum absolute  */
/*                       error.                                  */
/*                                                               */
/*   Returns:  The root MSE.                                     */
/*****************************************************************/

/*****************************************************************/
int LevelLex(size_t n, const size_t *nLevels, size_t *Level);
/*****************************************************************/
/*   Purpose:  Generate the next level combination of the n      */
/*             levels in Level (lexicographic order), where      */
/*             variable j has levels 0,..., nLevels[j].          */
/*                                                               */
/*   Returns:  ALL_DONE if no more combinations remain;          */
/*             OK       otherwise.                               */
/*****************************************************************/

/*****************************************************************/
size_t Combinations(size_t n, size_t m);
/*****************************************************************/
/*   Purpose:  Return n choose m.                                */
/*****************************************************************/

/*****************************************************************/
real Pythag(real a, real b);
/*****************************************************************/
/*   Purpose:  Return sqrt(a^2 + b^2).                           */
/*****************************************************************/


/* libout.c: */

#define SEV_WARNING      0
#define SEV_ERROR        1
#define SEVERITY_STRS    {"Warning", "Error"}

void      SetLogFile(FILE *fp);
FILE      *GetLogFile(void);

/*****************************************************************/
void FileOutput(FILE *fp, const string Format, ...);
/*****************************************************************/
/*   Purpose:  Send output in variable argument list to fp.      */
/*             If fp == stdout, then also send output to the log */
/*             file.                                             */
/*****************************************************************/

void      Output(const string Format, ...);
void      Error(const string Format, ...);
void      Fatal(const string Format, ...);
void      Incompatibility(const string Format, ...);

/*****************************************************************/
void OutputVA(const string Format, va_list Args);
/*****************************************************************/
/*   Purpose:  Send output in variable argument list to stdout   */
/*             and to the log file.                              */
/*****************************************************************/

/*****************************************************************/
void ErrorToMat(const string Severity, const string Format,
          va_list Args);
/*****************************************************************/
/*   Purpose:  Save error, warning, etc. messages in a matrix.   */
/*****************************************************************/

/*****************************************************************/
void ErrorMatOut(void);
/*****************************************************************/
/*   Purpose:  Output matrix of error, warning, etc. messages.   */
/*****************************************************************/

/*****************************************************************/
void OutputTemp(const string Format, ...);
/*****************************************************************/
/*   Purpose:  Output temporary message to stdout, which will be */
/*             overwritten by next temporary message.            */
/*****************************************************************/

/*****************************************************************/
void OutTime(void);
/*****************************************************************/
/* Purpose:    Output time and date.                             */
/*****************************************************************/

void      OutDataRow(FILE *OutFile, const real *x, size_t n);

/* #define CodeBug(FuncName)     Fatal("Code bug in %s.\n", FuncName)*/
/* #define IllegalType(FuncName) \
          (Fatal("Codebug: Illegal type in %s.\n", FuncName)) */

#define CodeBug(problem) \
{ \
     Output("Code bug detected: %s, file %s, line %d\n", \
               problem, __FILE__, __LINE__); \
     exit(1); \
}

#define CodeCheck(expression) \
{ \
     if (!(expression)) \
     { \
          Output("Code check failed: %s, file %s, line %d\n", \
                    #expression, __FILE__, __LINE__); \
          exit(1); \
     } \
}

/* libperm.c: */

/*****************************************************************/
void PermRand(size_t n, size_t *Perm);
/*****************************************************************/
/*   Purpose:  Generate a random permutation of the n            */
/*             nonnegative integers (not necessarily 1,...,n)    */
/*             in Perm.                                          */
/*****************************************************************/

/*****************************************************************/
int PermLex(size_t n, size_t *Perm);
/*****************************************************************/
/*   Purpose:  Generate the next permutation (lexicographic      */
/*             order) of the n nonnegative integers (not         */
/*             necessarily 1,...,n) in Perm.                     */
/*                                                               */
/*   Returns:  ALL_DONE if no more permutations remain;          */
/*             OK       otherwise.                               */
/*****************************************************************/


/* libprob.c: */

/*****************************************************************/
real CDFNorm(real z, boolean Lower);
/*****************************************************************/
/*   Purpose:  Returns tail area of the standard normal          */
/*             distribution, lower tail if Lower is TRUE, upper  */
/*             tail otherwise.                                   */
/*****************************************************************/

/*****************************************************************/
real CDFInvNorm(real p);
/*****************************************************************/
/*   Purpose:  Returns inverse of the normal distribution for    */
/*             probability p; i.e., the integral from -infinity  */
/*             to the returned value is p.                       */
/*                                                               */
/*   Returns:  Normal quantile if 0 <= p <= 1;                   */
/*             -REAL_MAX       otherwise.                        */
/*****************************************************************/

/*****************************************************************/
real PDFNorm(real z);
/*****************************************************************/
/*   Purpose:  Returns p.d.f. of standard normal distribution.   */
/*****************************************************************/

real      Cor(real *x, real *y, size_t n);


/* librandn.c: */

int       RandInit(int _xcomp, int _ycomp, int _zcomp);
real      RandUnif(void);


/* libreg.c: */

/* Column names and types for region matrices. */
#define REG_COL_NAMES    {VARIABLE, SUPPORT, MIN, MAX, NUM_LEVELS, \
                          NUM_CATS, DISTRIBUTION, INCLUSIVE, \
                          CAND_GROUP, STEP, "CandMatIndex", WEIGHT, \
                          ANALYZE, TRANSFORMATION}
#define REG_COL_TYPES    {STRING, SIZE_T, REAL, REAL, SIZE_T, \
                          SIZE_T, SIZE_T, INTEGER, \
                          SIZE_T, REAL, SIZE_T, REAL, \
                          INTEGER, SIZE_T}

#define RegNumVars(R)              MatNumRows(R)
#define RegVar(R, i)               MatStrElem(   R, i,  0)
#define RegSupport(R, i)           MatSize_tElem(R, i,  1)
#define RegMin(R, i)               MatElem(      R, i,  2)
#define RegMax(R, i)               MatElem(      R, i,  3)
#define RegNumLevels(R, i)         MatSize_tElem(R, i,  4)
#define RegNumCats(R, i)           MatSize_tElem(R, i,  5)
#define RegDistrib(R, i)           MatSize_tElem(R, i,  6)
#define RegInclusive(R, i)         MatIntElem(   R, i,  7)
#define RegCandGroup(R, i)         MatSize_tElem(R, i,  8)
#define RegStep(R, i)              MatElem(      R, i,  9)
#define RegCandMatIndex(R, i)      MatSize_tElem(R, i, 10)
#define RegWt(R, i)                MatElem(      R, i, 11)
#define RegAnalyze(R, i)           MatIntElem(   R, i, 12)
#define RegTransformation(R, i)    MatSize_tElem(R, i, 13)

#define RegPutVar(R, i, s)         MatPutStrElem(   R, i,  0, s)
#define RegPutSupport(R, i, z)     MatPutSize_tElem(R, i,  1, z)
#define RegPutMin(R, i, r)         MatPutElem(      R, i,  2, r)
#define RegPutMax(R, i, r)         MatPutElem(      R, i,  3, r)
#define RegPutNumLevels(R, i, z)   MatPutSize_tElem(R, i,  4, z)
#define RegPutNumCats(R, i, z) \
                                   MatPutSize_tElem(R, i,  5, z)
#define RegPutDistrib(R, i, z)     MatPutSize_tElem(R, i,  6, z)
#define RegPutInclusive(R, i, b)   MatPutIntElem(   R, i,  7, b)
#define RegPutCandGroup(R, i, z)   MatPutSize_tElem(R, i,  8, z)
#define RegPutStep(R, i, r)        MatPutElem(      R, i,  9, r)
#define RegPutCandMatIndex(R, i, z) \
                                   MatPutSize_tElem(R, i, 10, z)
#define RegPutWt(R, i, r)          MatPutElem(      R, i, 11, r)
#define RegPutAnalyze(R, i, b)     MatPutIntElem(   R, i, 12, b)
#define RegPutTransformation(R, i, z) \
                                   MatPutSize_tElem(R, i, 13, z)

#define RegRand(Reg, i)  RegTransform(RandUnif(), Reg, i)

/* Implemented support types. */
#define SUPPORT_NAMES    {FIXED_STR, CONTINUOUS_STR, GRID_STR, \
                          DISCRETE_STR}
#define FIXED       0
#define CONTINUOUS  1
#define GRID        2
#define DISCRETE    3

/* Implemented distributions. */
#define DISTRIB_NAMES    {ARCSIN_STR, ARCTAN_STR, LOGUNIFORM_STR, \
                          NORMAL_STR, UNIFORM_STR}
#define ARCSIN      0
#define ARCTAN      1
#define LOGUNIFORM  2
#define NORMAL      3
#define UNIFORM     4

/* Implemented transformations. */
#define TRANSFORM_NAMES  {NONE_STR, LOG_STR}
#define NONE        0
#define LOG         1

/*****************************************************************/
void RegAlloc(size_t nVars, Matrix *Reg);
/*****************************************************************/
/*   Purpose:  Allocate space for a region matrix with nVars     */
/*             variables.                                        */
/*****************************************************************/

/*****************************************************************/
int RegExtract(const Matrix *XDescrip, const string XDescripName,
     const string ColExt, Matrix *Reg);
/*****************************************************************/
/*   Purpose:  Extract region Reg from XDescrip.                 */
/*                                                               */
/*   Returns:  INPUT_ERR    if user's matrix is illegal;         */
/*             OK           otherwise.                           */
/*****************************************************************/

/*****************************************************************/
int RegCandCompat(const Matrix *Cand, Matrix *Reg);
/*****************************************************************/
/*   Purpose:  Update Reg for DISCRETE variables.                */
/*                                                               */
/*   Returns:  INCOMPAT_ERR if there is an incompatibility;      */
/*             OK           otherwise.                           */
/*****************************************************************/

/*****************************************************************/
void RegRandPt(const Matrix *Reg, real *x);
/*****************************************************************/
/*   Purpose:  Generate a random point x in region Reg.          */
/*                                                               */
/*   Comment:  Calling routine must allocate space for x.        */
/*****************************************************************/

/*****************************************************************/
real RegTransform(real u, const Matrix *Reg, size_t j);
/*****************************************************************/
/*   Purpose:  Transform Uniform[0, 1] u into valid value for    */
/*             variable j of XReg.                               */
/*                                                               */
/*   Returns:  Transformed value.                                */
/*****************************************************************/

/*****************************************************************/
real RegTransformCont(real u, real a, real b, size_t DistribNum);
/*****************************************************************/
/*   Purpose:  Transform Uniform[0, 1] u into the continuous     */
/*             region [a, b].                                    */
/*                                                               */
/*   Returns:  Transformed value.                                */
/*****************************************************************/

/*****************************************************************/
real RegLevel(const Matrix *Reg, size_t j, size_t LevelIndex);
/*****************************************************************/
/*   Purpose:  Return level LevelIndex of variable j.            */
/*****************************************************************/

/*****************************************************************/
real RegLevelWt(const Matrix *Reg, size_t j, size_t LevelIndex);
/*****************************************************************/
/*   Purpose:  Return (integration) weight for level LevelIndex  */
/*             of variable j.                                    */
/*****************************************************************/

/*****************************************************************/
size_t RegGroupIndices(const Matrix *Reg, size_t j,
     size_t *Index);
/*****************************************************************/
/*   Purpose:  On return, Index contains the indices of the      */
/*             variables in the same group as variable j         */
/*             (including j).                                    */
/*                                                               */
/*   Returns:  Number of variables in the group                  */
/*             (1 for ungrouped).                                */
/*                                                               */
/*   Comments: Calling routine must allocate space for Index.    */
/*****************************************************************/

/*****************************************************************/
size_t RegGroupings(const Matrix *Reg, size_t **GroupSize,
          Matrix *Index);
/*****************************************************************/
/*   Purpose:  Returns the number of groups.  On return,         */
/*             *GroupSize[i] is the number of variables in group */
/*             i, and column i of Index contains the indices of  */
/*             the variables in group i.                         */
/*                                                               */
/*   Returns:  Number of groups.                                 */
/*                                                               */
/*   Comments: GroupSize and Index allocated here.               */
/*****************************************************************/

/*****************************************************************/
void RegLevelsGroup(const Matrix *Reg, size_t GroupSize,
          const size_t *Index, size_t LevelIndex, real *x);
/*****************************************************************/
/*   Purpose:  Load level LevelIndex of the GroupSize variables  */
/*             with indices Index into x.                        */
/*                                                               */
/*   Comments: Calling routine must allocate space for x.        */
/*****************************************************************/

/*****************************************************************/
boolean RegIsCand(const Matrix *Reg);
/*****************************************************************/
/*   Purpose:  Is Reg a candidate set (i.e., all discrete        */
/*             variables in the same group)?                     */
/*****************************************************************/


/* libsort.c: */

void QuickIndex(const real *x, size_t n, size_t *Index);
int  CompIndex();
void QuickRank(const real *x, size_t n, size_t *Rank);

/*****************************************************************/
void QuickReal(size_t n, const real *x);
/*****************************************************************/
/*   Purpose:  Sort x[0], ..., x[n-1], smallest to largest.      */
/*                                                               */
/*   Comment:  Could be adapted to find ranks for any type of    */
/*             item using (void *), etc.                         */
/*****************************************************************/

/*****************************************************************/
int CompReal(/* real *x1, real *x2 */);
/*****************************************************************/
/*   Purpose:  Compare *x1 and *x2 for qsort called by QuickReal.*/
/*****************************************************************/

/* libstr.c: */

string    StrCatAlloc(const string s1, const string s2);
string    StrPaste(size_t n, ...);
int       stricmp(const char *s, const char *t);
int       strnicmp(const char *s, const char *t, size_t n);
size_t    StrVecCmp(const string *s, const string *t, size_t n);
string    StrDup(const string s);

/*****************************************************************/
string StrReplace(const string NewStr, string Target);
/*****************************************************************/
/*   Purpose:  Replace Target with NewStr.                       */
/*                                                               */
/*   Returns:  New pointer to Target.                            */
/*****************************************************************/

/*****************************************************************/
int StrToInt(const string s, int *i);
/*****************************************************************/
/*   Purpose:  Convert a string to an integer.                   */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*****************************************************************/

/*****************************************************************/
int StrToReal(const string s, real *r);
/*****************************************************************/
/*   Purpose:  Convert a string to a real.                       */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*****************************************************************/

/*****************************************************************/
int StrToSize_t(const string s, size_t *z);
/*****************************************************************/
/*   Purpose:  Convert a string to a size_t.                     */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*****************************************************************/

/*****************************************************************/
string StrFromInt(int i);
/*****************************************************************/
/*   Purpose:  Return a string corresponding to an integer.      */
/*                                                               */
/*   Comment:  Calling routine should duplicate the string       */
/*             before the next call of StrFromInt,               */
/*             StrFromSize_t, or StrFromReal.                    */
/*****************************************************************/

/*****************************************************************/
string StrFromReal(real r, const string Flags, int Precision,
          char Conversion);
/*****************************************************************/
/*   Purpose:  Return a string corresponding to a real.          */
/*             Flags are as defined for sprintf.                 */
/*             Precision is as defined for sprintf ( < 0 for     */
/*             default).                                         */
/*             Conversion can be 'g', 'e', or 'f'.               */
/*                                                               */
/*   Comment:  Calling routine should duplicate the string       */
/*             before the next call of StrFromInt,               */
/*             StrFromSize_t, or StrFromReal.                    */
/*****************************************************************/

/*****************************************************************/
string StrFromSize_t(size_t z);
/*****************************************************************/
/*   Purpose:  Return a string corresponding to a size_t.        */
/*                                                               */
/*   Comment:  Calling routine should duplicate the string       */
/*             before the next call of StrFromInt,               */
/*             StrFromSize_t, or StrFromReal.                    */
/*****************************************************************/

size_t    StrIndex(const string Target, const string *StrArray,
               size_t NumStrings);

/*****************************************************************/
size_t StrNumberOf(const string *StrArray);
/*****************************************************************/
/*   Purpose:  Return number of strings in StrArray, not         */
/*             counting the terminating NULL string.             */
/*****************************************************************/

int       StrBrackets(string Token, string *Sub, string *NextToken);

#define NumStr(s)    (sizeof(s) / sizeof(string *))
#define StrPlural(n) (((n) != 1) ? "s" : "")


/* libtempl.c: */

/*****************************************************************/
List *TemplAlloc(const string Name, int Type, List *L);
/*****************************************************************/
/*   Purpose:  Allocate a template of the specified name and     */
/*             type by adding it to the list L.                  */
/*                                                               */
/*   Returns:  Pointer to the updated list.                      */
/*****************************************************************/

/*****************************************************************/
template *TemplPtr(const string Name, const List *L);
/*****************************************************************/
/*   Purpose:  Return template corresponding to Name.            */
/*****************************************************************/

/*****************************************************************/
List *TemplIntAlloc(const string Name, int Min, int Max, List *L);
/*****************************************************************/
/*   Purpose:  Allocate an INTEGER template.                     */
/*                                                               */
/*   Comment:  TemplRealAlloc, TemplSize_tAlloc, and             */
/*             TemplStrAlloc are similar.                        */
/*                                                               */
/*   Returns:  Pointer to the updated list.                      */
/*****************************************************************/

List *TemplRealAlloc(const string Name, real Min, real Max,
          List *L);
List *TemplSize_tAlloc(const string Name, size_t Min, size_t Max,
          List *L);
List *TemplStrAlloc(const string Name, const string *LegalStr,
          size_t NumLegalStr, List *L);

/*****************************************************************/
size_t TemplIntCheck(const template *T, const int *i, size_t n);
/*****************************************************************/
/*   Purpose:  Check an array of values against a template.      */
/*                                                               */
/*   Returns:  j        if the value with index j is illegal;    */
/*             INDEX_OK otherwise.                               */
/*                                                               */
/*   Comment:  TemplRealCheck, TemplSize_tCheck, and             */
/*             TemplStrCheck are similar.                        */
/*****************************************************************/

size_t TemplRealCheck(const template *T, const real *r, size_t n);
size_t TemplSize_tCheck(const template *T, const size_t *z,
          size_t n);
size_t TemplStrCheck(const template *T, const string *s, size_t n);

/*****************************************************************/
void TemplError(const template *T);
/*****************************************************************/
/*   Purpose:  Output an error message.                          */
/*****************************************************************/


/* libvec.c */

/*****************************************************************/
size_t VecTempl(const vector *v, const template *T);
/*****************************************************************/
/*   Purpose:  Check a vector against a template.                */
/*                                                               */
/*   Returns:  j        if element j of v is illegal;            */
/*             INDEX_OK otherwise.                               */
/*****************************************************************/

/*****************************************************************/
void VecToStr(const vector *v, string *s);
/*****************************************************************/
/*   Purpose:  Convert a vector to an array of strings.          */
/*                                                               */
/*   Comment:  Calling routine should allocate s to be of        */
/*             sufficient length, but the strings themselves are */
/*             StrDup'd here.                                    */
/*****************************************************************/

/*****************************************************************/
size_t VecFromStr(const string *s, vector *v);
/*****************************************************************/
/*   Purpose:  Convert an array of strings to a vector.          */
/*                                                               */
/*   Returns:  j        if s[j] is illegal for the vector type;  */
/*             INDEX_OK otherwise.                               */
/*****************************************************************/
