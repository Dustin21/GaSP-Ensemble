/*****************************************************************/
/*   Copyright (c) William J. Welch 1990--99.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*   Version: 1999.06.17                                         */
/*****************************************************************/

#define MAT_DEFINED

/* Note: All matrix types are stored by column. */

/* Matrix shapes: */
#define   RECT      0    /* Rectangular.      */
#define   SYM       1    /* Symmetric.        */
#define   UP_TRIANG 2    /* Upper triangular. */

/* An UP-TRIANG matrix can also be viewed as lower-triangular, */
/* stored by rows.                                             */

typedef struct MatrixStruct
{
     size_t    NumRows;
     size_t    NumCols;
     int       Shape;         /* RECT,... as above. */
     int       Type;          /* INTEGER, MIXED, REAL, SIZE_T, */
                              /* or STRING.                    */
     int       *ColType;      /* Column types (INTEGER, REAL, */
                              /* SIZE_T or STRING).           */
     int       **IntElem;
     real      **Elem;
     size_t    **Size_tElem;
     string    **StrElem;
     boolean   Labelled;
     string    Text;
     string    *RowName;
     string    *ColName;
     boolean   Initialized;
     struct MatrixStruct *Next;
} matrix;

/* Temporary until Matrix replaced by matrix everywhere. */
typedef matrix Matrix;

/* These macros allow user routines to access members of the  */
/* Matrix structure, while hiding the internal data           */
/* representations.                                           */
/* Note, they all work on a pointer to a Matrix.              */

#define MatNumRows(M)         ((M)->NumRows)
#define MatPutNumRows(M,n)    ((M)->NumRows = (n))     /* Care! */
#define MatNumCols(M)         ((M)->NumCols)
#define MatPutNumCols(M,n)    ((M)->NumCols = (n))     /* Care! */

#define MatShape(M)           ((M)->Shape)
#define MatPutShape(M,NewShape) \
                              ((M)->Shape = NewShape)
#define MatColLen(M,j)        ((M)->Shape == RECT ? \
                                   (M)->NumRows : j + 1)

#define MatType(M)            ((M)->Type)
#define MatPutType(M,NewType) ((M)->Type = NewType)
#define MatColType(M,j)       ((M)->ColType[j])
#define MatPutColType(M,j,Type) \
                              ((M)->ColType[j] = (Type))
#define MatColTypes(M)        ((M)->ColType)

/* Return the column index corresponding to a column name. */
#define MatColIndex(M,Name) \
     ( (MatLabelled(M)) ? \
     StrIndex(Name,(M)->ColName, (M)->NumCols) : INDEX_ERR)

/* Assign column pointers. */
#define MatPutIntCol(M,j,Ptr) ((M)->IntElem[j] = (Ptr))
#define MatPutCol(M,j,Ptr)    ((M)->Elem[j] = (Ptr))
#define MatPutSize_tCol(M,j,Ptr) \
                              ((M)->Size_tElem[j] = (Ptr))
#define MatPutStrCol(M,j,Ptr) ((M)->StrElem[j] = (Ptr))

/* Return an individual element. */
#define MatIntElem(M,i,j)     ((M)->IntElem[j][i])
#define MatElem(M,i,j)        ((M)->Elem[j][i])
#define MatSize_tElem(M,i,j)  ((M)->Size_tElem[j][i])
#define MatStrElem(M,i,j)     ((M)->StrElem[j][i])

/* Put an individual element into a matrix. */
#define MatPutIntElem(M,i,j,e)     ((M)->IntElem[j][i] = e)
#define MatPutElem(M,i,j,e)        ((M)->Elem[j][i] = e)
#define MatPutSize_tElem(M,i,j,e)  ((M)->Size_tElem[j][i] = e)
#define MatPutStrElem(M,i,j,e) \
          ((M)->StrElem[j][i] = StrReplace(e, (M)->StrElem[j][i]))

#define MatLabelled(M)        ((M)->Labelled)
#define MatText(M)            ((M)->Text)
#define MatPutText(M,NewText) \
          ((M)->Text = StrReplace(NewText, (M)->Text))

#define MatRowNames(M)        ((M)->RowName)    /* All row names. */
/* Particular row name (default is row number). */
#define MatRowName(M,i) \
          (((M)->RowName != NULL && (M)->RowName[i] != NULL) ? \
          (M)->RowName[i] : StrFromSize_t((i) + 1))
#define MatPutRowName(M,i,Name) \
          ((M)->RowName[i] = StrReplace(Name, (M)->RowName[i]))
#define MatColNames(M)        ((M)->ColName)    /* All col names. */
/* Particular column name.*/
#define MatColName(M,j)       ((M)->ColName != NULL ? \
                                   (M)->ColName[j] : NULL)
#define MatPutColName(M,j,Name) \
          ((M)->ColName[j] = StrReplace(Name, (M)->ColName[j]))

#define MatPutInit(M, b)      ((M)->Initialized = b)
#define MatInitialized(M)     ((M)->Initialized)

#define MatEmpty(M)           ((M)->NumRows == 0 && \
                                        (M)->NumCols == 0)

/* matalloc.c: */

/*****************************************************************/
void MatInit(int Shape, int Type, boolean Labelled, Matrix *M);
/*****************************************************************/
/*   Purpose:  Initialize a matrix as 0 x 0.                     */
/*****************************************************************/

void MatReAllocate(size_t NewNumRows, size_t NewNumCols,
          const int *NewColType, Matrix *M);
void MatFree(Matrix *M);
void MatColReAlloc(size_t NewLen, size_t j, Matrix *M);
size_t MatColumnAdd(const string ColName, int NewColType, Matrix *M);

#define MatReAlloc(NewNumRows, NewNumCols, M) \
          MatReAllocate(NewNumRows, NewNumCols, NULL, M)
#define MatAllocate(NumRows,NumCols,Shape,Type,ColType,Labelled,M)\
          {MatInit(Shape, Type, Labelled, M); \
           MatReAllocate(NumRows, NumCols, ColType, M);}
#define MatAlloc(NumRows,NumCols,Shape,M) \
          MatAllocate(NumRows, NumCols, Shape, REAL, NULL, NO, M)

#define MatIntColAdd(ColName, M) \
          (MatIntCol(M, MatColumnAdd(ColName, INTEGER, M)))
#define MatColAdd(ColName, M) \
          (MatCol(M, MatColumnAdd(ColName, REAL, M)))
#define MatSize_tColAdd(ColName, M) \
          (MatSize_tCol(M, MatColumnAdd(ColName, SIZE_T, M)))
#define MatStrColAdd(ColName, M) \
          (MatStrCol(M, MatColumnAdd(ColName, STRING, M)))

/*****************************************************************/
void MatRowAdd(size_t nVars, const string *VarName,
          const real *row, matrix *M);
/*****************************************************************/
/*   Purpose:  Put the nVars variables in VarNames into a new    */
/*             row of M.                                         */
/*****************************************************************/


/* matblas.c: */

/*****************************************************************/
real VecDotProd(size_t n, const real *a, const real *b);
/*****************************************************************/
/* Purpose:    Return a[0] * b[0] + ... + a[n-1] * b[n-1].       */
/*                                                               */
/* Comment:    Use VecDotProdRange for dot product of elements i */
/*             through j.                                        */
/*****************************************************************/

/* Old parameter order: REMOVE. */
#define DotProd(a, b, n)      VecDotProd(n, a, b)

/*****************************************************************/
real VecDotProdRange(size_t i, size_t j, const real *a,
     const real *b);
/*****************************************************************/
/* Purpose:    Return a[i] * b[i] + ... + a[j] * b[j].           */
/*                                                               */
/* Comment:    Use VecDotProd for dot product of first n         */
/*             elements.                                         */
/*****************************************************************/

real      VecSS(const real *a, size_t n);    /* REORDER PARS. */
real      VecSSCor(const real *a, size_t n); /* REORDER PARS. */
void      VecInit(real s, size_t n, real *a);
void      VecStrInit(const string s, size_t n, string *a);
real      VecSum(const real *a, size_t n);   /* REORDER PARS. */
real      VecSumAbs(size_t n, const real *a);
void      VecAddScalar(real s, size_t n, real *a);
void      VecMultScalar(real s, size_t n, real *a);

/* REORDER PARS. (n, s, a, b) */
void      VecAddVec(real s, const real *a, size_t n, real *b);

#define   VecAddVecRange(i, j, s, a, b) \
               VecAddVec(s, (a) + (i), \
               (i <= j) ? (j) - (i) + 1 : 0, (b) + (i))
                    
void      VecAddVecNew(size_t n, const real *a, const real *b,
               real *c);
void      VecMultVec(real *a, size_t n, real *b); /* REORDER PARS. */
real      VecMax(const real *a, size_t n);        /* REORDER PARS. */
real      VecMin(const real *a, size_t n);        /* REORDER PARS. */

/*****************************************************************/
boolean VecHasNA(size_t n, real *a);
/*****************************************************************/
/*   Purpose:  Does a[0],..., a[n-1] contain at least one NA?    */
/*****************************************************************/

void      GivRot(real *a, real *b, real *c, real *s);
void      VecIntCopy(const int *a, size_t n, int *b);  /* REORDER PARS. */
void      VecCopy(const real *a, size_t n, real *b);    /* REORDER PARS. */
void      VecSize_tCopy(const size_t *a, size_t n, size_t *b);/* REORDER PARS. */
void      VecStrCopy(const string *a, size_t n, string *b);/* REORDER PARS. */

/*****************************************************************/
void VecCopyStride(size_t n, size_t aStride, const real *a,
          size_t bStride, real *b);
/*****************************************************************/
/*   Purpose:  Copy a[i*aStride] to b[i*bStride] for             */
/*             i = 0,..., n - 1.                                 */
/*****************************************************************/

void      VecCopyStride(size_t n, size_t aStride, const real *a,
               size_t bStride, real *b);
void      VecCopyIndex(size_t n, const size_t *aIndex, const real *a,
               const size_t *bIndex, real *b);
/* REORDER PARS. */
string    *VecIntStr(const int *a, size_t n, string *s);
/* REORDER PARS. */
string    *VecRealStr(const real *a, size_t n, string *s);
/* REORDER PARS. */
string    *VecSize_tStr(const size_t *a, size_t n, string *s);
/* REORDER PARS. */
size_t    VecStrInt(const string *s, size_t n, int *a);
/* REORDER PARS. */
size_t    VecStrReal(const string *s, size_t n, real *a);
/* REORDER PARS. */
size_t    VecStrSize_t(const string *s, size_t n, size_t *a);

/*****************************************************************/
size_t VecSize_tIndex(size_t Target, size_t n, const size_t *a);
/*****************************************************************/
/*   Purpose:  Find first i such that Target matches a[i].       */
/*                                                               */
/*   Returns:  i          if a[i] matches Target;                */
/*             INDEX_ERR  otherwise.                             */
/*                                                               */
/*   Comments: Combine with StrIndex in libstr.c, etc.           */
/*****************************************************************/

/* Shouldn't this be MatNumRows?  Where is this used? */
/* #define MatColInit(M, j, r) \
          VecInit(r, MatNumCols(M), MatCol(M, j)) */

/* matcopy.c: */

/*****************************************************************/
void MatCopySub(size_t m, size_t n, size_t SrcRowOffset,
          size_t SrcColOffset, const Matrix *Src,
          size_t DestRowOffset, size_t DestColOffset, Matrix *Dest);
/*****************************************************************/
/*   Purpose:  Copy m x n submatrix from Src, starting at        */
/*             (SrcRowOffset, SrcColOffset), to Dest, starting   */
/*             at (DestRowOffset, DestColOffset).                */
/*                                                               */
/*   Comments: Calling routine must allocate space for Dest.     */
/*             Labels are not copied.                            */
/*             MatCopy() is a macro with a simplified            */
/*             parameter list for copying an entire matrix.      */
/*****************************************************************/

/*****************************************************************/
void MatDup(const Matrix *Src, Matrix *Dest);
/*****************************************************************/
/*   Purpose:  Duplicates matrix Src to matrix Dest.             */
/*                                                               */
/*   Comments: Space is allocated here for Dest.                 */
/*****************************************************************/

/*****************************************************************/
void MatDupIndex(size_t nIn, const size_t *In, const Matrix *M,
          Matrix *NewM);
/*****************************************************************/
/*   Purpose:  Select rows indexed by In from M and put them in  */
/*             NewM.                                             */
/*                                                               */
/*   Comment:  NewM is allocated here.                           */
/*****************************************************************/

/*****************************************************************/
void MatCopyColSub(size_t m, size_t j, size_t SrcOffset,
          const Matrix *Src, size_t k, size_t DestOffset,
          Matrix *Dest);
/*****************************************************************/
/*   Purpose:  Copy m elements of column j of Src, starting at   */
/*             SrcOffset, to column k of Dest, starting at       */
/*             DestOffset.                                       */
/*                                                               */
/*   Comment:  MatCopyCol() is a macro with a simplified         */
/*             parameter list for copying an entire column.      */
/*****************************************************************/

/*****************************************************************/
void MatCopyRow(size_t i, const Matrix *Src, size_t k,
          Matrix *Dest);
/*****************************************************************/
/*   Purpose:  Copy row i of Src to row k of Dest.               */
/*             If both matrices are labelled, then the row label */
/*             is also copied.                                   */
/*                                                               */
/*   Comments: Calling routine must allocate space for Dest.     */
/*             Both Src and Dest must be RECT matrices.          */
/*****************************************************************/

/* Copy column j of Src to column k of Dest. */
#define MatCopyCol(j, Src, k, Dest) \
     MatCopyColSub(MatNumRows(Src), j, 0, Src, k, 0, Dest)

/* Copy Src to Dest, returns Dest. */
#define MatCopy(Src, Dest)    MatCopySub(MatNumRows(Src), \
          MatNumCols(Src), 0, 0, Src, 0, 0, Dest)


/* mateig.c: */

/*****************************************************************/
int MatEig(boolean SortValues, matrix *S, real *eVal, matrix *V);
/*****************************************************************/
/* Purpose:    Compute eigenvalue decomposition, V Diag(e) V',   */
/*             of symmetric matrix S, where the columns of V are */
/*             the eigenvectors and eVal contains the            */
/*             eigenvalues.                                      */
/*             If SortValues is TRUE, then the eigenvalues (and  */
/*             corresponding eigenvectors) are sorted largest to */
/*             smallest.                                         */
/*                                                               */
/* Returns:    NUMERIC_ERR if MatEigTriDiag fails to converge;   */
/*             OK          otherwise.                            */
/*                                                               */
/* Comments:   Calling routine must allocate space for eVal.     */
/*             V must be allocated RECT by the calling routine.  */
/*             If S and V are distinct, then S is unaltered.     */
/*             If V = S, then S is overwritten; in this case S   */
/*             must be allocated RECT.                           */
/*****************************************************************/

/*****************************************************************/
void MatTriDiag(matrix *S, real *d, real *e, matrix *Z);
/*****************************************************************/
/* Purpose:    Reduce a symmetric matrix S to Z T Z', where T is */
/*             a symmetric tridiagonal matrix, and Z is the      */
/*             orthogonal transformation matrix.                 */
/*             On exit, d contains the diagonals of T, and e     */
/*             contains the subdiagonals (with e[0] = 0.0).      */ 
/*                                                               */
/* Comments:   Calling routine must allocate space for d and e.  */
/*             Z must be allocated RECT by the calling routine.  */
/*             If S and Z are distinct, then S is unaltered.     */
/*             If Z = S, then S is overwritten; in this case S   */
/*             must be allocated RECT.                           */
/*                                                               */
/*             Translation of tred2, Handbook for Automatic      */
/*             Computation, Volume II, pp. 212--226.             */
/*****************************************************************/

/*****************************************************************/
int MatEigTriDiag(boolean SortValues, real *d, real *e, matrix *Z);
/*****************************************************************/
/* Purpose:    Compute the eigenvalues and eigenvectors of the   */
/*             symmetric tridiagonal matrix with diagonal        */
/*             elements in d and subdiagonal elements in e (e[0] */
/*             is arbitrary).                                    */
/*             On exit, d contains the eigenvalues in ascending  */
/*             order, e has been destroyed, and the columns of Z */
/*             contain the eigenvectors.                         */
/*             If SortValues is TRUE, then the eigenvalues (and  */
/*             corresponding eigenvectors) are sorted largest to */
/*             smallest.                                         */
/*             Usually, on entry d, e, and Z are from            */
/*             MatTriDiag; this gives the eigen decomposition of */
/*             a symmetric matrix.  To obtain the eigenvectors   */
/*             of a tridiagonal matrix, Z should be the identity */
/*             on entry.                                         */
/*                                                               */
/* Returns:    NUMERIC_ERR if the algorithm failed to converge;  */
/*             OK          otherwise.                            */
/*                                                               */
/* Comments:   Z must be allocated RECT by the calling routine.  */
/*                                                               */
/*             Translation of imtql2, Handbook for Automatic     */
/*             Computation, Volume II, 241--248.                 */
/*****************************************************************/


/* matio.c: */

/*****************************************************************/
int MatRead(FILE *InpFile, int Type, Matrix *M);
/*****************************************************************/
/*   Purpose:  Read a matrix from a file.                        */
/*                                                               */
/*   Returns:  OK or an error number.                            */
/*                                                               */
/*   Comment:  Only RECT, Labelled, REAL or STRING matrices can  */
/*             be read.                                          */
/*****************************************************************/

int       MatReadABlock(FILE *InpFile, int Type, Matrix *Block,
               boolean *Finished);

/*****************************************************************/
void MatWrite(Matrix *M, FILE *OutFile);
/*****************************************************************/
/*   Purpose:  Write a matrix to a file.                         */
/*                                                               */
/*   Comment:  Only RECT, Labelled matrices can be written.      */
/*             A matrix row may be written on several output     */
/*             lines.  Use MatWriteBlock() for "nice" output.    */
/*****************************************************************/

/*****************************************************************/
void MatWriteBlock(Matrix *M, boolean CaseLabels, FILE *OutFile);
/*****************************************************************/
/*   Purpose:  Write a matrix to a file in blocks of columns     */
/*             that fit on a single output line.                 */
/*                                                               */
/*   Comment:  Only RECT, Labelled matrices can be written.      */
/*****************************************************************/

/*****************************************************************/
size_t MatColWidth(const Matrix *M, size_t j, int *Precision,
               char *Conversion);
/*****************************************************************/
/*   Purpose:  Return the column width necessary for outputting  */
/*             column j of M.                                    */
/*             If the column is of type real, then, on return,   */
/*             *Precision will be the minimum precision that     */
/*             does not lose accuracy, and *Conversion will be   */
/*             'e' or 'f'.                                       */
/*****************************************************************/

/*****************************************************************/
size_t MatCaseWidth(const Matrix *M, boolean *RightAdj);
/*****************************************************************/
/*   Purpose:  Return the column width necessary for outputting  */
/*             the case column of M.                             */
/*****************************************************************/

/*****************************************************************/
string MatElemToStr(const Matrix *M, size_t i, size_t j,
          int Precision, char Conversion);
/*****************************************************************/
/*   Purpose:  Return a string representing element i, j of M.   */
/*             Precision and Conversion are only used if the     */
/*             element is real.                                  */
/*****************************************************************/

/*****************************************************************/
void MatPrint(const Matrix *M);
/*****************************************************************/
/*   Purpose:  Print a (small) real matrix of any shape to       */
/*             stdin.  Intended for debugging.                   */
/*****************************************************************/


/* matqr.c: */

size_t    QRLS(Matrix *F, real *y, Matrix *Q, Matrix *R, real *c,
               real *res);


/* matsym.c: */

/*****************************************************************/
void MatSymCol(const Matrix *S, size_t ColIndex, real *col);
/*****************************************************************/
/*   Purpose:  Load column ColIndex of S into col.               */
/*                                                               */
/*   Comment:  Calling routine must allocate space for col.      */
/*****************************************************************/

/*****************************************************************/
void MatSymUpdate(real w, const real *v, Matrix *S);
/*****************************************************************/
/*   Purpose:  Add w * v * v' to S.                              */
/*****************************************************************/

/*****************************************************************/
real MatSymQuadForm(const real *v, const Matrix *S);
/*****************************************************************/
/*   Purpose:  Return a' S a.                                    */
/*****************************************************************/


/* mattri.c: */

/*****************************************************************/
int TriForSolve(const Matrix *R, const real *b, size_t StartOff,
          real *x);
/*****************************************************************/
/*   Purpose:  Forward solve triangular system of equations      */
/*             R'x = b for x.                                    */
/*             Calling with TriForSolve(R, b, StartOff, b) will  */
/*             overwrite b with x.                               */
/*             Solving starts at element with offset StartOff    */
/*             (useful if new columns of R have been appended).  */
/*                                                               */
/*   Returns:  NONUNIQ_ERR if more than one solution is possible */
/*                         (zero is taken);                      */
/*             NUMERIC_ERR if equations cannot be solved because */
/*                         of a zero diagonal element of R;      */
/*             OK          otherwise.                            */
/*                                                               */
/*   Comment:  Calling routine must allocate space for x         */
/*             (unless x = b).                                   */
/*                                                               */
/*   Version:  1991 May 20                                       */
/*****************************************************************/

int       TriBackSolve(const Matrix *R, const real *b, real *x);
void      TriDet(const Matrix *R, real *d1, int *d2);

/*****************************************************************/
real TriCond(const Matrix *R);
/*****************************************************************/
/* Purpose:    Return estimate of condition number of R.         */
/*                                                               */
/* Comment:    From LINPACK, STRCO, p. C.86.                     */
/*****************************************************************/

/*****************************************************************/
size_t TriCholesky(const Matrix *S, size_t FirstOff, Matrix *R);
/*****************************************************************/
/*   Purpose:  Cholesky decomposition R'R of symmetric matrix S. */
/*             The decomposition starts at the column with       */
/*             offset FirstOff (useful if appending new columns).*/
/*                                                               */
/*   Returns:  j         if the leading j x j minor is not       */
/*                       positive definite;                      */
/*             OK        otherwise.                              */
/*                                                               */
/*   Comment:  See SPOFA, LINPACK, p. C.28.                      */
/*             The calling routine must allocate space for R.    */
/*             Adapted to complete decomposition if rank         */
/*             deficient.                                        */
/*             Calling with TriCholesky(S,..., S) will overwrite */
/*             S.                                                */
/*****************************************************************/

void      TriRect(const Matrix *X, Matrix *R);
void      TriUpdate(const real *xrow, real wt, Matrix *R, real *c,
               real *s);
int       TriDownDate(const real *xrow, real wt, Matrix *R,
               real *c, real *s);
void TriPerm(size_t FirstOff, size_t LastOff, Matrix *R, real *c,
               real *s);

/*****************************************************************/
void TriVec(const Matrix *R, const real *b, real *c);
/*****************************************************************/
/*   Purpose:  Compute c = R * b.                                */
/*                                                               */
/*   Comment:  Calling routine must allocate space for c.        */
/*             TriVec(R, b, b) overwrites b with R * b.          */
/*****************************************************************/


/* matutil.c: */

int       *MatIntCol(const Matrix *M, size_t j);
real      *MatCol(const Matrix *M, size_t j);
size_t    *MatSize_tCol(const Matrix *M, size_t j);
string    *MatStrCol(const Matrix *M, size_t j);
int       *MatIntColFind(const Matrix *M, const string Name,
               boolean HardFail);
real      *MatColFind(const Matrix *M, const string Name,
               boolean HardFail);
size_t    *MatSize_tColFind(const Matrix *M, const string Name,
               boolean HardFail);
string    *MatStrColFind(const Matrix *M, const string Name,
               boolean HardFail);
void      *MatVoidCol(const Matrix *M, size_t j);

/*****************************************************************/
void MatInitValue(real s, Matrix *M);
/*****************************************************************/
/*   Purpose:  Initialize all elements of the REAL matrix M to   */
/*             the value s.                                      */
/*****************************************************************/

/*****************************************************************/
void MatRow(const Matrix *M, size_t RowIndex, real *row);
/*****************************************************************/
/*   Purpose:  Load row RowIndex of M into row.                  */
/*                                                               */
/*   Comment:  Calling routine must allocate space for row.      */
/*                                                               */
/*   96.01.10  TriRow subsumed.                                  */
/*****************************************************************/

/*****************************************************************/
void MatRowPut(const real *row, size_t RowIndex, Matrix *M);
/*****************************************************************/
/*   Purpose:  Load row into row RowIndex of M.                  */
/*                                                               */
/*   Comment:  Replaced TriRowPut.                               */
/*****************************************************************/

/*****************************************************************/
void MatMultElemWise(const Matrix *A, Matrix *B);
/*****************************************************************/
/*   Purpose:  Replace B_ij by A_ij * B_ij.                      */
/*****************************************************************/

/*****************************************************************/
void MatVec(const Matrix *M, const real *x, real *y);
/*****************************************************************/
/*   Purpose:  Compute y = M'x.                                  */
/*                                                               */
/*   Comment:  Calling routine must allocate space for y.        */
/*****************************************************************/

void      MatStack(const Matrix *M, boolean ByCols, real *v);
void      MatUnStack(const real *v, boolean ByCols, Matrix *M);
int       MatMerge(Matrix *M1, Matrix *M2);
size_t    MatColConvert(size_t j, int NewColType, Matrix *M);

