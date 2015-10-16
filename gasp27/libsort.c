/*****************************************************************/
/*   ROUTINES FOR SORTING, RANKS, ETC.                           */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--94.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

/*******************************+++*******************************/
/*                                                               */
/*   void      QuickIndex(const real *x, size_t n,               */
/*                  size_t *Index)                               */
/*                                                               */
/*   Purpose:  Build index table.                                */
/*                                                               */
/*   Args:     x         Items to index.                         */
/*             n         Number of items.                        */
/*             Index     Output: x[Index[0]] is the smallest     */
/*                       element of x, etc.                      */
/*                                                               */
/*   Comment:  Could be adapted to find index table for any type */
/*             of item using (void *), etc.                      */
/*                                                               */
/*   Version:  1991 June 18                                      */
/*                                                               */
/*****************************************************************/

static const real *xx;   /* Used to communicate with CompIndex. */

void QuickIndex(const real *x, size_t n, size_t *Index)
{
     size_t    i;

     xx = x;

     for (i = 0; i < n; i++)
          Index[i] = i;

     qsort((void *) Index, n, sizeof(size_t), CompIndex);

     return;
}

/*******************************+++*******************************/
/*                                                               */
/*   int       CompIndex(size_t *i1, size_t *i2)                 */
/*                                                               */
/*   Purpose:  Compare xx[*i1] and xx[*i2] (i.e., x[*i1] and     */
/*             x[*i2] in QuickIndex) for qsort called by         */
/*             QuickIndex.                                       */
/*                                                               */
/*   Version:  1991 June 10                                      */
/*                                                               */
/*****************************************************************/

int CompIndex(size_t *i1, size_t *i2)
{
     if (xx[*i1] < xx[*i2])
          return -1;
     else if (xx[*i1] > xx[*i2])
          return 1;
     else
          return 0;
}

/*******************************+++*******************************/
/*                                                               */
/*   void      QuickRank(const real *x, size_t n, size_t *Rank)  */
/*                                                               */
/*   Purpose:  Rank x[0], ..., x[n-1].                           */
/*                                                               */
/*   Args:     x         Items to rank.                          */
/*             n         Number of items.                        */
/*             Rank      Output: Rank[i] is the rank of x[i].    */
/*                                                               */
/*   Comment:  Could be adapted to find ranks for any type of    */
/*             item using (void *), etc.                         */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*                                                               */
/*****************************************************************/

void QuickRank(const real *x, size_t n, size_t *Rank)
{
     size_t    *Index = NULL;
     size_t    i;

     /* Build index table. */
     Index = AllocSize_t(n, NULL);
     QuickIndex(x, n, Index);

     /* Convert indices to ranks. */
     for (i = 0; i < n; i++)
          Rank[Index[i]] = i;

     AllocFree(Index);
}

/*******************************+++*******************************/
void QuickReal(size_t n, const real *x)
/*****************************************************************/
/*   Purpose:  Sort x[0], ..., x[n-1], smallest to largest.      */
/*                                                               */
/*   Comment:  Could be adapted to find ranks for any type of    */
/*             item using (void *), etc.                         */
/*                                                               */
/*   Version:  1994 February 27                                  */
/*****************************************************************/
{
     qsort(x, n, sizeof(real), CompReal);
}

/*******************************+++*******************************/
int CompReal(real *x1, real *x2)
/*****************************************************************/
/*   Purpose:  Compare *x1 and *x2 for qsort called by QuickReal.*/
/*                                                               */
/*   Version:  1994 February 27                                  */
/*****************************************************************/
{
     if (*x1 < *x2)
          return -1;
     else if (*x1 > *x2)
          return 1;
     else
          return 0;
}

