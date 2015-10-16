/*****************************************************************/
/*                                                               */
/*   ROUTINES FOR GENERATING PERMUTATIONS                        */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--94.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

/*******************************+++*******************************/
void PermRand(size_t n, size_t *Perm)
/*****************************************************************/
/*   Purpose:  Generate a random permutation of the n            */
/*             nonnegative integers (not necessarily 1,...,n)    */
/*             in Perm.                                          */
/*                                                               */
/*   Version:  1994 January 28                                   */
/*****************************************************************/
{
     size_t    i, j, Temp;

     if (n == 0)
          return;

     for (i = n - 1; i > 0; i--)
     {
          j = (size_t) ((i + 1) * RandUnif());
          Temp    = Perm[i];
          Perm[i] = Perm[j];
          Perm[j] = Temp;
     }
}

/*******************************+++*******************************/
int PermLex(size_t n, size_t *Perm)
/*****************************************************************/
/*   Purpose:  Generate the next permutation (lexicographic      */
/*             order) of the n nonnegative integers (not         */
/*             necessarily 1,...,n) in Perm.                     */
/*                                                               */
/*   Returns:  ALL_DONE if no more permutations remain;          */
/*             OK       otherwise.                               */
/*                                                               */
/*   Version:  1994.01.28                                        */
/*****************************************************************/
{
     size_t i, j, Temp;

     if (n <= 1)
          return ALL_DONE;

     i = n - 1;
     while (i > 0 && Perm[i-1] >= Perm[i])
          i--;
     if (i == 0)
          return(ALL_DONE);

     i--;

     j = n - 1;
     while (Perm[j] <= Perm[i])
          j--;

     /* Interchange Perm[i] and Perm[j]. */
     Temp = Perm[i];
     Perm[i] = Perm[j];
     Perm[j] = Temp;

     if (i < n - 2)
          for (i++, j = n - 1; i < j; i++, j--)
          {
               Temp = Perm[i];
               Perm[i] = Perm[j];
               Perm[j] = Temp;
          }

     return OK;
}


