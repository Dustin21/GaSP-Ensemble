/*****************************************************************/
/*   MEMORY-ALLOCATION ROUTINES                                  */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--95.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "implem.h"
#include "define.h"
#include "lib.h"

/* Pointers that have been allocated by AllocGeneric()  */
/* and not yet freed by AllocFree().  These are kept to */
/* facilitate debugging.                                */
size_t nPointers = 0;
void   **Pointer = NULL;

/*******************************+++*******************************/
/*                                                               */
/*   char      *AllocChar(size_t n, char *p) etc.                */
/*                                                               */
/*   Purpose:  Allocate or reallocate a vector of n items of     */
/*             the following types:                              */
/*                                                               */
/*             Type      Routine                                 */
/*                                                               */
/*             char      AllocChar                               */
/*             string    AllocStr                                */
/*             string *  AllocPtrStr                             */
/*             real      AllocReal                               */
/*             real *    AllocPtrReal                            */
/*             int       AllocInt                                */
/*             int *     AllocPtrInt                             */
/*             size_t    AllocSize_t                             */
/*             size_t *  AllocPtrSize_t                          */
/*                                                               */
/*   Args:     n         Input     Number of items.              */
/*             p         Input     Pointer to one of the above   */
/*                                 types (NULL for first         */
/*                                 allocation).                  */
/*                       Output    Points to the n items.        */
/*                                                               */
/*   Return:   Pointer to the n items.                           */
/*                                                               */
/*   Comment:  exit(1) is called if there is insufficient memory.*/
/*                                                               */
/*   Version:  1991 December 3                                   */
/*                                                               */
/*****************************************************************/

char *AllocChar(size_t n, char *p)
{
     return (char *) AllocGeneric(n, sizeof(char), p);
}

string *AllocStr(size_t n, string *p)
{
     return (string *) AllocGeneric(n, sizeof(string), p);
}

string **AllocPtrStr(size_t n, string **p)
{
     return (string **) AllocGeneric(n, sizeof(string *), p);
}

real *AllocReal(size_t n, real *p)
{
     return (real *) AllocGeneric(n, sizeof(real), p);
}

real **AllocPtrReal(size_t n, real **p)
{
     return (real **) AllocGeneric(n, sizeof(real *), p);
}

int *AllocInt(size_t n, int *p)
{
     return (int *) AllocGeneric(n, sizeof(int), p);
}

int **AllocPtrInt(size_t n, int **p)
{
     return (int **) AllocGeneric(n, sizeof(int *), p);
}

size_t *AllocSize_t(size_t n, size_t *p)
{
     return (size_t *) AllocGeneric(n, sizeof(size_t), p);
}

size_t **AllocPtrSize_t(size_t n, size_t **p)
{
     return (size_t **) AllocGeneric(n, sizeof(size_t *), p);
}

/*******************************+++*******************************/
void *AllocGeneric(size_t n, size_t Size, void *p)
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
/*   Return:   Pointer to the n items (NULL if n == 0).          */
/*                                                               */
/*   Comment:  exit(1) is called if there is insufficient memory.*/
/*                                                               */
/*   Version:  1995 September 22                                 */
/*****************************************************************/
{
     size_t    i;

     if (n > 0 && p == NULL)
     {
          p = calloc(n, Size);
          Pointer = (void **) realloc(Pointer,
                    (++nPointers) * sizeof(void *));
          if (Pointer != NULL)
               Pointer[nPointers-1] = p;
     }
     else if (n > 0 && p != NULL)
     {
          i = AllocFindPtr(p);
          p = realloc(p, n * Size);
          Pointer[i] = p;
     }
     else if (n == 0 && p != NULL)
     {
          AllocFree(p);
          p = NULL;
     }

     /* No action required if n == 0 && p == NULL. */

     if ( (p == NULL && n > 0) ||
               (nPointers > 0 && Pointer == NULL) )
     {
          Fatal("Insufficient memory.\n");
          exit(1);
     }

     return p;
}

/*******************************+++*******************************/
size_t AllocFindPtr(void *p)
/*****************************************************************/
/*   Purpose:  Return i such that Pointer[i] = p.                */
/*                                                               */
/*   Version:  1995 September 22                                 */
/*****************************************************************/
{
     size_t    i, ii;

     /* Start search from the end of Pointer.       */
     /* (A frequently reallocated pointer is likely */
     /* to be near the end.)                        */
     for (ii = 0; ii < nPointers; ii++)
     {
          i = nPointers - 1 - ii;
          if (Pointer[i] == p)
               break;
     }

     CodeCheck(ii < nPointers);

     return i;
}

/*******************************+++*******************************/
void AllocFree(void *p)
/*****************************************************************/
/*   Purpose:  Free p.                                           */
/*                                                               */
/*   Version:  1995 September 22                                 */
/*****************************************************************/
{
     size_t    i, j;

     if (p != NULL)
     {
          i = AllocFindPtr(p);

          for (j = i; j < nPointers - 1; j++)
               Pointer[j] = Pointer[j+1];

          nPointers--;

          free(p);
     }
}

/*******************************+++*******************************/
/*                                                               */
/*   string    *AllocStrFree(size_t OldLen, size_t NewLen,       */
/*                  string *s)                                   */
/*                                                               */
/*   Purpose:  Re-allocate a vector of strings, freeing any      */
/*             "excess" strings, and initializing new strings to */
/*             NULL.                                             */
/*                                                               */
/*   Args:     OldLen    Old length.                             */
/*             NewLen    New length.                             */
/*             s         Input:  Pointer to the strings.         */
/*                       Output: Re-allocated pointer.           */
/*                                                               */
/*   Returns:  Re-allocated pointer.                             */
/*                                                               */
/*   Version:  1992 April 16                                     */
/*                                                               */
/*****************************************************************/

string *AllocStrFree(size_t OldLen, size_t NewLen, string *s)
{
     size_t    i;

     for (i = NewLen; i < OldLen; i++)
          AllocFree(s[i]);

     s = AllocStr(NewLen, s);

     for (i = OldLen; i < NewLen; i++)
          s[i] = NULL;

     return s;
}

/*******************************+++*******************************/
size_t AllocMax(size_t Size)
/*****************************************************************/
/*   Purpose:  Return maximum number of objects of a given size  */
/*             that could be allocated.                          */
/*                                                               */
/*   Version:  1992 March 16                                     */
/*****************************************************************/
{
     void      *p;
     size_t    Failure, Try, Success;

     if ( (p = calloc(SIZE_T_MAX, Size)) != NULL)
     {
          Success = SIZE_T_MAX;
          free(p);
     }

     else
     {
          Success = 0;
          Failure = SIZE_T_MAX;

          while (Failure - Success > 1)
          {
               Try = Success + (Failure - Success) / 2;
               if ( (p = calloc(Try, Size)) != NULL)
               {
                    Success = Try;
                    free(p);
               }
               else
                    Failure = Try;
          }
     }

     return Success;
}
