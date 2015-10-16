/*****************************************************************/
/*   ROUTINES TO MANAGE DATA TEMPLATES                           */
/*                                                               */
/*   Copyright (c) William J. Welch 1992--95.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

/*******************************+++*******************************/
List *TemplAlloc(const string Name, int Type, List *L)
/*****************************************************************/
/*   Purpose:  Allocate a template of the specified name and     */
/*             type by adding it to the list L.                  */
/*                                                               */
/*   Returns:  Pointer to the updated list.                      */
/*                                                               */
/*   1995 February 21                                            */
/*   2006.02.18: Duplicate of Name changed, not Name (bug fix)   */
/*****************************************************************/
{
     char      *DotPtr;
     template  *T;

     if ( (T = TemplPtr(Name, L)) == NULL)
     {
          /* New template. */

          T = (template *) AllocGeneric(1, sizeof(template), NULL);

          T->Name = StrDup(Name);
          /* If Name includes an extension, remove extension. */
          if ( (DotPtr = strrchr(T->Name, '.')) != NULL)
               *DotPtr = NULL;

          T->NumLegalStr = 0;
          T->LegalStr    = NULL;

          L = ListAdd(T->Name, UNKNOWN, T, L);
     }

     T->Type = Type;

     return L;
}

/******************************+++********************************/
template *TemplPtr(const string Name, const List *L)
/*****************************************************************/
/*   Purpose:  Return template corresponding to Name or NULL.    */
/*                                                               */
/*   1994 December 27                                            */
/*   2006.02.18: Duplicate of Name changed, not Name (bug fix)   */
/*****************************************************************/
{
     char      *DotPtr;
     string    NameDup;
     template  *T;

     NameDup = StrDup(Name);

     /* If Name includes an extension, remove extension. */
     if ( (DotPtr = strrchr(NameDup, '.')) != NULL)
          *DotPtr = NULL;

     T = (template *) ListData(NameDup, UNKNOWN, L);

     AllocFree(NameDup);

     return T;
}

/******************************+++********************************/
List *TemplIntAlloc(const string Name, int Min, int Max, List *L)
/*****************************************************************/
/*   Purpose:  Allocate an INTEGER template.                     */
/*                                                               */
/*   Returns:  Pointer to the updated list.                      */
/*                                                               */
/*   Comment:  TemplRealAlloc, TemplSize_tAlloc, and             */
/*             TemplStrAlloc are similar.                        */
/*                                                               */
/*   96.01.18: CodeBug replaced by CodeCheck.                    */
/*                                                               */
/*   Version:  1996 January 18                                   */
/*****************************************************************/
{
     template  *T;

     L = TemplAlloc(Name, INTEGER, L);

     T = TemplPtr(Name, L);
     CodeCheck(T != NULL);

     T->MinInt = Min;
     T->MaxInt = Max;

     return L;
}

List *TemplRealAlloc(const string Name, real Min, real Max,
          List *L)
{
     template  *T;

     L = TemplAlloc(Name, REAL, L);

     T = TemplPtr(Name, L);
     CodeCheck(T != NULL);

     T->MinReal = Min;
     T->MaxReal = Max;

     return L;
}

List *TemplSize_tAlloc(const string Name, size_t Min, size_t Max,
          List *L)
{
     template  *T;

     L = TemplAlloc(Name, SIZE_T, L);

     T = TemplPtr(Name, L);
     CodeCheck(T != NULL);

     T->MinSize_t = Min;
     T->MaxSize_t = Max;

     return L;
}

List *TemplStrAlloc(const string Name, const string *LegalStr,
          size_t NumLegalStr, List *L)
{
     size_t    i;
     template  *T;

     L = TemplAlloc(Name, STRING, L);

     T = TemplPtr(Name, L);
     CodeCheck(T != NULL);

     if (NumLegalStr > 0)
     {
          T->LegalStr = AllocStr(NumLegalStr, NULL);
          for (i = 0; i < NumLegalStr; i++)
               T->LegalStr[i] = StrDup(LegalStr[i]);
     }

     T->NumLegalStr = NumLegalStr;

     return L;
}

/*******************************+++*******************************/
size_t TemplIntCheck(const template *T, const int *i, size_t n)
/*****************************************************************/
/*   Purpose:  Check an array of values against a template.      */
/*                                                               */
/*   Returns:  j        if the value with index j is illegal;    */
/*             INDEX_OK otherwise.                               */
/*                                                               */
/*   Comment:  TemplRealCheck, TemplSize_tCheck, and             */
/*             TemplStrCheck are similar.                        */
/*                                                               */
/*   Version:  1994 September 8                                  */
/*****************************************************************/
{
     size_t    j;

     if (T == NULL)
          return INDEX_OK;

     for (j = 0; j < n; j++)
          if (i[j] < T->MinInt || i[j] > T->MaxInt)
               return j;

     return INDEX_OK;
}

size_t TemplRealCheck(const template *T, const real *r, size_t n)
{
     size_t    j;

     if (T == NULL)
          return INDEX_OK;

     for (j = 0; j < n; j++)
          if (r[j] < T->MinReal || r[j] > T->MaxReal)
               return j;

     return INDEX_OK;
}

size_t TemplSize_tCheck(const template *T, const size_t *z,
          size_t n)
{
     size_t    j;

     if (T == NULL)
          return INDEX_OK;

     for (j = 0; j < n; j++)
          if (z[j] < T->MinSize_t || z[j] > T->MaxSize_t)
               return j;

     return INDEX_OK;
}

size_t TemplStrCheck(const template *T, const string *s, size_t n)
{
     size_t    j;

     if (T == NULL || T->LegalStr == NULL)
          return INDEX_OK;

     for (j = 0; j < n; j++)
          if (TemplStrIndex(s[j], T) == INDEX_ERR)
               return j;

     return INDEX_OK;
}

/*******************************+++*******************************/
void TemplError(const template *T)
/*****************************************************************/
/*   Purpose:  Output an error message.                          */
/*                                                               */
/*   Version:  1992 February 14                                  */
/*****************************************************************/
{
     size_t    j;

     if (T->Type == INTEGER)
          Error("%s should take integer values between "
                    "%d and %d.\n", T->Name, T->MinInt, T->MaxInt);

     else if (T->Type == REAL)
          Error("%s should take numerical values between "
                    "%g and %g.\n", T->Name, T->MinReal,
                    T->MaxReal);

     else if (T->Type == SIZE_T)
          Error("%s should take integer values between "
                    "%u and %u.\n", T->Name, T->MinSize_t,
                    T->MaxSize_t);

     else if (T->Type == STRING)
     {
          Error("%s should be equal to one of the strings: \n",
                    T->Name);
          for (j = 0; j < T->NumLegalStr - 1; j++)
               Output("%s, ", T->LegalStr[j]);
          Output("%s.\n", T->LegalStr[j]);
     }
}
