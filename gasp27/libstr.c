/*****************************************************************/
/*   STRING-MANIPULATION ROUTINES                                */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--93.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

extern size_t  nPointers;

static char Buf[MAXTOK];

/*******************************+++*******************************/
/*                                                               */
/*   string    StrCatAlloc(const string s1, const string s2)     */
/*                                                               */
/*   Purpose:  Allocate space and concatenate strings s1 and s2. */
/*                                                               */
/*   Returns:  Pointer to the concatenated string.               */
/*                                                               */
/*   Version:  1991 July 5                                       */
/*                                                               */
/*****************************************************************/

string StrCatAlloc(const string s1, const string s2)
{
     string    CatStr = NULL;
     size_t    Len;

     Len = ((s1 != NULL) ? strlen(s1) : 0)
               + ((s2 != NULL) ? strlen(s2) : 0) + 1;
     CatStr = AllocChar(Len, NULL);
     *CatStr = NULL;

     if (s1 != NULL)
          strcpy(CatStr, s1);
     if (s2 != NULL)
          strcat(CatStr, s2);

     return CatStr;
}

/*******************************+++*******************************/
string StrPaste(size_t n, ...)
/*****************************************************************/
/*   Purpose:  Allocate space and paste together n strings;      */
/*             e.g., StrPaste(3, "x", ".", "y") produces "x.y".  */
/*                                                               */
/*   Returns:  Pointer to the concatenated string.               */
/*                                                               */
/*   Version:  1995 June 7                                       */
/*****************************************************************/
{
     string    CatStr, s;
     size_t    i, Len;
     va_list   Arg;

     va_start(Arg, n);

     for (Len = 0, CatStr = NULL, i = 0; i < n; i++)
     {
          if ( (s = va_arg(Arg, string)) != NULL)
          {
               Len += strlen(s);
               CatStr = AllocChar(Len + 1, CatStr);
               strcat(CatStr, s);
          }
     }

     va_end(Arg);

     return CatStr;
}

/*******************************+++*******************************/
int stricmp(const char *s, const char *t)
/*****************************************************************/
/*   Purpose:  Case insensitive (hence "i") version of strcmp.   */
/*                                                               */
/*   Returns:  < 0 if s < t;                                     */
/*               0 if s == t;                                    */
/*             > 0 if s > t;                                     */
/*                                                               */
/*   Version:  1992 March 2                                      */
/*****************************************************************/
{
     if (s != NULL && t != NULL)
     {
          for ( ; tolower(*s) == tolower(*t); s++, t++)
               if (*s == NULL)
                    return 0;

          return tolower(*s) - tolower(*t);
     }

     else if (s == NULL && t == NULL)
          return 0;

     else if (s == NULL)
          return -1;

     else
          return 1;
}

/*******************************+++*******************************/
/*                                                               */
/*   int strnicmp(const char *s, const char *t, size_t n)        */
/*                                                               */
/*   Purpose:  Case insensitive (hence "i") version of strncmp.  */
/*             Only first n characters of s and t compared.      */
/*                                                               */
/*   Returns:  < 0 if s < t;                                     */
/*               0 if s == t;                                    */
/*             > 0 if s > t;                                     */
/*                                                               */
/*   Version:  1991 May 8                                        */
/*                                                               */
/*****************************************************************/

int strnicmp(const char *s, const char *t, size_t n)
{
     size_t i;

     for (i = 0; i < n && tolower(*s) == tolower(*t); i++, s++, t++)
          if (*s == NULL)
               return 0;

     return (i < n) ? tolower(*s) - tolower(*t) : 0;
}

/*******************************+++*******************************/
/*                                                               */
/*   size_t    StrVecCmp(const string *s, const string *t,       */
/*                  size_t n)                                    */
/*                                                               */
/*   Purpose:  Pairwise comparison of s[0] with t[0],...,        */
/*             s[n-1] with t[n-1].                               */
/*                                                               */
/*   Returns:  INDEX_OK  if the two arrays are consistent;       */
/*             i         if strings s[i] and t[i] are the first  */
/*                       inconsistent pair.                      */
/*                                                               */
/*   Version:  1991 December 5                                   */
/*                                                               */
/*****************************************************************/

size_t StrVecCmp(const string *s, const string *t, size_t n)
{
     size_t    i;

     if (s == NULL || t == NULL)
          return INDEX_OK;

     for (i = 0; i < n; i++)
          if (s[i] != NULL && t[i] != NULL &&
                    stricmp(s[i], t[i]) != 0)
               return i;

     return INDEX_OK;
}

/*******************************+++*******************************/
/*                                                               */
/*   string    StrDup(const string s)                            */
/*                                                               */
/*   Purpose:  Allocate space, and duplicate string.             */
/*                                                               */
/*   Returns:  Pointer to duplicate.                             */
/*                                                               */
/*   Version:  1991 June 11                                      */
/*                                                               */
/*****************************************************************/

string StrDup(const string s)
{
     string DupStr = NULL;

     if (s == NULL)
          return NULL;

     else
     {
          DupStr = AllocChar(strlen(s) + 1, NULL);
          strcpy(DupStr, s);

          return DupStr;
     }
}

/*******************************+++*******************************/
string StrReplace(const string NewStr, string Target)
/*****************************************************************/
/*   Purpose:  Replace Target with NewStr.                       */
/*                                                               */
/*   Returns:  New pointer to Target.                            */
/*                                                               */
/*   Version:  1994 February 16                                  */
/*****************************************************************/
{
     if (Target != NULL)
          AllocFree(Target);
     Target = StrDup(NewStr);

     return Target;
}

/*******************************+++*******************************/
int StrToInt(const string s, int *i)
/*****************************************************************/
/*   Purpose:  Convert a string to an integer.                   */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*                                                               */
/*   Version:  1992 March 27                                     */
/*****************************************************************/
{
     char *EndPtr = NULL;
     int  ErrNum;
     long l;

     l = strtol(s, &EndPtr, 10);

     if (*EndPtr == NULL && l >= (long) INT_MIN &&
               l <= (long) INT_MAX)
     {
          *i = (int) l;
          ErrNum = OK;
     }

     else if (stricmp(s, NOT_AVAIL) == 0)
     {
          *i = NA_INT;
          ErrNum = OK;
     }

     else
     {
          *i = 0;
          ErrNum = INPUT_ERR;
     }

     return ErrNum;
}

/*******************************+++*******************************/
int StrToReal(const string s, real *r)
/*****************************************************************/
/*   Purpose:  Convert a string to a real.                       */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*                                                               */
/*   Version:  1992 March 4                                      */
/*****************************************************************/
{
     char *EndPtr = NULL;
     int  ErrNum;

     *r = strtod(s, &EndPtr);

     if (*EndPtr == NULL)
          ErrNum = OK;

     else if (stricmp(s, NOT_AVAIL) == 0)
     {
          *r = NA_REAL;
          ErrNum = OK;
     }

     else
     {
          *r = 0.0;
          ErrNum = INPUT_ERR;
     }

     return ErrNum;
}

/*******************************+++*******************************/
int StrToSize_t(const string s, size_t *z)
/*****************************************************************/
/*   Purpose:  Convert a string to a size_t.                     */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*                                                               */
/*   Version:  1992 March 9                                      */
/*****************************************************************/
{
     char           *EndPtr = NULL;
     int            ErrNum;
     unsigned long  u;

     u = strtoul(s, &EndPtr, 10);

     if (*EndPtr == NULL && u <= UINT_MAX)
     {
          *z = (size_t) u;
          ErrNum = OK;
     }

     else if (stricmp(s, NOT_AVAIL) == 0)
     {
          *z = NA_SIZE_T;
          ErrNum = OK;
     }

     else
     {
          *z = 0;
          ErrNum = INPUT_ERR;
     }

     return ErrNum;
}

/*******************************+++*******************************/
string StrFromInt(int i)
/*****************************************************************/
/*   Purpose:  Return a string corresponding to an integer.      */
/*                                                               */
/*   Comment:  Calling routine should duplicate the string       */
/*             before the next call of StrFromInt,               */
/*             StrFromSize_t, or StrFromReal.                    */
/*                                                               */
/*   Version:  1992 March 26                                     */
/*****************************************************************/
{
     if (i == NA_INT)
          strcpy(Buf, NOT_AVAIL);
     else
          sprintf(Buf, "%d", i);

     return Buf;
}

/*******************************+++*******************************/
string StrFromReal(real r, const string Flags, int Precision,
          char Conversion)
/*****************************************************************/
/* Purpose:    Return a string corresponding to a real.          */
/*             Flags are as defined for sprintf.                 */
/*             Precision is as defined for sprintf ( < 0 for     */
/*             default).                                         */
/*             Conversion can be 'g', 'e', or 'f'.               */
/*                                                               */
/* Comment:    Calling routine should duplicate the string       */
/*             before the next call of StrFromInt,               */
/*             StrFromSize_t, or StrFromReal.                    */
/*                                                               */
/* 1993.11.14:                                                   */
/* 1999.07.01: REAL_MAX gives INFINITY_TXT.                      */
/*****************************************************************/
{
     char Format[10];

     if (Precision < 0)
          Precision = PRECISION;   /* Default. */

     if (r == NA_REAL)
          strcpy(Buf, NOT_AVAIL);
     else if (r == REAL_MAX)     
          strcpy(Buf, INFINITY_TXT);
     else
     {
          /* Format = "%" + Flags + ".*" + Conversion character. */
          strcpy(Format, "%");
          strcat(Format, Flags);
          strcat(Format, ".*");
          strncat(Format, &Conversion, 1);

          sprintf(Buf, Format, Precision, r);
     }

     return Buf;
}

/*******************************+++*******************************/
string StrFromSize_t(size_t z)
/*****************************************************************/
/*   Purpose:  Return a string corresponding to a size_t.        */
/*                                                               */
/*   Comment:  Calling routine should duplicate the string       */
/*             before the next call of StrFromInt,               */
/*             StrFromSize_t, or StrFromReal.                    */
/*                                                               */
/*   Version:  1992 March 26                                     */
/*****************************************************************/
{
     if (z == NA_SIZE_T)
          strcpy(Buf, NOT_AVAIL);
     else
          sprintf(Buf, "%u", z);

     return Buf;
}

/*******************************+++*******************************/
/*                                                               */
/*   size_t    StrIndex(const string Target,                     */
/*                  const string *StrArray, size_t NumStrings)   */
/*                                                               */
/*   Purpose:  Find first i such that Target matches StrArray[i].*/
/*                                                               */
/*   Returns:  i          if StrArray[i] matches Target;         */
/*             INDEX_ERR  otherwise.                             */
/*                                                               */
/*   Version:  1991 June 19                                      */
/*                                                               */
/*****************************************************************/

size_t StrIndex(const string Target, const string *StrArray,
          size_t NumStrings)
{
     size_t    i;

     for (i = 0; i < NumStrings; i++)
          if (stricmp(Target, StrArray[i]) == 0)
               return i;

     return INDEX_ERR;
}

/*******************************+++*******************************/
size_t StrNumberOf(const string *StrArray)
/*****************************************************************/
/*   Purpose:  Return number of strings in StrArray, not         */
/*             counting the terminating NULL string.             */
/*                                                               */
/*   Version:  1992 May 28                                       */
/*****************************************************************/
{
     size_t    Num;

     Num = 0;

     if (StrArray != NULL)
          while (StrArray[Num] != NULL)
               Num++;

     return Num;
}

/*******************************+++*******************************/
/*                                                               */
/*   int       StrBrackets(string Token, string *Sub,            */
/*                  string *NextToken)                           */
/*                                                               */
/*   Purpose:  Parse a token for the first pair of opening and   */
/*             closing brackets.                                 */
/*                                                               */
/*   Args:     Token     Input:  the token to be parsed;         */
/*                       Output: the token will be terminated    */
/*                               to exclude the first '['.       */
/*             Sub       Output: a null-terminated string, the   */
/*                               subscript inside '[' and ']'.   */
/*             NextToken Output: the string after the ']'.       */
/*                                                               */
/*   Returns:  OK or INPUT_ERR.                                  */
/*                                                               */
/*   Version:  1993 June 4                                       */
/*                                                               */
/*****************************************************************/

int StrBrackets(string Token, string *Sub, string *NextToken)
{
     char *Close, *Open;

     Open  = strchr(Token, '[');
     Close = strchr(Token, ']');

     *Sub = NULL;
     *NextToken = NULL;

     if (Open == NULL && Close == NULL)
          /* No brackets. */
          return OK;
     else if (Open != NULL && Close != NULL &&  Open + 1 < Close)
     {
          /* Legal brackets. */
          *Open = NULL;
          *Sub = Open + 1;
          *Close = NULL;
          *NextToken = Close + 1;
          return OK;
     }
     else
          return INPUT_ERR;
}

