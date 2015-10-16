/*****************************************************************/
/*   ROUTINES TO MANAGE VECTOR OBJECTS                           */
/*                                                               */
/*   Copyright (c) William J. Welch 1994.                        */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

/*******************************+++*******************************/
size_t VecTempl(const vector *v, const template *T)
/*****************************************************************/
/*   Purpose:  Check a vector against a template.                */
/*                                                               */
/*   Returns:  j        if element j of v is illegal;            */
/*             INDEX_OK otherwise.                               */
/*                                                               */
/*   96.01.18: CodeBug parameters changed or CodeBug replaced by */
/*             CodeCheck.                                        */
/*   96.01.22: CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     size_t    BadIndex, Len;

     BadIndex = OK;

     CodeCheck(VecType(v) == TemplType(T));

     Len = VecLength(v);

     switch (VecType(v))
     {
          case INTEGER:
               BadIndex = TemplIntCheck(T, VecInts(v), Len);
               break;

          case REAL:
               BadIndex = TemplRealCheck(T, VecReals(v), Len);
               break;

          case SIZE_T:
               BadIndex = TemplSize_tCheck(T, VecSize_ts(v), Len);
               break;

          case STRING:
               BadIndex = TemplStrCheck(T, VecStrs(v), Len);
               break;

          default:
               CodeBug("Illegal type");
     }

     return BadIndex;
}

/*******************************+++*******************************/
void VecToStr(const vector *v, string *s)
/*****************************************************************/
/*   Purpose:  Convert a vector to an array of strings.          */
/*                                                               */
/*   Comment:  Calling routine should allocate s to be of        */
/*             sufficient length, but the strings themselves are */
/*             StrDup'd here.                                    */
/*                                                               */
/*   96.01.18: CodeBug parameters changed.                       */
/*   96.01.22: CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     size_t    j, Len;

     Len = VecLength(v);

     switch (VecType(v))
     {
          case INTEGER:
               for (j = 0; j < Len; j++)
                    s[j] = StrDup(StrFromInt(VecInt(v, j)));
               break;

          case REAL:
               for (j = 0; j < Len; j++)
                    s[j] = StrDup(StrFromReal(VecReal(v, j), "",
                              PRECISION, 'g'));
               break;

          case SIZE_T:
               for (j = 0; j < Len; j++)
                    s[j] = StrDup(StrFromSize_t(VecSize_t(v, j)));
               break;

          case STRING:
               for (j = 0; j < Len; j++)
                    s[j] = StrDup(VecStr(v, j));
               break;

          default:
               CodeBug("Illegal type");
     }
}

/*******************************+++*******************************/
size_t VecFromStr(const string *s, vector *v)
/*****************************************************************/
/*   Purpose:  Convert an array of strings to a vector.          */
/*                                                               */
/*   Returns:  j        if s[j] is illegal for the vector type;  */
/*             INDEX_OK otherwise.                               */
/*                                                               */
/*   96.01.18: CodeBug parameters changed.                       */
/*   96.01.22: CodeBug parameters changed.                       */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     size_t    j, Len;

     Len = VecLength(v);

     switch (VecType(v))
     {
          case INTEGER:
               for (j = 0; j < Len; j++)
                    if (StrToInt(s[j], &VecInt(v, j)) != OK)
                         break;
               break;

          case REAL:
               for (j = 0; j < Len; j++)
                    if (StrToReal(s[j], &VecReal(v, j)) != OK)
                         break;
               break;

          case SIZE_T:
               for (j = 0; j < Len; j++)
                    if (StrToSize_t(s[j], &VecSize_t(v, j)) != OK)
                         break;
               break;

          case STRING:
               for (j = 0; j < Len; j++)
                    VecPutStr(v, j, s[j]);
               break;

          default:
               CodeBug("Illegal type");
     }

     return (j == Len) ? INDEX_OK : j;
}
