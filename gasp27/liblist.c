/*****************************************************************/
/*                                                               */
/*   ROUTINES TO MAINTAIN A SINGLY-LINKED LIST                   */
/*                                                               */
/*   Copyright (c) William J. Welch 1991--92.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

/*******************************+++*******************************/
List *ListAdd(const string Name, int Type, void *Data, List *L)
/*****************************************************************/
/*   Purpose:  Add a new member to the beginning of a list.      */
/*                                                               */
/*   Returns:  Pointer to the new list head.                     */
/*                                                               */
/*   Comment:  The calling routine must allocate void *Data for  */
/*             the new member.                                   */
/*                                                               */
/*   Version:  1992 April 24                                     */
/*****************************************************************/
{
     List *NewMember;

     NewMember = AllocGeneric(1, sizeof(List), NULL);

     NewMember->Name = StrDup(Name);
     NewMember->Type = Type;
     NewMember->Data = Data;

     /* Put at the beginning of the list. */
     NewMember->Next = L;

     return NewMember;
}

/*******************************+++*******************************/
void *ListData(const string Name, int Type, const List *L)
/*****************************************************************/
/*   Purpose:  Return a pointer to the data for a named member   */
/*             of a list.                                        */
/*                                                               */
/*   Returns:  Pointer if Name is found and Type is correct;     */
/*             NULL    otherwise.                                */
/*                                                               */
/*   Version:  1992 April 29                                     */
/*****************************************************************/
{
     List *Member;

     /* Look for Name in the list. */
     for (Member = L; Member != NULL; Member = Member->Next)
          if (stricmp(Name, Member->Name) == 0)
               break;   /* Found it. */

     if (Member != NULL &&
               (Member->Type == Type || Type == UNKNOWN))
          return Member->Data;
     else
          return NULL;
}
