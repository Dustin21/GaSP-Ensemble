/*****************************************************************/
/*   ROUTINES TO MANAGE LINEAR MODEL FUNCTIONS (E.G., X^2)       */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--91.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"
#include "model.h"

/*******************************+++*******************************/
int ModFnParse(string Comp, int *Fn)
/*****************************************************************/
/*   Purpose:  Parse Comp for a linear-model function.           */
/*             On return, *Fn is the function found, and Comp    */
/*             has all function-related characters removed.      */
/*                                                               */
/*   Returns:  INPUT_ERR or OK.                                  */
/*                                                               */
/*   Version:  1992 March 2                                      */
/*****************************************************************/
{
     char *Carat;
     int  ErrNum;

     ErrNum = OK;
     Carat = strchr(Comp, '^');
     if (Carat == NULL)
          *Fn = 0;
     else if (StrToInt(Carat + 1, Fn) != OK || *Fn < 2)
     {
          *Fn = 0;
          Error("The exponent in model component \"%s\" should "
                    "be an integer >= 2.\n", Comp);
          ErrNum = INPUT_ERR;
     }
     else
          *Carat = NULL;

     return ErrNum;
}

/*******************************+++*******************************/
real ModFn(real x, int fn)
/*****************************************************************/
/*   Purpose:  Return linear-model component when function fn is */
/*             applied to x.                                     */
/*                                                               */
/*   Version:  1991 May 22                                       */
/*****************************************************************/
{
     if (fn >= 2)
          return pow(x, (double) fn);
     else
          return x;
}

