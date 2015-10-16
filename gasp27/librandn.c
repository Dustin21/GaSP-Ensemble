/*****************************************************************/
/*                                                               */
/*   PSEUDO-RANDOM-NUMBER GENERATORS                             */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--91.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

/* Default initial seeds for the 3 components of the  */
/* congruential generator (can be re-set by RandInit) */
static int xcomp = 1, ycomp = 15000, zcomp = 30000;

/*******************************+++*******************************/
/*                                                               */
/*   int RandInit(int _xcomp, int _ycomp, int _zcomp)            */
/*                                                               */
/*   Purpose:  Initialise sequence of pseudo-random numbers      */
/*             by setting the 3 components of the generator to   */
/*             _xcomp, _ycomp, and _zcomp (all 3 must be in the  */
/*             range 1 to 30000).                                */
/*                                                               */
/*   Returns:  RANGE_ERR if _xcomp, _ycomp, or _zcomp is out of  */
/*                       range;                                  */
/*             OK        otherwise.                              */
/*                                                               */
/*   Version:  1991 May 22                                       */
/*                                                               */
/*****************************************************************/

int RandInit(int _xcomp, int _ycomp, int _zcomp)
{
     if (_xcomp < 1 || _xcomp > 30000 ||
              _ycomp < 1 || _ycomp > 30000 ||
              _zcomp < 1 || _zcomp > 30000)
          return RANGE_ERR;
     else
     {
          xcomp = _xcomp;
          ycomp = _ycomp;
          zcomp = _zcomp;
          return OK;
     }
}
/*******************************+++*******************************/
/*                                                               */
/*   real      RandUnif(void)                                    */
/*                                                               */
/*   Purpose:  Return a uniform (0, 1) random number.            */
/*                                                               */
/*   Comments: Adapted from Algorithm AS 183, Applied            */
/*             Statistics 1982.                                  */
/*             Use RandInit to initialise (otherwise, same       */
/*             numbers generated every run).                     */
/*             Integer arithmetic up to 30323 is required.       */
/*                                                               */
/*   Version:  1991 May 22                                       */
/*                                                               */
/*****************************************************************/
 
real RandUnif(void)
{
     real u;
 
     xcomp = 171 * (xcomp % 177) -  2 * (xcomp / 177);
     ycomp = 172 * (ycomp % 176) - 35 * (ycomp / 176);
     zcomp = 170 * (zcomp % 178) - 63 * (zcomp / 178);
 
     if (xcomp < 0)
          xcomp = xcomp + 30269;
     if (ycomp < 0)
          ycomp = ycomp + 30307;
     if (zcomp < 0)
          zcomp = zcomp + 30323;
 
     u = xcomp / 30269.0 + ycomp / 30307.0 + zcomp / 30323.0;
     while (u > 1.0)
          u = u - 1.0;
     return u;
}
