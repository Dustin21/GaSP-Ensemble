/*****************************************************************/
/*   LOW-LEVEL ROUTINES FOR INPUT.                               */
/*                                                               */
/*   Copyright (c) William J. Welch 1996.                        */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

/*******************************+++*******************************/
int InputSize_t(FILE *FilePtr, size_t *z)
/*****************************************************************/
/*   Purpose:  Read size_t z from file FilePtr.                  */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*                                                               */
/*   Version:  1996.04.04                                        */
/*****************************************************************/
{
     if (fscanf(FilePtr, "%d", z) != 1)
     {
          Error(CANNOT_READ_SIZE_T);
          return INPUT_ERR;
     }
     else
          return OK;
}

/*******************************+++*******************************/
int InputReal(FILE *FilePtr, real *r)
/*****************************************************************/
/*   Purpose:  Read real r from file FilePtr.                    */
/*                                                               */
/*   Returns:  OK or an error condition.                         */
/*                                                               */
/*   Version:  1996.04.04                                        */
/*****************************************************************/
{
     if (fscanf(FilePtr, "%lg", r) != 1)
     {
          Error(CANNOT_READ_REAL);
          return INPUT_ERR;
     }
     else
          return OK;
}

