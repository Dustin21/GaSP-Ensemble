/*****************************************************************/
/*   ROUTINES FOR PROBABILITY DISTRIBUTIONS, ETC.                */
/*                                                               */
/*   Copyright (c) William J. Welch 1994--96.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "define.h"
#include "implem.h"
#include "lib.h"


#define AA1     0.398942280444
#define AA2     0.399903438504
#define AA3     5.75885480458
#define AA4    29.8213557808
#define AA5     2.62433121679
#define AA6    48.6959930692
#define AA7     5.92885724438

#define BB1     0.398942280385
#define BB2     3.8052e-8
#define BB3     1.00000615302
#define BB4     3.98064794e-4
#define BB5     1.98615381364
#define BB6     0.151679116635
#define BB7     5.29330324926
#define BB8     4.8385912808
#define BB9    15.1508972451
#define BB10    0.742380924027
#define BB11   30.789933034
#define BB12    3.99019417011

/* z at which lower tail area becomes 1.0 */
/* = (decimal digits + 9) / 3.            */
#define LT1     8.0

/* z at which upper tail area becomes 0        */
/* = sqrt(-2 * (ln(smallest real) + 1)) - 0.3. */
#define UT0    37.32

/*******************************+++*******************************/
real CDFNorm(real z, boolean Lower)
/*****************************************************************/
/*   Purpose:  Returns tail area of the standard normal          */
/*             distribution, lower tail if Lower is TRUE, upper  */
/*             tail otherwise.                                   */
/*                                                               */
/*   Version:  1996 January 15                                   */
/*****************************************************************/
{
     real p, y;

     /* Always work with positive z. */
     if (z < 0.0)
     {
          Lower = !Lower;
          z = -z;
     }

     /* Compute upper tail probability. */
     if (z <= LT1 || (!Lower && z <= UT0))
     {
          y = 0.5 * z * z;

          if (z <= 1.28)
          {
               p = 0.5 - z * (AA1 - AA2 * y
                         / ( y + AA3 - AA4 / (y + AA5 + AA6
                         / (y + AA7))));
          }
          else
               p = BB1 * exp(-y)  / (z - BB2
                           + BB3  / (z + BB4
                           + BB5  / (z - BB6
                           + BB7  / (z + BB8
                           - BB9  / (z + BB10
                           + BB11 / (z + BB12))))));

     }
     else
          p = 0.0;

     if (Lower)
          p = 1.0 - p;

     return p;
}


#define A0      2.50662823884
#define A1    -18.61500062529
#define A2     41.39119773534
#define A3    -25.44106049637
#define B1     -8.47351093090
#define B2     23.08336743743
#define B3    -21.06224101826
#define B4      3.13082909833

/* Sum:       143.70383558076 */

#define C0     -2.78718931138
#define C1     -2.29796479134
#define C2      4.85014127135
#define C3      2.32121276858
#define D1      3.54388924762
#define D2      1.63706781897

/* Sum:        17.43746520924 */

#define SPLIT 0.42

/*******************************+++*******************************/
real CDFInvNorm(real p)
/*****************************************************************/
/*   Purpose:  Returns inverse of the normal distribution for    */
/*             probability p; i.e., the integral from -infinity  */
/*             to the returned value is p.                       */
/*                                                               */
/*   Returns:  Normal quantile if 0 <= p <= 1;                   */
/*             -REAL_MAX       otherwise.                        */
/*                                                               */
/*   Version:  1996 January 12                                   */
/*****************************************************************/
{
     real InvN, q, r;

     if (fabs(q = p - 0.5) <= SPLIT)
     {
          r = q * q;
          InvN = q * (((A3 * r + A2) * r + A1) * r + A0) /
                    ((((B4 * r + B3) * r + B2) * r + B1) * r + 1.0);
     }
     else
     {
          r = (q > 0.0) ? 1.0 - p : p;

          if (r > 0.0)
          {
               r = sqrt(-log(r));
               InvN = (((C3 * r + C2) * r + C1) * r + C0) /
                       ((D2 * r + D1) * r + 1.0);
               if (q < 0.0)
                    InvN = -InvN;
           }
           else
               InvN = -REAL_MAX;
     }

     return InvN;
}

#define RECIP_ROOT_2_PI  0.3989422804014327

/*******************************+++*******************************/
real PDFNorm(real z)
/*****************************************************************/
/*   Purpose:  Returns p.d.f. of standard normal distribution.   */
/*                                                               */
/*   Version:  1996 January 15                                   */
/*****************************************************************/
{
     return RECIP_ROOT_2_PI * exp(-0.5 * z * z);
}

/*******************************+++*******************************/
real Cor(real *x, real *y, size_t n)
/*****************************************************************/
/*   Purpose:  Returns correlation between x[0], ..., x[n-1] and */
/*             y[0], ..., y[n-1].                                */
/*                                                               */
/*   Version:  1994 March 2                                      */
/*****************************************************************/
{
     real c, MeanX, MeanY, sX, sY;

     if (n == 0)
          return 0.0;

     /* Correct x and y for their means. */
     /* (Can this be avoided?)           */
     MeanX = VecSum(x, n) / n;
     MeanY = VecSum(y, n) / n;
     VecAddScalar(-MeanX, n, x);
     VecAddScalar(-MeanY, n, y);

     sX = sqrt(VecSS(x, n));
     sY = sqrt(VecSS(y, n));
     if (sX > 0.0 && sY > 0.0)
          c = DotProd(x, y, n) / sX  / sY;
     else
          c = 0.0;

     /* Restore x and y. */
     VecAddScalar(MeanX, n, x);
     VecAddScalar(MeanY, n, y);

     return c;
}

