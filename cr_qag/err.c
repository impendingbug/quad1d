/**
 *  \file  err.c
 *  \brief  An almost-exact copy of integration/err.c from GSL source code.
 *
 *  \author  Sandy Pratama
 *
 *  \internal
 *       Created:  04-08-14
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 2014 Sandy Pratama
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#include <math.h>
#include <float.h>

static double rescale_error (double err, const double result_abs, const double result_asc) ;

static double
rescale_error (double err, const double result_abs, const double result_asc)
{
  err = fabs(err) ;

  if (result_asc != 0 && err != 0)
      {
        double scale = pow((200 * err / result_asc), 1.5) ;

        if (scale < 1)
          {
            err = result_asc * scale ;
          }
        else
          {
            err = result_asc ;
          }
      }
  if (result_abs > DBL_MIN / (50 * DBL_EPSILON))
    {
      double min_err = 50 * DBL_EPSILON * result_abs ;

      if (min_err > err)
        {
          err = min_err ;
        }
    }

  return err ;
}

