/**
 *  \file  cr_qk.c
 *  \brief  Extension of integration/qk.c from GSL source code, to accomodate complex-valued functions with
 *          real-valued free parameters.
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
 *  Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#include "cr_quad1d.h"
#include "cr_function.h"

#include <math.h>
#include <complex.h>

#include "err.c"

void
cr_gsl_integration_qk(const int n, const double xgk[], const double wg[], const double wgk[],
                      double complex fv1[], double complex fv2[],
                      const cr_function *f, double a, double b,
                      double complex *result, double complex *abserr,
                      double complex *resabs, double complex *resasc)
{

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs(half_length);
  const double complex f_center = FN_EVAL(f, center);

  double complex result_gauss;

  if (n % 2)
	  result_gauss = 0.;
  else
	  result_gauss = f_center * wg[n / 2 - 1];

  double complex result_kronrod = f_center * wgk[n - 1];
  double complex result_abs = fabs(creal(result_kronrod)) + I * fabs(cimag(result_kronrod));

  for (int j = 0; j < (n - 1) / 2; j++)
	{
		const int jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
		const double abscissa = half_length * xgk[jtw];
		const double complex fval1 = FN_EVAL(f, center - abscissa);
		const double complex fval2 = FN_EVAL(f, center + abscissa);
		const double complex fsum = fval1 + fval2;
		fv1[jtw] = fval1;
		fv2[jtw] = fval2;
		result_gauss += wg[j] * fsum;
		result_kronrod += wgk[jtw] * fsum;
		result_abs += wgk[jtw] * ((fabs(creal(fval1)) + fabs(creal(fval2))) + I * (fabs(cimag(fval1)) + fabs(cimag(fval2))));
	}

  for (int j = 0; j < n / 2; j++)
	{
		const int jtwm1 = j * 2;
		const double abscissa = half_length * xgk[jtwm1];
		const double complex fval1 = FN_EVAL(f, center - abscissa);
		const double complex fval2 = FN_EVAL(f, center + abscissa);
		fv1[jtwm1] = fval1;
		fv2[jtwm1] = fval2;
		result_kronrod += wgk[jtwm1] * (fval1 + fval2);
		result_abs += wgk[jtwm1] * ((fabs(creal(fval1)) + fabs(creal(fval2))) + I * (fabs(cimag(fval1)) + fabs(cimag(fval2))));
	}

  double complex mean = result_kronrod * 0.5;

  const double complex temp = f_center - mean;
  double complex result_asc = wgk[n - 1] * (fabs(creal(temp)) + I * fabs(cimag(temp)));

  for (int j = 0; j < n - 1; j++)
	{
		const double complex temp1 = fv1[j] - mean;
		const double complex temp2 = fv2[j] - mean;
		result_asc += wgk[j] * ((fabs(creal(temp1)) + fabs(creal(temp2))) + I * (fabs(cimag(temp1)) + fabs(cimag(temp2))));
	}

  /* scale by the width of the integration region */

  double complex err = (result_kronrod - result_gauss) * half_length;

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  *result = result_kronrod;
  *resabs = result_abs;
  *resasc = result_asc;
  *abserr = rescale_error(creal(err), creal(result_abs), creal(result_asc)) +
                          I * rescale_error(cimag(err), cimag(result_abs), cimag(result_asc));

}
