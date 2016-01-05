/**
 *  \file  cr_qag.c
 *  \brief  Implementation of a 1D adaptive quadrature for complex-valued functions with real-valued free parameters.
 *
 * Based on QAG adaptive quadrature of GSL library.
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

#include "cr_quad1d.h"
#include "cr_function.h"
#include <gsl/gsl_errno.h>

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>

#include "initialise.c"
#include "set_initial.c"
#include "qpsrt.c"
#include "cr_qpsrt.c"
#include "cr_util.c"

static int
cr_qag(const cr_function *f,
       const double a, const double b,
       const double complex epsabs, const double complex epsrel,
       const size_t limit,
       qag_workspace *workspace_re,
       qag_workspace *workspace_im,
       double complex *result, double complex *abserr,
       cr_gsl_integration_rule *q);

static int
cr_qag_ext_re(const cr_function *f,
              const double epsabs, const double epsrel,
              const size_t limit, 
              qag_workspace *workspace_re,
              qag_workspace *workspace_im,
              size_t iteration_re, size_t iteration, 
              int roundoff_type1, int roundoff_type2, 
              double complex area, double complex *errsum,
              cr_gsl_integration_rule *q);

static int
cr_qag_ext_im(const cr_function *f,
              const double epsabs, const double epsrel,
              const size_t limit, 
              qag_workspace *workspace_re,
              qag_workspace *workspace_im,
              size_t iteration_im, size_t iteration, 
              int roundoff_type1, int roundoff_type2, 
              double complex area, double complex *errsum,
              cr_gsl_integration_rule *q);

int
cr_gsl_integration_qag(const cr_function *f,
                       double a, double b,
                       double complex epsabs, double complex epsrel, size_t limit,
                       int key,
                       qag_workspace *workspace_re,
                       qag_workspace *workspace_im,
                       double complex *result, double complex *abserr)
{
    cr_gsl_integration_rule *integration_rule;

    switch (key) 
    {
        case INTEG_GAUSS15:
            integration_rule = cr_gsl_integration_qk15;
            break;
        case INTEG_GAUSS21:
            integration_rule = cr_gsl_integration_qk21;
            break;
        case INTEG_GAUSS31:
            integration_rule = cr_gsl_integration_qk31; 
            break;
        case INTEG_GAUSS41:
            integration_rule = cr_gsl_integration_qk41;
            break;      
        case INTEG_GAUSS51:
            integration_rule = cr_gsl_integration_qk51;
            break;      
        case INTEG_GAUSS61:
            integration_rule = cr_gsl_integration_qk61;
            break;      
        default:
            GSL_ERROR("value of key does not specify a known integration rule", 
                                GSL_EINVAL);
}

    const int status = cr_qag(f, a, b, epsabs, epsrel, limit,
                              workspace_re, workspace_im, 
                              result, abserr, 
                              integration_rule);

    return status;
}

//**************************************************************************************

struct cr_iu_params
{
    const double a;
    const cr_function *const f ;
};

static double complex
cr_iu_transform(double t, void *params)
{
    const struct cr_iu_params *p = params;
    double x = p->a + (1 - t) / t;
    double complex y = FN_EVAL(p->f, x);
    return (y / t) / t;
}

int
cr_gsl_integration_qagiu(const cr_function *f, double a,
                         double complex epsabs, double complex epsrel, size_t limit,
                         int key,
                         qag_workspace *workspace_re,
                         qag_workspace *workspace_im,
                         double complex *result, double complex *abserr)
{
	cr_gsl_integration_rule *integration_rule;

	switch (key)
	{
		case INTEG_GAUSS15:
			integration_rule = cr_gsl_integration_qk15;
			break;
		case INTEG_GAUSS21:
			integration_rule = cr_gsl_integration_qk21;
			break;
		case INTEG_GAUSS31:
			integration_rule = cr_gsl_integration_qk31;
			break;
		case INTEG_GAUSS41:
			integration_rule = cr_gsl_integration_qk41;
			break;
		case INTEG_GAUSS51:
			integration_rule = cr_gsl_integration_qk51;
			break;
		case INTEG_GAUSS61:
			integration_rule = cr_gsl_integration_qk61;
			break;
		default:
			GSL_ERROR("value of key does not specify a known integration rule",
								GSL_EINVAL);
	}

  	struct cr_iu_params transform_params = {a, f};
	const cr_function f_transform = {&cr_iu_transform, &transform_params};

  	return (cr_qag(&f_transform, 0., 1., epsabs, epsrel, limit,
	               workspace_re, workspace_im,
	               result, abserr,
	               integration_rule));
}

struct cr_il_params
{
    const double b;
    const cr_function *const f ;
};

static double complex
cr_il_transform(double t, void *params)
{
    const struct cr_il_params *p = params;
    double x = p->b - (1 - t) / t;
    double complex y = FN_EVAL(p->f, x);
    return (y / t) / t;
}

int
cr_gsl_integration_qagil(const cr_function *f, double b,
                         double complex epsabs, double complex epsrel, size_t limit,
                         int key,
                         qag_workspace *workspace_re,
                         qag_workspace *workspace_im,
                         double complex *result, double complex *abserr)
{
	cr_gsl_integration_rule *integration_rule;

	switch (key)
	{
		case INTEG_GAUSS15:
			integration_rule = cr_gsl_integration_qk15;
			break;
		case INTEG_GAUSS21:
			integration_rule = cr_gsl_integration_qk21;
			break;
		case INTEG_GAUSS31:
			integration_rule = cr_gsl_integration_qk31;
			break;
		case INTEG_GAUSS41:
			integration_rule = cr_gsl_integration_qk41;
			break;
		case INTEG_GAUSS51:
			integration_rule = cr_gsl_integration_qk51;
			break;
		case INTEG_GAUSS61:
			integration_rule = cr_gsl_integration_qk61;
			break;
		default:
			GSL_ERROR("value of key does not specify a known integration rule",
								GSL_EINVAL);
	}

  	struct cr_il_params transform_params = {b, f};
	const cr_function f_transform = {&cr_il_transform, &transform_params};

  	return (cr_qag(&f_transform, 0., 1., epsabs, epsrel, limit,
	               workspace_re, workspace_im,
	               result, abserr,
	               integration_rule));
}

struct cr_i_params
{
    const cr_function *const f;
};

static double complex 
cr_i_transform(double t, void *params)
{
  const struct cr_i_params *p = params;
  double x = (1 - t) / t;
  double complex y = FN_EVAL(p->f, x) + FN_EVAL(p->f, -x);
  return (y / t) / t;
}

int
cr_gsl_integration_qagi(const cr_function *f,
                        double complex epsabs, double complex epsrel, size_t limit,
                        int key,
                        qag_workspace *workspace_re,
                        qag_workspace *workspace_im,
                        double complex *result, double complex *abserr)
{
	cr_gsl_integration_rule *integration_rule;

	switch (key)
	{
		case INTEG_GAUSS15:
			integration_rule = cr_gsl_integration_qk15;
			break;
		case INTEG_GAUSS21:
			integration_rule = cr_gsl_integration_qk21;
			break;
		case INTEG_GAUSS31:
			integration_rule = cr_gsl_integration_qk31;
			break;
		case INTEG_GAUSS41:
			integration_rule = cr_gsl_integration_qk41;
			break;
		case INTEG_GAUSS51:
			integration_rule = cr_gsl_integration_qk51;
			break;
		case INTEG_GAUSS61:
			integration_rule = cr_gsl_integration_qk61;
			break;
		default:
			GSL_ERROR("value of key does not specify a known integration rule",
								GSL_EINVAL);
	}

  	struct cr_i_params transform_params = {f};
	const cr_function f_transform = {&cr_i_transform, &transform_params};

  	return (cr_qag(&f_transform, 0., 1., epsabs, epsrel, limit,
	               workspace_re, workspace_im,
	               result, abserr,
	               integration_rule));
}

//**************************************************************************************

static int
cr_qag(const cr_function *f,
       const double a, const double b,
       const double complex epsabs, const double complex epsrel,
       const size_t limit,
       qag_workspace *workspace_re,
       qag_workspace *workspace_im,
       double complex *result, double complex *abserr,
       cr_gsl_integration_rule *q)
{
  /* Initialize results */
	initialise(workspace_re, a, b);
	initialise(workspace_im, a, b);

	*result = 0.;
	*abserr = 0.;

	if (limit > workspace_re->limit)
	{
		GSL_ERROR("iteration limit exceeds available workspace", GSL_EINVAL) ;
	}

    const double epsabs_re = creal(epsabs);
    const double epsabs_im = cimag(epsabs);
    const double epsrel_re = creal(epsrel);
    const double epsrel_im = cimag(epsrel);
	if (epsabs_re <= 0 && (epsrel_re < 50 * DBL_EPSILON || epsrel_re < 0.5e-28))
	{
		GSL_ERROR("tolerance for real part cannot be achieved with given epsabs and epsrel",
							 GSL_EBADTOL);
	}
	if (epsabs_im <= 0 && (epsrel_im < 50 * DBL_EPSILON || epsrel_im < 0.5e-28))
	{
		GSL_ERROR("tolerance for imaginary part cannot be achieved with given epsabs and epsrel",
							 GSL_EBADTOL);
	}

	/* perform the first integration */
	double complex result0, abserr0, resabs0, resasc0;
	q(f, a, b, &result0, &abserr0, &resabs0, &resasc0);

	double tolerance_re, tolerance_im; 
	int extension = 0;

	set_initial_result(workspace_re, creal(result0), creal(abserr0));
	set_initial_result(workspace_im, cimag(result0), cimag(abserr0));
	/* Test on accuracy */
	tolerance_re = fmax(epsabs_re, epsrel_re * fabs(creal(result0)));
	tolerance_im = fmax(epsabs_im, epsrel_im * fabs(cimag(result0)));

	/* need IEEE rounding here to match original quadpack behavior */
	//const double round_off_re = GSL_COERCE_DBL(50. * DBL_EPSILON * creal(resabs0));
	//const double round_off_im = GSL_COERCE_DBL(50. * DBL_EPSILON * cimag(resabs0));
	const double round_off_re = 50. * DBL_EPSILON * creal(resabs0);
	const double round_off_im = 50. * DBL_EPSILON * cimag(resabs0);


	if ((creal(abserr0) <= round_off_re && creal(abserr0) > tolerance_re)	||
            (cimag(abserr0) <= round_off_im && cimag(abserr0) > tolerance_im))
	{
		*result = result0;
		*abserr = abserr0;

		GSL_ERROR("cannot reach tolerance because of roundoff error "
							 "on first attempt", GSL_EROUND);
	}
	else if ((creal(abserr0) <= tolerance_re && creal(abserr0) != creal(resasc0)) || 
                 creal(abserr0) == 0.0) 
	{
		if ((cimag(abserr0) <= tolerance_im && cimag(abserr0) != cimag(resasc0)) || cimag(abserr0) == 0.0)
		{
			*result = result0;
			*abserr = abserr0;

			return GSL_SUCCESS;
		}
		extension = 2;
	}
	else if ((cimag(abserr0) <= tolerance_im && cimag(abserr0) != creal(resasc0)) || 
                 cimag(abserr0) == 0.0) 
	{
		extension = 1;
	}
	else if (limit == 1)
	{
		*result = result0;
		*abserr = abserr0;

		GSL_ERROR("a maximum of one iteration was insufficient", GSL_EMAXITER);
	}

	double complex area = result0;
	double complex errsum = abserr0;

	size_t iteration = 1, iteration_re = 1, iteration_im = 1;

	int roundoff_type1_re = 0, roundoff_type2_re = 0, roundoff_type1_im = 0, roundoff_type2_im = 0;
	int error_type = 1;			//if not changed after iteration, limit is too low 

	if (!extension)
	{
		do
		{
			/* Bisect the subinterval with the largest real error estimate */
			double a_i, b_i;
			double complex r_i, e_i;
			cr_retrieve_re(workspace_re, workspace_im, &a_i, &b_i, &r_i, &e_i);

			double a1 = a_i; 
			double b1 = 0.5 * (a_i + b_i);
			double a2 = b1;
			double b2 = b_i;

			double complex area1, area2, error1, error2, resabs1, resabs2, resasc1, resasc2;
			q(f, a1, b1, &area1, &error1, &resabs1, &resasc1);
			q(f, a2, b2, &area2, &error2, &resabs2, &resasc2);

			double complex area12 = area1 + area2;
			double complex error12 = error1 + error2;

			errsum += (error12 - e_i);
			area += area12 - r_i;

			if (creal(resasc1) != creal(error1) && creal(resasc2) != creal(error2))
			{
				const double delta_re = creal(r_i) - creal(area12);

				if (fabs(delta_re) <= 1.0e-5 * fabs(creal(area12)) && creal(error12) >= 0.99 * creal(e_i))
					roundoff_type1_re++;

				if (iteration_re >= 10 && creal(error12) > creal(e_i))
					roundoff_type2_re++;
			}

			cr_update_re(workspace_re, workspace_im, a1, b1, area1, error1, a2, b2, area2, error2);
			iteration++;
			iteration_re++;

			tolerance_re = fmax(epsabs_re, epsrel_re * fabs(creal(area)));
			tolerance_im = fmax(epsabs_im, epsrel_im * fabs(cimag(area)));

			if (creal(errsum) > tolerance_re) 
			{ 
				if (roundoff_type1_re >= 6 || roundoff_type2_re >= 20)
				{
					error_type = 2; /* round off error */ break;
				}	
				/* set error flag in the case of bad integrand behaviour at
					 a point of the integration range */

				if (subinterval_too_small(a1, a2, b2))
				{
					error_type = 3;  break;
				}

				if (cimag(errsum) <= tolerance_im)
				{
					extension = 1;  break;
				}
			}
			else
			{
				if (cimag(errsum) > tolerance_im)
				{
					extension = 2;  break;
				}
				else
				{
					extension = 0;  error_type = 0;  break;
				}
			}

			if (!(iteration < limit))
				break;

			//Continue by bisecting the subinterval with the largest imaginary error
			cr_retrieve_im(workspace_re, workspace_im, &a_i, &b_i, &r_i, &e_i);

			a1 = a_i; 
			b1 = 0.5 * (a_i + b_i);
			a2 = b1;
			b2 = b_i;

			q(f, a1, b1, &area1, &error1, &resabs1, &resasc1);
			q(f, a2, b2, &area2, &error2, &resabs2, &resasc2);

			area12 = area1 + area2;
			error12 = error1 + error2;

			errsum += (error12 - e_i);
			area += area12 - r_i;

			if (cimag(resasc1) != cimag(error1) && cimag(resasc2) != cimag(error2))
			{
				const double delta_im = cimag(r_i) - cimag(area12);

				if (fabs(delta_im) <= 1.0e-5 * fabs(cimag(area12)) && cimag(error12) >= 0.99 * cimag(e_i))
					roundoff_type1_im++;

				if (iteration_im >= 10 && cimag(error12) > cimag(e_i))
					roundoff_type2_im++;
			}

			cr_update_im(workspace_re, workspace_im, a1, b1, area1, error1, a2, b2, area2, error2);
			iteration++;
			iteration_im++;

			tolerance_re = fmax(epsabs_re, epsrel_re * fabs(creal(area)));
			tolerance_im = fmax(epsabs_im, epsrel_im * fabs(cimag(area)));

			if (cimag(errsum) > tolerance_im) 
			{ 
				if (roundoff_type1_im >= 6 || roundoff_type2_im >= 20)
				{
					error_type = 2; /* round off error */ break;
				}	
				/* set error flag in the case of bad integrand behaviour at
					 a point of the integration range */

				if (subinterval_too_small(a1, a2, b2))
				{
					error_type = 3;  break;
				}

				if (creal(errsum) <= tolerance_re)
				{
					extension = 2;  break;
				}
			}
			else
			{
				if (creal(errsum) > tolerance_re)
				{
					extension = 1;  break;
				}
				else
				{
					extension = 0;  error_type = 0;  break;
				}
			}
		}
		while (iteration < limit);
	}

	//test if extension is needed
	if (iteration < limit)
	{
		switch (extension)
		{
			case 1 :
				error_type = cr_qag_ext_re(f, epsabs_re, epsrel_re, limit, workspace_re,
                                           workspace_im, iteration_re, iteration,
                                           roundoff_type1_re, roundoff_type2_re,
                                           area, &errsum, q);
                break;
			case 2 :
				error_type = cr_qag_ext_im(f, epsabs_im, epsrel_im, limit, workspace_re, 
                                           workspace_im, iteration_im, iteration,
                                           roundoff_type1_im, roundoff_type2_im,
                                           area, &errsum, q);
                break;
		}
	}
    *result = sum_results(workspace_re) + I * sum_results(workspace_im);
    *abserr = errsum;

	switch (error_type)
	{	
		case 0 :
			return GSL_SUCCESS;	
		case 2 :
			GSL_ERROR ("roundoff error prevents tolerance from being achieved", GSL_EROUND);
		case 3 :
			GSL_ERROR ("bad integrand behavior found in the integration interval", GSL_ESING);
		case 1 :
			GSL_ERROR ("maximum number of subdivisions reached", GSL_EMAXITER);
		default :
			GSL_ERROR ("could not integrate function", GSL_EFAILED);
	}
}


static int
cr_qag_ext_re(const cr_function *f,
              const double epsabs, const double epsrel,
              const size_t limit, 
              qag_workspace *workspace_re,
              qag_workspace *workspace_im,
              size_t iteration_re, size_t iteration, 
              int roundoff_type1, int roundoff_type2, 
              double complex area, double complex *errsum,
              cr_gsl_integration_rule *q)
{
	double complex errsum_loc = *errsum;
	int error_type = 1;		//if not changed, subdivision limit is not enough

	do
	{
		double a_i, b_i;
		double complex r_i, e_i;
		cr_retrieve_re(workspace_re, workspace_im, &a_i, &b_i, &r_i, &e_i);

		double a1 = a_i; 
		double b1 = 0.5 * (a_i + b_i);
		double a2 = b1;
		double b2 = b_i;

		double complex area1, area2, error1, error2, resabs1, resabs2, resasc1, resasc2;
		q(f, a1, b1, &area1, &error1, &resabs1, &resasc1);
		q(f, a2, b2, &area2, &error2, &resabs2, &resasc2);

		double complex area12 = area1 + area2;
		double complex error12 = error1 + error2;

		errsum_loc += (error12 - e_i);
		area += area12 - r_i;

		if (creal(resasc1) != creal(error1) && creal(resasc2) != creal(error2))
		{
			const double delta_re = creal(r_i) - creal(area12);

			if (fabs(delta_re) <= 1.0e-5 * fabs(creal(area12)) && creal(error12) >= 0.99 * creal(e_i))
				roundoff_type1++;

			if (iteration_re >= 10 && creal(error12) > creal(e_i))
				roundoff_type2++;
		}

		cr_update_1st(workspace_re, workspace_im, a1, b1, area1, error1, a2, b2, area2, error2);
		iteration++;
		iteration_re++;

		double tolerance_re = fmax(epsabs, epsrel * fabs(creal(area)));

        if (creal(errsum_loc) > tolerance_re) 
        { 
            if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
            {
                error_type = 2; /* round off error */ break;
            }	
            /* set error flag in the case of bad integrand behaviour at
                 a point of the integration range */

            if (subinterval_too_small(a1, a2, b2))
            {
                error_type = 3;  break;
            }
        }
        else
        {
            error_type = 0;  break;
        }
	}
	while (iteration < limit); 

	*errsum = errsum_loc;
	return error_type;
}


static int
cr_qag_ext_im(const cr_function *f,
              const double epsabs, const double epsrel,
              const size_t limit, 
              qag_workspace *workspace_re,
              qag_workspace *workspace_im,
              size_t iteration_im, size_t iteration, 
              int roundoff_type1, int roundoff_type2, 
              double complex area, double complex *errsum,
              cr_gsl_integration_rule *q)
{
    double complex errsum_loc = *errsum;
    int error_type = 1;		//if not changed, subdivision limit is not enough

    do
    {
        double a_i, b_i;
        double complex r_i, e_i;
        cr_retrieve_im(workspace_re, workspace_im, &a_i, &b_i, &r_i, &e_i);

        double a1 = a_i; 
        double b1 = 0.5 * (a_i + b_i);
        double a2 = b1;
        double b2 = b_i;

        double complex area1, area2, error1, error2, resabs1, resabs2, resasc1, resasc2;
        q(f, a1, b1, &area1, &error1, &resabs1, &resasc1);
        q(f, a2, b2, &area2, &error2, &resabs2, &resasc2);

        double complex area12 = area1 + area2;
        double complex error12 = error1 + error2;

        errsum_loc += (error12 - e_i);
        area += area12 - r_i;

        if (cimag(resasc1) != cimag(error1) && cimag(resasc2) != cimag(error2))
        {
            const double delta_im = cimag(r_i) - cimag(area12);

            if (fabs(delta_im) <= 1.0e-5 * fabs(cimag(area12)) && cimag(error12) >= 0.99 * cimag(e_i))
                roundoff_type1++;

            if (iteration_im >= 10 && cimag(error12) > cimag(e_i))
                roundoff_type2++;
        }

        cr_update_2nd(workspace_re, workspace_im, a1, b1, area1, error1, a2, b2, area2, error2);
        iteration++;
        iteration_im++;

        double tolerance_im = fmax(epsabs, epsrel * fabs(cimag(area)));

        if (cimag(errsum_loc) > tolerance_im) 
        { 
            if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
            {
                error_type = 2; /* round off error */ break;
            }	
            /* set error flag in the case of bad integrand behaviour at
                 a point of the integration range */

            if (subinterval_too_small(a1, a2, b2))
            {
                error_type = 3;  break;
            }
        }
        else
        {
            error_type = 0;  break;
        }
    }
    while (iteration < limit); 

	*errsum = errsum_loc;
	return error_type;
}
