/**
 *  \file  cr_quad1d.h
 *  \brief  1D quadratures for complex-valued functions
 *
 *  <+DETAILED+>
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

#ifndef cr_quad1d_INCLUDED
#define cr_quad1d_INCLUDED

#include "quad1d/cr_function.h"

#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>


#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

//**************************************************************************************
/* Workspace for adaptive integrators */
typedef struct qag_workspace
{
	size_t limit;
	size_t size;
	size_t nrmax;
	size_t i;
	size_t maximum_level;
	double *alist;
	double *blist;
	double *rlist;
	double *elist;
	size_t *order;
	size_t *level;
} qag_workspace;

qag_workspace *
qag_workspace_alloc(const size_t n);

void
qag_workspace_free(const qag_workspace *w);

/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them.  */

enum
{
    INTEG_GAUSS15 = 1,      /* 15 point Gauss-Kronrod rule */
    INTEG_GAUSS21 = 2,      /* 21 point Gauss-Kronrod rule */
    INTEG_GAUSS31 = 3,      /* 31 point Gauss-Kronrod rule */
    INTEG_GAUSS41 = 4,      /* 41 point Gauss-Kronrod rule */
    INTEG_GAUSS51 = 5,      /* 51 point Gauss-Kronrod rule */
    INTEG_GAUSS61 = 6       /* 61 point Gauss-Kronrod rule */
};

#ifdef __cplusplus
int
cr_gsl_integration_qag(const cr_function *f,
                       double a, double b,
                       std::complex<double> epsabs, std::complex<double> epsrel, size_t limit,
                       int key,
                       qag_workspace *workspace_real,
                       qag_workspace *workspace_imag,
                       std::complex<double> *result, std::complex<double> *abserr);

int
cr_gsl_integration_qagiu(const cr_function *f, double a,
                         std::complex<double> epsabs, std::complex<double> epsrel, size_t limit,
                         int key,
                         qag_workspace *workspace_re,
                         qag_workspace *workspace_im,
                         std::complex<double> *result, std::complex<double> *abserr);

int
cr_gsl_integration_qagil(const cr_function *f, double b,
                         std::complex<double> epsabs, std::complex<double> epsrel, size_t limit,
                         int key,
                         qag_workspace *workspace_re,
                         qag_workspace *workspace_im,
                         std::complex<double> *result, std::complex<double> *abserr);
int
cr_gsl_integration_qagi(const cr_function *f,
                        std::complex<double> epsabs, std::complex<double> epsrel, size_t limit,
                        int key,
                        qag_workspace *workspace_re,
                        qag_workspace *workspace_im,
                        std::complex<double> *result, std::complex<double> *abserr);
#else
int
cr_gsl_integration_qag(const cr_function *f,
                       double a, double b,
                       double complex epsabs, double complex epsrel, size_t limit,
                       int key,
                       qag_workspace *workspace_real,
                       qag_workspace *workspace_imag,
                       double complex *result, double complex *abserr);
int
cr_gsl_integration_qagiu(const cr_function *f, double a,
                         double complex epsabs, double complex epsrel, size_t limit,
                         int key,
                         qag_workspace *workspace_re,
                         qag_workspace *workspace_im,
                         double complex *result, double complex *abserr);
int
cr_gsl_integration_qagil(const cr_function *f, double b,
                         double complex epsabs, double complex epsrel, size_t limit,
                         int key,
                         qag_workspace *workspace_re,
                         qag_workspace *workspace_im,
                         double complex *result, double complex *abserr);
int
cr_gsl_integration_qagi(const cr_function *f,
                        double complex epsabs, double complex epsrel, size_t limit,
                        int key,
                        qag_workspace *workspace_re,
                        qag_workspace *workspace_im,
                        double complex *result, double complex *abserr);
#endif

//**************************************************************************************
#ifndef __cplusplus
typedef void cr_gsl_integration_rule(const cr_function *f, double a, double b,
                                     double complex *result, double complex *abserr,
                                     double complex *resabs, double complex *resasc);

void cr_gsl_integration_qk15(const cr_function *f, double a, double b,
                             double complex *result, double complex *abserr,
                             double complex *resabs, double complex *resasc);

void cr_gsl_integration_qk21(const cr_function *f, double a, double b,
                             double complex *result, double complex *abserr,
                             double complex *resabs, double complex *resasc);

void cr_gsl_integration_qk31(const cr_function *f, double a, double b,
                             double complex *result, double complex *abserr,
                             double complex *resabs, double complex *resasc);

void cr_gsl_integration_qk41(const cr_function *f, double a, double b,
                             double complex *result, double complex *abserr,
                             double complex *resabs, double complex *resasc);

void cr_gsl_integration_qk51(const cr_function *f, double a, double b,
                             double complex *result, double complex *abserr,
                             double complex *resabs, double complex *resasc);

void cr_gsl_integration_qk61(const cr_function *f, double a, double b,
                             double complex *result, double complex *abserr,
                             double complex *resabs, double complex *resasc);

void cr_gsl_integration_qk(const int n, const double xgk[],
                           const double wg[], const double wgk[],
                           double complex fv1[], double complex fv2[],
                           const cr_function *f, double a, double b,
                           double complex *result, double complex *abserr,
                           double complex *resabs, double complex *resasc);
#endif //ifndef C++

//**************************************************************************************

/* Workspace for fixed-order Gauss-Legendre integration */

typedef struct glfixed_table
{
    size_t n;         /* number of points */
    double *x;        /* Gauss abscissae/points */
    double *w;        /* Gauss weights for each abscissae */
    int precomputed;  /* high precision abscissae/weights precomputed? */
} glfixed_table;

glfixed_table *
glfixed_table_alloc(size_t n);

void
glfixed_table_free(glfixed_table * t);

#ifdef __cplusplus
std::complex<double>
cr_gsl_integration_glfixed(const cr_function *f,
                           double a, double b,
                           const glfixed_table *t);
#else
double complex
cr_gsl_integration_glfixed(const cr_function *f,
                           double a, double b,
                           const glfixed_table *t);
#endif
//******************************************************************************************
__END_DECLS

#endif
