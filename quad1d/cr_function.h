/**
 *  \file  cr_function.h
 *  \brief  Definition of an arbitrary complex-valued function with one real-valued free parameter.
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

#ifndef cr_function_INCLUDED
#define cr_function_INCLUDED

#include <complex.h>

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

#define FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary complex-valued function with double real parameters */

typedef struct cr_function
{
#ifdef __cplusplus
    std::complex<double> (*function)(double x, void *params);
#else
    double complex (*function)(double x, void *params);
#endif
    void *params;
} cr_function;

__END_DECLS

#endif
