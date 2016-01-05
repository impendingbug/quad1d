/**
 *  \file  test_cag_inf.cpp
 *  \brief  
 *
 *  <+DETAILED+>
 *
 *  \author  Sandy Pratama
 *
 *  \internal
 *       Created:  05-01-16
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 2014 Sandy Pratama
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#include <boost/math/constants/constants.hpp>

#include "quad1d.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>

#define SPIT(n) (std::cout << #n " = " << n << std::endl)

using namespace std;
using cx_double = complex<double>;

cx_double f_expected(cx_double p_sqr, double q) {
    if (p_sqr.real() <= 0.) throw logic_error("f_expected() : re(p_sqr) must be > 0)");
    return exp(0.25 * q * q / p_sqr) * sqrt(M_PI / p_sqr);
}

cx_double f_integrand(double x, cx_double p_sqr, double q) {
    return exp(-p_sqr * x * x + q * x);
}

struct ftor {
    ftor(cx_double p_sqr, double q): p_sqr(p_sqr), q(q) { }
    cx_double operator()(double x) const { ++counter; return f_integrand(x, p_sqr, q); }
    cx_double p_sqr;
    double q;
    mutable size_t counter {0};
};


int main() {
    cout.precision(17);
    cout << scientific;

    cx_double abserr;

    using namespace quad1d;
    auto p_sqr = 1. + cx_double(0.,1.) * 2.;
    auto q = 2.;
    ftor func(p_sqr, q);

    Cag<decltype(func)> quad;
    Cag_gsl<decltype(func)> quad_gsl;
    
    SPIT(f_expected(p_sqr, q));
    auto res_gsl = quad_gsl.integrate(func, quad_gsl.neg_inf(), quad_gsl.pos_inf(), abserr);
    SPIT(res_gsl);
    SPIT(abserr);
    SPIT(func.counter);

    func.counter = 0;

    auto res = quad.integrate(func, quad.neg_inf(), quad.pos_inf(), abserr);
    SPIT(res);
    SPIT(abserr);
    SPIT(func.counter);
}
