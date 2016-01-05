/**
 *  \file  test_cag.cpp
 *  \brief  
 *
 *  <+DETAILED+>
 *
 *  \author  Sandy Pratama
 *
 *  \internal
 *       Created:  13-12-15
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 2014 Sandy Pratama
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#include "quad1d.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <functional>

#define SPIT(n) (std::cout << #n " = " << n << std::endl)

using namespace std;
using cx_double = complex<double>;

cx_double f_expected(cx_double mu) {
    auto temp = pow(2., mu);
    return -log(2.) / (temp * mu) + (temp - 1.) / (temp * mu * mu);
}

cx_double f_integrand(double x, cx_double mu) {
    return log(1. + x) / pow(1. + x, mu + 1.);
}

struct ftor {
    explicit ftor(cx_double mu): mu(mu) { }
    cx_double operator()(double x) const { ++counter; return f_integrand(x, mu); }
    cx_double mu;
    mutable size_t counter {0};
};


int main() {
    cout.precision(17);
    cout << scientific;

    cx_double abserr;

    using namespace quad1d;
    auto mu = 10. + cx_double(0.,1.) * 100.;

    ftor func(mu);

    Cag<ftor> quad;
    Cag_gsl<ftor> quad_gsl;
    
    SPIT(f_expected(mu));
    auto res_gsl = quad_gsl.integrate(func, 0., 1., abserr);
    SPIT(res_gsl);
    SPIT(abserr);
    SPIT(func.counter);

    func.counter = 0;

    auto res = quad.integrate(func, 0., 1., abserr);
    SPIT(res);
    SPIT(abserr);
    SPIT(func.counter);
}
