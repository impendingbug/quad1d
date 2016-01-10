/**
 *  \file  test_perform.cpp
 *  \brief  
 *
 *  <+DETAILED+>
 *
 *  \author  Sandy Pratama
 *
 *  \internal
 *       Created:  09-01-16
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 2014 Sandy Pratama
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */
#include "quad1d/quad1d.hpp"
#include <cstddef>
#include <complex>
#include <iostream>

using namespace std;

complex<double> f_expected(const complex<double>& mu) {  //the closed-form solution
    auto temp = pow(2., mu);
    return -log(2.) / (temp * mu) + (temp - 1.) / (temp * mu * mu);
}
complex<double> f_integrand(double x, const complex<double>& mu) {  //the integrand
    return log(1. + x) / pow(1. + x, mu + 1.);
}
struct Ftor {  //not thread-safe
    explicit Ftor(const complex<double>& mu): mu(mu) { }
    complex<double> operator()(double x) const { ++counter; return f_integrand(x, mu); }
    complex<double> mu;
    mutable size_t counter {0};
};

int main() {
    cout.precision(17); cout << scientific;
    const auto mu = 40. + complex<double>(0., 1.) * 500.;  //a high value of mu = a rapidly-varying integrand
    cout << "Expected value = " << f_expected(mu) << endl;

    Ftor integrand(mu);
    quad1d::Cag<Ftor> quad;
    complex<double> abserr;
    auto res = quad.integrate(integrand, 0., 1., abserr);
    cout << "Result = " << res << endl;
    cout << "Estimated error = " << abserr << endl;
    cout << "Function calls = " << integrand.counter << endl;

    integrand.counter = 0;

    quad1d::Cag_gsl<Ftor> quad_gsl;
    res = quad_gsl.integrate(integrand, 0., 1., abserr);
    cout << "Result (GSL) = " << res << endl;
    cout << "Estimated error (GSL) = " << abserr << endl;
    cout << "Function calls (GSL) = " << integrand.counter << endl;
}
