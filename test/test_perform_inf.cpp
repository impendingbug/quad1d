/**
 *  \file  test_perform_inf.cpp
 *  \brief  
 *
 *  <+DETAILED+>
 *
 *  \author  Sandy Pratama
 *
 *  \internal
 *       Created:  10-01-16
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
#include <complex>
#include <stdexcept>
#include <iostream>

using namespace std;

complex<double> f_expected(const complex<double>& p, double q) {  //the closed-form solution
    if (p.real() <= 0.) throw logic_error("f_expected() : re(p) must be > 0)");
    return exp(0.25 * q * q / p) * sqrt(M_PI / p);
}
complex<double> f_integrand(double x, const complex<double>& p, double q) {  //the integrand
    return exp(-p * x * x + q * x);
}
struct Ftor {  //not thread-safe
    Ftor(complex<double> p, double q): p(p), q(q) { }
    complex<double> operator()(double x) const { ++counter; return f_integrand(x, p, q); }
    complex<double> p;
    double q;
    mutable size_t counter {0};
};

int main() {
    cout.precision(17); cout << scientific;
    const auto p = 1. + complex<double>(0.,1.) * 2.;
    const auto q = 2.;
    cout << "Expected value = " << f_expected(p, q) << endl;

    Ftor integrand(p, q);
    quad1d::Cag<Ftor> quad;
    complex<double> abserr;
    auto res = quad.integrate(integrand, quad.neg_inf(), quad.pos_inf(), abserr);
    cout << "Result = " << res << endl;
    cout << "Estimated error = " << abserr << endl;
    cout << "Function calls = " << integrand.counter << endl;

    integrand.counter = 0;

    quad1d::Cag_gsl<Ftor> quad_gsl;
    res = quad_gsl.integrate(integrand, quad.neg_inf(), quad.pos_inf(), abserr);
    cout << "Result (GSL) = " << res << endl;
    cout << "Estimated error (GSL) = " << abserr << endl;
    cout << "Function calls (GSL) = " << integrand.counter << endl;
}

