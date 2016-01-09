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
struct Ftor {  //a wrapper functor
    Ftor(complex<double> p, double q): p(p), q(q) { }
    complex<double> operator()(double x) const { return f_integrand(x, p, q); }
    complex<double> p;
    double q;
};

int main() {
    cout.precision(17); cout << scientific;
    const auto p = 1. + complex<double>(0.,1.) * 2.;
    const auto q = 2.;
    cout << "Expected value = " << f_expected(p, q) << endl;

    quad1d::Cag<Ftor> quad;
    complex<double> abserr;
    auto res = quad.integrate(Ftor(p, q), quad.neg_inf(), quad.pos_inf(), abserr);  //infinite interval
    // For semi-infinite interval, use (lower, upper) = (quad.neg_inf(), b) or (a, quad.pos_inf()) 
    cout << "Result = " << res << endl;
    cout << "Estimated error = " << abserr << endl;
}
