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
#include "quad1d/quad1d.hpp"
#include <cmath>
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
struct Ftor {  //a wrapper functor
    explicit Ftor(const complex<double>& mu): mu(mu) { }
    complex<double> operator()(double x) const { return f_integrand(x, mu); }
    complex<double> mu;
};
complex<double> mu_global;  //not recommended, only for illustration
complex<double> dumb_integrand(double x) { return f_integrand(x, mu_global); }  //a dumb wrapper function

int main() {
    cout.precision(17); cout << scientific;
    const auto mu = 40. + complex<double>(0., 1.) * 500.;  //a high value of mu = a rapidly-varying integrand
    cout << "Expected value = " << f_expected(mu) << endl;

    Ftor ftor_integrand(mu);
    quad1d::Cag<Ftor> quad;  //does "proper" adaptive quadrature, taking Ftor objects as integrands
    //      ^ explicit Cag<Ftor>(complex<double> epsrel = {1.E-10, 1.E-10}, complex<double> epsabs = {0., 0.},
    //                           std::size_t max_limit = 2000, GK_rule rule = GK_rule::GK31)
    // The real (imaginary) part of epsrel and epsabs are used in integration of the real (imaginary) part
    complex<double> abserr;
    auto res = quad.integrate(ftor_integrand, 0., 1., abserr);
    //         ^ complex<double> Cag<Ftor>::integrate(const Ftor& integrand, double a, double b,
    //                                                complex<double>& abserr, std::size_t limit = 2000)
    cout << "Result = " << res << endl;
    cout << "Estimated error = " << abserr << endl;

    // Other callables are also fine, provided you give quad1d::Cag the correct template type argument, e.g.,
    // A raw function:
    quad1d::Cag<> r_quad;  //uses the default type argument
    mu_global = mu;
    res = r_quad.integrate(dumb_integrand, 0., 1., abserr);
    // A lambda expression:
    auto l_integrand = [mu](double x){ return f_integrand(x, mu); };
    quad1d::Cag<decltype(l_integrand)> l_quad;
    res = l_quad.integrate(l_integrand, 0., 1., abserr);
    // Or a std::function, if you have to:
    function<complex<double>(double)> sf_integrand = l_integrand;
    quad1d::Cag<decltype(sf_integrand)> sf_quad;
    res = sf_quad.integrate(sf_integrand, 0., 1., abserr);
}
