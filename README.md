#quad1d
A collection of C routines (requires C99 or above) and C++ wrappers (requires C++11 or above) for efficient 1D numerical integrations of complex-valued functions. Features an adaptive quadrature based on Gauss-Kronrod rules and a fixed-order Gauss-Legendre quadrature, both of which perform "**proper**" integrations of complex-valued functions, i.e., *automatic* decomposition into two real-domain integrations **without** wasted integrand evaluations (see below for more details). These routines thus require approximately **half** the number of function evaluations, when compared with the commonly-deployed, *manual* decomposition into real domains. Also includes a unified C++ interface to some of the commonly-used quadratures from [GSL](http://www.gnu.org/software/gsl/). The latter serve as invaluable fallbacks when the "proper" adaptive quadrature chokes on difficult integrands.

##Dependencies
[GSL](http://www.gnu.org/software/gsl/) and [Boost](http://www.boost.org/)

##Why?
Most, if not all, freely-available C/C++ routines for 1D integrations (including wrappers for [QUADPACK](http://www.netlib.org/quadpack/)) only deal with real-valued integrands. This situation has (understandably) persisted for long since any integration of a complex-valued function can always be reformulated, or at least implemented, as two real-domain integrations. What is often overlooked, however, is the price one has to pay in doing so.

Let me elaborate a bit on that. In some cases I have encountered, a clean reformulation into two real-valued integrands is not feasible, either because the integrand contains complex-valued special functions (Bessel's or Hankel's with complex-valued arguments, etc), because it has to be evaluated recursively, or because its mathematical form is parameter-dependent. And even when it can be done cleanly, it often requires an unpleasant amount of algebraic manipulations, added programming complexity due to the resulting hideous expressions, and, worst of all, it often results in far-from-satisfying performance gain (relative to *manual* decomposition of the integrand). To be clear, what I mean by manually-decomposing is simply integrating two real-valued wrapper functions, both of which evaluate (call) the full complex-valued function, but return only its real or imaginary part. And, like most people, I always opted for this manual approach, thinking "*Man, this wrapper function always throws away the imaginary (or real) part. What a waste.*"

Still I was a happy camper, until the primary bottleneck in my calculation was exactly these wasted function evaluations. So I decided to trim some of my weekends and face it like a gentleman. Since then, I have used this "proper" adaptive quadrature extensively and have found no replacement for it. And although I am aware that my problem with conventional quadratures is an edge case which is unlikely to be of concern to most of their users, I figure it is still worth sharing considering the significant speedup it offers and the far-outweighed memory overhead it incurs. Finally, even if your integrands are too difficult for it, I am sure you will still benefit from the convenient C++ interface to the more-specialized quadratures from GSL.

##Usage Examples
Let us consider the integral `$\int_0^1 \frac{ \log(1+x) }{ (1+x)^{\mu+1} }$` with `$\mu$` being a complex number, which has a simple closed-form solution,
```C++
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
complex<double> mu_global;  //not recommended, only for illustration purpose (see below)
complex<double> dumb_integrand(double x) { return f_integrand(x, mu_global); }  //a dumb wrapper function

int main() {
    cout.precision(17); cout << scientific;
    const auto mu = 40. + complex<double>(0., 1.) * 500.;  //a high value of mu = a rapidly-varying integrand
    cout << "Expected value = " << f_expected(mu) << endl;

    quad1d::Cag<Ftor> quad;  //does "proper" adaptive quadrature, taking Ftor objects as integrands
    //      ^ explicit Cag<Ftor>(complex<double> epsrel = {1.E-10, 1.E-10}, complex<double> epsabs = {0., 0.},
    //                           std::size_t max_limit = 2000, GK_rule rule = GK_rule::GK31)
    // The real (imaginary) part of epsrel and epsabs are used in integration of the real (imaginary) part
    complex<double> abserr;
    auto res = quad.integrate(Ftor(mu), 0., 1., abserr);
    //         ^ complex<double> Cag<Ftor>::integrate(const Ftor& integrand, double a, double b,
    //                                                complex<double>& abserr, std::size_t limit = 2000)
    cout << "Result = " << res << endl;
    cout << "Estimated error = " << abserr << endl;

}
```
On my PC, this outputs:
```sh
Expected value = (-3.92401191889583165e-06,-6.31885976678998859e-07)
Result = (-3.92401191889603409e-06,-6.31885976679864209e-07)
Estimated error = (2.30461538976541977e-17,3.88190305613418434e-17)
```
Note that `quad1d::Cag` can also be used for other callables, provided you give it the correct template type argument, e.g.,
```C++
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
```    
To compute integrals over **infinite** or **semi-infinite** intervals, you only need to explicitly set the upper or lower limit (or both) of the integral to negative or positive infinity. For this purpose, `quad1d::Cag` defines its own infinity types, so that it can enforce safe template instantiations. That is, depending on whether the upper or lower limit is finite, the compiler will instantiate the right overload of `quad1d::Cag::integrate` to do the job, while ensuring that only those pairs of integration limits that are supported will compile. 
