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
#include <cstddef>
#include <iostream>

using namespace std;

complex<double> f_expected(complex<double> mu) {  //the closed-form solution
    auto temp = pow(2., mu);
    return -log(2.) / (temp * mu) + (temp - 1.) / (temp * mu * mu);
}
complex<double> f_integrand(double x, complex<double> mu) {  //the integrand
    return log(1. + x) / pow(1. + x, mu + 1.);
}
struct Ftor {  //a wrapper functor
    explicit Ftor(complex<double> mu): mu(mu) { }
    complex<double> operator()(double x) const { ++counter; return f_integrand(x, mu); }
    complex<double> mu;
    mutable size_t counter {0};
};

int main() {
    cout.precision(17); cout << scientific;

    auto mu = 40. + complex<double>(0.,1.) * 500.;  //a high value of mu yields a rapidly-varying integrand
    Ftor integrand(mu);
    cout << "Expected value = " << f_expected(mu) << endl;
    
    quad1d::Cag<Ftor> quad;
    complex<double> abserr;
    auto res = quad.integrate(integrand, 0., 1., abserr);
    cout << "Result = " << res << endl;
    cout << "Estimated error = " << abserr << endl;
    cout << "Number of function evaluations = " << integrand.counter << endl;
}
```
