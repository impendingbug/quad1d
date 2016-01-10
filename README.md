#quad1d
A collection of C routines (requires C99 or above) and C++ wrappers (requires C++11 or above) for efficient 1D numerical integrations of complex-valued functions. Features an adaptive quadrature based on Gauss-Kronrod rules and a fixed-order Gauss-Legendre quadrature, both of which perform "**proper**" integrations of complex-valued functions, i.e., *automatic* decomposition into two real-domain integrations **without** wasted integrand evaluations (see **Performance**). These routines thus require approximately **half** the number of function evaluations, when compared with the commonly-deployed, *manual* decomposition into real domains. Also includes a unified C++ interface to some of the commonly-used quadratures from [GSL](http://www.gnu.org/software/gsl/). The latter serve as invaluable fallbacks when the "proper" adaptive quadrature chokes on difficult integrands.

##Dependencies
[GSL](http://www.gnu.org/software/gsl/) and [Boost](http://www.boost.org/)

##Why?
Most, if not all, freely-available C/C++ routines for 1D integrations (including wrappers for [QUADPACK](http://www.netlib.org/quadpack/)) only deal with real-valued integrands. This situation has (understandably) persisted for long since any integration of a complex-valued function can always be reformulated, or at least implemented, as two real-domain integrations. What is often overlooked, however, is the price one has to pay in doing so.

Let me elaborate a bit on that. In some cases I have encountered, a clean reformulation into two real-valued integrands is not feasible, either because the integrand contains complex-valued special functions (Bessel's or Hankel's with complex-valued arguments, etc), because it has to be evaluated recursively, or because its mathematical form is parameter-dependent. And even when it can be done cleanly, it often requires an unpleasant amount of algebraic manipulations, added programming complexity due to the resulting hideous expressions, and, worst of all, it often results in far-from-satisfying performance gain (relative to *manual* decomposition of the integrand). To be clear, what I mean by manually-decomposing is simply integrating two real-valued wrapper functions, both of which evaluate (call) the full complex-valued function, but return only its real or imaginary part. And, like most people, I always opted for this manual approach, thinking "*Man, this wrapper function always throws away the imaginary (or real) part. What a waste.*"

Still I was a happy camper, until the primary bottleneck in my calculation was exactly these wasted function evaluations. So I decided to trim some of my weekends and face it like a gentleman. Since then, I have used this "proper" adaptive quadrature extensively and have found no replacement for it. And although I am aware that my problem with conventional quadratures is an edge case which is unlikely to be of concern to most of their users, I figure it is still worth sharing considering the significant speedup it offers and the far-outweighed memory overhead it incurs. Finally, even if your integrands are too difficult for it, I am sure you will still benefit from the convenient C++ interface to the more-specialized quadratures from [GSL](http://www.gnu.org/software/gsl/).

##Usage Examples
####Example 1
Let us consider the integral `$\int_0^1 \frac{ \log(1+x) }{ (1+x)^{\mu+1} }{\rm d}x$` with `$\mu$` being a complex number, which has a simple closed-form solution,
```C++
#include "quad1d/quad1d.hpp"
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
complex<double> mu_global;  //not recommended, only for illustration purpose
complex<double> dumb_integrand(double x) { return f_integrand(x, mu_global); }  //a dumb wrapper function
...
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
####Example 2
To compute integrals over **infinite** or **semi-infinite** intervals, you only need to explicitly set the upper or lower limit (or both) of the integral to negative or positive infinity. For this purpose, `quad1d::Cag` defines its own infinity types, so that it can enforce safe template instantiations. That is, depending on whether the upper or lower limit is finite, the compiler will instantiate the right overload of `quad1d::Cag::integrate` for the job, and only those pairs of integration limits that are supported will compile. As an example, let us compute $\int_{-\infty}^\infty \exp(-px^2 + qx) {\rm d}x$ for $\Re(p) > 0$,
```C++
#include "quad1d/quad1d.hpp"
#include <complex>
#include <stdexcept>
#include <iostream>

using namespace std;

complex<double> f_expected(const complex<double>& p, double q) {  //the closed-form solution
    if (p.real() <= 0.) throw logic_error("f_expected() : re(p) must be > 0)");
    return exp(0.25 * q * q / p) * sqrt(M_PI / p);
}
complex<double> f_integrand(double x, const complex<double>& p, double q) { return exp(-p * x * x + q * x); }
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
```
which outputs
```sh
Expected value = (8.37912751875862560e-01,-1.18061875415731632e+00)
Result = (8.37912751875862005e-01,-1.18061875415731632e+00)
Estimated error = (6.85431097035958387e-14,1.17011564912306070e-13)
```
*Usage notes:*
* The callables are passed to the `integrate` method by a reference-to-`const` (be careful detaching your thread). A functor object passed to it thus must have a `const` overload of its `operator()`. But it need not be copy- nor move-constructible. I always find this more convenient. An added benefit is that a single `quad1d::Cag` instance can handle many functor objects of different types by making use of polymorphism.
* `quad1d::Cag` derives its interval-bisection procedure and error-estimation mechanism from [`gsl_integration_qag`](https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html). It is thus not intended for integrands containing singularities, for which `quad1d::Cag_gsl` is more suitable (see **Performance**).
* Integrations over (semi-)infinite intervals involve variable transformations which might introduce integrable singularities in the integrands. In such cases, `quad1d::Cag_gsl` is the one to reach for as it internally uses [`gsl_integration_qags`](https://www.gnu.org/software/gsl/manual/html_node/QAGS-adaptive-integration-with-singularities.html#QAGS-adaptive-integration-with-singularities) when handling (semi-)infinite intervals.
* An instance of `std::bad_alloc` is thrown if the constructor fails to allocate memory for integration workspace.
* The `integrate` method throws `std::logic_error` or `std::runtime_error`. In the latter case, the integration result can still be accessed if one uses the version of `integrate` that returns `void` (see `quad1d/quad1d.hpp`).
* The C routines wrapped by `quad1d::Cag` are declared in `quad1d/cr_quad1d.h`. Its API closely follows the [quadratures from GSL](https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html#Numerical-Integration).

##Performance
In the following, the performance of `quad1d::Cag` is compared against `quad1d::Cag_gsl`. The latter is a "smart" wrapper for (real-domain) adaptive quadratures from GSL. It internally performs two independent real-domain integrations, each of which discards the real or imaginary component of the integrand. This is in contrast with `quad1d::Cag` which integrates both components "in parallel" such that every integrand evaluation can be fully utilized (see below for more details). For the comparison I reuse the integrand in **Example 1** and modify the functor to include a call counter,
```C++
struct Ftor {  //not thread-safe
    explicit Ftor(const complex<double>& mu): mu(mu) { }
    complex<double> operator()(double x) const { ++counter; return f_integrand(x, mu); }
    complex<double> mu;
    mutable size_t counter {0};
};
```
which I then use as follows,
```C++
...
    Ftor integrand(mu);
    quad1d::Cag<Ftor> quad;
    complex<double> abserr;
    auto res = quad.integrate(integrand, 0., 1., abserr);
    cout << "Result = " << res << endl;
    cout << "Estimated error = " << abserr << endl;
    cout << "Function calls = " << integrand.counter << endl;

    integrand.counter = 0;

    quad1d::Cag_gsl<Ftor> quad_gsl;  //A wrapper class for (real-domain) adaptive quadratures from GSL
    res = quad_gsl.integrate(integrand, 0., 1., abserr);
    cout << "Result (GSL) = " << res << endl;
    cout << "Estimated error (GSL) = " << abserr << endl;
    cout << "Function calls (GSL) = " << integrand.counter << endl;
```
to give
```sh
Expected value = (-3.92401191889583165e-06,-6.31885976678998859e-07)
Result = (-3.92401191889603409e-06,-6.31885976679864209e-07)
Estimated error = (2.30461538976541977e-17,3.88190305613418434e-17)
Function calls = 1147
Result (GSL) = (-3.92401191889603239e-06,-6.31885976679851080e-07)
Estimated error (GSL) = (2.89859102927190174e-17,3.88301820442891792e-17)
Function calls (GSL) = 2232
```
Similarly, reusing the integrand in **Example 2** gives
```sh
Expected value = (8.37912751875862560e-01,-1.18061875415731632e+00)
Result = (8.37912751875862005e-01,-1.18061875415731632e+00)
Estimated error = (6.85431097035958387e-14,1.17011564912306070e-13)
Function calls = 806
Result (GSL) = (8.37912751875862449e-01,-1.18061875415731610e+00)
Estimated error (GSL) = (4.75964513994949376e-11,1.55688200563444364e-11)
Function calls (GSL) = 2220
```
in which `Cag` is almost three (instead of two) times more efficient. This is expected since, for (semi-)infinite intervals, `Cag_gsl` uses [`gsl_integration_qags`](https://www.gnu.org/software/gsl/manual/html_node/QAGS-adaptive-integration-with-singularities.html#QAGS-adaptive-integration-with-singularities) which inherently does more work than [`gsl_integration_qag`](https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html). Also note that in these two comparisons, `Cag` and `Cag_gsl` both integrate with the same requested accuracy and Gauss-Kronrod rule.

*More details:*
`quad1d::Cag` is as reliable as [`gsl_integration_qag`](https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html) as they both share the same underlying algorithm (from [QUADPACK](http://www.netlib.org/quadpack/)). This algorithm keeps track of the estimated error of each integration subinterval, and at each iteration, it bisects the subinterval with the largest error. The iteration stops when either the requested accuracy has been achieved, the maximum number of iterations has been reached, or it is decided that further iterations will no longer improve accuracy due to round-off errors. `Cag` simply ensures that *every* integrand evaluation is used to improve integration accuracy of both the real and imaginary components. To achieve this, the two real-domain integrations are forced to share the same set of subintervals. The right to pick which subinterval to bisect is passed back and forth between the two domains, thus retaining the desired capability of isolating "bad" regions in both domains. If, for instance, the real-part integration achieves its requested accuracy before the imaginary part, all upcoming bisections are chosen by the imaginary-part integration, while the real-part result still benefits from each of these bisections. 

##Compiling and Linking
To compile, `cd` into the outer `quad1d` directory and run `make`. This will create the static library `libquad1d.a` which you can link to by, e.g.,
```sh
g++ -std=c++11 -Wall -I./quad1d -c example.cpp
g++ example.o -L./quad1d -lquad1d -lgsl -lgslcblas -lm -o example
```
assuming that `quad1d` directory is in `./`, while the GSL library and the Boost headers are installed in their default locations. The relevant header files can be found in `quad1d/quad1d` subdirectory. Thus `example.cpp` needs either
```C++
#include "quad1d/quad1d.hpp"  //use C++ interface
```
or
```C
#include "quad1d/cr_quad1d.h"  //use C interface
```
(see also **Usage Examples**).

