#quad1d
A collection of C routines (requires C99 or above) and C++ wrappers (requires C++11 or above) for efficient 1D numerical integrations of complex-valued functions. Features an adaptive quadrature based on Gauss-Kronrod rules and a fixed-order Gauss-Legendre quadrature, both of which perform "**proper**" integrations of complex-valued functions, i.e., *automatic* decomposition into two real-domain integrations **without** wasted integrand evaluations (see below for more details). These routines thus require approximately **half** (or even less for difficult cases) the number of function evaluations, when compared to the commonly-deployed, *manual* decomposition into real domains. Also includes a unified C++ interface to some of the commonly-used quadratures from [GSL](http://www.gnu.org/software/gsl/). The latter serve as invaluable fallbacks when the "proper" adaptive quadrature chokes on difficult integrands.

##Dependencies
[GSL](http://www.gnu.org/software/gsl/) and [Boost](http://www.boost.org/)

##Why?
Most (if not all) 
