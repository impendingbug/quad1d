##quad1d

A collection of C routines (requires C99 or above) and C++ wrappers (requires C++11 or above) for 1D numerical integrations of complex-valued functions. Features an adaptive quadrature based on Gauss-Kronrod rules and a fixed-order Gauss-Legendre quadrature, both of which perform "*proper*" integrations of complex-valued functions, i.e., automatic decomposition into two real-domain integrations **without** wasted integrand evaluations.
