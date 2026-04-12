# Compressible Spline-Based Hyperelasticity

This repository provides an implementation of spline-based strain-energy functions for modeling compressible material behavior.

The stored-energy function is formulated as a truncated expansion in a selected set of kinematic invariants, here ${\bar{I}_1,\bar{I}_2,J}$. The model retains single-variable contributions and pairwise interaction terms, while neglecting higher-order couplings.

The strain-energy density is expressed as an additive composition of separable and non-separable contributions:

$$
\Psi^{(\bar{I}_1,\bar{I}_2,J)} =
\Psi^{(\bar{I}_1)} +
\Psi^{(\bar{I}_2)} +
\Psi^{(J)} +
g^{(J)}  h^{(\bar{I}1)} +
i^{(J)}  j^{(\bar{I}2)} +
k^{(\bar{I}_1)}  l^{(\bar{I}2)}
$$

This representation enables the modeling of both separable and coupled material responses.

The first three terms are univariate spline functions, whereas the remaining terms introduce coupling through multiplicative combinations of univariate splines. This construction can be interpreted as a low-rank approximation of a fully multivariate spline, retaining expressiveness while controlling the number of unknowns.

Training of the spline-based strain-energy function is performed via alternating optimization, which leads to a linear dependence of the stress on the spline coefficients. The coefficients correspond directly to interpolation values of the strain-energy function, providing a transparent and physically interpretable link between the model parameters and the energy landscape in invariant space.

Quick Start

## Quick start

Navigate to the main folder and run: `run_fit_eq()`.

This will generate an output file: `fit_result_eq.mat`.

To visualize the results, run: `postprocess_eq('fit_result_eq.mat')`.

In `run_fit_eq.m`, the variables `opts.use_XX` can be set to true or false to activate or deactivate specific terms in the energy function.
