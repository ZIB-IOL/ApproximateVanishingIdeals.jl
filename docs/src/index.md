# AVI.jl

[![][docs-img]][docs-url]

This package is a toolbox for polynomial feature extraction and transformation using the Oracle Approximate Vanishing Ideal Algorithm.

## Overview

The Oracle Approximate Vanishing Ideal ($\texttt{OAVI}$) algorithm was designed to compute the vanishing ideal of a set of points. Instead of adopting the then common approach of
using singular value decomposition, $\texttt{OAVI}$ finds vanishing polynomials by solving a convex optimization problem of the form
```math
\min_{x \in C} f(x),
```
where $f$ is a differentiable convex function and $C$ is a compact and convex set. Usually $f$ will be of the form
```math
f(x) = \|Ax + b\|_2^2.
```

## Installation
The most recent release is available via:
```julia
Pkg.add(url="https://github.com/ZIB-IOL/AVI.jl")
```

## Getting started
We provide built-in oracles that construct the objective function and feasible region and solve the optimization problem with a version of the Frank-Wolfe (conditional gradients) algorithm implemented in [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl/tree/master). Obtaining a basic feature transformation for some random data $\texttt{X}$ is as simple as:

```julia
using AVI
using random

X = rand(10000, 5)
X_transformed, sets = fit_oavi(X);
```
`X_transformed` holds the transformed data and `sets` keeps track of important sets. It is recommended to adjust some keyword arguments for better results. See the examples section for more information.

## Documentation and Examples
To explore the contents of the package and see examples with a more detailed look at the different keyword arguments, go to the [documentation](https://github.com/ZIB-IOL/AVI).


[docs-img]: https://img.shields.io/badge/docs-latest%20release-blue.svg
[docs-url]: https://github.com/ZIB-IOL/AVI


