# RigorousInvariantMeasures

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/RigorousInvariantMeasures.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliadynamics.github.io/RigorousInvariantMeasures.jl/dev/)
[![Build Status](https://github.com/JuliaDynamics/RigorousInvariantMeasures.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaDynamics/RigorousInvariantMeasures.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaDynamics/RigorousInvariantMeasures.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/RigorousInvariantMeasures.jl)

This Package provides methods for the rigorous approximation of Absolutely Continuous Invariant Measures for one dimensional dynamical systems,
using the results in [[2]](#2) and subsequent [[1]](#1).

## Mathematical background
By Birkhoff Ergodic Theorem we know that if a dynamical system T admits an ergodic invariant measure μ, for μ-almost every point x the frequence of the visits to a set E corresponds to the measure of the set with respect to μ.

Therefore, being able to approximate invariant measures with a large basin is interesting to investigate the statistical properties of the dynamical system T. [[5]](#4)

The existence of absolutely continuous invariant measures for one dimensional maps is a delicate topic; with this package we present approximation schemes for the invariant measures of system that satisfy a Lasota-Yorke inequality through the use of a coarse-fine scheme and a posteriori estimates on the mixing time, i.e., this means that the mixing time is estimated by our algorithm and we do not need an a priori estimate, which is usually difficult to obtain.

The Ulam approximation schemes works under relatively weak hypothesis on the dynamics and was used in [[4]](#4) to approximate the invariant measure for the geometric Lorenz 1-dimensional map.

We are currently working on the implementation of the Ulam scheme for system with additive uniform noise, as the one used in
[[2]](#2)

## Basic Usage
Examples of usage are present in the directory examples.

```julia
using InvariantMeasures
D = Mod1Dynamic(x -> 4x + 0.01InvariantMeasures.sinpi(8x))
B = Ulam(1024)
Q = DiscretizedOperator(B, D)
```

The code snippet above defines a dynamic obtained by reducing f(x) = 4x+0.01 sin(8πx) modulo 1, a basis B associated to the Ulam discretization on a partition of 1024 homogenous intervals, and computes the discretized operator Q, a Markov chain whose entries are P[T(x)∈ Iᵢ | x ∈ Iⱼ].

Note the usage of `InvariantMeasures.sinpi(8*x)` rather than `Base.sinpi` or `Base.sin(8\pi*x)`. This detail is required to ensure that f(1) == 4 exactly.

```julia
norms = norms = powernormbounds(B, D; Q=Q)
```

This function computes the L¹ norm of Q^k, for k = 1,2,...,k_max (up to a sufficiently large number of powers to observe decay) when restricted to the space U of average-0 vectors. This gives us the a posteriori estimate for the mixing time of the Markov chain and is used in our rigorous estimate.

```julia
w = invariant_vector(B, Q)
distance_from_invariant(B, D, Q, w, norms)
```
This computes a (non-rigorous) approximation of the invariant measure of D; then  `distance_from_invariant` computes an upper bound for the L¹ distance between w and the density of the absolutely continuous invariant measure of the system.

Inside the examples it is showed how to use the coarse-fine scheme to obtain better L¹ bounds and reduce the computational time.

### Caveat
The function `sinpi` in the interval arithmetic package that we are using relies on the `CRlibm.jl` package, which currently works only under Linux. So the examples that use trigonometric functions only work on this OS.

## References
<a id="1">[1]</a> Galatolo S., Monge M., Nisoli I., Poloni F. A general framework for the rigorous computation of invariant densities and the coarse-fine strategy
https://doi.org/10.48550/arXiv.2212.05017

<a id="2">[2]</a>
Galatolo S., Monge M., Nisoli I., Existence of noise induced order, a computer aided proof Nonlinearity 33 (9), 4237 (2020)

<a id="3">[3]</a>
Galatolo S., Nisoli I., An elementary approach to rigorous approximation of invariant
measures SIAM J. Appl Dyn Sys.13 pp. 958-985 (2014)

<a id="4">[4]</a> Galatolo S., Nisoli I. Rigorous computation of invariant measures and fractal dimension for maps with contracting fibers: 2D Lorenz-like maps  
Ergodic Theory and Dynamical Systems 36 (6), 1865-1891 (2016)

<a id="5">[5]</a> Viana M., Olivera K. Foundations of Ergodic Theory
Cambridge studies in advanced mathematics, Cambridge University Press 2016
