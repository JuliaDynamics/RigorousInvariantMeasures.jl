# InvariantMeasures

[![Build Status](https://travis-ci.com/orkolorko/InvariantMeasures.jl.svg?branch=master)](https://travis-ci.com/orkolorko/InvariantMeasures.jl)

This Package provides methods for the rigorous approximation of Absolutely Continuous Invariant Measures for one dimensional dynamical systems,
using the results in [[2]](#2) and subsequent works.

## Mathematical background
By Birkhoff Ergodic Theorem we know that if a dynamical system T admits an ergodic invariant measure μ, for μ-almost every point x the frequence of the visits to a set E corresponds to the measure of the set with respect to μ.

Therefore, being able to approximate invariant measures with a large basin is interesting to investigate the statistical properties of the dynamical system T. [[4]](#4)

The existence of absolutely continuous invariant measures for one dimensional maps is a delicate topic; with this package we present approximation schemes for the invariant measures of system that satisfy a Lasota-Yorke inequality through the use of a coarse-fine scheme and a posteriori estimates on the mixing time, i.e., this means that the mixing time is estimated by our algorithm and we do not need an a priori estimate, which is usually difficult to obtain.

The Ulam approximation schemes works under relatively weak hypothesis on the dynamics and was used in [[3]](#3) to approximate the invariant measure for the geometric Lorenz 1-dimensional map.

We are currently working on the implementation of the Ulam scheme for system with additive uniform noise, as the one used in
[[1]](#1)

## Caveat
The function sinpi is currently only working using the CRlibm.jl
package, so, be careful when defining the dynamics. The example below, indeed, needs this function to work.

## Basic Usage
Examples of usage are present in the directory examples.

```julia
using InvariantMeasures
D = Mod1Dynamic(x -> 4*x + 0.01*InvariantMeasures.sinpi(8*x))
B = Ulam(1024)
Q = DiscretizedOperator(B, D)
```

The code snippet above defines a dynamic D(x) = 4x+0.01 sin(8πx),
a basis B associated to the Ulam discretization on a partition of 1024 homogenous interval and computes the discretized operator Q, a Markov chain whose entries are P(T(x)∈ Iᵢ | x ∈ Iⱼ).

Note the usage of `InvariantMeasures.sinpi(8*x)` rather than `Base.sinpi` or `Base.sin(8\pi*x)`. This detail is required to ensure that D(8) = 4 exactly.

```julia
norms = norms_of_powers(weak_norm(B), m, Q, integral_covector(B))
```

the function norm of powers computes the L¹ norm of Q when restricted to the space of average 0 vectors. This gives us the a posteriori estimate for the mixing time of the Markov chain and is used in our
rigorous estimate.

```julia
w = invariant_vector(B, Q)
distance_from_invariant(B, D, Q, w, norms)
```
the computed approximation of the invariant measure of D is stored in w, and `distance_from_invariant` computes an upper bound for the L¹ distance between w and the density of the absolutely continuous invariant measure of the system.

Inside the examples it is showed how to use the coarse-fine scheme to obtain better L¹ bounds and reduce the computational time.

## References
<a id="1">[1]</a>
Galatolo S., Monge M., Nisoli I., Existence of noise induced order, a computer aided proof Nonlinearity 33 (9), 4237 (2020)

<a id="2">[2]</a>
Galatolo S., Nisoli I., An elementary approach to rigorous approximation of invariant
measures SIAM J. Appl Dyn Sys.13 pp. 958-985 (2014)

<a id="3">[3]</a> Galatolo S., Nisoli I. Rigorous computation of invariant measures and fractal dimension for maps with contracting fibers: 2D Lorenz-like maps  
Ergodic Theory and Dynamical Systems 36 (6), 1865-1891 (2016)

<a id="4">[4]</a> Viana M., Olivera K. Foundations of Ergodic Theory
Cambridge studies in advanced mathematics, Cambridge University Press 2016
