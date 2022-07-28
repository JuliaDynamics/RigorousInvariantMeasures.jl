# Userguide

The main objects involved in the approximation are the following:
1. A dynamic object
2. A basis object 

```jldoctest userguide
julia> using RigorousInvariantMeasures

julia> D0 = mod1_dynamic(x->4*x+0.5*x*(1-x), full_branch = true)
Piecewise-defined dynamic with 4 branches

julia> D = D0∘D0∘D0
RigorousInvariantMeasures.ComposedDynamic((RigorousInvariantMeasures.ComposedDynamic((Piecewise-defined dynamic with 4 branches, Piecewise-defined dynamic with 4 branches), Piecewise-defined dynamic with 16 branches), Piecewise-defined dynamic with 4 branches), Piecewise-defined dynamic with 64 branches)

julia> B = Hat(1024)
Hat{LinRange{Float64, Int64}}(range(0.0, stop=1.0, length=1025))
```

## Building the discretized operator
Once a basis is chosen we call [RigorousInvariantMeasures.DiscretizedOperator](@ref)
to compute the discretized operator.

Remark that a discretized operator can be of two types:
- [RigorousInvariantMeasures.IntegralPreservingDiscretizedOperator](@ref)
- [RigorousInvariantMeasures.NonIntegralPreservingDiscretizedOperator](@ref)

The type of the discretized operator is prescribed by the basis;
while an integral preserving operator is stored simply as a matrix,
a non integral preserving operator is stored as a triple
``(L, e, w)``; the operator ``L`` corresponds to the matrix,
while ``e, w`` are used to guarantee that the operator ``Q = L + e*w``
preserves the integral.

This is fundamental in our theory, since the rigorous estimate depends on the fact that the 
discretized operator preserves the space of average ``0`` functions.

```jldoctest userguide
julia> Q = DiscretizedOperator(B, D);
```


## Bounding the norms of the discretized operator

To compute our rigorous error bound we need to compute rigorously
upper bounds for the norms of the discretized operator restricted 
``||Q^k|_{U_0}||`` to the space of average ``0`` functions.

This is done through the use of [powernormbounds](@ref)

```jldoctest userguide
julia> norms = powernormbounds(B, D; Q=Q);
```

This function infers the strong and weak norm from the basis and uses a priori 
and a posteriori information to estimate these norms.

## Computing the invariant vector and the rigorous error
To compute the approximation and the error we use the functions 
[invariant_vector](@ref) and [distance_from_invariant](@ref).

Distance from invariant returns us an upper bound between 
``w`` and the density of the invariant a.c.i.m. with respect to the 
weak norm of the basis.

```jldoctest userguide
julia> w = invariant_vector(B, Q);

julia> error = distance_from_invariant(B, D, Q, w, norms)
0.000593413243533947

julia> weak_norm(B)
Linf
```

## Using the coarse-fine estimates
We can use now the coarse fine bounds to bound the norm of the powers
of a finer discretization.

The function [finepowernormbounds](@ref) uses the computed norms 
from the coarse discretization, the coefficients of the Doeblin-Fortet-Lasota-Yorke
inequality and a computed error bound on the norm of ``Q_f``. 

We first define a finer basis and compute the operator
```jldoctest userguide; filter = r".*" 
julia> B_fine = Hat(16384);

julia> Q_fine = DiscretizedOperator(B_fine, D);
```

Then, we call the coarse fine routines
```jldoctest userguide
julia> normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
1.0433672005775962

julia> norms_fine = finepowernormbounds(B, B_fine, D, norms; normQ_fine=normQ_fine);

julia> w_fine = invariant_vector(B_fine, Q_fine);
     
julia> error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)
7.574963167252225e-5
```
