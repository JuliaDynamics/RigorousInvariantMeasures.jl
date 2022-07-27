# Userguide

The main objects involved in the approximation are the following:
1. A dynamic object
2. A basis object 

```julia
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

```julia
julia> Q = DiscretizedOperator(B, D);

julia> Q.e
1024-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 ⋮
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0

julia> Q.w
1×1024 adjoint(::Vector{Interval{Float64}}) with eltype Interval{Float64}:
 [6.01662e-06, 6.01663e-06]  …  [-2.08295e-08, -2.08294e-08]
```


## Bounding the norms of the discretized operator

To compute our rigorous error bound we need to compute rigorously
upper bounds for the norms of the discretized operator restricted 
``||Q^k|_{U_0}||`` to the space of average ``0`` functions.

This is done through the use of [powernormbounds](@ref)

```julia
julia> norms = powernormbounds(B, D; Q=Q)
16-element Vector{Float64}:
 1.0433224088954731
 0.18200605351394591
 0.0017499771147111234
 2.0992469139922695e-5
 1.4111467136923114e-7
 4.001217780270299e-9
 1.2902427238537547e-9
 2.469474454461355e-10
 2.9623453839089177e-12
 1.9913350475646103e-14
 5.646305321395677e-16
 1.6009743725151178e-17
 5.162542127548024e-18
 4.1803039533246575e-19
 2.8100659082311227e-21
 7.96776519879092e-23
```

This function infers the strong and weak norm from the basis and uses a priori 
and a posteriori information to estimate these norms.

## Computing the invariant vector and the rigorous error
To compute the approximation and the error we use the functions 
[invariant_vector](@ref) and [distance_from_invariant](@ref).

Distance from invariant returns us an upper bound between 
``w`` and the density of the invariant a.c.i.m. with respect to the 
weak norm of the basis.

```jldoctests userguide
julia> w = invariant_vector(B, Q)
1024-element Vector{Float64}:
 0.9589595620450252
 0.9590346071917383
 0.9591096666147971
 0.9591847403127435
 ⋮
 1.0420701699276038
 1.0421170085206535
 1.042163852287063

julia> error = distance_from_invariant(B, D, Q, w, norms)
0.0010164506833733997

julia> weak_norm(B)
Linf
```

## Using the coarse-fine estimates
We can use now the coarse fine bounds to bound the norm of the powers
of a finer discretization.

The function [finepowernormbounds](@ref) uses the computed norms 
from the coarse discretization, the coefficients of the Doeblin-Fortet-Lasota-Yorke
inequality and a computed error bound on the norm of ``Q_f``. 

```@meta
DocTestFilters = r"Computing preimages and derivatives...\(.*\)"
```


```julia
julia> B_fine = Hat(16384);

julia> Q_fine = DiscretizedOperator(B_fine, D);
Computing preimages and derivatives... 100%|████████████████████████████████████| Time: 0:00:03

julia> normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
1.0433672005775962

julia> norms_fine = finepowernormbounds(B, B_fine, D, norms; normQ_fine=normQ_fine)
16-element Vector{Float64}:
   1.0433672005775962
   1.08861511524113
   1.1358253052955953
   0.29514500875596494
   0.010963157005008568
   0.0012921514468915017
   0.0009827177377044122
   0.0009728796288068277
   0.0009725396338062495
   0.00012019081151646843
   1.4166059186520524e-5
   1.6696553617038016e-6
   1.2698201466606996e-6
   9.65734151998878e-7
   9.560660678797541e-7
   1.5530473100406846e-7

julia> w_fine = invariant_vector(B_fine, Q_fine)
  16384-element Vector{Float64}:
   0.958867661434263
   0.9588723523398917
   0.9588770433012951
   0.9588817343184616
   0.9588864253913828
   0.9588911165200567
   0.9588958077045152
   0.9589004989447679
   0.9589051902407743
   0.9589098815925523
   ⋮
   1.0424681996544376
   1.0424286590509482
   1.0423891175997508
   1.042349575300802
   1.0423100321540628
   1.0422704881595075
   1.042230943317113
   1.0421913976268897
   1.0421518510887897
   
julia> error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)
0.0001305956645536314
```

```@meta
DocTestFilters = nothing
```