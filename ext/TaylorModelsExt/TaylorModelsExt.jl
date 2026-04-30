module TaylorModelsExt

export discretizationlogder

using RigorousInvariantMeasures
using IntervalArithmetic
using IntervalOptimisation: maximise

import RigorousInvariantMeasures: Observable, ProjectedFunction, integrateobservable

import TaylorModels

function integrate(f, I; steps = 1024, degree = 6)
    lo = inf(I)
    hi = sup(I)
    l = 2 * radius(I)
    int_center = interval(0.0)
    int_error = interval(0.0)
    for i = 1:steps
        left = lo + (i - 1) * (l / steps)
        right = lo + i * l / steps
        r = (1 / 2) * (l / steps)
        J = interval(left, right)
        TM = TaylorModels.TaylorModel1(degree, J)
        FM = f(TM)
        #@info FM
        for i = 0:Int64(floor(degree / 2))
            int_center += 2 * (FM.pol[2*i] * r^(2 * i + 1)) / (2 * i + 1)
            int_error += 2 * FM.rem * r
        end
    end
    return int_center + int_error
end

function adaptive_integration(f, I::Interval; tol = 2^-10, steps = 8, degree = 6) # tol 2^-10, steps = 8 are default values
    lo = inf(I)
    hi = sup(I)
    l = 2 * radius(I)
    int_value = interval(0.0)
    for i = 1:steps
        left = lo + (i - 1) * (l / steps)
        right = lo + i * l / steps
        Istep = interval(left, right)
        val = integrate(f, Istep)
        if radius(val) < tol
            int_value += val
        else
            I₁, I₂ = bisect(I)
            val₁ = adaptive_integration(
                f,
                I₁;
                tol = tol / 2,
                steps = steps,
                degree = degree + 2,
            )
            val₂ = adaptive_integration(
                f,
                I₂;
                tol = tol / 2,
                steps = steps,
                degree = degree + 2,
            )
            int_value += val₁ + val₂
        end
    end
    return int_value
end

# `Observable` and `ProjectedFunction` structs are defined in the main package
# (src/Observables.jl). This extension supplies the Ulam-specific constructors
# below.

### TODO: Actually some assumptions are made, as the fact that
# the Ulam base is equispaced

import TaylorSeries
"""
    discretizationlogder(B::Ulam, D::PwMap; degree = 7)

Compute the discretization of the logarithm of the derivative 
of the dynamics 'D' on the Ulam basis 'B', using a Taylor series 
expansion of degree 'degree' 
"""
function discretizationlogder(B::Ulam, D::PwMap; degree = 7)
    v = zeros(Interval{Float64}, length(B))
    # we first compute the indexes in the Ulam approximation 
    # of the endpoints of the branches 
    infbound = emptyinterval()

    for (i, br) in enumerate(D.branches)
        ind_X1 = Int64(floor(length(B) * inf(br.X[1]))) + 1
        ind_X2 = min(Int64(floor(length(B) * inf(br.X[2]))) + 1, length(B))
        dom = hull(br.X[1], br.X[2])
        fprime = derivative(br.f)
        for i = ind_X1:ind_X2
            I = intersect_interval(interval(B.p[i], B.p[i+1]), dom)
            r = interval(radius(I))
            Tmid = TaylorSeries.Taylor1([interval(mid(I)), interval(1)], degree)
            Tint = TaylorSeries.Taylor1([I, interval(1)], degree)
            Fmid = log(abs(fprime(Tmid)))
            Fint = log(abs(fprime(Tint)))
            infbound = hull(infbound, abs(Fint[0]))
            ϵ = mag(Fint[degree] - Fmid[degree])

            for k = 0:Int64(floor(Float64(degree) / 2))
                v[i] += 2 * (Fmid[2*k] * r^(2 * k + 1)) / (2 * k + 1)
            end
            v[i] += interval(-ϵ, ϵ) * r^(degree + 1) / (degree + 1)
        end
        v[ind_X1] += 2 * interval(radius(br.X[1])) * abs(log(fprime(br.X[1])))
        v[ind_X2] += 2 * interval(radius(br.X[2])) * abs(log(fprime(br.X[2])))
    end

    v *= length(B)
    return ProjectedFunction(B, v, infbound, nothing)
end

function discretizationlogder_fast(B, D::PwMap)
    v = zeros(Interval{Float64}, length(B))

    @error "Not implemented yet!"
end

"""
    VariationBound(f; steps = 1024)

Rigorous variation bound for a C¹ function on `[0,1]`, computed as
``\\int_0^1 |f'(x)|\\,dx`` via a Riemann upper bound over a uniform
partition with `steps` subintervals. `f'` is obtained from Taylor-series
automatic differentiation, so `f` must be callable on a
`TaylorSeries.Taylor1{Interval{Float64}}`.

Equivalent to `WklSeminorm(f; k = 1, l = 1, steps = steps)`.
"""
function VariationBound(f; steps = 1024)
    total = interval(0.0)
    h = interval(1.0) / steps
    for i = 1:steps
        I = interval((i - 1) / steps, i / steps)
        Tx = TaylorSeries.Taylor1([I, interval(1.0)], 1)
        fprime_I = f(Tx)[1]
        total += abs(fprime_I) * h
    end
    return total
end

@doc raw"""
    WklSeminorm(f; k = 1, l = 1, steps = 1024)

Rigorous bound on the Sobolev seminorm ``\|f^{(k)}\|_{L^l}`` for a
function ``f: [0,1] \to ℝ``, computed as a Riemann-style upper bound
over a uniform partition with `steps` panels. For `k = 1, l = 1` this
returns the total variation and matches [`VariationBound`](@ref).

Algorithm: on each panel ``I_i`` of width ``h = 1/\text{steps}``, expand
`f` as an interval Taylor series of order `k`; the k-th Taylor
coefficient times ``k!`` is an interval enclosure of ``f^{(k)}`` over
``I_i``. The Riemann-box upper bound is

```math
\|f^{(k)}\|_{L^l}^l \;\leq\; \sum_i |f^{(k)}(I_i)|^l \cdot h,
```

then take the ``l``-th root. Result is an interval enclosure of an
upper bound on the seminorm.

`f` must be callable on `TaylorSeries.Taylor1{Interval{Float64}}`.
Typical use: pass `Wk1_seminorm = WklSeminorm(f; k = 2, l = 1)` to the
Fourier `ProjectedFunction` constructor for a `W^{2,1}` basis.
"""
function WklSeminorm(f; k::Integer = 1, l::Integer = 1, steps::Integer = 1024)
    k ≥ 1 || throw(ArgumentError("WklSeminorm requires k ≥ 1; got k = $k"))
    l ≥ 1 || throw(ArgumentError("WklSeminorm requires l ≥ 1; got l = $l"))
    fact_k = interval(Float64(factorial(k)))
    total = interval(0.0)
    h = interval(1.0) / steps
    for i = 1:steps
        I = interval((i - 1) / steps, i / steps)
        Tx = TaylorSeries.Taylor1([I, interval(1.0)], k)
        f_kth = f(Tx)[k] * fact_k
        contrib = l == 1 ? abs(f_kth) : abs(f_kth)^l
        total += contrib * h
    end
    return l == 1 ? total : total^(interval(1.0) / interval(Float64(l)))
end


@doc raw"""
    ProjectedFunction(B::Ulam, f::Function;
                      tol = 2^-10,
                      var_bound = VariationBound(f),
                      weak_dual_bound = maximise(x -> abs(f(x)), interval(0,1))[1])

Discretize `f` on the Ulam basis. Computes both:

- `weak_dual_bound`: a Taylor-model-driven bound on ``\|f\|_{L^∞}``
  (the dual of the weak `L¹` norm). Auto-computed via
  `IntervalOptimisation.maximise`; can be overridden if a tighter
  user-known bound is available.
- `proj_error = var_bound / length(B)`: the L¹ projection error
  ``\|f - π_N f\|_{L^1} \leq \mathrm{Var}(f)/N``. `var_bound` defaults
  to the Riemann TV bound from `VariationBound`.

The discrete coefficient vector `v[i] = N · ∫_{I_i} f dx` matches the
old separate `Observable`/`ProjectedFunction` constructors.
"""
function ProjectedFunction(
    B::Ulam,
    f::Function;
    tol = 2^-10,
    var_bound = VariationBound(f),
    weak_dual_bound = maximise(x -> abs(f(x)), interval(0, 1))[1],
)
    v = zeros(Interval{Float64}, length(B))
    for i = 1:length(B)
        I = interval(B.p[i], B.p[i+1])
        v[i] = adaptive_integration(f, I; tol = tol, steps = 1, degree = 2) * length(B)
    end
    proj_error = var_bound / length(B)
    return ProjectedFunction(B, v, weak_dual_bound, proj_error)
end

RigorousInvariantMeasures.projection(B::Ulam, f::Function; kwargs...) =
    ProjectedFunction(B, f; kwargs...)




#= function discretizationlogder(B, D::PwMap; degree = 7)
    v = zeros(Interval{Float64}, length(B))
    infbound  = emptyinterval()
    endpoints = D.endpoints
    delicate_indexes = Int64.(floor.(length(B)*[inf(x) for x in endpoints])).+1
    for i in 1:(length(delicate_indexes)-1)
        for j in delicate_indexes[i]:delicate_indexes[i+1]-1
            I = interval(B.p[j], B.p[j+1])
            r = interval(radius(I))
            Tmid = TaylorSeries.Taylor1([interval(mid(I)), interval(1)], degree)
            Tint = TaylorSeries.Taylor1([I, interval(1)], degree)
            Fmid = log(TaylorModels.derivative(D.Ts[i](Tmid)))
            Fint = log(TaylorModels.derivative(D.Ts[i](Tint)))   
            infbound = hull(infbound, abs(Fint[0]))
            ϵ = mag(Fint[degree]-Fmid[degree]) 

            for k in 0:Int64(floor(Float64(degree)/2))
                v[j]+=2*(Fmid[2*k]*r^(2*k+1))/(2*k+1)            
            end
            v[j]+=interval(-ϵ, ϵ)*r^(degree+1)/(degree+1)
        end
    end

    #correction since the endpoints may be wide intervals
    for i in 2:length(endpoints)-1
        x = endpoints[i]
        Tx = TaylorSeries.Taylor1([x, interval(1)], 1)

        corr = 2*interval(radius(x))*(abs(log(TaylorModels.derivative(D.Ts[i-1](Tx))))
                +abs(log(TaylorModels.derivative(D.Ts[i](Tx)))))[0]
        v[delicate_indexes[i]]+=corr
    end
    v*=length(B)
    return Observable(B, v, infbound)
end =#

function integrateobservable(B::Ulam, ϕ::ProjectedFunction, f::Vector, error)
    val = (ϕ.v)' * f
    return val / length(B) + sup(ϕ.weak_dual_bound) * interval(-error, error)
end



end
