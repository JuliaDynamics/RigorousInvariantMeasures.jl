abstract type NormKind end
struct L1 <: NormKind end
struct L2 <: NormKind end
struct Linf <: NormKind end
struct Lipschitz <: NormKind end
struct TotalVariation <: NormKind end
struct ℓ1 <: NormKind end
struct ℓinf <: NormKind end

# Sobolev W^{k,l} norm
struct C1 <: NormKind end
struct W{k,l} <: NormKind end
order(::Type{W{k,l}}) where {k,l} = k
regularity(::Type{W{k,l}}) where {k,l} = l

@doc raw"""
    L2μ <: NormKind

``L^2`` against the arcsine measure ``dμ = dx/(π\sqrt{x(1-x)})`` on ``[0,1]``,
as opposed to [`L2`](@ref), which is ``L^2(dx)``.

The distinction matters for the Chebyshev basis: the ``T_m`` are orthogonal in
``L^2(μ)`` and not in ``L^2(dx)``, so under `L2μ` Parseval holds on the
coefficients and the norm interface reduces to the same one-liners as the
Fourier bases, whereas under `L2` the (dense) Gram matrix enters. The two are
related on polynomials of degree < n by
[`l2_measure_conversion_bounds`](@ref).
"""
struct L2μ <: NormKind end

# Analytic strip norm
struct Aη <: NormKind
    η::Float64
end

# Adjoint analytic norm
struct Cω <: NormKind end



"""
Functions to deal with various types of norms and seminorms
"""

using FastRounding
using IntervalArithmetic, IntervalOptimisation
using TaylorSeries: Taylor1
using SparseArrays: getcolptr


"""
'Absolute value' definition that returns mag(I) for an interval and abs(x) for a real
"""
abs_or_mag(x::Number) = Float64(abs(x), RoundUp)
abs_or_mag(x::Interval) = Float64(mag(x), RoundUp)

"""
Computes a rigorous upper bound for z*z'
"""
z_times_conjz(z::Complex) =
    square_round(abs_or_mag(real(z)), RoundUp) ⊕₊ square_round(abs_or_mag(imag(z)), RoundUp)
abs_or_mag(z::Complex) = sqrt_round(z_times_conjz(z), RoundUp)
