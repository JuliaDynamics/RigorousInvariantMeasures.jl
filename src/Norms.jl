abstract type NormKind end
struct L1 <: NormKind end
struct L2 <: NormKind end
struct Linf <: NormKind end
struct Lipschitz <: NormKind end
struct TotalVariation <: NormKind end
struct ℓ1 <: NormKind end
struct ℓinf <: NormKind end



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
z_times_conjz(z::Complex) = square_round(abs_or_mag(real(z)), RoundUp) ⊕₊ square_round(abs_or_mag(imag(z)), RoundUp)
abs_or_mag(z::Complex) = sqrt_round(z_times_conjz(z), RoundUp)