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

# Analytic strip norm
struct Aη <: NormKind
    η::Float64
end

@doc raw"""
    Eρ(ρ) <: NormKind

Analytic norm on the **Bernstein ellipse** ``E_ρ`` — the image of the circle
``|w| = ρ`` (``ρ > 1``) under the Joukowski map ``z = (w + w^{-1})/2``, i.e. the
ellipse with foci ``\pm 1`` and semi-axes ``(ρ \pm ρ^{-1})/2``.

``\|f\|_{E_ρ} = \sup_{E_ρ} |f|``. This is the Chebyshev counterpart of [`Aη`](@ref),
which measures analyticity on a strip for the Fourier bases: a strip is the
natural neighbourhood of the circle, an ellipse the natural neighbourhood of
``[-1,1]``.

See [`bernstein_parameter`](@ref) and [`bernstein_expansion`](@ref) for the
rigorous test that a map expands such an ellipse.
"""
struct Eρ <: NormKind
    ρ::Float64
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
