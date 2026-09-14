abstract type NormKind end

@doc raw"""
    L1 <: NormKind

The ``L^1`` norm against Lebesgue measure ``dx`` on the domain. This is the
weak norm of the Ulam and hat bases. See [`L1μ`](@ref) for the counterpart
against the arcsine measure, which is what the Chebyshev bases use.
"""
struct L1 <: NormKind end

@doc raw"""
    L2 <: NormKind

The ``L^2`` norm against Lebesgue measure ``dx`` on the domain, and the
induced ``\ell^2`` operator norm on matrices. See [`L2μ`](@ref) for the
counterpart against the arcsine measure.
"""
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
@doc raw"""
    Aη(η) <: NormKind

Analytic norm on the **strip** of half-width ``η`` around the real axis, i.e.
``\|f\|_{A_η} = \sum_k |\hat f_k| e^{η|k|}`` on the Fourier coefficients.

Since ``|e^{2\pi i k z}| \le e^{2\pi |k| \operatorname{Im} z|}``, this weighted
``\ell^1`` norm dominates ``\sup |f|`` on the strip, so bounds stated for the
sup norm remain valid. It is the strong norm of [`FourierAnalytic`](@ref).

[`Eρ`](@ref) is the Chebyshev counterpart: a strip is the natural
neighbourhood of the circle, a Bernstein ellipse the natural neighbourhood of
``[-1,1]``.
"""
struct Aη <: NormKind
    η::Float64
end

@doc raw"""
    L2μ <: NormKind
    L1μ <: NormKind

``L^2`` and ``L^1`` against the arcsine measure ``dμ = dx/(π\sqrt{x(1-x)})`` on
``[0,1]``, as opposed to [`L2`](@ref) and [`L1`](@ref), which are against ``dx``.

These are the natural weak and auxiliary norms for the Chebyshev bases, because
the two facts the analysis rests on both live against ``μ``:

- the ``T_m`` are orthogonal in ``L^2(μ)``, so Parseval holds on the
  coefficients and the norm interface reduces to the same one-liners as the
  Fourier bases;
- the coefficients obey ``|\hat b_k| = |2\int f T_k\,dμ| \le 2\|f\|_{L^1(μ)}``,
  uniformly in ``k``. There is no such bound against ``L^1(dx)``: a width-``ε``
  sliver at an endpoint has ``\|f\|_{L^1(dx)} \sim ε`` but
  ``\|f\|_{L^1(μ)} \sim \sqrt{2ε}/π``.

Keeping both against ``μ`` also makes ``\|v\|_{L^1(μ)} \le \|v\|_{L^2(μ)}`` a
plain Cauchy–Schwarz with constant 1, since ``μ`` is a probability measure;
pairing ``L^1(dx)`` with ``L^2(dμ)`` instead costs a factor ``π/(2\sqrt2)``.

Use [`l2_measure_conversion_bounds`](@ref) to return to ``L^2(dx)``.
"""
struct L2μ <: NormKind end
struct L1μ <: NormKind end

# The docstring above covers both measures; attach it to this binding too, so
# that `[`L1μ`](@ref)` resolves.
@doc (@doc L2μ) L1μ

@doc raw"""
    Eρ(ρ) <: NormKind

Analytic norm on the **Bernstein ellipse** ``E_ρ`` — the image of the circle
``|w| = ρ`` (``ρ > 1``) under the Joukowski map ``z = (w + w^{-1})/2``, i.e. the
ellipse with foci ``\pm 1`` and semi-axes ``(ρ \pm ρ^{-1})/2``.

``\|f\|_{E_ρ} = \sum_k |\hat b_k| ρ^k`` on the Chebyshev coefficients, matching the
convention of [`Aη`](@ref), which is likewise a weighted ``\ell^1`` norm on the
Fourier coefficients. Since ``|T_k| \le (ρ^k + ρ^{-k})/2 \le ρ^k`` on ``E_ρ``,
this dominates ``\sup_{E_ρ}|f|``, so bounds stated for the sup norm (such as the
Trefethen projection estimates) remain valid. This is the Chebyshev counterpart of [`Aη`](@ref),
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
