using .RigorousInvariantMeasures: Dual
using .RigorousInvariantMeasures: NormKind, L1, L2, Aη, W, Cω, VarInterval, dfly, WeightedSobolev, circle_smoothness, strong_kind
import .RigorousInvariantMeasures: opnormbound, normbound, restrict_to_average_zero

export Fourier, evalFourier, FourierPoints, assemble_common, eval_on_dual, dual_nodes

using IntervalArithmetic

abstract type Fourier <: Basis end


FourierPoints(n, T) = [Interval{T}(i) / (n) for i = 0:n-1]

Base.lastindex(B::Fourier) = length(B)



function evalFourier(coeff, x)
    @assert length(coeff) % 2 == 1

    K = length(coeff)
    N = (K - 1) ÷ 2
    pos_coeff = coeff[1:N+1]
    neg_coeff = [typeof(x)(0); reverse(coeff[N+2:K])]

    return evalpoly(exp(2 * pi * im * x), pos_coeff) +
           evalpoly(exp(-2 * pi * im * x), neg_coeff)
end

"""
Make so that B[j] returns the basis function of coordinate j
"""
function Base.getindex(B::Fourier, i::Int)

    K = length(B)
    @boundscheck 1 <= i <= K || throw(BoundsError(B, i))

    N = (K - 1) ÷ 2

    if i == 1
        return x -> typeof(x)(1.0)
    end

    if 1 < i <= N + 1
        return x -> exp(2 * pi * im * (i - 1) * x)
    end

    if N + 2 <= i <= K
        L = i + (-2 * N - 2)
        return x -> exp(2 * pi * im * L * x)
    end
end

@doc raw"""
    is_refinement(Bf::Fourier, Bc::Fourier)

Whether `Bf` refines `Bc`. The Fourier spaces are nested frequency truncations,
``\mathcal U_{N_c} \subseteq \mathcal U_{N_f}`` iff ``N_c \le N_f``.

Note the argument order: the contract of [`is_refinement`](@ref) is
`(fine, coarse)`, as for `Ulam` and `HatNP`. This method used to be written
`(Bc, Bf)`, so a correct `is_refinement(fine, coarse)` call returned `false`
and `norms_of_powers_from_coarser_grid` logged "The fine basis is not a
refinement of the coarse basis" on every coarse–fine run. It is also now `≤`,
matching `Ulam`, so a basis refines itself.
"""
is_refinement(Bf::Fourier, Bc::Fourier) = length(Bc) <= length(Bf)
integral_covector(B::Fourier; T = Float64) = [Interval{T}(1); zeros(length(B) - 1)]'
one_vector(B::Fourier) = [1; zeros(length(B) - 1)]

Base.length(S::AverageZero{T}) where {T<:Fourier} = length(S.basis) - 1

function Base.iterate(S::AverageZero{T}, state = 1) where {T<:Fourier}
    B = S.basis
    i = state
    if i == length(B)
        return nothing
    end
    v = zeros(length(B))
    v[i+1] = 1
    return v, state + 1
end

###############################################################################
# Shared norm interface for all Fourier bases (dispatch on Fourier)
###############################################################################

# --- aux_weak_bound: ||v||_{L¹} ≤ M₂ · ||v||_{L²} ---
# On [0,1]: ||v||_{L¹} ≤ ||v||_{L²} by Cauchy-Schwarz → M₂ = 1
aux_weak_bound(B::Fourier) = 1.0

# --- weak_by_strong_and_aux_bound: ||v||_{L²} ≤ S₁·||v||_s + S₂·||v||_{L¹} ---
# Both Aη and W^{k,1}: ||v||_{L²} ≤ ||v||_{L∞} ≤ ||v||_{Aη} or ||v||_{W^{1,1}} ≤ ||v||_{W^{k,1}}
# So (S₁, S₂) = (1, 0); for WeightedSobolev see the method below.

@doc raw"""
    _fourier_sobolev_projection_error(B::Fourier) -> Float64

Bound on ``\|f - πf\|_{L^2}`` per unit of ``\|f^{(s)}\|_{L^1}``, for ``f`` in ``W^{s,1}``
of the circle, ``s`` the order of the strong norm (`W{s,l}`, or [`WeightedSobolev`](@ref)
with ``a_s = 1``), ``K`` = `B.k` the frequency cutoff and ``n`` = `length(B.p)` the number
of points.  It bounds the ``L^1`` error as well, ``dx`` being a probability measure.

The projection is the one the assembler applies: the DFT of the samples at ``j/n``,
truncated to ``|m| \le K``, so
``πf = \sum_{|m|\le K} \hat c^{(n)}_m e_m`` with ``\hat c^{(n)}_m = \sum_{l\in\mathbb Z}
\hat c_{m+ln}``.  Integrating by parts ``s`` times on the circle,
``|\hat c_m| \le V/(2π|m|)^s`` with ``V = \|f^{(s)}\|_{L^1}``.  Then, by Parseval,

```math
\|f-πf\|_{L^2}^2 = \sum_{|m|>K} |\hat c_m|^2 + \sum_{|m|\le K}
\Big|\sum_{l\ne0} \hat c_{m+ln}\Big|^2 .
```

The tail is at most ``2V^2(2π)^{-2s}\sum_{m>K} m^{-2s} \le 2V^2(2π)^{-2s}
K^{1-2s}/(2s-1)``.  For ``|m| \le K`` and ``l \ne 0``, ``|m+ln| \ge |l|n - K \ge
|l|(n-K)``, so the aliased sum is at most ``2ζ(s)V(2π)^{-s}(n-K)^{-s}`` with
``ζ(s) \le s/(s-1)``.  Hence

```math
\|f-πf\|_{L^2} \le \frac{V}{(2π)^s}\Big(\frac{2K^{1-2s}}{2s-1}
+ (2K+1)\Big(\frac{2s}{s-1}\Big)^2 (n-K)^{-2s}\Big)^{1/2},
```

evaluated in interval arithmetic.  Requires ``s \ge 2``, ``K \ge 1`` and ``n \ge 2K+1``.
"""
function _fourier_sobolev_projection_error(B::Fourier)
    s = _fourier_sobolev_order(B.strong)
    K, n = B.k, length(B.p)
    s >= 2 || throw(ArgumentError("the Fourier W^{s,1} projection error needs s ≥ 2, got $s"))
    K >= 1 || throw(ArgumentError("the Fourier projection error needs K ≥ 1"))
    n >= 2K + 1 || throw(ArgumentError("the Fourier projection error needs n ≥ 2K + 1"))
    si, Ki, ni = interval(s), interval(K), interval(n)
    tail = 2 * Ki^(1 - 2s) / (2si - 1)
    alias = (2Ki + 1) * (2si / (si - 1))^2 * (ni - Ki)^(-2s)
    return sup(sqrt(tail + alias) / (2 * interval(pi))^s)
end
_fourier_sobolev_order(::W{s,l}) where {s,l} = s
_fourier_sobolev_order(n::WeightedSobolev) = n.k

# WeightedSobolev strong norm ||f||_a = Σ_{j>=1} a_j ||f^(j)||_1 (a_s = 1). The seminorm has
# no ||f||_1 term, so it goes to the auxiliary norm:
#   ||v||_{L²} <= ||v||_∞ <= ||v||_1 + ||v'||_1 <= ||v||_1 + ||v||_a / a_1,
#   ||v^(j)||_{L¹} <= ||v^(j)||_{L²} <= (2πK)^j ||v||_{L²} for trigonometric polynomials
#   of degree K, so ||v||_a <= Σ_{j>=1} a_j (2πK)^j ||v||_{L²}.
function weak_by_strong_and_aux_bound(B::Fourier)
    B.strong isa WeightedSobolev || return (1.0, 0.0)
    a = B.strong.a
    a[1] > 0 || throw(ArgumentError("the weight a_1 must be positive, got $(a[1])"))
    return (1.0 ⊘₊ a[1], 1.0)
end

# --- bound_weak_norm_from_linalg_norm: ||v||_{L²} ≤ W₁·||v̂||_{ℓ¹} + W₂·||v̂||_{ℓ∞} ---
# Parseval: ||v||_{L²} = ||v̂||_{ℓ²} ≤ ||v̂||_{ℓ¹} → (1, 0)
bound_weak_norm_from_linalg_norm(B::Fourier) = (1.0, 0.0)

# --- bound_linalg_norm_L1_from_weak: ||v̂||_{ℓ¹} ≤ A · ||v||_{L²} ---
# Cauchy-Schwarz: ||v̂||_{ℓ¹} ≤ √n · ||v̂||_{ℓ²} = √n · ||v||_{L²}
bound_linalg_norm_L1_from_weak(B::Fourier) =
    sqrt_round(Float64(length(B), RoundUp), RoundUp)

# --- bound_linalg_norm_L∞_from_weak: ||v̂||_{ℓ∞} ≤ A · ||v||_{L²} ---
# ||v̂||_{ℓ∞} ≤ ||v̂||_{ℓ²} = ||v||_{L²} → A = 1
bound_linalg_norm_L∞_from_weak(B::Fourier) = 1.0

# --- opnormbound and normbound: delegate to NormBounds.jl L2 implementations ---
# The 2-argument opnormbound(::Type{L2}, M) is defined in NormBounds.jl
opnormbound(B::Fourier, ::Type{L2}, M::AbstractVecOrMat{S}) where {S} =
    opnormbound(L2, M)
normbound(B::Fourier, ::Type{L2}, v) = normbound(L2, v)

function _fourier_circle_check(B::Fourier, D)
    B.strong isa Union{W,WeightedSobolev} || return nothing
    s = _fourier_sobolev_order(B.strong)
    s >= 2 || return nothing
    r = circle_smoothness(D)
    (r !== nothing && r >= s) || throw(
        ArgumentError(
            "a Fourier basis of Sobolev order $s needs a map asserted C^$s on the circle " *
            "(CircleMap{r}, r ≥ $s); circle_smoothness(D) = $r",
        ),
    )
    return nothing
end

# --- invariant_measure_strong_norm_bound: standard DFLY pattern B/(1-A) ---
#
# For a strong norm of order s >= 2 (W{s,l} or WeightedSobolev) the higher-order DFLY of
# the Symbolics extension bounds ||(Lf)^(j)||_1 on (0,1); on the circle Lf must also be in
# W^{s,1} of the circle, i.e. Lf, ..., (Lf)^(s-1) must match across 0 ~ 1. That holds when
# T is C^s on the circle, which the caller asserts with a CircleMap{r}, r >= s. The
# coefficients must also be stated in the basis's own norm.
function invariant_measure_strong_norm_bound(
    B::Fourier,
    D::Dynamic;
    dfly_coefficients = nothing,
)
    _fourier_circle_check(B, D)
    dfly_coefficients === nothing &&
        (dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D))
    c = strong_kind(dfly_coefficients)
    (B.strong isa WeightedSobolev && c != B.strong) && throw(
        ArgumentError("the DFLY coefficients are stated in $c, not in the basis's $(B.strong)"),
    )
    (B.strong isa W && c isa WeightedSobolev) && throw(
        ArgumentError(
            "the DFLY coefficients are stated in $c, not in $(B.strong); build the basis " *
            "on that norm",
        ),
    )
    A, Bcoeff = dfly_coefficients
    @assert A < 1.0
    return Bcoeff ⊘₊ (1.0 ⊖₋ A)
end

# --- bound_weak_norm_abstract: bound ||L||_{L²→L²} ---
# For Markov transfer operators: ||L||_{L¹} ≤ 1, ||L||_{L∞} ≤ 1
# By Riesz-Thorin: ||L||_{L²} ≤ 1
# For general maps: use dfly_coefficients[2] ⊕₊ 1.0 as Hat does
function bound_weak_norm_abstract(B::Fourier, D = nothing; dfly_coefficients = nothing)
    if dfly_coefficients === nothing
        D isa Dynamic && _fourier_circle_check(B, D)
        dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D)
    end
    return dfly_coefficients[2] ⊕₊ 1.0
end

# --- restrict_to_average_zero for Fourier: U^0 = span{e₂,...,eₙ} ---
# For Fourier, integral_covector = [1, 0, ..., 0], so restriction is just BM[2:end, 2:end]
restrict_to_average_zero(B::Fourier, BM, f) = BM[2:end, 2:end]

# --- weak_dual_norm_bound and integral_pairing ---
#
# Weak norm of the Fourier bases is L²; dual is L² (Hilbert). Parseval gives
# ‖φ_v‖_{L²} = ‖v‖_{ℓ²} = √(Σ |v_n|²).

import ..RigorousInvariantMeasures: weak_dual_norm_bound, integral_pairing,
    Observable

function weak_dual_norm_bound(B::Fourier, v::AbstractVector)
    return sup(sqrt(sum(abs2(c) for c in v)))
end

@doc raw"""
    integral_pairing(ϕ::Observable{<:Fourier}, ρ, ρ_w_error;
                     ρ_dual_weak_bound = weak_dual_norm_bound(ϕ.B, ρ))

For Fourier bases (weak `L²`), the pairing
``\int_0^1 ϕ_N(x)\,ρ_N(x)\,dx`` equals ``\sum_n \hat ϕ_n\,\overline{\hat ρ_n}``.
The result is real-valued for real signals; we return the real part of the
sum.
"""
function integral_pairing(
    ϕ::Observable{<:Fourier},
    ρ::AbstractVector,
    ρ_w_error;
    ρ_dual_weak_bound = weak_dual_norm_bound(ϕ.B, ρ),
)
    @assert length(ϕ.v) == length(ρ)
    # Lift ρ to interval arithmetic so the inner product encloses the
    # floating-point rounding from the sum and the per-term multiplications.
    ρi = _lift_to_interval(ρ)
    pairing = sum(ϕ.v[i] * conj(ρi[i]) for i in eachindex(ϕ.v))
    val = real(pairing)
    err_proj =
        ϕ.proj_error === nothing ? 0.0 :
        sup(ϕ.proj_error * ρ_dual_weak_bound)
    err_density = sup(ϕ.weak_dual_bound) * ρ_w_error
    err = err_proj + err_density
    return val + interval(-err, err)
end

# --- Fourier multiplication (convolution) ---

# Direct convolution of two N-mode Fourier coefficient vectors in FFT layout
# (DC at index 1, then positive freqs 1..k_freq, then negative freqs -k_freq..-1).
# Output: length 2N-1 vector in NATURAL ordering, freqs -2k_freq..2k_freq, where
# index l + 2k_freq + 1 holds frequency l.
function _fourier_direct_convolution(v1::AbstractVector, v2::AbstractVector, k_freq::Integer)
    N = 2 * k_freq + 1
    @assert length(v1) == length(v2) == N
    out = zeros(eltype(v1), 2N - 1)
    fft_to_freq(i) = (i ≤ k_freq + 1) ? (i - 1) : (i - 1 - N)
    @inbounds for i_m = 1:N
        m = fft_to_freq(i_m)
        for i_n = 1:N
            n = fft_to_freq(i_n)
            l = m + n
            out[l+2*k_freq+1] += v1[i_m] * v2[i_n]
        end
    end
    return out
end

# Reorder a natural-ordered length-(2k_freq+1) vector (freqs -k_freq..k_freq) into
# FFT layout (freqs 0,1,...,k_freq,-k_freq,...,-1).
function _natural_to_fft_order(v_natural::AbstractVector, k_freq::Integer)
    N = 2 * k_freq + 1
    @assert length(v_natural) == N
    v_fft = similar(v_natural, N)
    v_fft[1] = v_natural[k_freq+1]
    for m = 1:k_freq
        v_fft[m+1] = v_natural[k_freq+1+m]
        v_fft[N-m+1] = v_natural[k_freq+1-m]
    end
    return v_fft
end

@doc raw"""
    p1 * p2  for two `ProjectedFunction{<:Fourier}` on the same basis

Multiply two Fourier projections by treating each operand as a
trigonometric polynomial (the discrete coefficient vector) plus an
opaque L²-error bound. All bounds are derived from the finite
coefficient vectors alone — provenance (the original W^{k,1} seminorm,
function class, …) is *not* used. This keeps multiplication composable:
the result is again a trigonometric polynomial with an L² error bound,
ready to be passed to `*` or `+` again.

The arithmetic:

1. Convolve the two N-mode coefficient vectors → 2N−1 modes.
2. Truncate the central N modes back into the basis layout.
3. Parseval gives ``\|\text{discarded modes}\|_{\ell^2}`` directly from
   the convolution output.

The bounds:

- ``\|fg\|_{L^2} \leq \|\hat f\|_{\ell^1} \cdot \|g\|_{L^2}`` (Young's
  inequality, with ``\|\hat f\|_{\ell^1}`` the surrogate for
  ``\|f\|_{L^\infty}``).
- ``\|fg - π_N(φ_f φ_g)\|_{L^2} \leq \|\hat f\|_{\ell^1}\,p_2.\text{proj\_error}
  + \|\hat g\|_{\ell^1}\,p_1.\text{proj\_error} + \|\text{discarded modes}\|_{\ell^2}``.

Each input's `proj_error` is the L² distance from the original function
to its trigonometric-polynomial approximant; the multiplication treats
that distance as the only thing it knows.
"""
function Base.:*(
    p1::ProjectedFunction{<:Fourier},
    p2::ProjectedFunction{<:Fourier},
)
    @assert length(p1.v) == length(p2.v) "Fourier multiplication requires same length"
    @assert p1.B.k == p2.B.k "Fourier multiplication requires same frequency cutoff"
    k_freq = p1.B.k
    N = 2 * k_freq + 1

    conv_natural = _fourier_direct_convolution(p1.v, p2.v, k_freq)
    # Kept: freqs in [-k_freq, k_freq] → natural indices [k_freq+1 : 3k_freq+1].
    kept_natural = conv_natural[k_freq+1:3*k_freq+1]
    # Discarded: low [1:k_freq] (freqs -2k_freq..-k_freq-1) and
    #            high [3k_freq+2:4k_freq+1] (freqs k_freq+1..2k_freq).
    disc_low = conv_natural[1:k_freq]
    disc_high = conv_natural[3*k_freq+2:4*k_freq+1]
    v_new = _natural_to_fft_order(kept_natural, k_freq)

    # ℓ² of discarded modes (Parseval) — fully rigorous via interval sum.
    l2_disc_sq = sum(abs2, disc_low) + sum(abs2, disc_high)
    L2_discarded = sup(sqrt(l2_disc_sq))

    # ℓ¹ of input coefficient vectors as a (loose) surrogate for ‖φ_·‖_{L^∞}.
    f_linf = sup(sum(abs, p1.v))
    g_linf = sup(sum(abs, p2.v))

    # weak_dual_bound: Young's gives ‖fg‖_{L²} ≤ ‖f̂‖_{ℓ¹} · ‖g‖_{L²}.
    weak_dual_new = _mul_bound(f_linf, p2.weak_dual_bound)

    # proj_error: Hölder + Parseval discarded.
    proj_err_new = _add_bound(
        _add_bound(
            _mul_bound(p1.proj_error, g_linf),
            _mul_bound(f_linf, p2.proj_error),
        ),
        L2_discarded,
    )

    return ProjectedFunction(p1.B, v_new, weak_dual_new, proj_err_new)
end

abstract type FourierDual <: Dual end

function eval_on_dual(B::Fourier, computed_dual::FourierDual, ϕ) end

@doc raw"""
    dual_nodes(B::Fourier, computed_dual) -> (x, labels, weights)

Flat description of a computed dual, used by the fast assembler.

`eval_on_dual(B, computed_dual, ϕ_m)` for the basis function
``ϕ_m(x) = e^{2πimx}`` is, for every Fourier basis, of the form

```math
w^{(m)}[\ell] = \sum_{j\,:\,\texttt{labels}[j] = \ell}
                \texttt{weights}[j]\; e^{2πi m x_j},
```

with **real** weights. Exposing the three vectors lets `assemble_common`
sweep all frequencies at once instead of re-evaluating a rigorous complex
exponential for every (basis function, node) pair; see the comment there.

Weights being real is what makes the ``m \mapsto -m`` conjugate symmetry
valid, so any new `Fourier` dual must preserve that.
"""
function dual_nodes end

# `assemble_common(::Fourier, D; …)` lives in the FFTWExt extension; loading
# `using FFTW` makes it available. Without FFTW loaded, callers will hit a
# MethodError on this name.
function assemble_common end
