###############################################################################
# Lasota-Yorke inequalities for the analytic norms
#
# For the analytic norms the contraction does not come from expansion of the
# map, as it does for BV or Sobolev, but from the transfer operator enlarging
# the domain of analyticity: L sends functions analytic on a strip (or ellipse)
# to functions analytic on a *strictly larger* one. That gain is invisible in
# L¹ or L², so it cannot be read off a real-variable estimate — it has to come
# from a complex evaluation, which is what `annulus_expansion` and
# `bernstein_expansion` provide.
#
# Two facts turn the gain into a two-norm inequality:
#
#   (i)  the coefficients of the *image* inherit the decay of the larger domain,
#            |ĉ_k(Lf)| ≤ ‖Lf‖_{larger} · (decay)^{|k|},
#   (ii) each coefficient is separately controlled by an L¹-type norm,
#            |ĉ_k(f)| ≤ const · ‖f‖_{L¹} ,
#        with no dependence on k.
#
# Splitting the series at |k| = K and using (i) above K, (ii) below, gives
#
#     ‖Lf‖_strong ≤ A(K)·‖f‖_strong + B(K)·‖f‖_aux ,
#
# where A(K) → 0 geometrically and B(K) grows geometrically. So any contraction
# factor can be bought, at a computable price in B. Letting K → ∞ formally is
# the degenerate case A = 0, the same shape as the noise DFLY
# (`dfly(VarInterval, L1, ::DiscretizedNoiseKernelUlam)` returns `(0.0, …)`),
# where the smoothing plays the role of the analytic gain.
###############################################################################

export analytic_dfly, analytic_dfly_choose_K, annulus_expansion, strip_expansion

# Sum of the tail weights above K: 2 Σ_{j>K} q^j = 2 q^{K+1}/(1-q), q < 1.
function _tail_factor(q, K::Integer)
    qi = interval(q)
    @assert sup(qi) < 1 "the analytic gain must be strictly positive (q < 1)"
    return 2 * qi^(K + 1) / (1 - qi)
end

# Sum of the head weights: 2 Σ_{j=1}^K w^j = 2w(w^K - 1)/(w - 1), w = e^{2πη} > 1.
_strip_head(w, K::Integer) = 2 * w * (w^K - 1) / (w - 1)

@doc raw"""
    analytic_dfly(strong::Aη, η′, C, K) -> (A, B)

Lasota–Yorke constants for the Fourier analytic norm
``\|f\|_{A_η} = \sum_k |\hat c_k| e^{2πη|k|}``, with auxiliary norm ``L^1(dx)``:

```math
\|Lf\|_{A_η} \le A\,\|f\|_{A_η} + B\,\|f\|_{L^1}.
```

`η′ > η` and `C` must give the coefficient bound
``|\hat c_k(Lf)| \le C\,e^{-2πη'|k|}\|f\|_{A_η}`` for every `k`. It holds with
``C \ge \sup_{S_{η'}}|Lf| / \sup_{S_η}|f|`` by Cauchy's estimate on the lines
``\operatorname{Im} x = \mp η'`` and ``\sup_{S_η}|f| \le \|f\|_{A_η}``; and with ``C = 1``
under the sign condition of [`annulus_coefficient_decay`](@ref), by moving the line of
integration of ``\int f e^{-2πikT}`` (Blumenthal–Nisoli–Taylor-Crush 2025, proof of
Lemma 3.7).  `C` is not ``\|L\|_{A_η \to A_{η'}}``, which the argument does not bound.
Splitting the Fourier series at `K`:

```math
A(K) = C\,\frac{2q^{K+1}}{1-q},\quad q = e^{-2π(η′-η)},
\qquad
B(K) = \sum_{|k| \le K} e^{2πη|k|} .
```

The `B` term uses ``|\hat c_k(Lf)| \le \|Lf\|_{L^1} \le \|f\|_{L^1}``: the
coefficient bound is uniform in `k`, and the transfer operator is an ``L^1``
contraction, so no constant is lost.
"""
function analytic_dfly(strong::Aη, η′::Real, C::Real, K::Integer)
    η = interval(strong.η)
    η′ > strong.η || error("need η′ > η: the strip must be strictly enlarged")
    q = exp(-2 * interval(π) * (interval(η′) - η))
    A = interval(C) * _tail_factor(q, K)

    r = exp(2 * interval(π) * η)                 # weight of one mode
    B = 1 + _strip_head(r, K)                    # Σ_{|k| ≤ K} e^{2πη|k|}
    return (sup(A), sup(B))
end

@doc raw"""
    analytic_dfly(strong::Eρ, ρ′, C, K; L1μ_bound = 1.0) -> (A, B)

Lasota–Yorke constants for the Chebyshev analytic norm
``\|f\|_{E_ρ} = \sum_k |\hat b_k| ρ^k``, with auxiliary norm ``L^1(μ)``, `μ` the
arcsine measure:

```math
\|Lf\|_{E_ρ} \le A\,\|f\|_{E_ρ} + B\,\|f\|_{L^1(μ)} ,
\qquad
A(K) = C\,\frac{2q^{K+1}}{1-q},\; q = ρ/ρ',
\qquad
B(K) = 2 D_μ \sum_{k \le K} ρ^k .
```

`ρ′ > ρ` and `C` must give ``|\hat b_k(Lf)| \le 2C\,ρ'^{-k}\|f\|_{E_ρ}`` for ``k \ge 1``,
which holds with ``C \ge \sup_{E_{ρ'}}|Lf| / \sup_{E_ρ}|f|`` by Cauchy's estimate on the
ellipse and ``\sup_{E_ρ}|f| \le \|f\|_{E_ρ}``; `C` is not ``\|L\|_{E_ρ \to E_{ρ'}}``.
`L1μ_bound` ``\ge D_μ`` is required: it has no default (see below).

!!! warning "The auxiliary norm is L¹(μ), not L¹(dx)"
    The Chebyshev coefficients are ``\hat b_k = 2\int f\,T_k\,dμ``, so
    ``|\hat b_k| \le 2\|f\|_{L^1(μ)}`` — but there is **no** uniform bound by
    ``\|f\|_{L^1(dx)}``: concentrating `f` on a width-`ε` sliver at an endpoint
    gives ``\|f\|_{L^1(dx)} \sim ε`` against ``\|f\|_{L^1(μ)} \sim \sqrt{2ε}/π``,
    a ratio diverging like ``ε^{-1/2}``.

    Consequently `L1μ_bound` = ``D_μ = \|L\|_{L^1(μ)}`` is **not** automatically
    1, unlike the Lebesgue case where the transfer operator is a contraction. It
    must be supplied, from a bound established independently; it used to default
    to 1, a placeholder with nothing behind it.
"""
function analytic_dfly(strong::Eρ, ρ′::Real, C::Real, K::Integer; L1μ_bound::Real)
    ρ = interval(strong.ρ)
    ρ′ > strong.ρ || error("need ρ′ > ρ: the ellipse must be strictly enlarged")
    q = ρ / interval(ρ′)
    A = interval(C) * _tail_factor(q, K)

    B = 2 * interval(L1μ_bound) * (ρ^(K + 1) - 1) / (ρ - 1)   # 2 D_μ Σ_{k ≤ K} ρ^k
    return (sup(A), sup(B))
end

@doc raw"""
    analytic_dfly_choose_K(strong, gain, C; target_A = 0.5, max_K = 4096, kwargs...)

Smallest `K` for which [`analytic_dfly`](@ref) achieves `A ≤ target_A`, returned
as `(K, A, B)`. Since `A(K)` decreases geometrically and `B(K)` grows
geometrically, this is the cheapest usable inequality: any smaller `K` fails the
target, any larger one inflates `B` for nothing.

Returns `nothing` if the target is not reachable within `max_K`.
"""
function analytic_dfly_choose_K(
    strong::NormKind,
    gain::Real,
    C::Real;
    target_A::Real = 0.5,
    max_K::Integer = 4096,
    kwargs...,
)
    for K = 0:max_K
        A, B = analytic_dfly(strong, gain, C, K; kwargs...)
        A <= target_A && return (K, A, B)
    end
    return nothing
end

@doc raw"""
    annulus_expansion(f, η; n = 1024) -> (η_in, η_out)

Certified enlargement of the strip ``|\mathrm{Im}\,x| \le η`` under a circle map
`f`, the Fourier counterpart of [`bernstein_expansion`](@ref).

Under ``x \mapsto z = e^{2πix}`` the strip becomes the annulus
``e^{-2πη} \le |z| \le e^{2πη}``, so the image of each boundary circle is
enclosed and its extreme modulus recorded:

```math
η_{out} = \frac{\log \min_{|z| = e^{2πη}} |f(z)|}{2π},
\qquad
η_{in} = \frac{-\log \max_{|z| = e^{-2πη}} |f(z)|}{2π}.
```

`f` acts on the `z` variable and must accept a `Complex{Interval}`. The strip is
enlarged when both exceed `η`; `min(η_in, η_out)` is then the certified `η′` for
[`analytic_dfly`](@ref).
"""
function annulus_expansion(f, η; n::Integer = 1024)
    ηi = η isa Interval ? η : interval(η)
    twoπ = 2 * interval(π)
    r_out, r_in = exp(twoπ * ηi), exp(-twoπ * ηi)

    lo_out, hi_in = Inf, -Inf
    for j = 1:n
        θ = interval((j - 1) / n, j / n)
        c, s = cos(twoπ * θ), sin(twoπ * θ)
        lo_out = min(lo_out, inf(_cabs(f(complex(r_out * c, r_out * s)))))
        hi_in = max(hi_in, sup(_cabs(f(complex(r_in * c, r_in * s)))))
    end
    η_out = inf(log(interval(lo_out)) / twoπ)
    η_in = inf(-log(interval(hi_in)) / twoπ)
    return (η_in, η_out)
end

"""
    strip_expansion(f, η; n = 1024) -> η′

`min` of the two values from [`annulus_expansion`](@ref): the largest strip
half-width certified to be reached from `η`.
"""
function strip_expansion(f, η; n::Integer = 1024)
    η_in, η_out = annulus_expansion(f, η; n = n)
    return min(η_in, η_out)
end

###############################################################################
# A priori operator bounds from the complex neighbourhood
#
# These must NOT be read off the assembled matrix: the DFLY is the a priori
# input that justifies the discretization, so deriving it from the discretized
# operator would be circular. Everything below is computed from the map and its
# derivative on the enclosed complex neighbourhood.
###############################################################################

export min_modulus_on_ellipse, min_modulus_on_circle, analytic_transfer_bound,
    analytic_dfly_degenerate

@doc raw"""
    min_modulus_on_ellipse(g, ρ; n = 1024)

Rigorous lower bound on ``\min_{∂E_ρ} |g|``, by the same boundary covering as
[`bernstein_expansion`](@ref).

If `g` has no zero in the closed ellipse — which is the case for `T'` of an
expanding map — the minimum modulus principle makes this a bound on the whole of
``E_ρ``, not just its boundary.
"""
function min_modulus_on_ellipse(g, ρ; n::Integer = 1024)
    ρi = ρ isa Interval ? ρ : interval(ρ)
    lo = Inf
    for j = 1:n
        θ = interval((j - 1) / n, j / n)
        lo = min(lo, inf(_cabs(g(bernstein_point(ρi, θ)))))
    end
    return lo
end

"""
    min_modulus_on_circle(g, η; n = 1024)

Lower bound on `|g|` over both circles `|z| = e^{±2πη}`, the strip counterpart
of [`min_modulus_on_ellipse`](@ref).
"""
function min_modulus_on_circle(g, η; n::Integer = 1024)
    ηi = η isa Interval ? η : interval(η)
    twoπ = 2 * interval(π)
    lo = Inf
    for r in (exp(twoπ * ηi), exp(-twoπ * ηi)), j = 1:n
        θ = interval((j - 1) / n, j / n)
        lo = min(lo, inf(_cabs(g(complex(r * cos(twoπ * θ), r * sin(twoπ * θ))))))
    end
    return lo
end

@doc raw"""
    analytic_transfer_bound(min_derivative, nbranches) -> C

``C`` with ``\|Lf\|_{∞,\,\text{nbhd}} \le C\,\|f\|_{∞,\,\text{nbhd}}``, for the
transfer operator ``Lf(z) = \sum_k f(g_k(z))\,g_k'(z)``.

Since ``|g_k'| = 1/|T'\circ g_k|`` and the inverse branches land in the
neighbourhood where `min_derivative` was certified,

```math
C \;\le\; \frac{\#\text{branches}}{\min |T'|} .
```

`min_derivative` comes from [`min_modulus_on_ellipse`](@ref) or
[`min_modulus_on_circle`](@ref) — i.e. from the complex neighbourhood, never
from the assembled matrix.
"""
function analytic_transfer_bound(min_derivative::Real, nbranches::Integer)
    min_derivative > 0 || error("T' must be bounded away from 0 on the neighbourhood")
    return Float64(nbranches, RoundUp) ⊘₊ min_derivative
end

@doc raw"""
    analytic_dfly_degenerate(strong, gain, C₂) -> (A, 0.0)

The degenerate Lasota–Yorke ``\|Lf\|_{A_η} \le A\,\|f\|_{A_η}``, with no
auxiliary term — that is, `A` is simply the **continuity constant** of `L` on
the analytic space, valid for every `f`, not on any subspace.

`C₂` bounds ``\sup_{S_{η'}}|Lf| \le C_2\sup_{S_η}|f|``, a bound between sup norms on the
two domains, not ``\|L\|_{A_η \to A_{η'}}``.  With ``M := \sup_{S_{η'}}|Lf|``, Cauchy's
estimate on the lines ``\operatorname{Im} x = \mp η'`` gives
``|\hat c_k(Lf)| \le Me^{-2πη'|k|}``, and ``\sup_{S_η}|f| \le \|f\|_{A_η}``, so, with
``δ = η' - η``,

```math
\|Lf\|_{A_η} = \sum_k e^{2πη|k|}|\hat c_k(Lf)|
   \;\le\; M\sum_k e^{-2πδ|k|}
   \;\le\; G(δ)\,C_2\,\|f\|_{A_η},
```

the geometric sum converging because the neighbourhood is enlarged. For `Eρ` the weight is
`ρ^k`, the ratio is `ρ/ρ'`, and Cauchy's estimate on the ellipse is ``|\hat b_0| \le M``,
``|\hat b_k| \le 2Mρ'^{-k}``, which gives the same ``G``.

This is the single-space setting of Nisoli, *Certified spectral approximation of
transfer operators and the Gauss map*, arXiv:2602.19435, where the analytic
(Hardy-space) case needs only ``\|L\|_{B \to B} \le C`` together with the
truncation bound coming from the domain enlargement — the strong–weak DFLY
scale of that paper's Appendix A being a separate setting. Compactness, not a
contraction factor, is what drives the spectral certification, so `A` here need
not be less than 1.
"""
function analytic_dfly_degenerate(strong::Aη, η′::Real, C₂::Real)
    η′ > strong.η || error("need η′ > η")
    q = exp(-2 * interval(π) * (interval(η′) - interval(strong.η)))
    G = 1 + 2 * q / (1 - q)                    # Σ_{k ∈ ℤ} q^{|k|}
    return (sup(interval(C₂) * G), 0.0)
end

function analytic_dfly_degenerate(strong::Eρ, ρ′::Real, C₂::Real)
    ρ′ > strong.ρ || error("need ρ′ > ρ")
    q = interval(strong.ρ) / interval(ρ′)
    # Cauchy estimate on the ellipse (Trefethen, ATAP, Thm 8.1) is |b̂₀| ≤ M and
    # |b̂_k| ≤ 2M ρ'^{-k} for k ≥ 1 — the factor 2 matters. Σ_k |b̂_k|ρ^k is then
    # bounded by M(1 + 2Σ_{k≥1} q^k), the same G as the Fourier case.
    G = 1 + 2 * q / (1 - q)
    return (sup(interval(C₂) * G), 0.0)
end

###############################################################################
# dfly methods, so the analytic norms drop into the existing RIM pipeline
###############################################################################

import .RigorousInvariantMeasures: dfly, branches, derivative

# Branch k, reparametrized so its domain is [-1,1]:  s ↦ x_k(s).
_branch_chart(br) = s -> br.X[1] + (br.X[2] - br.X[1]) * (s + 1) / 2

@doc raw"""
    dfly(norm::Aη, ::Type{L1}, D::PwMap; C₂, n = 1024) -> (A, 0.0)
    dfly(norm::Eρ, ::Type{L1}, D::PwMap; C₂, n = 1024) -> (A, 0.0)

Degenerate Lasota–Yorke for the analytic norms: `A` is the continuity constant
of `L` on the analytic space, and the auxiliary constant is `0`.

Same call shape as every other `dfly`, so `powernormbounds`,
`refine_norms_of_powers` and the coarse–fine workflow run unchanged once the
basis fixes its strong and weak norms. Because `B = 0` the auxiliary norm is
multiplied by zero and plays no part; the weak norm is free to be whatever the
approximation error is measured in, which for these bases is `L2`.

The geometry needs only **forward enclosures** — the branches of a `PwMap`
evaluate at `Complex{Interval}` directly, so the enlarged neighbourhood is
certified by [`strip_expansion`](@ref) / [`bernstein_expansion`](@ref) applied to
each branch reparametrized onto its own domain. No complex derivative is
involved.

`C₂` bounds `L` between sup norms on the two domains,
``\sup_{S_{η'}}|Lf| \le C_2\sup_{S_η}|f|``; the factor `G(δ)` from the enlargement then
gives `A = G(δ)C₂` (see [`analytic_dfly_degenerate`](@ref)).

`C₂` is the one ingredient the forward enclosure does not determine on its own:
the transfer operator carries the weight ``|g_k'| = 1/|T_k'\circ g_k|``, so by
default it is computed as ``\sum_k 1/\min|T_k'|`` with the minimum taken over the
enclosed neighbourhood — valid on the closed neighbourhood by the minimum
modulus principle, ``T_k'`` being zero-free for an expanding map. Pass `C₂`
explicitly to override it with a sharper estimate for your operator.

!!! note "A ≥ 1 is expected"
    `A` is a continuity constant, not a contraction factor: compactness is what
    drives the certification. Helpers assuming the classical shape — notably
    `invariant_measure_strong_norm_bound`, which forms `B/(1-A)` — do not apply
    and say so rather than return a meaningless number.
"""
function dfly(norm::Aη, ::Type{L1}, D::PwMap; C₂::Union{Real,Nothing} = nothing,
              n::Integer = 1024)
    η = interval(norm.η)
    η′, C = Inf, 0.0
    for br in branches(D)
        lo_der = Inf
        for j = 1:n, σ in (1, -1)
            t = br.X[1] + (br.X[2] - br.X[1]) * interval((j - 1) / n, j / n)
            x = complex(t, σ * η)
            η′ = min(η′, inf(abs(imag(br.f(x)))))
            lo_der = min(lo_der, inf(_cabs(derivative(br.f, x))))
        end
        C = C ⊕₊ (1.0 ⊘₊ lo_der)
    end
    C₂ = C₂ === nothing ? C : C₂
    η′ > norm.η ||
        error("the strip is not enlarged (η′ = $η′ ≤ η = $(norm.η)); try a smaller η")
    return DFLYCoefficients(norm, L1, analytic_dfly_degenerate(norm, η′, C₂)...)
end

function dfly(norm::Eρ, ::Type{L1}, D::PwMap; C₂::Union{Real,Nothing} = nothing,
              n::Integer = 1024)
    ρ = interval(norm.ρ)
    ρ′, C = Inf, 0.0
    for br in branches(D)
        # The Bernstein ellipse lives around the GLOBAL interval, so the branch
        # must be read in the global chart t = 2x - 1 (`to_symmetric_interval`),
        # under which F_k maps the sub-arc [t₁,t₂] ONTO [-1,1].
        #
        # This used to use `_branch_chart`, rescaling each branch's own domain
        # to [-1,1]. That makes F_k a bijection of [-1,1] with derivative
        # f′(x)(X₂-X₁) — for the Lanford map 1.0961 → 1.0000 → 0.9039 across
        # branch 1 — i.e. a near-isometry, never uniformly expanding, so the
        # test failed for every ρ (ρ′/ρ ≈ 0.92, flat in ρ). The rescaling
        # normalizes away exactly the contraction the argument needs: g_k
        # carries all of [0,1] onto a subinterval, a factor ≈2 in the global
        # coordinate. With the global chart the Lanford map expands for every
        # ρ ≥ 1.2, with q = ρ/ρ′ minimal (0.8605) near ρ = 3.
        F = to_symmetric_interval(br.f)
        ρ′ = min(ρ′, bernstein_expansion(F, ρ; n = n))
        lo_der = Inf
        for j = 1:n
            # global chart: t ∈ ∂E_ρ  ↦  x = (t+1)/2
            x = (bernstein_point(ρ, interval((j - 1) / n, j / n)) + 1) / 2
            lo_der = min(lo_der, inf(_cabs(derivative(br.f, x))))
        end
        C = C ⊕₊ (1.0 ⊘₊ lo_der)
    end
    C₂ = C₂ === nothing ? C : C₂
    ρ′ > norm.ρ ||
        error("the ellipse is not enlarged (ρ′ = $ρ′ ≤ ρ = $(norm.ρ)); try a smaller ρ")
    return DFLYCoefficients(norm, L1, analytic_dfly_degenerate(norm, ρ′, C₂)...)
end

# checked: invariant_measure_strong_norm_bound -- the methods below are the existing Aη and
# Eρ methods of the generic function, corrected in place (they returned A); the B/(1-A) of
# the other bases is the same formula, without the aux-norm check the analytic pairs need.
@doc raw"""
    invariant_measure_strong_norm_bound(B, D; dfly_coefficients)

Strong-norm bound ``B/(1-A)`` on the invariant density ``h``, ``\ell(h) = 1``, from
``\|Lf\|_s \le A\|f\|_s + B\|f\|_{L^1}`` with ``A < 1``: ``\|h\|_s = \|Lh\|_s \le A\|h\|_s + B``.

The degenerate pair of `dfly(::Aη, L1, D)` and `dfly(::Eρ, L1, D)`, with ``B = 0`` and `A` a
continuity constant ``\ge 1``, gives no bound, and is refused with an `ArgumentError`.  (These
methods used to return `A` in that case, on the argument
``\|h\|_s = \|Lh\|_s \le A\|h\|_s``, which bounds nothing.)  For the Fourier strip, the pair of
[`annulus_dfly_l1`](@ref) in its full-norm form, or of [`analytic_dfly`](@ref), has ``A < 1``
and an ``L^1`` second term; pass it as `dfly_coefficients`.  For `Eρ` the second term of
[`analytic_dfly`](@ref) is in ``L^1(μ)``, where ``\|h\|_{L^1(μ)} \ne 1``, so only a pair
stated with the Lebesgue ``L^1`` norm is accepted.
"""
function invariant_measure_strong_norm_bound(
    B::FourierAnalytic{Aη},
    D::Dynamic;
    dfly_coefficients = basis_dfly(B, D),
)
    check_dfly_coefficients(B, D, dfly_coefficients)
    return _analytic_fixed_point_bound(dfly_coefficients)
end

function invariant_measure_strong_norm_bound(
    B::Chebyshev{Eρ},
    D::Dynamic;
    dfly_coefficients = basis_dfly(B, D),
)
    check_dfly_coefficients(B, D, dfly_coefficients)
    return _analytic_fixed_point_bound(dfly_coefficients)
end

# B/(1 - A) for an analytic pair, refused unless A < 1 and the second term is Lebesgue L¹.
function _analytic_fixed_point_bound(c)
    aux_kind(c) === L1 ||
        throw(ArgumentError("the second term must be the Lebesgue L¹ norm, got $(aux_kind(c))"))
    A, Bc = c
    A < 1 || throw(ArgumentError(
        "A = $A ≥ 1: a continuity constant bounds no fixed point; pass a pair with A < 1"))
    return Bc ⊘₊ (1.0 ⊖₋ A)
end

###############################################################################
# Single-space certification: L is compact on the analytic space, so the
# resolvent of the abstract operator follows from the resolvent of the
# discretization plus ONE truncation number. No Lasota-Yorke enters.
#
# This is the setting of Nisoli, "Certified spectral approximation of transfer
# operators and the Gauss map", arXiv:2602.19435, Section 2 (Assumption 2.1,
# Corollary 2.3), as opposed to the strong-weak (DFLY) framework of its
# Appendix A. For the analytic norms the single-space route is the right one and
# is enormously sharper: on the Lanford map it gives R(1,L) <= 2.46, against
# 1.3e52 from the Appendix A lifting, which has to buy `a < 1` by splitting the
# coefficient series and pays b ~ rho^K in the numerator.
###############################################################################

export bernstein_geometry, strip_geometry, analytic_truncation_error,
       single_space_resolvent_lift

@doc raw"""
    bernstein_geometry(strong::Eρ, D::PwMap; n = 1024) -> (ρ′, C₂)

The two numbers the analytic theory rests on, computed once from the map:

* `ρ′ > ρ` — the enlarged Bernstein parameter, certified by
  [`bernstein_expansion`](@ref) applied to each branch in the GLOBAL chart
  ``t = 2x-1``. It means the inverse branches contract: ``g_k(E_{ρ'}) ⊆ E_ρ``.
* `C₂ = Σ_k 1/\min_{∂E_ρ}|T_k'|`, rounded up, the minimum enclosed on `n` arcs of
  ``∂E_ρ``.  When every ``T_k'`` is zero-free on ``E_ρ`` the minimum modulus principle
  puts ``\min_{E_ρ}|T_k'|`` on ``∂E_ρ``, and `C₂` bounds ``L`` between sup norms,
  ``\sup_{E_{ρ'}}|Lf| \le C_2\sup_{E_ρ}|f|``; this routine does not check that hypothesis
  ([`certify_compatibility`](@ref) does).  `C₂` is not an operator norm between the
  coefficient spaces: in ``ℓ^1``, Cauchy's estimate gives ``\|L\|_{E_ρ\to E_ρ} \le GC_2``
  ([`analytic_dfly_degenerate`](@ref)) and ``\sum_{k>K}|b_k(Lf)|ρ^k \le 2C_2q^{K+1}/(1-q)\|f\|_{E_ρ}``;
  in ``ℓ^2``, ``\|L\|_{E_{ρ,2}\to E_{ρ',2}} \le 2C_2/(1 - (r/ρ)^2)^{1/2}`` when ``C_2`` is the
  constant of a smaller ellipse ``E_r \to E_{ρ'}``.

These are exactly the quantities `dfly(::Eρ, ::Type{L1}, ::PwMap)` computes
internally and then discards; exposing them avoids recomputing.  Throws if
`ρ′ ≤ ρ`.
"""
function bernstein_geometry(strong::Eρ, D::PwMap; n::Integer = 1024)
    ρ′, C = _bernstein_geometry(strong, D; n = n)
    ρ′ > strong.ρ ||
        error("the ellipse is not enlarged (ρ′ = $ρ′ ≤ ρ = $(strong.ρ))")
    return (ρ′, C)
end

# The two numbers of `bernstein_geometry`, without the check ρ′ > ρ.
function _bernstein_geometry(strong::Eρ, D::PwMap; n::Integer = 1024)
    ρ = interval(strong.ρ)
    ρ′, C = Inf, 0.0
    for br in branches(D)
        ρ′ = min(ρ′, bernstein_expansion(to_symmetric_interval(br.f), ρ; n = n))
        lo_der = Inf
        for j = 1:n
            x = (bernstein_point(ρ, interval((j - 1) / n, j / n)) + 1) / 2
            lo_der = min(lo_der, inf(_cabs(derivative(br.f, x))))
        end
        C = C ⊕₊ (1.0 ⊘₊ lo_der)
    end
    return (ρ′, C)
end

@doc raw"""
    strip_geometry(strong::Aη, D::PwMap; n = 1024) -> (η′, C₂)

Fourier counterpart of [`bernstein_geometry`](@ref): `η′` is the minimum of the certified
lower bounds of ``|\operatorname{Im} T|`` on `n` arcs of the boundary lines
``\operatorname{Im} x = \pm η`` of each branch, and ``C_2 = \sum_k 1/\min|T_k'|`` over the
same arcs, rounded up.  Given the inclusion ``τ_k(S_{η'}) \subset S_η`` of the inverse
branches, which this routine does not prove (`certify_compatibility(::AnnulusGeometry, T;
degree)` does, together with the check that ``T'`` does not vanish inside the strip), `C₂`
bounds ``L`` between sup norms, ``\sup_{S_{η'}}|Lf| \le C_2\sup_{S_η}|f|``; it is not
``\|L\|_{A_η → A_{η'}}``.
"""
function strip_geometry(strong::Aη, D::PwMap; n::Integer = 1024)
    η′, C = Inf, 0.0
    for br in branches(D)
        lo_der = Inf
        for j = 1:n, σ in (1, -1)
            t = br.X[1] + (br.X[2] - br.X[1]) * interval((j - 1) / n, j / n)
            x = complex(t, σ * interval(strong.η))
            η′ = min(η′, inf(abs(imag(br.f(x)))))
            lo_der = min(lo_der, inf(_cabs(derivative(br.f, x))))
        end
        C = C ⊕₊ (1.0 ⊘₊ lo_der)
    end
    η′ > strong.η || error("the strip is not enlarged (η′ = $η′ ≤ η = $(strong.η))")
    return (η′, C)
end

@doc raw"""
    analytic_truncation_error(strong, gain, C₂, K) -> ε_K

``ε_K \ge \|(I - Π_K)L\|_{\mathcal B \to \mathcal B}`` on the analytic space
``\mathcal B``, i.e. the truncation number of Assumption 2.1(3).

One line, from the domain gain alone. `C₂` must bound ``\|L\|_{E_ρ \to E_{ρ'}}`` in the
coefficient norm in use, so that ``\|Lf\|_{E_{ρ'}} \le C_2`` for `f` in the unit ball of
``E_ρ``.  The `C₂` of [`bernstein_geometry`](@ref) and [`strip_geometry`](@ref) bounds ``L``
between sup norms and is not that constant (see the docstring of `bernstein_geometry`).
Writing
``Lf = \sum_j b_j T_j``,

```math
\|(I-Π_K)Lf\|_{E_ρ} = \sum_{j>K}|b_j|ρ^j
 = \sum_{j>K}|b_j|ρ'^j\Big(\frac{ρ}{ρ'}\Big)^{j} \le C_2\,q^{K+1},
\qquad q = ρ/ρ' < 1 .
```

The Fourier case is identical with ``q = e^{-2π(η'-η)}``.

!!! note "Where the compactness lives"
    Not in `L` itself but in the INCLUSION ``A(E_{ρ'}) \hookrightarrow A(E_ρ)``,
    whose approximation numbers decay like ``q^K``. Factoring
    ``L = ι ∘ \tilde L`` with ``\tilde L : A(E_ρ) → A(E_{ρ'})`` bounded by `C₂`
    gives the bound above immediately. There is no separate
    ``\|L(I-Π_K)\|`` to estimate as long as the discretization is taken
    one-sided, ``L_K := Π_K L``, which has the same nonzero spectrum as
    ``Π_K L Π_K`` (since ``σ(AB)\setminus\{0\} = σ(BA)\setminus\{0\}``) and is
    therefore represented by the same assembled matrix.
"""
function analytic_truncation_error(strong::Eρ, ρ′::Real, C₂::Real, K::Integer)
    q = interval(strong.ρ) / interval(ρ′)
    return sup(interval(C₂) * q^(K + 1))
end

function analytic_truncation_error(strong::Aη, η′::Real, C₂::Real, K::Integer)
    q = exp(-2 * interval(π) * (interval(η′) - interval(strong.η)))
    return sup(interval(C₂) * q^(K + 1))
end

@doc raw"""
    single_space_resolvent_lift(R_K, ε_K) -> R

Corollary 2.3 of the Gauss-map paper: if ``α := ε_K R(z,L_K) < 1`` then
``z ∈ ρ(L)`` and

```math
R(z, L) \;\le\; \frac{R(z, L_K)}{1 - ε_K R(z, L_K)} .
```

Returns `Inf` when `α ≥ 1`.

Both arguments must be measured in the SAME norm — for the analytic bases, the
``ℓ^2`` Bernstein norm of [`bernstein_l2_resolvent_bound`](@ref), which is a
diagonal reweighting and so is directly accessible to the verified SVD.

The result is a certificate on the ABSTRACT operator `L`. It is therefore a
property of the map, not of the discretization: computed once at a small `K` it
may be reused at any basis size and any working precision, and never needs
recomputing.
"""
function single_space_resolvent_lift(R_K::Real, ε_K::Real)
    α = ε_K ⊗₊ R_K
    α < 1 || return Inf
    return R_K ⊘₊ (1.0 ⊖₋ α)
end
