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
# (`dfly(TotalVariation, L1, ::DiscretizedNoiseKernelUlam)` returns `(0.0, …)`),
# where the smoothing plays the role of the analytic gain.
###############################################################################

export analytic_dfly, analytic_dfly_choose_K, annulus_expansion, strip_expansion

# Sum of the tail weights above K: 2 Σ_{j>K} q^j = 2 q^{K+1}/(1-q), q < 1.
function _tail_factor(q, K::Integer)
    qi = interval(q)
    @assert sup(qi) < 1 "the analytic gain must be strictly positive (q < 1)"
    return 2 * qi^(K + 1) / (1 - qi)
end

@doc raw"""
    analytic_dfly(strong::Aη, η′, C, K) -> (A, B)

Lasota–Yorke constants for the Fourier analytic norm
``\|f\|_{A_η} = \sum_k |\hat c_k| e^{2πη|k|}``, with auxiliary norm ``L^1(dx)``:

```math
\|Lf\|_{A_η} \le A\,\|f\|_{A_η} + B\,\|f\|_{L^1}.
```

`η′ > η` is the enlarged strip certified by [`annulus_expansion`](@ref) and `C`
bounds ``\|L\|_{A_η \to A_{η′}}``. Splitting the Fourier series at `K`:

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
    B = 1 + 2 * r * (r^K - 1) / (r - 1)          # Σ_{|k| ≤ K} e^{2πη|k|}
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

`ρ′ > ρ` is the enlarged ellipse certified by [`bernstein_expansion`](@ref) and
`C` bounds ``\|L\|_{E_ρ \to E_{ρ'}}``.

!!! warning "The auxiliary norm is L¹(μ), not L¹(dx)"
    The Chebyshev coefficients are ``\hat b_k = 2\int f\,T_k\,dμ``, so
    ``|\hat b_k| \le 2\|f\|_{L^1(μ)}`` — but there is **no** uniform bound by
    ``\|f\|_{L^1(dx)}``: concentrating `f` on a width-`ε` sliver at an endpoint
    gives ``\|f\|_{L^1(dx)} \sim ε`` against ``\|f\|_{L^1(μ)} \sim \sqrt{2ε}/π``,
    a ratio diverging like ``ε^{-1/2}``.

    Consequently `L1μ_bound` = ``D_μ = \|L\|_{L^1(μ)}`` is **not** automatically
    1, unlike the Lebesgue case where the transfer operator is a contraction. It
    must be supplied; the default of 1 is a placeholder and is only valid if you
    have established it independently.
"""
function analytic_dfly(strong::Eρ, ρ′::Real, C::Real, K::Integer; L1μ_bound::Real = 1.0)
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

`C₂` is the bound between the two domains,
``\|L\|_{A_η \to A_{η'}}``, and the enlargement is converted into a constant by
splitting the weight at the gain `δ`:

```math
\|Lf\|_{A_η} = \sum_k e^{-2πkδ}\,e^{2πkη'}|\hat c_k(Lf)|
   \;\le\; \Big(\sum_k e^{-2πkδ}\Big)\,\|Lf\|_{A_{η'}}
   \;\le\; G(δ)\,C_2\,\|f\|_{A_η},
```

the geometric sum converging exactly because the neighbourhood is enlarged. For
`Eρ` the weight is `ρ^k` and the ratio is `ρ/ρ'`.

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
    G = 1 / (1 - q)                            # Σ_{k ≥ 0} q^k
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

`C₂` bounds `L` between the two domains, ``\|L\|_{A_η \to A_{η'}}``; the Hölder
factor `G(δ)` from the enlargement then gives `A = G(δ)C₂` (see
[`analytic_dfly_degenerate`](@ref)).

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
    return analytic_dfly_degenerate(norm, η′, C₂)
end

function dfly(norm::Eρ, ::Type{L1}, D::PwMap; C₂::Union{Real,Nothing} = nothing,
              n::Integer = 1024)
    ρ = interval(norm.ρ)
    ρ′, C = Inf, 0.0
    for br in branches(D)
        chart = _branch_chart(br)
        # F_k maps [-1,1] onto [-1,1]; its ellipse image is what must expand
        ρ′ = min(ρ′, bernstein_expansion(s -> 2 * br.f(chart(s)) - 1, ρ; n = n))
        lo_der = Inf
        for j = 1:n
            x = chart(bernstein_point(ρ, interval((j - 1) / n, j / n)))
            lo_der = min(lo_der, inf(_cabs(derivative(br.f, x))))
        end
        C = C ⊕₊ (1.0 ⊘₊ lo_der)
    end
    C₂ = C₂ === nothing ? C : C₂
    ρ′ > norm.ρ ||
        error("the ellipse is not enlarged (ρ′ = $ρ′ ≤ ρ = $(norm.ρ)); try a smaller ρ")
    return analytic_dfly_degenerate(norm, ρ′, C₂)
end

@doc raw"""
    invariant_measure_strong_norm_bound(B, D; dfly_coefficients)

Strong-norm bound on the invariant density.

The classical DFLY route is ``B/(1-A)``, which needs `A < 1` and `B > 0`. In the
degenerate analytic case `B = 0` and `A` is a continuity constant, and the bound
is simply **`A`** — the density lies in the image of `L`, so the continuity
constant already controls it.
"""
function invariant_measure_strong_norm_bound(
    B::FourierAnalytic{Aη},
    D::Dynamic;
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
)
    return first(dfly_coefficients)
end

function invariant_measure_strong_norm_bound(
    B::Chebyshev{Eρ},
    D::Dynamic;
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
)
    return first(dfly_coefficients)
end
