###############################################################################
# Strip (annulus) geometry of the Fourier analytic setting, and its certification
# for an analytic circle map: the counterpart of EllipseGeometry.jl.
#
# A function on the circle analytic on the strip S_η = {|Im x| < η} is, in the variable
# z = e^{2πix}, analytic on the annulus e^{-2πη} < |z| < e^{2πη}.  The geometry is a
# triple 0 < r < ρ < R of half-widths, chosen by the user.  `certify_compatibility(G, T;
# degree)` proves, for an analytic lift T of the circle map, that the inverse branches
# extend analytically to S_R with τ_j(S_R) ⊂ S_r, and returns a certified
# W ≥ sup_{S_R} Σ_j |τ_j'|, so that sup_{S_R}|Lf| ≤ W sup_{S_r}|f|.
#
# checked: certify_compatibility -- the EllipseGeometry method checks each branch of a
# PwMap on E_r; the strip needs one lift T on a period rectangle, so this file adds a
# method of the same function for AnnulusGeometry.  The search loop is the shared
# _search_geometry of EllipseGeometry.jl.
###############################################################################

export AnnulusGeometry, AnnulusCompatibility, search_annulus_geometry

@doc raw"""
    AnnulusGeometry(r, ρ, R)

The strip half-widths ``0 < r < ρ < R`` of a Fourier analytic setting: functions are
bounded on ``S_r = \{|\operatorname{Im} x| < r\}``, the strong norm is that of ``S_ρ``, and
the transfer operator is to map functions analytic on ``S_r`` to functions analytic on
``S_R``.  In the variable ``z = e^{2πix}`` the strip ``S_η`` is the annulus
``e^{-2πη} < |z| < e^{2πη}``.  The values are stored as the `Float64`s passed, exactly; as
for [`EllipseGeometry`](@ref), two geometries built from the same numbers are different
objects.  Whether a map is compatible with the geometry is proved by
[`certify_compatibility`](@ref).
"""
mutable struct AnnulusGeometry
    # mutable, with constant fields, so that `===` is identity (see EllipseGeometry)
    const r::Float64
    const ρ::Float64
    const R::Float64
    function AnnulusGeometry(r::Real, ρ::Real, R::Real)
        r, ρ, R = Float64(r), Float64(ρ), Float64(R)
        0 < r < ρ < R || throw(ArgumentError("need 0 < r < ρ < R, got $r, $ρ, $R"))
        return new(r, ρ, R)
    end
end

@doc raw"""
    AnnulusCompatibility

What [`certify_compatibility`](@ref) proves for an [`AnnulusGeometry`](@ref) `G` and a lift
`lift` of degree `degree`: every inverse branch maps ``S_R`` into ``S_r``, and
``\sup_{S_R}\sum_j|τ_j'| \le W``.  `R_image` is the certified lower bound on
``|\operatorname{Im} T|`` on the boundary lines of ``S_r``; it exceeds `G.R`.
"""
struct AnnulusCompatibility{FT}
    G::AnnulusGeometry
    lift::FT
    degree::Int
    W::Float64
    R_image::Float64
end

@doc raw"""
    certify_compatibility(G::AnnulusGeometry, T; degree, n = 1024, cover = 64)
        -> AnnulusCompatibility

Proves that the circle map with lift `T` is compatible with `G`, or throws an
`ArgumentError`.  `T` must accept `Complex{Interval}` and `Taylor1` of them, and the
statement assumes what no finite computation can check: `T` is analytic on a neighbourhood
of the closed strip ``\bar S_r`` and ``T(x + 1) = T(x) + d`` with ``d = `` `degree`
``\ne 0``.  The routine checks

1. ``T(1) - T(0)`` encloses ``d`` (a consistency check of the assumption, not a proof);
2. on the cover of the period rectangle ``[0,1] \times [-r, r]`` by `cover × cover` boxes,
   ``\operatorname{Re}T'`` has the sign of ``d``;
3. on `n` arcs of each boundary line ``\operatorname{Im} x = \pm r``,
   ``0 \le \operatorname{Re} x \le 1``, the enclosure of ``|\operatorname{Im} T|`` is
   ``> R``; the minimum of its lower bounds is `R_image`.

Since ``T'`` is 1-periodic, 2 holds on all of ``\bar S_r``, and ``T`` is injective there,
``\bar S_r`` being convex (Noshiro–Warschawski).  ``T(x) - dx`` is 1-periodic, hence
bounded on ``\bar S_r``, so ``\operatorname{Re} T(z) \to \pm\infty`` as
``\operatorname{Re} z \to \pm\infty`` uniformly in the strip.  Then ``T(S_r) \cap S_R`` is
open in ``S_R`` (``T`` is open), closed in ``S_R`` (a limit ``T(z_k) \to y \in S_R`` has
``\operatorname{Re} z_k`` bounded, so a subsequence ``z_k \to z \in \bar S_r`` with
``T(z) = y``, and ``z \notin \partial S_r`` by 3), and contains the real line.  As ``S_R`` is
connected, ``S_R \subset T(S_r)``, and the ``|d|`` inverse branches
``τ_j(y) = T^{-1}(y + j)``, ``j = 0, \dots, |d| - 1``, map ``S_R`` into ``S_r``.  By 2,
``1/T'`` is holomorphic, bounded and 1-periodic on ``\bar S_r``, so in the variable
``z = e^{2πix}`` it is holomorphic on a closed annulus and its maximum modulus is on the
boundary lines, where it is enclosed on the same `n` arcs:
``|τ_j'(y)| = 1/|T'(τ_j(y))| \le 1/\min_{\partial S_r}|T'|``, and
``W = |d| / \min_{\partial S_r}|T'|``, rounded up.

Remark.  Check 3 needs `T` on the boundary path only, not on a box of the plane.  If `T`
extends to complex intervals and ``T'`` is computed with `TaylorSeries`, the image of the
arc ``x + i r``, ``x \in I_i = [x_i, x_{i+1}]``, is contained in
```math
T(x_i + i r) + T'(I_i + i r)\,[0,\, x_{i+1} - x_i],
```
since ``T(x + ir) - T(x_i + ir) = \int_{x_i}^{x} T'(t + ir)\,dt`` and the average of
``T'`` on the segment lies in the box ``T'(I_i + ir)``, which is convex; the same holds on
``\operatorname{Im} x = -r``.  The routine uses the direct evaluation of `T` on
``I_i \pm ir`` instead.
"""
function certify_compatibility(G::AnnulusGeometry, T; degree::Integer,
                               n::Integer = 1024, cover::Integer = 64)
    d = _checked_degree(T, degree)
    fail = _lift_derivative_failure(T, d, G.r, cover)
    fail === nothing || throw(ArgumentError(
        "Re T' may vanish or have the wrong sign on S_$(G.r) (box $(fail.box), T' ∈ $(fail.d))"))
    R_image, W = _strip_boundary_bounds(T, d, G.r, n)
    R_image > G.R || throw(ArgumentError(
        "the images of ∂S_r are certified only outside S_$R_image, not outside S_R, R = $(G.R)"))
    isfinite(W) || throw(ArgumentError("|T'| is not bounded away from 0 on ∂S_r"))
    return AnnulusCompatibility(G, T, d, W, R_image)
end

@doc raw"""
    search_annulus_geometry(T, rs; degree, cost = (r, R_image, W) -> R_image / r,
                            n = 1024, cover = 64) -> NamedTuple

The counterpart of [`search_ellipse_geometry`](@ref) for the lift `T` of degree `degree`:
tries each half-width ``r > 0`` in `rs`, keeps those that pass checks 2 and 3 of
`certify_compatibility(::AnnulusGeometry, T; degree)` with `R_image > r`, and returns the
one that maximises `cost(r, R_image, W)`, with the same fields and `table` as
`search_ellipse_geometry`.  Throws an `ArgumentError` if no ``r`` passes, or if check 1
fails.

The search chooses ``r`` and certifies nothing by itself.  Any ``ρ, R`` with
``r < ρ < R < R_\mathrm{image}`` then pass
`certify_compatibility(AnnulusGeometry(r, ρ, R), T; degree)` with the same `n` and `cover`,
which is the certificate.
"""
function search_annulus_geometry(T, rs; degree::Integer,
                                 cost = (r, R_image, W) -> R_image / r,
                                 n::Integer = 1024, cover::Integer = 64)
    d = _checked_degree(T, degree)
    return _search_geometry(rs, cost, r -> begin
        r > 0 || throw(ArgumentError("need r > 0, got $r"))
        _lift_derivative_failure(T, d, r, cover) === nothing || return nothing
        return _strip_boundary_bounds(T, d, r, n)
    end)
end

# Check 1: the degree is nonzero and T(1) - T(0) encloses it.
function _checked_degree(T, degree::Integer)
    d = Int(degree)
    d != 0 || throw(ArgumentError("the degree must be nonzero"))
    in_interval(d, T(interval(1)) - T(interval(0))) ||
        throw(ArgumentError("T(1) - T(0) does not enclose the degree $d"))
    return d
end

# The interval [(i-1)/m, i/m], rounded outward.
_mesh_cell(i::Integer, m::Integer) = hull(interval(i - 1) / interval(m), interval(i) / interval(m))

# Check 2: `nothing` if Re T' has the sign of d on every box of the cover of
# [0,1] × [-r, r], otherwise the first box where that is not certified and the enclosure
# of T' on it.  The grid lines are shared by neighbouring boxes and the outer ones are
# exactly ±r, so the boxes cover the rectangle without gaps.
function _lift_derivative_failure(T, d::Integer, r::Float64, m::Integer)
    ys = [l == 0 ? -r : l == m ? r : -r + 2r * l / m for l = 0:m]
    for i = 1:m, l = 1:m
        box = complex(_mesh_cell(i, m), interval(ys[l], ys[l+1]))
        D = derivative(T, box)
        ok = d > 0 ? inf(real(D)) > 0 : sup(real(D)) < 0
        ok || return (box = box, d = D)
    end
    return nothing
end

# Check 3 and W: on n arcs of each line Im x = ±r, 0 ≤ Re x ≤ 1, the minimum of the
# lower bounds of |Im T|, and |d| / min |T'| rounded up.
function _strip_boundary_bounds(T, d::Integer, r::Float64, n::Integer)
    R_image, lo_der = Inf, Inf
    for j = 1:n, σ in (1, -1)
        x = complex(_mesh_cell(j, n), interval(σ * r))
        R_image = min(R_image, inf(abs(imag(T(x)))))
        lo_der = min(lo_der, inf(_cabs(derivative(T, x))))
    end
    return (R_image, Float64(abs(d)) ⊘₊ lo_der)
end

###############################################################################
# The Lasota–Yorke inequality with the L¹ norm in the second term, and the direct
# Galerkin defect, for the weighted ℓ¹ norm of the strip S_ρ of a certified geometry.
# Statements: cor:strip-LY-L1 and lem:strip-deltaN of the Diffusion paper
# (sections/appendix-strip.tex, draft); the coefficient bound is the one in the proof
# of Blumenthal–Nisoli–Taylor-Crush 2025, Lemma 3.7.
###############################################################################

export annulus_coefficient_decay, annulus_dfly_l1, annulus_dfly_l1_smallest_M,
       annulus_galerkin_defect

@doc raw"""
    annulus_coefficient_decay(C::AnnulusCompatibility; T = Float64) -> Interval

The rate ``e^{-2πR}`` of the coefficient bound
```math
|\widehat{Lf}_j| \le \min\bigl\{\|f\|_{L^1},\ e^{-2πR|j|}\sup_{S_r}|f|\bigr\},\qquad j \in \mathbb Z,
```
for the transfer operator ``L`` of the circle map whose lift `C.lift` is certified on the
geometry ``(r, ρ, R)`` of `C.G`.  Since ``T`` has integer degree, ``e^{-2πijT}`` is a function
on the circle of modulus 1 on the real line and
``\widehat{Lf}_j = \int_0^1 f(y)e^{-2πijT(y)}\,dy``, which gives the first bound.  For the
second, the line of integration is moved to ``\operatorname{Im} y = -\operatorname{sign}(jd)\,r``, where
``|e^{-2πijT}| = e^{2πj\operatorname{Im}T} \le e^{-2πR|j|}``: by check 2 of
[`certify_compatibility`](@ref), ``t \mapsto \operatorname{Im}T(x - it)`` has derivative
``-\operatorname{Re}T'(x - it)``, of the sign of ``-d``, on ``[0, r]``, and vanishes at ``t = 0``, so
``\operatorname{Im}T(x \mp ir)`` has the sign of ``\mp d``; by check 3 its modulus exceeds ``R``.
This is the estimate in the proof of [Blumenthal, Nisoli, Taylor-Crush 2025, Lemma 3.7]
read in the variable ``y`` of the strip.  For ``f`` in the weighted ``\ell^1`` space of
``S_ρ``, ``\sup_{S_r}|f| \le \|f\|_{A_ρ}`` since ``r < ρ``.
"""
annulus_coefficient_decay(C::AnnulusCompatibility; T = Float64) =
    exp(-2 * interval(T, π) * as_interval(T, C.G.R))

@doc raw"""
    annulus_dfly_l1(C::AnnulusCompatibility, M; T = Float64) -> NamedTuple

Constants of the Lasota–Yorke inequality of cor:strip-LY-L1 for the seminorm
``V(f) = \sum_{j\ne0}e^{2πρ|j|}|\hat f_j|`` of the strip ``S_ρ`` of `C.G`, with the ``L^1``
norm in the second term,
```math
V(Lf) \le A_MV(f) + (A_M + B'_M)\|f\|_{L^1},\qquad
A_M = \frac{2q^{M+1}}{1 - q},\quad q = e^{-2π(R - ρ)},\quad B'_M = 2\sum_{j=1}^Me^{2πρj},
```
as intervals of type `T`.  Split ``V(Lf)`` at ``|j| = M``: the head is bounded by
``B'_M\|f\|_{L^1}`` and the tail by ``A_M\|f\|_{A_ρ} \le A_M(V(f) + \|f\|_{L^1})``, both by
[`annulus_coefficient_decay`](@ref).  The same split gives the full-norm form
``\|Lf\|_{A_ρ} \le A_M\|f\|_{A_ρ} + (1 + B'_M)\|f\|_{L^1}``, which is
[`analytic_dfly`](@ref)`(Aη(ρ), R, 1, M)` (with ``C = 1`` justified by the coefficient bound).
"""
function annulus_dfly_l1(C::AnnulusCompatibility, M::Integer; T = Float64)
    M >= 0 || throw(ArgumentError("need M ≥ 0"))
    G = C.G
    ρ, R = as_interval(T, G.ρ), as_interval(T, G.R)
    q = exp(-2 * interval(T, π) * (R - ρ))
    w = exp(2 * interval(T, π) * ρ)
    return (; A_M = _tail_factor(q, M), Bp_M = _strip_head(w, M), q = q)
end

@doc raw"""
    annulus_dfly_l1_smallest_M(C::AnnulusCompatibility; Mmax = 2000, T = Float64)
        -> (M, constants)

The smallest ``M \le`` `Mmax` with ``\sup A_M < 1`` in [`annulus_dfly_l1`](@ref), and the
constants at that ``M``; throws if there is none.
"""
function annulus_dfly_l1_smallest_M(C::AnnulusCompatibility; Mmax::Integer = 2000, T = Float64)
    for M = 0:Mmax
        c = annulus_dfly_l1(C, M; T)
        sup(c.A_M) < 1 && return (M, c)
    end
    throw(ArgumentError("A_M ≥ 1 for every M ≤ $Mmax"))
end

@doc raw"""
    annulus_galerkin_defect(C::AnnulusCompatibility, N; M_w = sqrt(C.W), T = Float64) -> Interval

The bound of lem:strip-deltaN on ``δ_N = \|L - Π_NLΠ_N\|_{\mathcal B_s\to L^2}``, ``Π_N`` the
Fourier truncation at ``|j| \le N``, ``\|f\|_s = \|f\|_{L^1} + V(f)`` with ``V`` the seminorm of
[`annulus_dfly_l1`](@ref):
```math
δ_N \le \frac{2e^{-2πR(N+1)}}{1 - e^{-2πR}} + M_we^{-2πρ(N+1)} .
```
The first term bounds ``\|(I - Π_N)Lf\|_{L^2} \le \sum_{|j|>N}|\widehat{Lf}_j|`` by
[`annulus_coefficient_decay`](@ref) and ``\|f\|_{A_ρ} \le \|f\|_s``; the second
``\|Π_NL(I - Π_N)f\|_{L^2} \le M_w\sum_{|j|>N}|\hat f_j| \le M_we^{-2πρ(N+1)}V(f)``.  `M_w` must
bound ``\|L\|_{L^2\to L^2}``; the default ``W^{1/2}`` does, since by Cauchy–Schwarz over the
inverse branches ``|Lf|^2 \le (\sum_jτ_j')(\sum_j|f\circτ_j|^2τ_j')`` and
``\sum_jτ_j' \le W`` on the real line.
"""
function annulus_galerkin_defect(C::AnnulusCompatibility, N::Integer; M_w = nothing, T = Float64)
    N >= 0 || throw(ArgumentError("need N ≥ 0"))
    G = C.G
    ρ, R = as_interval(T, G.ρ), as_interval(T, G.R)
    Mw = M_w === nothing ? sqrt(as_interval(T, C.W)) : as_interval(T, M_w)
    twoπ = 2 * interval(T, π)
    return 2 * exp(-twoπ * R * (N + 1)) / (1 - exp(-twoπ * R)) + Mw * exp(-twoπ * ρ * (N + 1))
end
