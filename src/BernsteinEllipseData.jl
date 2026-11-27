###############################################################################
# Bernstein-ellipse data of a full-branch analytic map, in the global chart t = 2x - 1.
#
# The data are those of Assumption ass:analytic of the Diffusion paper's single-space
# appendix (appendix-single-space.tex): 1 < r < R and W with sup_{E_R}|Lf| ≤ W sup_{E_r}|f|,
# which holds when every inverse branch τ_j extends analytically to E_R with
# τ_j(E_R) ⊂ E_r and W ≥ sup_{E_R} Σ_j |τ_j'|.  The paper's ρ, with r < ρ < R, is stored
# beside them.
###############################################################################

import .RigorousInvariantMeasures: bernstein_geometry, bernstein_parameter, derivative

@doc raw"""
    BernsteinEllipseData

The ellipse data `(r, ρ, R, W)` of a full-branch map ``T`` of ``[0,1]`` in the chart
``t = 2x - 1``, certified by [`bernstein_ellipse_data`](@ref): every inverse branch
``τ_j`` extends analytically to ``E_R`` with ``τ_j(E_R) \subset E_r``, and
``\sup_{E_R}\sum_j|τ_j'| \le W``; and ``r < ρ < R``.  `R` is a lower bound and `W` an
upper bound, both exactly representable `Float64`s.
"""
struct BernsteinEllipseData
    r::Float64
    ρ::Float64
    R::Float64
    W::Float64
end

@doc raw"""
    bernstein_ellipse_data(D::PwMap; r, ρ, n = 1024, cover = 64) -> BernsteinEllipseData

Certify the ellipse data of [`BernsteinEllipseData`](@ref) for the full-branch map `D`.

`R` and `W` are the two numbers of [`bernstein_geometry`](@ref)`(Eρ(r), D; n)`: `R` is a
lower bound on the Bernstein parameter of ``F_j(\partial E_r)``, ``F_j`` the branch in the
global chart, and `W` is ``\sum_j 1/\min_{\partial E_r}|T_j'|``.  That routine looks at
``\partial E_r`` only.  Before its output is used this function checks, on a cover of the
closed ellipse ``E_r`` by `cover × cover` boxes (those whose Bernstein parameter can be
``\le r``), that ``\operatorname{Re}T_j'`` has the sign of the branch on every box, and
throws an `ArgumentError` otherwise.  On the convex set ``E_r`` this gives:

1. ``T_j' \ne 0`` on ``E_r``, so by the minimum modulus principle
   ``\min_{E_r}|T_j'| = \min_{\partial E_r}|T_j'|`` and `W` bounds
   ``\sum_j 1/|T_j'\circ τ_j| = \sum_j|τ_j'|`` wherever ``τ_j`` takes values in ``E_r``;
2. ``F_j`` is injective on ``E_r`` (``\operatorname{Re}F_j' > 0`` or ``< 0`` on a convex
   domain, Noshiro–Warschawski), so ``F_j(E_r)`` is the domain bounded by the Jordan curve
   ``F_j(\partial E_r)``, which lies outside the open ellipse of parameter `R`; that ellipse
   is connected and meets ``F_j(E_r)`` (it contains ``[-1,1] = F_j([t_j, t_{j+1}])``, the
   branch being full), so it lies in ``F_j(E_r)`` and ``τ_j = F_j^{-1}`` maps it into ``E_r``.

Throws unless `D` is full branch and ``r < ρ < R``.
"""
function bernstein_ellipse_data(D::PwMap; r::Real, ρ::Real, n::Integer = 1024,
                                cover::Integer = 64)
    is_full_branch(D) || throw(ArgumentError("the ellipse data need a full-branch map"))
    1 < r < ρ || throw(ArgumentError("need 1 < r < ρ, got r = $r, ρ = $ρ"))
    _check_branch_derivatives_on_ellipse(D, Float64(r), cover)
    R, W = bernstein_geometry(Eρ(Float64(r)), D; n = n)
    ρ < R || throw(ArgumentError("need ρ < R; the certified R is $R for r = $r"))
    return BernsteinEllipseData(Float64(r), Float64(ρ), Float64(R), Float64(W))
end

# Re T_j' has the sign of the branch on every box of a cover of the closed ellipse E_r
# (chart t = 2x - 1, x = (t + 1)/2).  Throws otherwise.
function _check_branch_derivatives_on_ellipse(D::PwMap, r::Float64, m::Integer)
    ri = interval(r)
    a = sup((ri + 1 / ri) / 2)     # semi-axes of E_r, rounded up
    b = sup((ri - 1 / ri) / 2)
    # Grid lines shared by neighbouring boxes, with the outer ones exactly ±a, ±b, so the
    # boxes cover the rectangle [-a,a]×[-b,b] ⊇ E_r without gaps.
    xs = [i == 0 ? -a : i == m ? a : -a + 2a * i / m for i = 0:m]
    ys = [l == 0 ? -b : l == m ? b : -b + 2b * l / m for l = 0:m]
    for (jb, br) in enumerate(branches(D)), i = 1:m, l = 1:m
        re = interval(xs[i], xs[i+1])
        im_ = interval(ys[l], ys[l+1])
        # Skip boxes certainly outside E_r.
        tbox = complex(re, im_)
        inf(bernstein_parameter(tbox)) > r && continue
        d = derivative(br.f, (tbox + 1) / 2)
        ok = is_increasing(br) ? inf(real(d)) > 0 : sup(real(d)) < 0
        ok || throw(ArgumentError(
            "Re T' of branch $jb may vanish or change sign on E_$r (box $tbox, " *
            "T' ∈ $d); the boundary values of bernstein_geometry do not bound it there"))
    end
    return nothing
end
