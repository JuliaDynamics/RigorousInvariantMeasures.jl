###############################################################################
# Bernstein-ellipse geometry of the analytic bases, and its certification for a map.
#
# The geometry is the triple 1 < r < ρ < R of Assumption ass:analytic of the Diffusion
# paper's single-space appendix (appendix-single-space.tex), chosen by the user: a basis
# built on it carries it and nothing about a dynamic.  `certify_compatibility(G, D)`
# then proves, for one map or for a family given by interval parameters, that every
# inverse branch τ_j extends analytically to E_R with τ_j(E_R) ⊂ E_r, and returns a
# certified W ≥ sup_{E_R} Σ_j |τ_j'|, so that sup_{E_R}|Lf| ≤ W sup_{E_r}|f|.
###############################################################################


@doc raw"""
    EllipseGeometry(r, ρ, R)

The ellipse parameters ``1 < r < ρ < R`` of an analytic basis, in the chart ``t = 2x - 1``:
functions are bounded on ``E_r``, the strong norm is that of ``E_ρ``, and the transfer
operator is to map functions analytic on ``E_r`` to functions analytic on ``E_R``.  The
values are stored as the `Float64`s passed, exactly.  Two geometries built from the same
numbers are different objects, and the bases and certificates built on one are not accepted
with the other.  Whether a map is compatible with the
geometry is proved by [`certify_compatibility`](@ref).
"""
mutable struct EllipseGeometry
    # mutable, with constant fields, so that `===` is identity: a basis and a certificate
    # are checked to share one geometry object, not merely three equal numbers
    const r::Float64
    const ρ::Float64
    const R::Float64
    function EllipseGeometry(r::Real, ρ::Real, R::Real)
        r, ρ, R = Float64(r), Float64(ρ), Float64(R)
        1 < r < ρ < R || throw(ArgumentError("need 1 < r < ρ < R, got $r, $ρ, $R"))
        return new(r, ρ, R)
    end
end

@doc raw"""
    EllipseCompatibility

What [`certify_compatibility`](@ref) proves: the map `dynamic` (a single map or a family
with interval parameters) is compatible with the geometry `G`, every inverse branch maps
``E_R`` into ``E_r``, and ``\sup_{E_R}\sum_j|τ_j'| \le W``.  `R_image` is the certified
lower bound on the Bernstein parameter of the images ``F_j(\partial E_r)``; it exceeds
`G.R`.
"""
struct EllipseCompatibility{DT}
    G::EllipseGeometry
    dynamic::DT
    W::Float64
    R_image::Float64
end

@doc raw"""
    certify_compatibility(G::EllipseGeometry, D::PwMap; n = 1024, cover = 64)
        -> EllipseCompatibility

Proves that the full-branch map `D` is compatible with `G`, branch by branch, or throws an
`ArgumentError`.  The branches may contain `Interval` parameters: every enclosure below
then holds for every map of the family, and one call certifies all of them.  With
``F_j(t) = 2T_j((t+1)/2) - 1`` the branch in the global chart:

1. `D` is full branch (`is_full_branch`);
2. on a cover of the closed ellipse ``E_r`` by `cover × cover` boxes (those whose Bernstein
   parameter can be ``\le r``), ``\operatorname{Re}T_j'`` has the sign of the branch;
3. the Bernstein parameter of ``F_j(\partial E_r)``, enclosed on `n` arcs by
   [`bernstein_expansion`](@ref), is ``> R``.

From 2, on the convex set ``E_r``, ``T_j' \ne 0`` and ``F_j`` is injective
(Noshiro–Warschawski), so ``F_j`` is a homeomorphism of the closed ``E_r`` onto a domain
bounded by the Jordan curve ``F_j(\partial E_r)``.  By 3 that curve lies outside ``E_R``;
``E_R`` is connected and contains ``[-1,1] = F_j([t_j, t_{j+1}])`` by 1, so it lies in
``F_j(E_r)``, and ``τ_j = F_j^{-1}`` maps it into ``E_r``.  Then
``|τ_j'(z)| = 1/|T_j'(τ_j(z))| \le 1/\min_{E_r}|T_j'|``, and since ``T_j'`` is zero-free
on ``E_r`` the minimum modulus principle puts that minimum on ``\partial E_r``, where it is
enclosed on `n` arcs: ``W = \sum_j 1/\min_{\partial E_r}|T_j'|``, rounded up.

The two boundary enclosures are those of `bernstein_geometry`, which looks at
``\partial E_r`` only; check 2 is what makes them valid inside.
"""
function certify_compatibility(G::EllipseGeometry, D::PwMap; n::Integer = 1024,
                               cover::Integer = 64)
    is_full_branch(D) || throw(ArgumentError("the certificate needs a full-branch map"))
    for (jb, br) in enumerate(branches(D))
        _check_branch_derivative_on_ellipse(br, jb, G.r, cover)
    end
    # min over the branches of the image parameter, and Σ_j 1/min_{∂E_r}|T_j'|
    R_image, W = bernstein_geometry(Eρ(G.r), D; n = n)
    R_image > G.R || throw(ArgumentError(
        "the images of ∂E_r are certified only outside E_$R_image, not outside E_R, R = $(G.R)"))
    isfinite(W) || throw(ArgumentError("|T'| is not bounded away from 0 on ∂E_r"))
    return EllipseCompatibility(G, D, W, R_image)
end

export search_ellipse_geometry

@doc raw"""
    search_ellipse_geometry(D::PwMap, rs; cost = (r, R_image, W) -> R_image / r,
                            n = 1024, cover = 64) -> NamedTuple

Tries each ``r > 1`` in `rs` for the full-branch map `D` and returns the one that
maximises `cost` among those that pass the two checks of [`certify_compatibility`](@ref)
on ``E_r``:

* ``\operatorname{Re}T_j'`` has the sign of the branch on the cover of ``E_r`` by
  `cover × cover` boxes, for every branch;
* `R_image > r`, where `R_image` is the certified lower bound on the Bernstein parameter of
  the images ``F_j(\partial E_r)`` enclosed on `n` arcs, minimised over the branches.

For these ``r``, `W` is the certified upper bound ``\sum_j 1/\min_{\partial E_r}|T_j'|``
of `certify_compatibility`, and `cost(r, R_image, W)` is called; the default is
``R_\mathrm{image}/r``.  The result has fields `r`, `R_image`, `W`, `cost` for the best
``r``, and `table`, one named tuple `(r, R_image, W, feasible, cost)` per entry of `rs`, in
order (`R_image = W = NaN` and `cost = -Inf` where the derivative check fails, `cost = -Inf`
where `R_image ≤ r`).  Throws an `ArgumentError` if no ``r`` passes.

The search chooses ``r`` and certifies nothing by itself.  Any ``ρ, R`` with
``r < ρ < R < R_\mathrm{image}`` then pass `certify_compatibility(EllipseGeometry(r, ρ, R), D)`
with the same `n` and `cover`, which is the certificate.
"""
function search_ellipse_geometry(D::PwMap, rs; cost = (r, R_image, W) -> R_image / r,
                                 n::Integer = 1024, cover::Integer = 64)
    is_full_branch(D) || throw(ArgumentError("the search needs a full-branch map"))
    table = NamedTuple{(:r, :R_image, :W, :feasible, :cost),
                       Tuple{Float64,Float64,Float64,Bool,Float64}}[]
    for r in rs
        r = Float64(r)
        r > 1 || throw(ArgumentError("need r > 1, got $r"))
        if any(br -> _branch_derivative_failure(br, r, cover) !== nothing, branches(D))
            push!(table, (r = r, R_image = NaN, W = NaN, feasible = false, cost = -Inf))
            continue
        end
        R_image, W = _bernstein_geometry(Eρ(r), D; n = n)
        feasible = R_image > r && isfinite(W)
        c = feasible ? Float64(cost(r, R_image, W)) : -Inf
        push!(table, (r = r, R_image = R_image, W = W, feasible = feasible, cost = c))
    end
    ok = findall(e -> e.feasible, table)
    isempty(ok) && throw(ArgumentError("no r in the list passes the two checks"))
    best = table[ok[argmax([table[i].cost for i in ok])]]
    return (r = best.r, R_image = best.R_image, W = best.W, cost = best.cost, table = table)
end

# Re T' has the sign of the branch on every box of a cover of the closed ellipse E_r
# (chart t = 2x - 1, x = (t + 1)/2).  Throws otherwise.
function _check_branch_derivative_on_ellipse(br, jb::Integer, r::Float64, m::Integer)
    fail = _branch_derivative_failure(br, r, m)
    fail === nothing || throw(ArgumentError(
        "Re T' of branch $jb may vanish or change sign on E_$r (box $(fail.box), T' ∈ $(fail.d))"))
    return nothing
end

# `nothing` if Re T' has the sign of the branch on every box of the cover of E_r, otherwise
# the first box where that is not certified and the enclosure of T' on it.
function _branch_derivative_failure(br, r::Float64, m::Integer)
    ri = interval(r)
    a = sup((ri + 1 / ri) / 2)     # semi-axes of E_r, rounded up
    b = sup((ri - 1 / ri) / 2)
    # Grid lines shared by neighbouring boxes, the outer ones exactly ±a, ±b, so the boxes
    # cover the rectangle [-a,a]×[-b,b] ⊇ E_r without gaps.
    xs = [i == 0 ? -a : i == m ? a : -a + 2a * i / m for i = 0:m]
    ys = [l == 0 ? -b : l == m ? b : -b + 2b * l / m for l = 0:m]
    for i = 1:m, l = 1:m
        tbox = complex(interval(xs[i], xs[i+1]), interval(ys[l], ys[l+1]))
        inf(bernstein_parameter(tbox)) > r && continue      # certainly outside E_r
        d = derivative(br.f, (tbox + 1) / 2)
        ok = is_increasing(br) ? inf(real(d)) > 0 : sup(real(d)) < 0
        ok || return (box = tbox, d = d)
    end
    return nothing
end

# A certificate must be for the geometry the basis carries, and for the map in use.
function _check_certificate(B, C::EllipseCompatibility, D = C.dynamic)
    B.geometry === nothing && throw(ArgumentError("the basis carries no EllipseGeometry"))
    B.geometry === C.G ||
        throw(ArgumentError("the certificate is for another EllipseGeometry than the basis's"))
    C.dynamic === D || throw(ArgumentError("the certificate is for another dynamic"))
    return nothing
end
