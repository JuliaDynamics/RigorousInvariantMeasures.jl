###############################################################################
# Local evaluators of the transfer operator
#
#     (Lf)(y) = Σ_j f(τ_j(y)) |τ_j'(y)|,
#
# the sum over the inverse branches τ_j of a PwMap whose image contains y:
#
# * `transfer_jet(D, f, Y, k)`: the k-jet of Lf at a point, or uniformly on an interval;
# * `transfer_taylor_model(D, f, I, k)`: a Taylor model of order k of Lf on an interval,
#   with its remainder (TaylorModelsExt, load with `using TaylorModels`);
# * `transfer_taylor_models(D, f, I, k)`: the same on the pieces of I cut at the endpoints
#   of the branch images, where Lf may jump.
#
# Both start from the preimage enclosure `preimage(Y, br)` of each branch.
#
# checked: transfer_jet -- moved here from the untracked src/MatrixFree.jl (never
# committed, no longer included), with _inverse_branch_jet and _pieces_mod1, and
# generalised from full-branch circle maps in Float64 to any PwMap in the interval type
# of Y; MatrixFree.jl is to call these when it comes back.
###############################################################################

import TaylorSeries
using TaylorSeries: Taylor1

export transfer_jet, transfer_taylor_model, transfer_taylor_models

# a thin interval of value v in the number type of X
_c(::Interval{T}, v) where {T} = interval(T, v)

# Split a point or interval enclosure, read modulo 1, into pieces contained in [0, 1].
function _pieces_mod1(c)
    Y = interval(c)
    s = floor(inf(Y))
    Yr = Y - _c(Y, s)
    sup(Yr) - inf(Yr) < 1 ||
        throw(ArgumentError("the enclosure $c is not shorter than the circle"))
    unit = hull(_c(Yr, 0), _c(Yr, 1))
    pieces = [intersect_interval(Yr, unit)]
    sup(Yr) > 1 && push!(pieces, intersect_interval(Yr - _c(Yr, 1), unit))
    return filter(!isempty_interval, pieces)
end

# The pieces of Y on which the branches are read: modulo 1 for a circle map, Y itself for
# an interval map.
_transfer_pieces(::PwMap{<:OnCircle}, Y) = _pieces_mod1(Y)
_transfer_pieces(::PwMap, Y) = [interval(Y)]

# Whether the branch contributes to Lf on all of Y (true) or on none of it (false); throws
# when neither is certified.
function _branch_covers(br::MonotonicBranch, Y::Interval)
    a, b = br.Y[1], br.Y[2]
    lo, hi = br.increasing ? (a, b) : (b, a)
    sup(lo) <= inf(Y) && sup(Y) <= inf(hi) && return true
    (sup(Y) < inf(lo) || inf(Y) > sup(hi)) && return false
    throw(ArgumentError(
        "the enclosure $Y is not certified to lie inside or outside the image " *
        "$(hull(a, b)) of a branch; split it at the image's endpoints"))
end

# σ of order k and sign·σ' of order k, where σ(t) = τ(y + t) − τ(y), from the branch jet
# at X ∋ τ(y).  S(s) = T(x0 + s) − T(x0) has coefficients enclosed by those of
# T(Taylor1([X, 1], k + 1)) with the constant term set to zero; its series inverse is σ.
function _inverse_branch_jet(br::MonotonicBranch, X::Interval, k::Integer)
    S = br.f(Taylor1([X, _c(X, 1)], k + 1))
    S[0] = _c(X, 0)
    σ = TaylorSeries.inverse(S)
    sgn = br.increasing ? 1 : -1
    dσ = Taylor1([_c(X, sgn * (r + 1)) * σ[r+1] for r = 0:k], k)
    return Taylor1(σ.coeffs[1:(k+1)], k), dσ
end

@doc raw"""
    transfer_jet(D::PwMap, f, Y, k; ϵ = 0.0, max_iter = 100) -> Taylor1

Certified Taylor jet of order `k` of ``Lf`` on the point or interval `Y`: for every
``y \in Y`` and ``m = 0, \dots, k`` the coefficient of ``t^m`` contains
``(Lf)^{(m)}(y)/m!``, where ``(Lf)(y) = \sum_j f(τ_j(y))\,|τ_j'(y)|``, the sum over the
branches ``j`` of `D` whose image contains ``y``.  So `factorial(m) * J[m]` encloses the
`m`-th derivative of ``Lf`` at ``y``.

For a circle map (`PwMap{OnCircle}`, see [`CircleMap`](@ref)) `Y` is read modulo 1, and an
interval that crosses ``0 \sim 1`` is split there, the jets of the pieces joined by hull.
A branch contributes when `Y` is certified to lie in its image and not when certified to
lie outside; otherwise an `ArgumentError` is thrown, and `Y` should be split at the
endpoints of the branch images.  `f` must be of class ``C^k`` on the preimages and accept
`Taylor1` with interval coefficients; each branch of `D` must be of class ``C^{k+1}`` and
accept the same.  The arithmetic is that of the interval type of `Y`; a real `Y` is read
as a thin `Interval{Float64}`.

The argument, for one branch and one piece of `Y`.  `X = preimage(Y, br)` contains
``x_0 = τ(y)`` for every ``y \in Y``.  Expanding the branch on `Taylor1([X, 1], k + 1)` gives
intervals that contain ``T^{(m)}(x_0)/m!``; with the constant term set to zero this encloses
the coefficients of ``S(s) = T(x_0 + s) - T(x_0)``.  Its series inverse
``σ(t) = τ(y + t) - τ(y)`` is computed by TaylorSeries' `inverse`, the Lagrange formula
``σ_m = [s^{m-1}](s/S)^m / m``: ring operations, a division by ``S'(0) = T'(x_0) \ne 0`` and
divisions by integers, so each coefficient of ``σ`` encloses the true one.  The jet of ``f``
on `Taylor1([X, 1], k)` encloses ``f^{(m)}(x_0)/m!``; composing it with ``σ`` truncated at
order ``k`` and multiplying by ``\pm σ'`` (the sign of the branch) gives the jet of
``f(τ(y+t))\,|τ'(y+t)|`` to order ``k``, using ``f`` to order ``k`` and ``T`` to order
``k + 1``.  These are polynomial operations on enclosures of quantities that, for each
``y``, are the true ones at ``x_0 = τ(y)``, so the result encloses the true coefficients for
every ``y \in Y``.
"""
function transfer_jet(D::PwMap, f, Y, k::Integer; ϵ = 0.0, max_iter = 100)
    k >= 0 || throw(ArgumentError("need k ≥ 0, got $k"))
    J = nothing
    for Yp in _transfer_pieces(D, Y)
        Jp = Taylor1(fill(_c(Yp, 0), k + 1), k)
        for br in branches(D)
            _branch_covers(br, Yp) || continue
            X = preimage(Yp, br; ϵ, max_iter)
            σ, dσ = _inverse_branch_jet(br, X, k)
            Jp += TaylorSeries.evaluate(f(Taylor1([X, _c(X, 1)], k)), σ) * dσ
        end
        J = J === nothing ? Jp.coeffs : hull.(J, Jp.coeffs)
    end
    return Taylor1(J, k)
end

@doc raw"""
    transfer_taylor_model(D::PwMap, f, I::Interval, k; ϵ = 0.0, max_iter = 100,
                          newton_steps = 20) -> TaylorModels.TaylorModel1

A Taylor model ``(P, Δ)`` of order `k` of ``Lf`` on the interval `I`, expanded at the
midpoint ``y_0`` of `I`: ``(Lf)(y) - P(y - y_0) \in Δ`` for every ``y \in I``.  Defined in
the `TaylorModels` extension; load it with `using TaylorModels`.  The hypotheses on `D`,
`f` and `I` are those of [`transfer_jet`](@ref) with `Y = I`, with `I` inside ``[0,1]``
also for a circle map; `f` and the branches must also accept real `TaylorModel1`
arguments (TaylorModels has `+ - * /` and `inv, sqrt, exp, log, sin, cos, tan, asin,
acos, atan, sinh, cosh, tanh` on them).

The argument, for one branch ``T`` with inverse ``τ`` on `I`.  Start from the model with
polynomial ``p_0 = \operatorname{mid}X`` and remainder ``Δ_0 = X - \operatorname{mid}X``,
`X = preimage(I, br)`, which contains ``τ``.  With a constant ``c`` (an enclosure of
``1/T'(\operatorname{mid}X)``, a free choice) and ``N(g) = g - c\,(T∘g - \mathrm{id})``, the
true ``τ`` satisfies ``τ = N(τ)``.  If ``τ - p \in Δ`` on `I`, the mean value theorem gives,
for each ``y``, ``N(τ)(y) - N(p)(y) = (1 - c\,T'(ξ))\,(τ(y) - p(y))`` with ``ξ`` between
``τ(y)`` and ``p(y)``, so in the hull ``R`` of the ranges of ``(p, Δ)`` and of ``p`` on `I`.
Hence ``τ - p' \in Δ_Q + K Δ``, where ``(p', Δ_Q)`` is the model of ``N(p)`` formed by
TaylorModels and ``K = 1 - c\,T'(R)``, ``T'(R)`` enclosed by an interval Taylor series.
This step is repeated `newton_steps` times, or until ``Δ`` stops shrinking; every iterate
encloses ``τ``.  Then ``|τ'| = \pm 1/T'(τ)``: the model of ``T'∘τ`` is TaylorModels' `rpa`
of ``T'`` on the model of ``τ``, with ``T'`` acting on a series ``u`` as
``(T∘\tilde u)'/\tilde u'``, ``\tilde u`` being ``u`` extended by one zero coefficient
(exact whenever ``u'(0) \ne 0``, as for the argument ``x_0 + s`` that `rpa` passes), and its
reciprocal is `inv` of that model.  The model of the branch's term is
``f(τ)\,(\pm 1/T'(τ))``, and the terms are summed.
"""
function transfer_taylor_model end

@doc raw"""
    transfer_taylor_models(D::PwMap, f, I::Interval, k; ϵ = 0.0, max_iter = 100,
                           newton_steps = 20, parts = 1)
        -> Vector{NamedTuple{(:piece, :cells)}}

Taylor models of ``Lf`` on `I`, organised by the cuts where ``Lf`` may jump.  Defined in the
`TaylorModels` extension.  `I` is cut at the endpoints of the branch images that lie inside
it; the result lists the pieces between cuts in order, covering `I`, and each piece ``K``
is covered by `cells`, a list of `(piece, model)` in order: on each cell ``C``,
``(Lf)(y) - P(y - y_C) \in Δ`` for every ``y`` in the interior of ``C``, ``y_C`` the
midpoint of ``C``.  Each ``K`` is partitioned into `parts` equal cells, whose endpoints are
free choices and not cuts: ``Lf`` is continuous across them when ``f`` is.

An endpoint of a branch image is an interval enclosure ``E``.  When ``E`` is thin it is a
cut point.  When it is not, the sliver ``E \cap I`` is a piece of its own with
`cells = nothing`: whether that branch contributes is not decided there, and ``Lf`` is not
modelled on it.  On every other piece the set of contributing branches is the same at every
interior point, and the model of a cell is the sum of the branch models of
[`transfer_taylor_model`](@ref) on it.  At a cut that is the endpoint of a branch image the
value of ``Lf`` is not claimed by either neighbouring piece.
"""
function transfer_taylor_models end

@doc raw"""
    transfer_taylor_models(D::PwMap, F::AbstractVector{<:NamedTuple{(:piece, :cells)}},
                           I::Interval, k; ϵ = 0.0, max_iter = 100, newton_steps = 20,
                           parts = 1)

The same for a piecewise ``f`` given in that form, for instance the output of a previous
call, so that the call can be iterated.  `F` is read as: ``f`` is continuous on each piece
with cells, and on each cell ``f(x) - p_q(x - x_q) \in Δ_q`` for every ``x`` in the interior
of the cell, ``x_q`` the expansion point of its model, of order `k`; on the pieces with
`cells = nothing` and at the endpoints of the pieces, ``f`` is not given.

The cuts are the branch image endpoints and, for each branch, the images ``T(c)`` of the
endpoints ``c`` of the pieces of `F` that lie in its domain.  On a piece ``K`` of the output
that is not a sliver no ``T(c)`` lies in the interior, so each contributing branch maps the
interior of ``K`` into the interior of the one piece of `F` that contains ``τ`` of the
midpoint of ``K`` (it is connected, contains that point and meets no ``c``).  If that piece
has no cells, ``K`` gets none; if no piece of `F` contains it, an `ArgumentError` is thrown.

The cells of ``K`` are bounded by the points of its partition into `parts` equal intervals
and by the images ``T(c')`` of the interior cell boundaries ``c'`` of those pieces of `F`.  On a cell ``C`` that contains no ``T(c')``,
``τ`` maps its interior into the interior of one input cell, found from ``τ`` of the
midpoint, and ``f∘τ`` is enclosed by that cell's polynomial evaluated on the model of
``τ - x_q``, plus ``Δ_q``.  When ``T(c')`` is a non-thin enclosure, the cell ``C`` it gives
is enclosed by a constant model: ``f∘τ`` on ``C`` lies in the hull of the ranges of the
input cells that ``τ(C)`` meets, ``f`` being continuous at ``c'``, and ``|τ'|`` in
``1/|T'|`` on the preimage of ``C``.
"""
transfer_taylor_models(::PwMap, ::AbstractVector, ::Interval, ::Integer)
