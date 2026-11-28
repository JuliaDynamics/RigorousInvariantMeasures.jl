# The Taylor model of Lf on an interval: `transfer_taylor_model`, documented in
# src/TransferJets.jl, where the argument is written out.

import RigorousInvariantMeasures: transfer_taylor_model, transfer_taylor_models
using RigorousInvariantMeasures: MonotonicBranch, PwMap, OnCircle, branches, preimage,
    derivative, _branch_covers, _c
using TaylorModels: TaylorModel1, polynomial, remainder, rpa, expansion_point, get_order
import TaylorSeries
using TaylorSeries: Taylor1

# T' acting on a series u with u'(0) ≠ 0: (T∘ũ)'/ũ', ũ being u extended by one zero
# coefficient, so that the derivative loses no order.  This is the series of T'∘u to the
# order of u.
function _derivative_series(T, u::Taylor1)
    n = TaylorSeries.get_order(u)
    ũ = Taylor1(vcat(u.coeffs, zero(u[0])), n + 1)
    return TaylorSeries.derivative(T(ũ)) / TaylorSeries.derivative(ũ)
end

# A Taylor model of the inverse branch τ on I, from the enclosure X of τ(I), by the
# mean-value form of the step g ↦ g − c(T∘g − id) (see the docstring of
# transfer_taylor_model).
function _inverse_branch_taylor_model(br::MonotonicBranch, X, x0, I, k, newton_steps)
    id = TaylorModel1(k, x0, I)
    cdom = I - x0
    zΔ = _c(X, 0)
    xm = _c(X, mid(X))
    c = _c(X, mid(inv(derivative(br.f, xm))))     # a free choice: a thin c ≈ 1/T'(mid X)
    p = Taylor1([xm; fill(zΔ, k)], k)
    Δ = X - xm
    for _ = 1:newton_steps
        ptm = TaylorModel1(p, zΔ, x0, I)
        Q = ptm - c * (br.f(ptm) - id)
        pI = p(cdom)
        R = hull(pI + Δ, pI)
        K = 1 - c * derivative(br.f, R)
        Δn = remainder(Q) + K * Δ
        diam(Δn) < diam(Δ) || break
        p, Δ = polynomial(Q), Δn
    end
    return TaylorModel1(p, Δ, x0, I)
end

# The sum over `brs` of the models f(τ)·|τ'| of the branches on I.
function _sum_branch_models(brs, f, I, k, ϵ, max_iter, newton_steps)
    x0 = _c(I, mid(I))
    S = TaylorModel1(Taylor1(fill(_c(I, 0), k + 1), k), _c(I, 0), x0, I)
    for br in brs
        X = preimage(I, br; ϵ, max_iter)
        P = _inverse_branch_taylor_model(br, X, x0, I, k, newton_steps)
        dT = rpa(u -> _derivative_series(br.f, u), P)
        sgn = br.increasing ? 1 : -1
        S += f(P) * (sgn * inv(dT))
    end
    return S
end

function _check_model_domain(D, I, k)
    k >= 0 || throw(ArgumentError("need k ≥ 0, got $k"))
    if D isa PwMap{<:OnCircle}
        0 <= inf(I) && sup(I) <= 1 ||
            throw(ArgumentError("for a circle map I must lie in [0, 1], got $I"))
    end
    return nothing
end

function transfer_taylor_model(D::PwMap, f, I::Interval, k::Integer; ϵ = 0.0,
                               max_iter = 100, newton_steps::Integer = 20)
    _check_model_domain(D, I, k)
    brs = [br for br in branches(D) if _branch_covers(br, I)]
    return _sum_branch_models(brs, f, I, k, ϵ, max_iter, newton_steps)
end

# Whether the branch contributes at every interior point of J (true) or at none (false),
# for a piece J whose interior meets no endpoint of the branch image.
function _branch_covers_interior(br::MonotonicBranch, J::Interval)
    lo, hi = br.increasing ? (br.Y[1], br.Y[2]) : (br.Y[2], br.Y[1])
    sup(lo) <= inf(J) && sup(J) <= inf(hi) && return true
    (sup(J) <= inf(lo) || inf(J) >= sup(hi)) && return false
    error("a branch image endpoint lies inside the piece $J; this should not happen")
end

# Cut I at the enclosures `encls`: returns the pieces in order, each with `true` when it
# is a sliver (the part of I inside a non-thin enclosure) and `false` otherwise.
function _cut_pieces(I::Interval, encls)
    a, b = inf(I), sup(I)
    cuts = Tuple{typeof(a),typeof(a)}[]
    for E in encls
        inf(E) < b && sup(E) > a || continue
        push!(cuts, (max(inf(E), a), min(sup(E), b)))
    end
    sort!(cuts)
    merged = Tuple{typeof(a),typeof(a)}[]
    for c in cuts
        if !isempty(merged) && c[1] <= merged[end][2]
            merged[end] = (merged[end][1], max(merged[end][2], c[2]))
        else
            push!(merged, c)
        end
    end
    pieces = Tuple{typeof(I),Bool}[]
    left = a
    for (lo, hi) in merged
        lo > left && push!(pieces, (interval(left, lo), false))
        hi > lo && push!(pieces, (interval(lo, hi), true))
        left = hi
    end
    left < b && push!(pieces, (interval(left, b), false))
    return pieces
end

_image_endpoints(D) = [E for br in branches(D) for E in br.Y]

# `parts` equal intervals covering J.  The partition points are free choices, shared by
# neighbouring intervals.
function _partition(J::Interval, parts::Integer)
    parts >= 1 || throw(ArgumentError("need parts ≥ 1, got $parts"))
    a, b = inf(J), sup(J)
    xs = [i == 0 ? a : i == parts ? b : a + (b - a) * i / parts for i = 0:parts]
    return [interval(xs[i], xs[i+1]) for i = 1:parts]
end

_zero_model(J, k) = TaylorModel1(Taylor1(fill(_c(J, 0), k + 1), k), _c(J, 0), _c(J, mid(J)), J)

# The model on J of the constant interval V: polynomial mid V, remainder V - mid V.
function _constant_model(V::Interval, J, k)
    m = _c(J, mid(V))
    return TaylorModel1(Taylor1([m; fill(_c(J, 0), k)], k), V - m, _c(J, mid(J)), J)
end

function transfer_taylor_models(D::PwMap, f, I::Interval, k::Integer; ϵ = 0.0,
                                max_iter = 100, newton_steps::Integer = 20,
                                parts::Integer = 1)
    _check_model_domain(D, I, k)
    return map(_cut_pieces(I, _image_endpoints(D))) do (K, sliver)
        sliver && return (piece = K, cells = nothing)
        brs = [br for br in branches(D) if _branch_covers_interior(br, K)]
        cells = [(piece = C,
                  model = _sum_branch_models(brs, f, C, k, ϵ, max_iter, newton_steps))
                 for C in _partition(K, parts)]
        return (piece = K, cells = cells)
    end
end

# The element of `list` (with a field `piece`) whose closed piece contains the
# enclosure X; `nothing` if there is none.
function _containing(list, X)
    for p in list
        inf(p.piece) <= inf(X) && sup(X) <= sup(p.piece) && return p
    end
    return nothing
end

# f∘τ·|τ'| on a cell C of the output, for one branch: by the model of the input cell
# containing τ(C), or, when C is a sliver around T(c') for an interior cell boundary c',
# by a constant model from the ranges of the input cells that τ(C) meets.
function _branch_term(br, Fc, C, csliver, k, ϵ, max_iter, newton_steps)
    sgn = br.increasing ? 1 : -1
    X = intersect_interval(preimage(C, br; ϵ, max_iter), Fc.piece)
    if csliver
        V = emptyinterval(X)
        for q in Fc.cells
            Y = intersect_interval(X, q.piece)
            isempty_interval(Y) && continue
            V = hull(V, q.model(Y - expansion_point(q.model)))
        end
        return _constant_model(V * (sgn * inv(derivative(br.f, X))), C, k)
    end
    q = _containing(Fc.cells, preimage(_c(C, mid(C)), br; ϵ, max_iter))
    q === nothing && error("no input cell contains the preimage of $C; this should not happen")
    X = intersect_interval(X, q.piece)
    P = _inverse_branch_taylor_model(br, X, _c(C, mid(C)), C, k, newton_steps)
    fP = q.model(P - expansion_point(q.model)) + remainder(q.model)
    dT = rpa(u -> _derivative_series(br.f, u), P)
    return fP * (sgn * inv(dT))
end

function transfer_taylor_models(D::PwMap, F::AbstractVector{<:NamedTuple{(:piece, :cells)}},
                                I::Interval, k::Integer; ϵ = 0.0, max_iter = 100,
                                newton_steps::Integer = 20, parts::Integer = 1)
    _check_model_domain(D, I, k)
    all(Fc -> Fc.cells === nothing || all(q -> get_order(q.model) == k, Fc.cells), F) ||
        throw(ArgumentError("the models of the piecewise f must have order k = $k"))
    indomain(br, c) = in_interval(c, hull(br.X[1], br.X[2]))
    # the cuts: the branch image endpoints, and the images under each branch of the
    # endpoints of the continuity pieces of F that lie in its domain
    encls = _image_endpoints(D)
    for br in branches(D), Fc in F, c in (inf(Fc.piece), sup(Fc.piece))
        indomain(br, c) && push!(encls, br.f(_c(I, c)))
    end
    return map(_cut_pieces(I, encls)) do (K, sliver)
        sliver && return (piece = K, cells = nothing)
        brs = [br for br in branches(D) if _branch_covers_interior(br, K)]
        # for each branch, the continuity piece of F containing τ(K)
        Fcs = map(brs) do br
            Fc = _containing(F, preimage(_c(K, mid(K)), br; ϵ, max_iter))
            Fc === nothing && throw(ArgumentError(
                "the piecewise f does not cover the preimage of $K under a branch"))
            Fc
        end
        any(Fc -> Fc.cells === nothing, Fcs) && return (piece = K, cells = nothing)
        # cell boundaries: the points of the partition of K into `parts`, and the images of
        # the interior cell boundaries of those pieces of F
        bounds = [_c(K, inf(C)) for C in _partition(K, parts)[2:end]]
        for (br, Fc) in zip(brs, Fcs), q in Fc.cells[2:end]
            c = inf(q.piece)
            indomain(br, c) && push!(bounds, br.f(_c(K, c)))
        end
        cells = NamedTuple{(:piece, :model)}[]
        for (C, csliver) in _cut_pieces(K, bounds)
            S = _zero_model(C, k)
            for (br, Fc) in zip(brs, Fcs)
                S += _branch_term(br, Fc, C, csliver, k, ϵ, max_iter, newton_steps)
            end
            push!(cells, (piece = C, model = S))
        end
        return (piece = K, cells = cells)
    end
end
