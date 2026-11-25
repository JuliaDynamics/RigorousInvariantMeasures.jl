using Symbolics
using SymbolicUtils
using IntervalArithmetic
using TaylorSeries: Taylor1

# This file adds methods to the package's own `dfly` and `WeightedSobolev`.
import RigorousInvariantMeasures: dfly, WeightedSobolev

@doc raw"""
    SymbL(n, f)

A symbolic transfer operator applied to a symbolic function,
``(L_n f)(x) = \sum_{y\in T^{-1}(x)} f(y)/|T'(y)|^n``.
"""
struct SymbL
    n::Any
    f::Num
end

@doc raw"""
    Diff(P::SymbL, ∂, Dist)

The derivative of ``L_n g`` as a sum of two symbolic transfer operators,

```math
\frac{d}{dx} L_n g = L_{n+1}\, g' + n\, L_n\big(g\,(1/T')'\big),
```

where `∂` is the derivative in the variable ``y`` of the preimage and `Dist` stands for
``(1/T')'(y)``.  The identity is exact on increasing branches; on a decreasing branch the
first term carries a factor ``\operatorname{sgn} T' = -1``, which the ``L^1`` bounds of
[`dfly_coefficient_matrix`](@ref) do not see.  Butterley, Kiamari and Liverani, *Locating
Ruelle-Pollicott resonances*, Eq. 3.2.
"""
Diff(P::SymbL, ∂::Differential, Dist) =
    [SymbL(P.n + 1, ∂(P.f)), SymbL(P.n, P.n * P.f * Dist)]

"""
    Diff(v::Vector{SymbL}, ∂, Dist)

The derivative of the sum of the operators in `v`, as the concatenation of the
derivatives of the summands.
"""
function Diff(v::Vector{SymbL}, ∂::Differential, Dist)
    w = SymbL[]
    for P in v
        append!(w, Diff(P, ∂, Dist))
    end
    return w
end

@doc raw"""
    compute_dfly_k_fi_DDi(k)

The ``k``-th derivative of ``Lf`` as a vector of symbolic transfer operators
``L_n(\cdot)``, obtained by applying [`Diff`](@ref) ``k`` times to ``L_1 f``.  In the
output `f[i]` stands for ``f^{(i-1)}`` and `DD[m]` for ``(1/T')^{(m)}``, both functions
of the preimage ``y``.  For ``k = 0`` the output is the single `SymbL(1, f[1])`.
"""
function compute_dfly_k_fi_DDi(k::Int)
    @variables x
    @variables (f(x))[1:k+1]
    @variables Dist(x)
    @variables (DD(x))[1:k+1]

    ∂ = Differential(x)

    der_dict = Dict(
        [
            [Dist => DD[1]]
            [∂(f[i]) => f[i+1] for i = 1:k]
            [∂(DD[i]) => DD[i+1] for i = 1:k]
        ],
    )

    v = SymbL(1, f[1])
    for _ = 1:k
        v = Diff(v, ∂, Dist)
        for j = 1:length(v)
            temporary = expand_derivatives(v[j].f, true)
            temporary = substitute(temporary, der_dict)
            temporary = simplify(temporary; expand = true)
            v[j] = SymbL(v[j].n, temporary)
        end
    end
    return v
end

###############################################################################
# From the symbolic expansion to integer monomials
###############################################################################

# The summands of an expanded expression, as (coefficient, product) pairs.
function _summands(e)
    e = Symbolics.unwrap(e)
    if SymbolicUtils.isadd(e)
        iszero(e.coeff) ||
            error("compute_dfly_k_fi_DDi produced a constant summand $(e.coeff)")
        return [(c, t) for (t, c) in e.dict]
    end
    return [(1, e)]
end

# The factors of a product, as (coefficient, Dict(factor => exponent)).
function _factors(h)
    SymbolicUtils.ismul(h) && return (h.coeff, Dict(h.dict))
    SymbolicUtils.ispow(h) && return (1, Dict(h.base => h.exp))
    return (1, Dict(h => 1))
end

function _exact_int(c, what)
    (c isa Number && isinteger(c)) ||
        error("non-integer coefficient $c in $what of the Faà di Bruno expansion")
    return Int(c)
end

@doc raw"""
    dfly_monomials(j) -> Vector{NamedTuple{(:n, :i, :coeff, :p)}}

The expansion of [`compute_dfly_k_fi_DDi`](@ref)`(j)` read as a list of monomials

```math
(Lf)^{(j)} = \sum_t \mathrm{coeff}_t\, L_{n_t}\Big(f^{(i_t)} \prod_{m=1}^{j}
\big((1/T')^{(m)}\big)^{p_{t,m}}\Big),
```

with integer `coeff`, `n`, `i` and exponents `p` (a vector of length `j`).  Every factor of
every summand must be recognised as one `f[i]` to the first power or a power of one
`DD[m]`; anything else throws, so that a change in the output format of Symbolics cannot
drop a term in silence.
"""
function dfly_monomials(j::Int)
    @variables x
    @variables (f(x))[1:j+1]
    @variables (DD(x))[1:j+1]
    fs = [Symbolics.unwrap(s) for s in collect(f)]
    DDs = [Symbolics.unwrap(s) for s in collect(DD)]

    v = compute_dfly_k_fi_DDi(j)
    v isa SymbL && (v = [v])

    out = NamedTuple{(:n, :i, :coeff, :p),Tuple{Int,Int,Int,Vector{Int}}}[]
    for P in v
        for (c, h) in _summands(P.f)
            cf, dict = _factors(h)
            coeff = _exact_int(c, "a summand") * _exact_int(cf, "a product")
            i = -1
            p = zeros(Int, j)
            for (fac, e) in dict
                e = _exact_int(e, "an exponent")
                fi = findfirst(s -> isequal(s, fac), fs)
                mi = findfirst(s -> isequal(s, fac), DDs)
                if fi !== nothing
                    (i == -1 && e == 1) ||
                        error("summand $h of the expansion is not linear in f")
                    i = fi - 1
                elseif mi !== nothing
                    mi <= j || error("summand $h contains DD[$mi] beyond order $j")
                    p[mi] += e
                else
                    error("unrecognised factor $fac in summand $h of the expansion")
                end
            end
            i >= 0 || error("summand $h of the expansion contains no derivative of f")
            push!(out, (n = Int(P.n), i = i, coeff = coeff, p = p))
        end
    end
    return out
end

###############################################################################
# Certified enclosures of the derivatives of 1/T' over boxes
###############################################################################

@doc raw"""
    inverse_derivative_enclosures(D::PwMap, k; boxes = 1024)

Interval enclosures, over a cover of the domain of every branch by `boxes` boxes per
branch, of ``(1/T')^{(m)}/m!`` for ``m = 0, \dots, k``.  Returns a vector with one entry per
box, each a vector `h` with `h[m+1]` enclosing ``(1/T')^{(m)}(y)/m!`` for every ``y`` in
the box.

On a box ``X`` the branch is evaluated on the interval Taylor series ``X + t`` of order
``k+1``; the coefficient of ``t^r`` of the result encloses ``T^{(r)}(y)/r!`` for every
``y \in X``, since each coefficient is the interval evaluation of an algebraic expression
in ``y``.  The series of ``1/T'`` is formed from that of ``T'`` by the division recurrence
of TaylorSeries, again in interval arithmetic, so its coefficients enclose
``(1/T')^{(m)}(y)/m!`` for every ``y \in X``.  The boxes are
``X_1 + (X_2 - X_1)\,[\frac{b-1}{N}, \frac{b}{N}]`` for the endpoint enclosures
``X_1, X_2`` of the branch domain; consecutive boxes share a floating-point endpoint, so
their union covers the domain.

The branch functions must accept `Taylor1{Interval{Float64}}`; a branch whose derivatives
are given only through `@define_with_derivatives` is not supported.  Throws if an
enclosure of ``T'`` contains zero or is unbounded.
"""
function inverse_derivative_enclosures(D::PwMap, k::Int; boxes::Int = 1024)
    out = Vector{Vector{Interval{Float64}}}()
    for br in branches(D)
        X1, X2 = br.X
        for b = 1:boxes
            X = X1 + (X2 - X1) * interval((b - 1) / boxes, b / boxes)
            t = Taylor1([X, interval(1.0)], k + 1)
            y = br.f(t)
            dT = Taylor1([(r + 1) * y[r+1] for r = 0:k], k)
            (all(isbounded, dT.coeffs) && !in_interval(0, dT[0])) || error(
                "the enclosure $(dT[0]) of T' on the box $X contains zero or is unbounded",
            )
            h = 1 / dT
            push!(out, [h[m] for m = 0:k])
        end
    end
    return out
end

# The enclosure, on one box, of |coeff| * |Π_m ((1/T')^{(m)})^{p_m}| * |1/T'|^{n-1}.
function _monomial_on_box(h, n, coeff, p)
    v = interval(abs(coeff)) * abs(h[1])^(n - 1)
    for (m, e) in enumerate(p)
        e == 0 && continue
        v *= abs(interval(factorial(m)) * h[m+1])^e
    end
    return v
end

@doc raw"""
    distortion_vals(D::PwMap, k; boxes = 1024) -> Matrix{Float64}

Upper bounds `vals[i+1, n]` of

```math
\sup_y \big|(1/T')^{(i)}(y)\big|\, \big|1/T'(y)\big|^{n-1},
\qquad i = 0,\dots,k,\quad n = 1,\dots,k+1,
```

over all branches of `D`, the supremum taken jointly on each box of
[`inverse_derivative_enclosures`](@ref) rather than as a product of separate suprema.
Each entry is the upper endpoint of an interval enclosure.
"""
function distortion_vals(D::PwMap, k::Int; boxes::Int = 1024)
    H = inverse_derivative_enclosures(D, k; boxes = boxes)
    vals = zeros(k + 1, k + 1)
    for i = 0:k, n = 1:k+1
        p = zeros(Int, k)
        i >= 1 && (p[i] = 1)
        # for i = 0 the factor (1/T')^{(0)} = 1/T' is one more power of |1/T'|
        np = i == 0 ? n + 1 : n
        vals[i+1, n] = maximum(sup(_monomial_on_box(h, np, 1, p)) for h in H)
    end
    return vals
end

@doc raw"""
    dfly_coefficient_matrix(D::PwMap, k; boxes = 1024) -> Matrix{Float64}

Upper bounds `C[j+1, i+1]` ``= c_{j,i}``, for ``0 \le i \le j \le k``, such that for every
``f \in W^{k,1}``

```math
\|(Lf)^{(j)}\|_{L^1} \le \sum_{i=0}^{j} c_{j,i}\, \|f^{(i)}\|_{L^1},
```

where ``L`` is the transfer operator of `D`; the entries with ``i > j`` are zero.

Statement and proof.  Let every branch of `D` be full and ``C^{k+1}`` with ``T' \ne 0``.
By [`dfly_monomials`](@ref), ``(Lf)^{(j)} = \sum_t \mathrm{coeff}_t L_{n_t}(f^{(i_t)}
M_t)`` with ``M_t = \prod_m ((1/T')^{(m)})^{p_{t,m}}``, up to a sign on decreasing
branches.  For each operator, changing variables on each branch,
``\|L_n g\|_{L^1} \le \int |g(y)|\,|T'(y)|^{1-n}\,dy``; hence

```math
\|(Lf)^{(j)}\|_{L^1} \le \sum_i \int |f^{(i)}(y)| \sum_{t:\,i_t = i}
|\mathrm{coeff}_t|\,|M_t(y)|\,|1/T'(y)|^{n_t-1}\,dy ,
```

and ``c_{j,i}`` is the supremum over ``y`` of the inner sum, bounded on each box of
[`inverse_derivative_enclosures`](@ref) by evaluating the sum on the box enclosures, and
then maximised over the boxes.  The diagonal ``c_{j,j}`` bounds ``\sup |1/T'|^{j}``, and
``c_{0,0} = 1``.  All arithmetic is interval arithmetic, and the entries are upper
endpoints.

Throws unless every branch is full: on a non-full branch ``Lf`` jumps at the image of an
endpoint and ``(Lf)^{(j)}`` is not a function.
"""
function dfly_coefficient_matrix(D::PwMap, k::Int; boxes::Int = 1024)
    is_full_branch(D) || throw(
        ArgumentError(
            "dfly_coefficient_matrix requires full branches: otherwise Lf jumps and " *
            "(Lf)^(j) is not in L1",
        ),
    )
    H = inverse_derivative_enclosures(D, k; boxes = boxes)
    C = zeros(k + 1, k + 1)
    for j = 0:k
        mons = dfly_monomials(j)
        for i = 0:j
            terms = [t for t in mons if t.i == i]
            isempty(terms) && continue
            C[j+1, i+1] = maximum(
                sup(sum(_monomial_on_box(h, t.n, t.coeff, t.p) for t in terms)) for h in H
            )
        end
    end
    return C
end

###############################################################################
# Weights of the strong norm and the constants (A, B)
###############################################################################

# The smallest weights a_1..a_k with a_k = 1 for which every row i >= 1 satisfies
# sum_{j>=i} a_j c_{j,i} <= A* a_i: a_i = sum_{j>i} a_j c_{j,i} / (A* - c_{i,i}), top down.
# Float64: the weights are a free choice, and (A, B) are recomputed for them in intervals.
function _dfly_weights(C::Matrix{Float64}, k::Int, Astar::Float64)
    a = zeros(k)
    a[k] = 1.0
    for i = k-1:-1:1
        a[i] = sum(a[j] * C[j+1, i+1] for j = i+1:k) / (Astar - C[i+1, i+1])
    end
    return a
end

# The target A*: A_target if given, after checking max_{i>=1} c_ii < A* < 1; otherwise the
# point of a Float64 grid in (max c_ii, 1) at which B(A*)/(1 - A*), in Float64, is least.
function _dfly_target(C::Matrix{Float64}, k::Int, A_target, grid::Int)
    cmax = maximum(C[i+1, i+1] for i = 1:k)
    if A_target === nothing
        Bover(As) = sum(_dfly_weights(C, k, As)[j] * C[j+1, 1] for j = 1:k) / (1 - As)
        cands = [cmax + (1 - cmax) * t / (grid + 1) for t = 1:grid]
        return cands[argmin(Bover.(cands))]
    end
    Astar = Float64(A_target)
    cmax < Astar < 1 ||
        throw(ArgumentError("A_target = $A_target must lie in (max_i c_ii, 1) = ($cmax, 1)"))
    return Astar
end

# Certified (A, B) for the weights a: the coefficient of ||f^(i)||_1 in ||Lf||_a is
# S_i = sum_{j>=i} a_j c_{j,i}; A = max_{i: a_i > 0} S_i / a_i, and a row with a_i = 0 must
# have S_i = 0. B = sum_{j>=1} a_j c_{j,0}.
function _dfly_AB(C::Matrix{Float64}, k::Int, a::Vector{Float64})
    A = 0.0
    for i = 1:k
        S = sum(interval(a[j]) * interval(C[j+1, i+1]) for j = i:k)
        if a[i] > 0
            A = max(A, sup(S / interval(a[i])))
        else
            sup(S) == 0 || error("the weight a_$i is zero but row $i of ||Lf||_a is not")
        end
    end
    B = sup(sum(interval(a[j]) * interval(C[j+1, 1]) for j = 1:k))
    return A, B
end

function _check_expanding(C::Matrix{Float64})
    C[2, 2] < 1 || error(
        "the map is not expanding: the bound on sup|1/T'| is $(C[2, 2]), not below 1",
    )
end

"""
    WeightedSobolev(D::PwMap, k; A_target = nothing, boxes = 1024, grid = 256)

The weighted seminorm whose weights `dfly(W{k,1}, L1, D; A_target)` would choose for `D`,
with the target ``A^*`` it was chosen for recorded as a label.  See
[`dfly`](@ref) for the weights and for the choice of ``A^*`` when `A_target` is omitted.
"""
function WeightedSobolev(
    D::PwMap,
    k::Integer;
    A_target = nothing,
    boxes::Int = 1024,
    grid::Int = 256,
)
    k >= 2 || throw(ArgumentError("WeightedSobolev(D, k) needs k ≥ 2, got $k"))
    C = dfly_coefficient_matrix(D, Int(k); boxes = boxes)
    _check_expanding(C)
    Astar = _dfly_target(C, Int(k), A_target, grid)
    return WeightedSobolev(k, _dfly_weights(C, Int(k), Astar); A_target = Astar)
end

###############################################################################
# dfly(W{k,1}, L1, PwMap)
###############################################################################

@doc raw"""
    dfly(::Type{W{k,1}}, ::Type{L1}, D::PwMap; A_target = nothing, weights = nothing,
         absorb = :weighted, boxes = 1024, grid = 256) -> DFLYCoefficients

Constants ``(A, B)`` for ``k \ge 2`` with

```math
\|Lf\|_a \le A\,\|f\|_a + B\,\|f\|_{L^1}, \qquad
\|f\|_a = \sum_{j=1}^{k} a_j \|f^{(j)}\|_{L^1},\quad a_k = 1,
```

returned as `DFLYCoefficients(WeightedSobolev(k, a; A_target), L1, A, B)`.  The strong norm
is the weighted seminorm [`WeightedSobolev`](@ref), which carries the weights; it is
neither the seminorm ``\|f'\|_{L^1}`` that `dfly(W{1,1}, L1, D)` uses nor the unweighted
sum.

With ``c_{j,i}`` from [`dfly_coefficient_matrix`](@ref),
``\|Lf\|_a \le \sum_{i\ge 1} S_i \|f^{(i)}\|_{L^1} + B\|f\|_{L^1}`` with
``S_i = \sum_{j\ge i} a_j c_{j,i}`` and ``B = \sum_{j\ge 1} a_j c_{j,0}``, so the inequality
holds with ``A = \max_{i \ge 1} S_i/a_i`` (a row with ``a_i = 0`` must have ``S_i = 0``).
For a target ``A^*`` with ``\max_{i\ge1} c_{i,i} < A^* < 1`` the weights are

```math
a_k = 1, \qquad a_i = \frac{\sum_{j>i} a_j c_{j,i}}{A^* - c_{i,i}}, \quad i = k-1,\dots,1,
```

the smallest that give ``S_i \le A^* a_i`` for every ``i``; since every ``c_{j,i} \ge 0``,
each ``a_i`` and hence ``B`` is the least possible for that ``A^*``, and ``B`` decreases as
``A^*`` grows.  The weights are computed in `Float64`, as a free choice; ``A`` and ``B`` are
then recomputed for those exact weights in interval arithmetic, so the returned ``A`` may
exceed ``A^*`` by rounding.

`A_target` sets ``A^*``.  When it is not given, ``A^*`` is the point of a `grid` of
`Float64` values in ``(\max_{i\ge1} c_{i,i}, 1)`` at which ``B(A^*)/(1-A^*)``, evaluated in
`Float64`, is least; the certified constants are then computed at that point.

`weights` fixes the norm instead, either a vector ``(a_1,\dots,a_k)`` with ``a_k = 1`` or a
`WeightedSobolev`, and skips the choice of ``A^*``: ``A`` and ``B`` are certified for those
weights, so that a family of maps can share one norm.  It cannot be combined with
`A_target`.

`absorb = :interpolation` is reserved for a route through an interpolation inequality
between ``\|f'\|_{L^1}``, ``\|f''\|_{L^1}`` and ``\|f\|_{L^1}``; it throws until that lemma
is approved.

Throws for ``l \ne 1``, unless the branches are full (see
[`dfly_coefficient_matrix`](@ref)), unless the bound on ``\sup|1/T'|`` is below 1, if
`A_target` is outside ``(\max_{i\ge1} c_{i,i}, 1)``, and if the certified ``A`` is not
below 1.
"""
function dfly(
    ::Type{W{k,l}},
    ::Type{L1},
    D::PwMap;
    A_target = nothing,
    weights = nothing,
    absorb::Symbol = :weighted,
    boxes::Int = 1024,
    grid::Int = 256,
) where {k,l}
    l == 1 || throw(ArgumentError("dfly(W{k,l}, L1, D) is implemented only for l = 1"))
    k >= 2 || throw(
        ArgumentError("dfly(W{$k,1}, L1, D) is not handled here; W{1,1} is in DFLY.jl"),
    )
    if absorb === :interpolation
        error("dfly(W{$k,1}, L1, D): pending approval of the interpolation lemma")
    elseif absorb !== :weighted
        throw(ArgumentError("absorb must be :weighted or :interpolation, got :$absorb"))
    end
    (weights === nothing || A_target === nothing) ||
        throw(ArgumentError("give either A_target or weights, not both"))

    C = dfly_coefficient_matrix(D, k; boxes = boxes)
    _check_expanding(C)

    if weights === nothing
        Astar = _dfly_target(C, k, A_target, grid)
        norm = WeightedSobolev(k, _dfly_weights(C, k, Astar); A_target = Astar)
    else
        norm = weights isa WeightedSobolev ? weights : WeightedSobolev(k, weights)
        norm.k == k ||
            throw(ArgumentError("the weights are for k = $(norm.k), not for k = $k"))
    end

    A, B = _dfly_AB(C, k, norm.a)
    A < 1 || error("dfly(W{$k,1}, L1, D): the certified A = $A is not below 1")
    return DFLYCoefficients(norm, L1, A, B)
end
