###############################################################################
# Legendre basis on [0,1]
#
# Weak space L²([0,1], dx).  The map x ↦ t = 2x - 1 is a unitary equivalence with
# L²([-1,1], dx/2), the space of the Diffusion paper's Legendre appendix
# (appendix-legendre.tex), so every matrix and constant below is the same in either
# convention: the basis element e_k(x) = (2k+1)^{1/2} P_k(2x-1) on [0,1] is the paper's
# P̃_k(t) = (2k+1)^{1/2} P_k(t) on [-1,1], and ℓ(f) = ½∫_{-1}^1 f dt = ∫_0^1 f dx.
###############################################################################

import .RigorousInvariantMeasures: evaluate, evaluate_integral, integral_covector
import .RigorousInvariantMeasures: one_vector, is_refinement, weak_norm
import .RigorousInvariantMeasures: bernstein_parameter, assemble, change_of_basis

@doc raw"""
    Legendre{T<:AbstractFloat, D, E} <: Basis

The orthonormal Legendre basis of ``L^2([0,1], dx)`` up to degree `N`,
``e_k(x) = (2k+1)^{1/2}P_k(2x-1)``, ``k = 0, …, N``, with ``P_k`` the Legendre
polynomials, ``\|P_k(2\cdot-1)\|_{L^2(dx)}^2 = 1/(2k+1)``.  The projection onto the span
is the ``L^2(dx)``-orthogonal one, ``\Pi_N``, so ``\|\Pi_N\| = 1`` and
``\int\Pi_Nf = \int f``.  A coefficient vector ``b`` stands for ``\sum_k b_k e_k`` and
its ``L^2(dx)`` norm is ``\|b\|_2``.  Under ``t = 2x - 1`` the basis is the orthonormal
basis ``\tilde P_k = (2k+1)^{1/2}P_k`` of ``L^2([-1,1], dx/2)``, so every matrix computed
here is also the matrix in that convention.

`T` is the floating-point type of the certified output (`Float64` or `BigFloat`).

Two constructors:

- `Legendre(N; T = Float64)`: the basis alone, `dynamic = nothing`, `ellipse = nothing`;
- `Legendre(N, D; r, ρ, T = Float64, ngeom = 1024, cover = 64)`: the basis built for the map
  `D`, storing `D` and its [`BernsteinEllipseData`](@ref) certified by
  [`bernstein_ellipse_data`](@ref).

[`galerkin_matrix`](@ref) needs only `D`; [`galerkin_matrix_from_chebyshev`](@ref)
needs the ellipse data.  The basis does not take part in `DiscretizedOperator`, whose
assembly goes through preimages and dual elements.
"""
struct Legendre{T<:AbstractFloat,DT,E} <: Basis
    N::Int
    dynamic::DT
    ellipse::E
    function Legendre{T}(N::Integer, dynamic = nothing, ellipse = nothing) where {T<:AbstractFloat}
        N >= 0 || throw(ArgumentError("the degree must be ≥ 0, got $N"))
        return new{T,typeof(dynamic),typeof(ellipse)}(Int(N), dynamic, ellipse)
    end
end
Legendre(N::Integer; T::Type{<:AbstractFloat} = Float64) = Legendre{T}(N)
function Legendre(N::Integer, D::PwMap; r::Real, ρ::Real, T::Type{<:AbstractFloat} = Float64,
                  ngeom::Integer = 1024, cover::Integer = 64)
    return Legendre{T}(N, D, bernstein_ellipse_data(D; r, ρ, n = ngeom, cover))
end

Base.length(B::Legendre) = B.N + 1
Base.show(io::IO, B::Legendre{T}) where {T} = print(io,
    "Legendre basis on [0,1], orthonormal in L²(dx), degree ≤ $(B.N) ($T)",
    B.dynamic === nothing ? "" : ", built for a dynamic")

weak_norm(::Legendre) = L2
is_refinement(Bf::Legendre, Bc::Legendre) = Bf.N >= Bc.N

# ∫_0^1 e_k dx = δ_{k0}: e_0 = 1, and e_k ⊥ e_0 for k ≥ 1.
evaluate_integral(B::Legendre, i, T = Float64) = i == 1 ? one(T) : zero(T)
integral_covector(B::Legendre; T = Float64) =
    permutedims([i == 1 ? one(T) : zero(T) for i = 1:length(B)])
one_vector(B::Legendre) = [1; zeros(length(B) - 1)]

@doc raw"""
    evaluate(B::Legendre, i, x)

``e_{i-1}(x) = (2i-1)^{1/2}P_{i-1}(2x-1)``, by the three-term recurrence
``(n+1)P_{n+1}(t) = (2n+1)tP_n(t) - nP_{n-1}(t)``.  For an `Interval` argument the
recurrence runs in interval arithmetic and the result encloses the value; interval
arithmetic does not see the cancellation in the recurrence, so the radius grows
geometrically with the degree and a high degree wants `Interval{BigFloat}`.
"""
function evaluate(B::Legendre, i::Integer, x)
    1 <= i <= length(B) || throw(BoundsError(B, i))
    n = i - 1
    return sqrt(_as_number(x, 2n + 1)) * _legendre_P(n, 2 * x - 1)
end
_as_number(::Interval{T}, v) where {T} = interval(T, v)
_as_number(x, v) = oftype(float(x), v)

function _legendre_P(n::Integer, t)
    p0 = one(t)
    n == 0 && return p0
    p1 = t
    for m = 1:n-1
        p0, p1 = p1, ((2m + 1) * t * p1 - m * p0) / (m + 1)
    end
    return p1
end

###############################################################################
# Exact changes of basis between Chebyshev and Legendre
###############################################################################

@doc raw"""
    legendre_to_chebyshev_rational(N) -> Matrix{Rational{BigInt}}

The exact matrix ``L`` with ``P_n = \sum_{k\le n} L_{k+1,n+1}T_k`` for
``0 \le n \le N``: column `n+1` holds the Chebyshev coefficients of ``P_n``.  It is
upper triangular, and ``L_{k+1,n+1} = 0`` unless ``n - k`` is even.

Computed by the Legendre recurrence ``(n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}``, with
``xT_0 = T_1`` and ``xT_k = (T_{k+1} + T_{k-1})/2`` for ``k \ge 1``, in exact
rational arithmetic.  The closed form, with ``\Lambda(z) = \Gamma(z+1/2)/\Gamma(z+1)``,
is ``L_{k+1,n+1} = \frac{2 - \delta_{k0}}{\pi}\Lambda(\frac{n-k}2)\Lambda(\frac{n+k}2)``
for ``n - k`` even (Alpert and Rokhlin 1991, SIAM J. Sci. Stat. Comput. 12, 158–179);
the test suite compares the two.  The same matrix relates the shifted families
``P_k(2x-1)`` and ``T_k(2x-1)`` on ``[0,1]``.  An internal helper: the basis-level
interface is [`legendre_basis`](@ref) and [`change_of_basis`](@ref).
"""
function legendre_to_chebyshev_rational(N::Integer)
    L = zeros(Rational{BigInt}, N + 1, N + 1)
    L[1, 1] = 1
    N >= 1 && (L[2, 2] = 1)
    for n = 1:N-1
        xP = _cheb_times_x(view(L, :, n + 1))
        L[:, n+2] = ((2n + 1) * xP - n * L[:, n]) / (n + 1)
    end
    return L
end

# Chebyshev coefficients of x·p, for p given by Chebyshev coefficients of degree < length - 1.
function _cheb_times_x(c::AbstractVector{Q}) where {Q}
    m = length(c)
    iszero(c[m]) || throw(ArgumentError("x·p would exceed the degree of the array"))
    out = zeros(Q, m)
    out[2] += c[1]
    for k = 1:m-2
        out[k+2] += c[k+1] / 2
        out[k] += c[k+1] / 2
    end
    return out
end

@doc raw"""
    chebyshev_to_legendre_rational(N) -> Matrix{Rational{BigInt}}

The exact matrix ``C`` with ``T_n = \sum_{k\le n} C_{k+1,n+1}P_k`` for
``0 \le n \le N``: column `n+1` holds the Legendre coefficients of ``T_n``.  It is
upper triangular and the inverse of [`legendre_to_chebyshev_rational`](@ref).

Computed by ``T_{n+1} = 2xT_n - T_{n-1}`` with
``xP_k = ((k+1)P_{k+1} + kP_{k-1})/(2k+1)``, in exact rational arithmetic.  The
closed form (Alpert and Rokhlin 1991) is ``C_{n+1,n+1} = \sqrt\pi/(2\Lambda(n))`` on the
diagonal for ``n \ge 1`` and, for ``k < n`` with ``n - k`` even,
``C_{k+1,n+1} = -n(k+\frac12)\Lambda(\frac{n-k-2}2)\Lambda(\frac{n+k-1}2)/((n+k+1)(n-k))``;
the test suite compares the two.  An internal helper, as the one above.
"""
function chebyshev_to_legendre_rational(N::Integer)
    C = zeros(Rational{BigInt}, N + 1, N + 1)
    C[1, 1] = 1
    N >= 1 && (C[2, 2] = 1)
    for n = 1:N-1
        xT = _legendre_times_x(view(C, :, n + 1))
        C[:, n+2] = 2 * xT - C[:, n]
    end
    return C
end

# Legendre coefficients of x·p, for p given by Legendre coefficients of degree < length - 1.
function _legendre_times_x(c::AbstractVector{Q}) where {Q}
    m = length(c)
    iszero(c[m]) || throw(ArgumentError("x·p would exceed the degree of the array"))
    out = zeros(Q, m)
    for k = 0:m-2
        out[k+2] += c[k+1] * (k + 1) // (2k + 1)
        k >= 1 && (out[k] += c[k+1] * k // (2k + 1))
    end
    return out
end

# A rational number enclosed in an interval of type T (BigFloat: current precision).
# Julia's conversion of a Rational with a RoundingMode is correctly rounded.
_rational_interval(::Type{T}, q::Rational) where {T<:AbstractFloat} =
    interval(T, T(q, RoundDown), T(q, RoundUp))

# (2k+1)^{1/2}, k = 0..N, enclosed.
_orthonormal_scale(::Type{T}, N::Integer) where {T} = [sqrt(interval(T, 2k + 1)) for k = 0:N]

@doc raw"""
    legendre_basis(Bc::Chebyshev; T = Float64) -> (Bl, cheb_to_leg, leg_to_cheb)

The Legendre basis of the same degree as `Bc`, for the same map and with the same
[`BernsteinEllipseData`](@ref), together with the matrices of
[`change_of_basis_matrices`](@ref).  `Bc` must have been built for a map by
`Chebyshev(n, D; r, ρ)`; the ellipse data are reused, not recomputed.
"""
function legendre_basis(Bc::Chebyshev; T::Type{<:AbstractFloat} = Float64)
    Bc.dynamic === nothing &&
        throw(ArgumentError("the Chebyshev basis was not built for a dynamic: use Chebyshev(n, D; r, ρ)"))
    Bl = Legendre{T}(length(Bc) - 1, Bc.dynamic, Bc.ellipse)
    return (Bl, change_of_basis_matrices(Bc, Bl)...)
end

@doc raw"""
    change_of_basis_matrices(Bc::Chebyshev, Bl::Legendre) -> (cheb_to_leg, leg_to_cheb)

`cheb_to_leg * c` encloses the orthonormal Legendre coefficients ``b`` of
``f = \sum_m c_mT_m(2x-1)``, and `leg_to_cheb * b` the Chebyshev coefficients of
``\sum_k b_ke_k``.  Both are `BallMatrix{T}`, upper triangular, with the entries of
[`chebyshev_to_legendre_rational`](@ref) and [`legendre_to_chebyshev_rational`](@ref)
scaled by ``(2k+1)^{\mp1/2}`` and enclosed by outward rounding.

Throws an `ArgumentError` unless both bases were built for a dynamic, for the same one
(the same object), with the same ellipse data, and have the same degree.

For an operator, `cheb_to_leg * Q * leg_to_cheb` is the matrix in the Legendre basis of
the operator that `Q` represents in the Chebyshev basis.  For the interpolatory Chebyshev
matrix that operator is ``\Pi^{\mathrm{Cheb}}_NL``, not the Galerkin ``\Pi_NL``; for the
latter see [`galerkin_matrix_from_chebyshev`](@ref).
"""
function change_of_basis_matrices(Bc::Chebyshev, Bl::Legendre{T}) where {T}
    (Bc.dynamic === nothing || Bl.dynamic === nothing) &&
        throw(ArgumentError("both bases must be built for a dynamic"))
    Bc.dynamic === Bl.dynamic ||
        throw(ArgumentError("the two bases were built for different dynamics"))
    Bc.ellipse == Bl.ellipse ||
        throw(ArgumentError("the two bases carry different ellipse data"))
    length(Bl) == length(Bc) ||
        throw(ArgumentError("the bases have degrees $(Bl.N) and $(length(Bc) - 1)"))
    N = Bl.N
    s = _orthonormal_scale(T, N)
    C = _rational_interval.(T, chebyshev_to_legendre_rational(N))
    L = _rational_interval.(T, legendre_to_chebyshev_rational(N))
    # b_j = Σ_m C[j,m] c_m / s_j;  c_m = Σ_k L[m,k] s_k b_k
    return (BallMatrix(C ./ s), BallMatrix(L .* permutedims(s)))
end

"""
    change_of_basis(Bc::Chebyshev, Bl::Legendre, v)
    change_of_basis(Bl::Legendre, Bc::Chebyshev, v)

The coefficients `v` in the first basis, converted to the second, with the matrices of
[`change_of_basis_matrices`](@ref), which throws unless the two bases share the dynamic;
the result is a `BallVector` that encloses them.
"""
change_of_basis(Bc::Chebyshev, Bl::Legendre, v::AbstractVector) =
    change_of_basis_matrices(Bc, Bl)[1] * v
change_of_basis(Bl::Legendre, Bc::Chebyshev, v::AbstractVector) =
    change_of_basis_matrices(Bc, Bl)[2] * v

###############################################################################
# The Galerkin matrix, route 1: Koopman duality and Taylor expansions with
# Cauchy-bounded remainders
###############################################################################

@doc raw"""
    galerkin_matrix(B::Legendre, D::PwMap; nsub = B.N + 8, order = max(48, 3B.N ÷ 4 + 16), κ = 4,
                    precision = max(256, 3B.N + 128)) -> BallMatrix

Certified enclosure of the matrix ``A_{j+1,k+1} = \langle Le_k, e_j\rangle_{L^2(dx)}``,
``0 \le j, k \le N``, of ``\Pi_NL\Pi_N`` in the orthonormal basis, with ``L`` the transfer
operator of `D` on ``L^2(dx)``.  By duality
``\langle Le_k, e_j\rangle = \int_0^1 e_k(x)\,e_j(T(x))\,dx
= \sum_b\int_{I_b}e_k(x)\,e_j(T_b(x))\,dx``,
the sum over the branches of `D`, so no inverse branch is needed.

Each branch domain ``I_b`` is split into `nsub` subintervals at the points
``a + (b-a)(1 - \cos(\pi i/\mathrm{nsub}))/2``, finer near the ends, where ``e_j\circ T``
grows fastest off the real line.  On a subinterval with centre ``c`` and half-length
``\le H``, in the variable ``s = (x - c)/H \in [-1, 1]``:

1. the Taylor coefficients of order ``\le d`` (`order`) of ``p_k(s) = e_k(c + Hs)`` and
   ``q_j(s) = e_j(T_b(c + Hs))`` are computed in `Interval{BigFloat}` at `precision` bits
   (the recurrence of `evaluate`, on Taylor polynomials; ``T_b`` evaluated on a
   `Taylor1`), and ``\int p_k^{\le d}q_j^{\le d}`` is integrated exactly against the moments
   ``\int s^m``;
2. the remainders are bounded by Cauchy's estimate on the disc ``|s| \le κ``: for ``g``
   analytic there with ``|g| \le M``, the Taylor coefficients satisfy ``|g_i| \le Mκ^{-i}``,
   so ``\sup_{|s|\le1}|g - g^{\le d}| \le Mκ^{-(d+1)}/(1 - 1/κ)``.  ``M`` comes from
   ``|P_n(u)| \le ρ(u)^n``, ``ρ(u)`` the Bernstein parameter of ``u``
   ([`bernstein_parameter`](@ref)), evaluated on a box containing the disc and on the
   image of that box under ``2T_b - 1``.  The bound ``|P_n| \le ρ^n`` on ``E_ρ``: the
   Chebyshev coefficients of ``P_n`` are nonnegative and sum to ``P_n(1) = 1``
   (Lemma lem:legendre-cheb of appendix-legendre.tex), and
   ``|T_m| \le (ρ^m + ρ^{-m})/2 \le ρ^m`` on ``E_ρ``;
3. ``|\int(PQ - pq)| \le |I|\,(R_k\sup|Q_j| + \sup|p_k|\,S_j)``, with ``R_k``, ``S_j``
   the two remainder bounds and ``|I|`` the length of the subinterval.

Where an endpoint of ``I_b`` is an interval of nonzero width (a branch point of a
`mod1_dynamic`), the integral over the uncertain part, of length at most that width, is
bounded by the width times ``\sup|e_k||e_j\circ T_b|`` over it.  The branches must accept
`Taylor1{Interval{BigFloat}}` and `Complex{Interval{BigFloat}}` arguments.

The subintervals run in parallel on `Threads.nthreads()` threads.  The output is a
`BallMatrix{T}` for `B::Legendre{T}`, rounded outward from the `BigFloat` enclosure.
Row 1 encloses ``(1, 0, …, 0)``, since ``\int Lf = \int f``.
"""
function galerkin_matrix(B::Legendre{T}, D::PwMap; nsub::Integer = B.N + 8,
                         order::Integer = max(48, 3B.N ÷ 4 + 16), κ::Real = 4,
                         precision::Integer = max(256, 3B.N + 128)) where {T}
    N, d = B.N, order
    κ > 1 || throw(ArgumentError("κ must exceed 1"))
    d >= 1 || throw(ArgumentError("order must be ≥ 1"))
    A = setprecision(BigFloat, precision) do
        TI = Interval{BigFloat}
        subs = Tuple{Any,BigFloat,BigFloat}[]
        slivers = Tuple{Any,TI}[]
        for br in branches(D)
            X1 = interval(BigFloat, inf(br.X[1]), sup(br.X[1]))
            X2 = interval(BigFloat, inf(br.X[2]), sup(br.X[2]))
            lo, hi = sup(X1), inf(X2)
            lo < hi || error("a branch domain is narrower than its endpoint enclosures")
            pts = [i == 0 ? lo : i == nsub ? hi :
                   lo + (hi - lo) * (1 - cospi(BigFloat(i) / nsub)) / 2 for i = 0:nsub]
            for i = 1:nsub
                push!(subs, (br, pts[i], pts[i+1]))
            end
            for X in (X1, X2)
                sup(X) > inf(X) && push!(slivers, (br, X))
            end
        end
        nt = max(1, min(Threads.nthreads(), length(subs)))
        chunks = [subs[i:nt:end] for i = 1:nt]
        tasks = map(chunks) do chunk
            Threads.@spawn setprecision(BigFloat, precision) do
                acc = fill(interval(BigFloat, 0), N + 1, N + 1)
                for (br, a, b) in chunk
                    acc .+= _legendre_subinterval_block(br, a, b, N, d, BigFloat(κ))
                end
                acc
            end
        end
        total = sum(fetch.(tasks))
        for (br, X) in slivers
            total .+= _legendre_sliver_block(br, X, N)
        end
        total
    end
    return _to_ballmatrix(T, A)
end

_to_ballmatrix(::Type{Float64}, A::AbstractMatrix{<:Interval}) = BallMatrix(
    map(x -> interval(Float64, Float64(inf(x), RoundDown), Float64(sup(x), RoundUp)), A))
_to_ballmatrix(::Type{BigFloat}, A::AbstractMatrix{<:Interval{BigFloat}}) = BallMatrix(A)

# Upper bounds of (2n+1)^{1/2} |P_n| over a set enclosed by u, n = 0..N.
function _legendre_sup_bounds(u, N::Integer)
    ρ = interval(sup(bernstein_parameter(u)))
    return [sup(sqrt(interval(BigFloat, 2n + 1)) * ρ^n) for n = 0:N]
end

function _legendre_subinterval_block(br, a::BigFloat, b::BigFloat, N::Int, d::Int, κ::BigFloat)
    TI = Interval{BigFloat}
    ai, bi = interval(a), interval(b)
    c = (a + b) / 2                                   # any point of [a,b]; exact as stored
    ci = interval(c)
    H = max(sup(bi - ci), sup(ci - ai))               # H ≥ max(b - c, c - a)
    Hi = interval(H)
    s = _orthonormal_scale(BigFloat, N)

    # p_k(s) = e_k(c + H s): u = 2x - 1 = (2c - 1) + 2H s is linear, so the recurrence
    # multiplies by a linear polynomial and the coefficients of order ≤ d are exact.
    u0, u1 = 2ci - 1, 2Hi
    Pt = Matrix{TI}(undef, N + 1, d + 1)
    pm = zeros(TI, d + 1); pm[1] = interval(BigFloat, 1)
    pc = zeros(TI, d + 1); pc[1] = u0; pc[2] = u1
    Pt[1, :] = pm
    N >= 1 && (Pt[2, :] = pc)
    for n = 1:N-1
        up = similar(pc)
        up[1] = u0 * pc[1]
        for i = 2:d+1
            up[i] = u0 * pc[i] + u1 * pc[i-1]
        end
        pn = ((2n + 1) * up - n * pm) / (n + 1)
        pm, pc = pc, pn
        Pt[n+2, :] = pc
    end
    Pt .*= s

    # q_j(s) = e_j(T_b(c + H s)): T_b on a Taylor polynomial of order d.
    xs = Taylor1([ci; Hi; zeros(TI, d - 1)], d)
    v = 2 * br.f(xs) - 1
    Qt = Matrix{TI}(undef, N + 1, d + 1)
    qm = Taylor1([interval(BigFloat, 1); zeros(TI, d)], d)
    qc = v
    Qt[1, :] = qm.coeffs
    N >= 1 && (Qt[2, :] = qc.coeffs)
    for n = 1:N-1
        qn = ((2n + 1) * v * qc - n * qm) / (n + 1)
        qm, qc = qc, qn
        Qt[n+2, :] = qc.coeffs
    end
    Qt .*= s

    # moments μ_m = ∫_a^b ((x - c)/H)^m dx = H (β^{m+1} - α^{m+1})/(m+1)
    α, β = (ai - ci) / Hi, (bi - ci) / Hi
    μ = [Hi * (β^(m + 1) - α^(m + 1)) / (m + 1) for m = 0:2d]
    Mm = [μ[i+l-1] for i = 1:d+1, l = 1:d+1]
    block = Qt * (Mm * permutedims(Pt))               # rows j, columns k

    # Cauchy bounds on the disc |x - c| ≤ κH, enclosed by a box.
    rad = interval(κ) * Hi
    Z = complex(ci + rad * interval(BigFloat, -1, 1), rad * interval(BigFloat, -1, 1))
    MP = _legendre_sup_bounds(2Z - 1, N)
    MQ = _legendre_sup_bounds(2 * br.f(Z) - 1, N)
    tailf = interval(κ)^(-(d + 1)) / (1 - 1 / interval(κ))
    R = [interval(m) * tailf for m in MP]
    S = [interval(m) * tailf for m in MQ]
    pn = [sum(abs.(Pt[k, :])) for k = 1:N+1]
    len = bi - ai
    for k = 1:N+1, j = 1:N+1
        e = sup(len * (R[k] * interval(MQ[j]) + pn[k] * S[j]))
        block[j, k] += interval(BigFloat, -e, e)
    end
    return block
end

# Contribution of the uncertain end of a branch domain: an integral over a part of X.
function _legendre_sliver_block(br, X::Interval{BigFloat}, N::Int)
    zero_im = interval(BigFloat, 0)
    MP = _legendre_sup_bounds(complex(2X - 1, zero_im), N)
    MQ = _legendre_sup_bounds(complex(2 * br.f(X) - 1, zero_im), N)
    w = interval(sup(X)) - interval(inf(X))
    return [(e = sup(w * interval(MP[k]) * interval(MQ[j])); interval(BigFloat, -e, e))
            for j = 1:N+1, k = 1:N+1]
end

###############################################################################
# The Galerkin matrix, route 2: RIM's Chebyshev interpolation, aliasing and tail
# from the ellipse data, then the exact change of basis
###############################################################################

@doc raw"""
    galerkin_matrix_from_chebyshev(B::Legendre, D::PwMap; Nprime) -> BallMatrix{Float64}

Certified enclosure of the same matrix as [`galerkin_matrix`](@ref), from the Chebyshev
interpolation coefficients of ``LT_m``, ``m \le N'``, that `assemble` returns for
`Chebyshev(Nprime, Eρ(r))` (an `FFTW` extension method: `using FFTW` first).  `Nprime`
must be a power of two with ``N' \ge N``, and `B` must have been built for `D` by
`Legendre(N, D; r, ρ)`, whose [`BernsteinEllipseData`](@ref) `(r, R, W)` give:

* ``\sup_{E_R}|Le_k| \le W\sup_{E_r}|e_k| \le W(2k+1)^{1/2}r^k =: M_k``
  (``|P_k| \le r^k`` on ``E_r``, as in [`galerkin_matrix`](@ref));
* the Chebyshev coefficients of ``g = Le_k`` satisfy ``|c_0| \le M_k``,
  ``|c_m| \le 2M_kR^{-m}`` (Trefethen, *ATAP*, Theorem 8.1);
* the interpolant in the ``N'+1`` Chebyshev points has coefficients
  ``\tilde c_m = c_m + \sum_{l\ge1}(c_{2lN'-m} + c_{2lN'+m})`` for ``0 < m < N'``,
  ``\tilde c_0 = c_0 + \sum_{l\ge1}c_{2lN'}`` and
  ``\tilde c_{N'} = c_{N'} + \sum_{l\ge1}c_{(2l+1)N'}`` (*ATAP*, Theorem 4.2), so for
  ``0 \le m \le N'``, ``|\tilde c_m - c_m| \le 2M_k(R^m + R^{-m})R^{-2N'}/(1 - R^{-2N'})``;
* the tail ``g - \sum_{m\le N'}c_mT_m`` has ``L^2(dx)`` norm at most its sup norm, at most
  ``\sum_{m>N'}|c_m| \le 2M_kR^{-(N'+1)}/(1 - R^{-1})``, which bounds its pairing with
  each ``e_j``.

Then ``\langle Le_k, e_j\rangle = \sum_{m\le N'}C_{j+1,m+1}c_m/(2j+1)^{1/2}
+ \langle\mathrm{tail}, e_j\rangle``, with ``C`` of
[`chebyshev_to_legendre_rational`](@ref).  The ``T_m`` with ``N < m \le N'`` are kept:
they have Legendre components of degree ``\le N``, so truncating the Chebyshev matrix at
``N`` would give ``\Pi^{\mathrm{Cheb}}_NL`` and not ``\Pi_NL``.
"""
function galerkin_matrix_from_chebyshev(B::Legendre{Float64,<:PwMap,BernsteinEllipseData},
                                        D::PwMap; Nprime::Integer)
    B.dynamic === D || throw(ArgumentError("the Legendre basis was not built for this dynamic"))
    N = B.N
    ispow2(Nprime) && Nprime >= N ||
        throw(ArgumentError("Nprime must be a power of two ≥ N, got $Nprime"))
    E = B.ellipse
    Bc = Chebyshev(Nprime, Eρ(E.r))
    Mt = assemble(Bc, D)                              # (N'+1)², interpolation coefficients
    size(Mt) == (Nprime + 1, Nprime + 1) || error("unexpected Chebyshev matrix size")
    s = _orthonormal_scale(Float64, N)
    Lc = _rational_interval.(Float64, legendre_to_chebyshev_rational(Nprime))[:, 1:N+1]
    C = _rational_interval.(Float64, chebyshev_to_legendre_rational(Nprime))[1:N+1, :]
    Ct = Matrix(Mt) * (Lc .* permutedims(s))          # interpolation coefficients of L e_k
    r, R, W = interval(E.r), interval(E.R), interval(E.W)
    Mk = [W * s[k+1] * r^k for k = 0:N]
    q = R^(-2Nprime) / (1 - R^(-2Nprime))
    for k = 0:N, m = 0:Nprime
        e = sup(2 * Mk[k+1] * (R^m + R^(-m)) * q)
        Ct[m+1, k+1] += interval(-e, e)
    end
    A = (C * Ct) ./ s
    for k = 0:N
        e = sup(2 * Mk[k+1] * R^(-(Nprime + 1)) / (1 - 1 / R))
        A[:, k+1] .+= interval(-e, e)
    end
    return BallMatrix(A)
end

###############################################################################
# Resolvent and norms in L²(dx)
###############################################################################

@doc raw"""
    l2_resolvent_bound_average_zero(A::BallMatrix, z = 1) -> (; bound, Rbar)

For `A` the Galerkin matrix of [`galerkin_matrix`](@ref), an upper bound `bound` on
``\|(zI - A_0)^{-1}\|_2``, where ``A_0 = A[2:\mathrm{end}, 2:\mathrm{end}]`` is the matrix of
``\Pi_NL\Pi_N`` on the mean-zero subspace ``V_0 = \mathrm{span}\{e_k : k \ge 1\}``, and
``\bar R = \max(1, \mathrm{bound})``.  The norm is that of ``L^2(dx)`` on ``V_0``.

``V_0`` is invariant: ``\langle Le_k, e_0\rangle = \int e_k = 0`` for ``k \ge 1``, which is
an identity, not a property of the computed enclosure.  The bound is
`svd_bound_L2_opnorm_inverse` of BallArithmetic on the ball ``zI - A_0``, valid for every
matrix in the ball; it is `Inf` when the certified smallest singular value is not positive.
"""
function l2_resolvent_bound_average_zero(A::BallMatrix{T}, z::Real = 1) where {T}
    n = size(A, 1) - 1
    A0 = BallMatrix(A.c[2:end, 2:end], A.r[2:end, 2:end])
    M = BallMatrix(Matrix{T}(T(z) * I, n, n)) - A0
    bound = svd_bound_L2_opnorm_inverse(M)
    return (; bound, Rbar = max(one(bound), bound))
end

@doc raw"""
    l2_transfer_norm_bound(D::PwMap; tol = 1e-3) -> Float64

An upper bound on ``\|L\|_{L^2(dx)\to L^2(dx)}`` for the transfer operator of `D`:
``\|L\|^2 \le \sup_y\sum_b|τ_b'(y)| \le \sum_b\sup_{I_b}1/|T'|``, since by Cauchy–Schwarz
``|Lf|^2 \le (\sum_b|τ_b'|)(\sum_b|f\circ τ_b|^2|τ_b'|)`` and the second factor integrates to
``\|f\|^2``.  Each supremum is `branch_inverse_derivative` (global optimisation in
interval arithmetic) and the result is ``(\sum_b\cdot)^{1/2}`` rounded up.
"""
function l2_transfer_norm_bound(D::PwMap; tol = 1e-3)
    W = 0.0
    for br in branches(D)
        W = W ⊕₊ Float64(sup(branch_inverse_derivative(br, tol)), RoundUp)
    end
    return sqrt_round(W, RoundUp)
end

###############################################################################
# Constants of the Legendre appendix (appendix-legendre.tex), certified
###############################################################################


@doc raw"""
    legendre_C_M(ρ, M; T = Float64) -> Interval

``C_M = (\sum_{j=0}^M(2j+1)ρ^{2j})^{1/2}``, the constant of Lemma lem:legendre-cheb:
``\|p\|_ρ \le C_M\|p\|`` for polynomials of degree ``\le M``.  A real argument is taken as
the exact value of the number passed.
"""
function legendre_C_M(ρ, M::Integer; T = Float64)
    ρi = as_interval(T, ρ)
    return sqrt(sum(interval(T, 2j + 1) * ρi^(2j) for j = 0:M))
end

@doc raw"""
    legendre_tau_M(ρ, R, M; T = Float64) -> Interval

``τ_M`` of Lemma lem:legendre-tail, ``\|(I - \Pi_M)g\|_ρ \le τ_M\|g\|_R`` for
``1 < ρ < R``, in the closed form
``τ_M^2 = \frac{q^{2(M+1)}}{1-q^2} + 2C_M\frac{(ρ/R^2)^{M+1}}{1-ρ/R^2}
+ C_M^2\frac{R^{-2(M+1)}}{1-R^{-2}}``, ``q = ρ/R``.
"""
function legendre_tau_M(ρ, R, M::Integer; T = Float64)
    ρi, Ri = as_interval(T, ρ), as_interval(T, R)
    1 < inf(ρi) && sup(ρi) < inf(Ri) || throw(ArgumentError("need 1 < ρ < R"))
    C = legendre_C_M(ρi, M; T)
    q = ρi / Ri
    p = ρi / Ri^2
    t = q^(2(M + 1)) / (1 - q^2) + 2C * p^(M + 1) / (1 - p) +
        C^2 * Ri^(-2(M + 1)) / (1 - Ri^(-2))
    return sqrt(t)
end

@doc raw"""
    legendre_l2_tail(ρ, N; T = Float64) -> Interval

``ρ^{-(N+1)}(1 - ρ^{-2})^{-1/2}``: by Lemma lem:legendre-tail,
``\|(I - \Pi_N)g\|_{L^2} \le`` this times ``\|g\|_ρ`` for ``ρ > 1``.
"""
function legendre_l2_tail(ρ, N::Integer; T = Float64)
    ρi = as_interval(T, ρ)
    inf(ρi) > 1 || throw(ArgumentError("need ρ > 1"))
    return ρi^(-(N + 1)) / sqrt(1 - ρi^(-2))
end

@doc raw"""
    legendre_l2_tail_geometric(K, θ, N; T = Float64) -> Interval

If the orthonormal Legendre coefficients of ``f`` satisfy ``|b_k| \le Kθ^k`` for
``k > N``, with ``0 \le θ < 1``, then
``\|(I - \Pi_N)f\|_{L^2} = (\sum_{k>N}|b_k|^2)^{1/2} \le Kθ^{N+1}(1 - θ^2)^{-1/2}``.
"""
function legendre_l2_tail_geometric(K, θ, N::Integer; T = Float64)
    θi = as_interval(T, θ)
    0 <= inf(θi) && sup(θi) < 1 || throw(ArgumentError("need 0 ≤ θ < 1"))
    return as_interval(T, K) * θi^(N + 1) / sqrt(1 - θi^2)
end

@doc raw"""
    legendre_dfly(ρ, R, r, W, M; T = Float64) -> NamedTuple
    legendre_dfly(E::BernsteinEllipseData, M; T = Float64)

The constants of Corollary cor:legendre-LY, ``\|Lf\|_ρ \le A_M\|f\|_ρ + B_M\|f\|``, as
intervals: `C_M` ([`legendre_C_M`](@ref)), `τ_M` ([`legendre_tau_M`](@ref)),
``C_2 = 2W/(1 - (r/ρ)^2)^{1/2}``, ``A_M = τ_MC_2``, ``B_M = C_MW^{1/2}`` and the bound
``M_w = W^{1/2}`` on ``\|L\|_{L^2}``.  Needs ``1 < r < ρ < R``.  Numbers are taken as
given; certified inputs are the fields of a [`BernsteinEllipseData`](@ref).
"""
function legendre_dfly(ρ, R, r, W, M::Integer; T = Float64)
    ρi, Ri, ri, Wi = as_interval(T, ρ), as_interval(T, R), as_interval(T, r), as_interval(T, W)
    1 < inf(ri) && sup(ri) < inf(ρi) || throw(ArgumentError("need 1 < r < ρ"))
    C_M = legendre_C_M(ρi, M; T)
    τ_M = legendre_tau_M(ρi, Ri, M; T)
    C₂ = 2Wi / sqrt(1 - (ri / ρi)^2)
    return (; C_M, τ_M, C₂, A_M = τ_M * C₂, B_M = C_M * sqrt(Wi), M_w = sqrt(Wi))
end
legendre_dfly(E::BernsteinEllipseData, M::Integer; T = Float64) =
    legendre_dfly(E.ρ, E.R, E.r, E.W, M; T)

@doc raw"""
    legendre_dfly_smallest_M(ρ, R, r, W; Mmax = 2000, T = Float64) -> (M, constants)
    legendre_dfly_smallest_M(E::BernsteinEllipseData; Mmax, T)

The smallest ``M \le`` `Mmax` with ``\sup A_M < 1`` in [`legendre_dfly`](@ref), and the
constants at that ``M``; throws if there is none.
"""
function legendre_dfly_smallest_M(ρ, R, r, W; Mmax::Integer = 2000, T = Float64)
    for M = 0:Mmax
        c = legendre_dfly(ρ, R, r, W, M; T)
        sup(c.A_M) < 1 && return (M, c)
    end
    throw(ArgumentError("A_M ≥ 1 for every M ≤ $Mmax"))
end
legendre_dfly_smallest_M(E::BernsteinEllipseData; Mmax::Integer = 2000, T = Float64) =
    legendre_dfly_smallest_M(E.ρ, E.R, E.r, E.W; Mmax, T)
