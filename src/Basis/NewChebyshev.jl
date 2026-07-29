@doc raw"""
    Chebyshev{S<:NormKind, WK<:NormKind, T<:AbstractVector} <: Basis

Chebyshev basis ``φ_i(x) = T_{i-1}(2x-1)`` on ``[0,1]``, carrying its strong and
weak norms as type parameters in the same style as
[`FourierAnalytic`](@ref).

Two weak norms are available:

- `C1` (the default, and what `Chebyshev(n, k)` builds) — the setting of
  Nisoli & Taylor-Crush, *Rigorous Computation of Linear Response for
  Intermittent Maps*, J. Stat. Phys. 190 (2023) 192, whose Theorems 3.13/3.14
  supply the projection errors used below.
- `L2`, selected by `Chebyshev(n, k, L2)`. **The `L2` here is
  ``L^2(μ)`` for the arcsine measure ``dμ = dx/(π\sqrt{x(1-x)})``**, the one that
  makes the ``T_m`` orthogonal, so Parseval holds on the coefficients and the
  norm interface collapses to the Fourier-like one-liners. Use
  [`l2_measure_conversion_bounds`](@ref) and
  [`gram_restrict_to_average_zero`](@ref) to move between this and
  ``L^2(dx)``.

The strong norm is ``W^{k,1}`` in both cases, so `dfly(W{k,1}, L1, D)` applies
unchanged — the Lasota–Yorke inequality is a statement about function spaces
and does not see the basis.
"""
struct Chebyshev{S<:NormKind,WK<:NormKind,T<:AbstractVector} <: Basis
    p::T
    k::Integer
    strong::S
    weak::WK
    # TODO: check in constructor that p is sorted, starts with 0 and ends with 1
end

ChebCouples(n, T) = hcat(
    [
        Interval{T}(pi)
        (reverse([Interval{T}(0.0); [j * Interval{T}(pi) / n for j = 1:n-1]]))
    ],
    [
        Interval{T}(0.0)
        (reverse([Interval{T}(1.0); [cos(j * Interval{T}(pi) / n) for j = 1:n-1]]) .+ 1) / 2
    ],
)

ChebPoints(n, T) = ChebCouples(n, T)[:, 2]

# Strong W^{k,1}; weak L². As for FourierAnalytic the weak norm is L², but for
# `L2` here is L²(dx), the ambient measure of the rest of the package; `C1` is
# the legacy Nisoli-Taylor-Crush norm.
function Chebyshev(n::Integer, k::Integer, ::Type{WK} = L2; T = Float64) where {WK<:NormKind}
    return Chebyshev(ChebPoints(n, T), k, W{k,1}(), WK())
end
Chebyshev(p::AbstractVector, k::Integer, ::Type{WK} = L2) where {WK<:NormKind} =
    Chebyshev(p, k, W{k,1}(), WK())

# Analytic (Bernstein-ellipse) strong norm, weak L2 — the Chebyshev counterpart
# of `FourierAnalytic(k, n; η)`. `k` is unused for this strong norm.
function Chebyshev(n::Integer, strong::Eρ; T = Float64)
    return Chebyshev(ChebPoints(n, T), 0, strong, L2μ())
end
Chebyshev(p::AbstractVector, strong::Eρ) = Chebyshev(p, 0, strong, L2μ())
Base.show(io::IO, B::Chebyshev) =
    print(io, "Chebyshev basis on $(length(B)) points, highest degree $(length(B)-1)")

"""
Return the size of the Chebyshev basis
"""
Base.length(B::Chebyshev) = length(B.p)

##########################################################
# These are the methods to rigorously enclose the image of a Chebyshev polynomial
# they need refactoring, but not now
#########################################################
"""
	_eval_T

Eval the Chebyshev polynomial up to degree n on an array of 
points in [-1, 1].

Not satisfactory, the intervals explode

"""
function _eval_T(n, x::Array{T}) where {T}
    k = length(x)
    M = zeros(T, k, n + 1)
    M[:, 1] = ones(k)
    M[:, 2] = x
    for i = 3:n+1
        M[:, i] = (2 * x) .* M[:, i-1] - M[:, i-2]
    end
    return M
end

""" 
	eval_Clenshaw_BackwardFirst

Eval a polynomial in Chebyshev basis, ClenshawBackward, using ball arithmetic
Following Viviane Ledoux, Guillaume Moroz 
"Evaluation of Chebyshev polynomials on intervals andapplication to root finding"
"""
function eval_Clenshaw_BackwardFirst(coeff::Vector{Interval{S}}, x::Interval{T}) where {S,T}
    coeff_a = mid.(coeff)
    coeff_r = radius.(coeff)
    a, r = midradius(x)
    m = length(coeff)
    u = zeros(Interval{T}, m + 1)
    ϵ = zeros(Interval{T}, m + 1)
    u[m] = coeff_a[m]
    for k in reverse(2:m-1)
        u_temp = 2 * a * u[k+1] - u[k+2] + Interval{T}(coeff_a[k])
        u[k], ϵ[k] = midradius(u_temp)
    end
    u_temp = a * u[2] - u[3] + Interval{T}(coeff_a[1])
    u[1], ϵ[1] = midradius(u_temp)

    e = zeros(Interval{T}, m + 1)
    e[m] = coeff_r[m]
    for k in reverse(2:m-1)
        e[k] = e[k+1] + 2 * r * abs(u[k+1]) + ϵ[k] + coeff_r[k]
    end
    e[1] = e[2] + r * abs(u[1]) + ϵ[1] + coeff_r[1]
    γ = sup(e[1])
    return u[1] + interval(-γ, γ)
end
eval_Clenshaw_BackwardFirst(coeff::Vector{Float64}, x::Interval) =
    eval_Clenshaw_BackwardFirst(Interval.(coeff), x)

function eval_Clenshaw_BackwardSecond(
    coeff::Vector{Interval{S}},
    x::Interval{T},
) where {S,T}
    coeff_a = mid.(coeff)
    coeff_r = radius.(coeff)
    a, r = midradius(x)
    m = length(coeff)
    u = zeros(Interval{T}, m + 1)
    ϵ = zeros(Interval{T}, m + 1)
    u[m] = coeff_a[m]
    for k in reverse(2:m-1)
        u_temp = 2 * a * u[k+1] - u[k+2] + Interval{T}(coeff_a[k])
        u[k], ϵ[k] = midradius(u_temp)
    end
    u_temp = 2 * a * u[2] - u[3] + Interval{T}(coeff_a[1])
    u[1], ϵ[1] = midradius(u_temp)

    e = zeros(Interval{T}, m + 1)
    e[m] = coeff_r[m]
    for k in reverse(2:m-1)
        e[k] = e[k+1] + (k + 1) * (2 * r * abs(u[k+1]) + ϵ[k] + coeff_r[k])
    end
    e[1] = e[2] + 2 * r * abs(u[1]) + ϵ[1] + coeff_r[1]
    γ = sup(e[1])
    return u[1] + interval(-γ, γ)
end


function Clenshaw(coeff, x)
    n = length(coeff)
    u = zeros(typeof(x), n + 1)
    u[n] = coeff[n]
    for k in reverse(2:n-1)
        u[k] = coeff[k] + 2 * x * u[k+1] - u[k+2]
    end
    u[1] = coeff[1] + x * u[2] - u[3]
    return u[1]
end

function ClenshawSecond(coeff, x::T) where {T<:Real}
    n = length(coeff)
    u = zeros(T, n + 1)
    u[n] = coeff[n]
    for k in reverse(2:n-1)
        u[k] = coeff[k] + 2 * x * u[k+1] - u[k+2]
    end
    u[1] = coeff[1] + 2 * x * u[2] - u[3]
    return u[1]
end


Clenshaw(coeff, x::Interval{T}) where {T} = eval_Clenshaw_BackwardFirst(coeff, x)
function ChebyshevDerivative(coeff, x::Interval{T}) where {T}
    n = length(coeff)
    coeff_der = [(i - 1) * Interval{T}(coeff[i]) for i = 2:n]
    #@info coeff
    #@info coeff_der
    return eval_Clenshaw_BackwardSecond(coeff_der, x)
end

evalChebyshev(coeff, x::Interval) =
    eval_Clenshaw_BackwardFirst(coeff, interval(mid.(2 * x - 1)))
evalChebyshevDerivative(coeff, x::Interval) = 2 * ChebyshevDerivative(coeff, 2 * x - 1)
function evalChebyschevCentered(coeff, x::Interval)
    m = interval(mid.(x))
    return evalChebyshev(coeff, m) + evalChebyshevDerivative(coeff, x) * (x - m)
end
###############################################################################
###############################################################################

function weak_projection_error(B::Chebyshev{S,C1}) where {S}
    n = Float64(length(B), RoundUp)
    ν = B.k
    νf = Float64(B.k, RoundUp)
    den = n ⊗₋ (νf ⊖₋ 2.0) ⊗₋ Float64(π, RoundDown) ⊗₋ reduce(⊗₋, [n - i for i = 2:ν-1])
    return (4.0 ⊗₊ (n + 1)) ⊘₊ den
end
function aux_normalized_projection_error(B::Chebyshev)
    n = Float64(length(B), RoundUp)
    ν = B.k
    νf = Float64(B.k, RoundUp)
    den = Float64(π, RoundDown) ⊗₋ νf ⊗₋ n ⊗₋ reduce(⊗₋, [n - i for i = 1:ν-1])
    return 2.0 ⊘₊ den
end

"""
	strong_weak_bound(B::Chebyshev)

V.A. Markov estimate from GRADIMIR MILOVANOVIC EXTREMAL PROBLEMS AND 
INEQUALITIES OF MARKOV-BERNSTEIN TYPE FOR POLYNOMIALS
"""
# TODO: Check the indexes

function strong_weak_bound(B::Chebyshev{S,C1}) where {S}
    n = length(B) - 1
    k = B.k - 1
    # we want to estimate the norm of f^(k) by the C1 norm of f, so we use Markov estimate
    # for derivative k-1
    den = reduce(⊗₋, [Float64(2k - 1 - 2 * i, RoundDown) for i = 0:k-1])
    num = reduce(⊗₊, [Float64(n^2 - i^2, RoundDown) for i = 0:k-1])
    return num ⊘₊ den ⊕₊ 1.0 # the 1.0 is to take into account the L1 norm of f
end
# aux_weak_bound: M₂ with ||v||_{L¹(dx)} ≤ M₂ ||v||_{weak}.
#
# For weak = L²(dx), Cauchy-Schwarz on the probability measure dx gives 1; for
# weak = C1, ||v||_{L¹} ≤ ||v||_∞ ≤ ||v||_{C¹} gives 1 as well.
aux_weak_bound(B::Chebyshev{S,WK}) where {S<:NormKind,WK<:NormKind} = 1.0

"""
    weak_by_strong_and_aux_bound(B::Chebyshev)

Returns `(S₁, S₂)` such that `||f||_{C1} ≤ S₁·||f||_{W^{k,1}} + S₂·||f||_{L1}`.

For k ≥ 2: Sobolev embedding gives ||f||_∞ ≤ ||f||_{L1} + ||f'||_{L1} and
||f'||_∞ ≤ ||f'||_{L1} + ||f''||_{L1}, so ||f||_{C1} ≤ 2·||f||_{W^{k,1}}.

For k = 1: uses Markov inequality ||p'||_∞ ≤ 2n²·||p||_∞ for polynomials of
degree n on [0,1], giving ||f||_{C1} ≤ (1 + 2n²)·||f||_{W^{1,1}}.
"""
function weak_by_strong_and_aux_bound(B::Chebyshev{S,C1}) where {S}
    if B.k >= 2
        return (2.0, 0.0)
    else
        n = Float64(length(B) - 1, RoundUp)
        return (1.0 ⊕₊ 2.0 ⊗₊ n ⊗₊ n, 0.0)
    end
end

"""
    bound_weak_norm_from_linalg_norm(B::Chebyshev)

Returns `(W₁, W₂)` such that `||f||_{C1} ≤ W₁·||ĉ||_{ℓ¹} + W₂·||ĉ||_{ℓ∞}`.

Since |Tⱼ(x)| ≤ 1: ||f||_∞ ≤ ||ĉ||_{ℓ¹}.
For the derivative: ||f'||_∞ ≤ Σ|cⱼ|·2j² ≤ 2(n-1)²·||ĉ||_{ℓ¹} where n = degree.
"""
function bound_weak_norm_from_linalg_norm(B::Chebyshev{S,C1}) where {S}
    n = Float64(length(B) - 1, RoundUp)
    nm1 = n ⊖₋ 1.0
    W₁ = 1.0 ⊕₊ 2.0 ⊗₊ nm1 ⊗₊ nm1
    return (W₁, 0.0)
end

"""
    bound_linalg_norm_L1_from_weak(B::Chebyshev)

Returns `A` such that `||ĉ||_{ℓ¹} ≤ A·||f||_{C1}`.

Chebyshev coefficients satisfy |cⱼ| ≤ 2·||f||_∞ for j ≥ 1 and |c₀| ≤ ||f||_∞,
so ||ĉ||_{ℓ¹} ≤ (2n-1)·||f||_∞ ≤ (2n-1)·||f||_{C1} where n = length(B).
"""
function bound_linalg_norm_L1_from_weak(B::Chebyshev{S,C1}) where {S}
    n = Float64(length(B), RoundUp)
    return 2.0 ⊗₊ n ⊖₋ 1.0
end

"""
    bound_linalg_norm_L∞_from_weak(B::Chebyshev)

Returns `A` such that `||ĉ||_{ℓ∞} ≤ A·||f||_{C1}`.

max_j |cⱼ| ≤ 2·||f||_∞ ≤ 2·||f||_{C1}.
"""
function bound_linalg_norm_L∞_from_weak(B::Chebyshev{S,C1}) where {S}
    return 2.0
end
weak_norm(B::Chebyshev) = typeof(B.weak)
aux_norm(B::Chebyshev) = L1
# Weak L²(μ) pairs with auxiliary L¹(μ): both facts the estimates rest on --
# Parseval and |b̂_k| ≤ 2‖f‖_{L¹(μ)} -- live against μ, and matching the two
# measures makes aux_weak_bound a plain Cauchy-Schwarz 1.
aux_norm(B::Chebyshev{S,L2μ}) where {S<:NormKind} = L1μ
strong_norm(B::Chebyshev) = typeof(B.strong)

"""
	Base.getindex(B::Chebyshev, i::Int)

Make so that B[j] returns a HatFunctionOnTorus with the j-th basis element
"""
function Base.getindex(B::Chebyshev, i::Int)
    n = length(B)
    v = zeros(n)
    v[i] = 1
    @boundscheck 1 <= i <= n || throw(BoundsError(B, i))
    return x -> evalChebyshev(v, x)
end

is_refinement(Bc::Chebyshev, Bf::Chebyshev) = length(Bc) < length(Bf)
integral_covector(B::Chebyshev; T = Float64) =
    [Interval{T}(1); 0; [0.5 * Interval{T}((-1)^n + 1) / (1 - n^2) for n = 2:length(B)-1]]'
one_vector(B::Chebyshev) = [1; zeros(length(B) - 1)]

Base.length(S::AverageZero{T}) where {T<:Chebyshev} = length(S.basis) - 1

function Base.iterate(S::AverageZero{T}, state = 1) where {T<:Chebyshev}
    B = S.basis
    i = state
    if i == length(B)
        return nothing
    end
    v = zeros(length(B))
    v[i+1] = 1
    v[1] = -mid.(integral_covector(B)[i+1])
    return v, state + 1
end

struct ChebyshevDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    x′::Vector{Interval}
end

function ChebDualBranch(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter)
    if is_increasing(br)
        endpoint_X = br.X[2]
        der = derivative(br.f)(endpoint_X)
        preim_der = preimages_and_derivatives(y, br, ylabel; ϵ, max_iter)
        return [preim_der[1]; endpoint_X],
        [preim_der[2]; length(preim_der[2]) + 1],
        [preim_der[3]; der]
    else
        endpoint_X = br.X[2]
        der = derivative(br.f)(endpoint_X)
        preim_der = preimages_and_derivatives(B.p, D, 1:length(B.p)-1; ϵ, max_iter)
        return [preim_der[1]; endpoint_X],
        [preim_der[2]; length(preim_with_der[2]) + 1],
        [preim_der[3]; der]
    end
end

function Dual(B::Chebyshev, D::PwMap; ϵ, max_iter)
    @assert is_full_branch(D)
    results =
        collect(ChebDualBranch(B.p, b, 1:length(B.p)-1; ϵ, max_iter) for b in branches(D))
    x = vcat((result[1] for result in results)...)
    xlabel = vcat((result[2] for result in results)...)
    x′ = vcat((result[3] for result in results)...)
    return x, xlabel, x′
end

Base.length(dual::ChebyshevDual) = length(dual.x)
Base.eltype(dual::ChebyshevDual) =
    Tuple{eltype(dual.xlabel),Tuple{eltype(dual.x),eltype(dual.x′)}}
function Base.iterate(dual::ChebyshevDual, state = 1)
    if state <= length(dual.x)
        return ((dual.xlabel[state], (dual.x[state], abs(dual.x′[state]))), state + 1)
    else
        return nothing
    end
end


# chebtransform and assemble(::Chebyshev, ::Dynamic; …) live in the FFTWExt
# extension. Loading `using FFTW` provides them.

using IntervalOptimisation

function infnormoffunction(B::Chebyshev, v)
    val = 0
    try
        val = maximize(x -> abs(evalChebyschevCentered(v, x)), interval(0, 1))[1]
    catch
        print("Refining grid")
        f(x) = abs(evalChebyshevCentered(v, x))
        ran = range_estimate(f, interval(0, 1), 5)
        Bval = hull(val, ran)
    end
    return val
end

function infnormofderivative(B::Chebyshev, v)
    val = interval(0)
    try
        val = maximize(x -> abs(evalChebyshevDerivative(v, x)), interval(0, 1))[1]
    catch
        print("Refining grid")
        f(x) = abs(evalChebyshevDerivative(v, x))
        ran = range_estimate(f, interval(0, 1), 5)
        val = hull(val, ran)
    end
    return val
end


is_integral_preserving(B::Chebyshev) = false

"""
    restrict_to_average_zero(B::Chebyshev, BM::BallMatrix, f; certification=nothing)

Restrict `BM` to the average-zero subspace U⁰ using a certified spectral projector.

For Chebyshev, `integral_covector(B) ≠ [1,0,...,0]`, so simple submatrix extraction
does not work. Instead, we compute the Riesz projector P₁ onto the eigenspace of
eigenvalue 1 (the invariant measure direction) via Schur decomposition, then return
`(I - P₁) * BM`.

If `certification` is provided (from `certify_spectral_gap`), the projector radii are
inflated by the certified projector error δ_P to account for the Schur factorization
error ‖ZTZ* - A‖₂.
"""
function restrict_to_average_zero(B::Chebyshev, BM::BallMatrix, f;
                                   certification = nothing)
    n = size(BM, 1)

    # Find eigenvalue closest to 1 via Schur decomposition
    F = LinearAlgebra.schur(Complex{Float64}.(BM.c))
    eigenvalues = diag(F.T)
    idx_1 = argmin(abs.(eigenvalues .- 1.0))

    # Compute Schur-based spectral projector for eigenvalue 1
    proj_result = compute_spectral_projector_schur(BM, [idx_1])
    P1 = proj_result.projector

    # Inflate projector radii by certified error if available
    if certification !== nothing
        δ_P = certification.projector_error
        P1 = BallMatrix(P1.c, P1.r .+ Float64(δ_P))
    end

    # (I - P₁) · BM acts on complement of eigenvalue 1
    I_n = BallMatrix(Matrix{eltype(BM.c)}(LinearAlgebra.I, n, n))
    return (I_n - P1) * BM
end

"""
    certify_spectral_gap(B::Chebyshev, Q::DiscretizedOperator;
                          samples=256, radius_factor=0.5)

Run the full spectral gap certification pipeline for a Chebyshev discretized operator.

Uses BallArithmetic's CertifScripts to:
1. Compute Schur decomposition with rigorous error bounds
2. Find the spectral gap (distance from eigenvalue 1 to next largest eigenvalue)
3. Certify the resolvent on a circle separating eigenvalue 1 from the rest
4. Compute projector error bound δ_P for `restrict_to_average_zero`

Returns a named tuple with fields:
- `spectral_gap`, `second_largest`, `ρ` (contour radius)
- `M_inf` (certified resolvent bound on contour)
- `projector_error` (δ_P for projector inflation)
- `certification`, `schur_data`, `eigenvalue_index`
"""
function certify_spectral_gap(B::Chebyshev, Q::DiscretizedOperator;
                               samples = 256, radius_factor = 0.5)
    BM = BallMatrix(Q.L)
    n = size(BM, 1)

    # Step 1: Schur decomposition with error bounds
    schur_data = CertifScripts.compute_schur_and_error(BM)
    S, errF, errT, norm_Z, norm_Z_inv = schur_data

    # Step 2: Find eigenvalue 1 and spectral gap
    eigenvalues = diag(S.T)
    idx_1 = argmin(abs.(eigenvalues .- 1.0))
    second_largest = maximum(abs(eigenvalues[i]) for i in 1:n if i != idx_1)
    spectral_gap = 1.0 - second_largest

    # Step 3: Choose contour Γ separating eigenvalue 1
    ρ = second_largest + radius_factor * spectral_gap
    circle = CertifScripts.CertificationCircle(0.0 + 0.0im, ρ; samples = samples)

    # Step 4: Certify resolvent on Γ (accounts for Schur→original error)
    cert = CertifScripts.run_certification(BM, circle; schur_data = schur_data)
    M_inf = cert.resolvent_original

    # Step 5: Projector error bound (Gauss.pdf Eq 16)
    # δ_P = ρ · M²_∞ · ε_K / (1 - ε_K · M_∞)
    ε_K = Float64(cert.errT)
    α = ε_K * M_inf
    if α >= 1.0
        @error "Small-gain condition fails: α = ε_K · M_∞ = $α ≥ 1"
    end
    δ_P = ρ * M_inf^2 * ε_K / (1.0 - α)

    return (spectral_gap = spectral_gap,
            second_largest = second_largest,
            ρ = ρ,
            M_inf = M_inf,
            projector_error = δ_P,
            α = α,
            certification = cert,
            schur_data = schur_data,
            eigenvalue_index = idx_1)
end

function opnormbound(B::Chebyshev, N::Type{C1}, v::Vector{S}) where {S}
    return normbound(B, N, v)
end

function opnormbound(B::Chebyshev, N::Type{C1}, w::LinearAlgebra.Adjoint)
    return normbound(B, N, w')
end

function opnormbound(B::Chebyshev, N::Type{C1}, A::Matrix{S}) where {S}
    n, m = size(A)
    norm = 0.0
    for i = 1:m
        norm = max(
            norm,
            Float64(opnormbound(B, N, A[:, i]), RoundUp),
        )
    end
    # D = log(N+1) bounds ||·||_{ℓ¹} ≤ D ||·||_{C¹} via Chebyshev coefficient decay
    # (Theorem 3.12 in Nisoli–Taylor-Crush 2023); columns must be unnormalized.
    return norm ⊗₊ Float64(log(m + 1), RoundUp)
end

normbound(B::Chebyshev, N::Type{C1}, v) =
    Float64(sup(infnormoffunction(B, v) + infnormofderivative(B, v)), RoundUp)

mutable struct NormCacherC1 <: NormCacher{C1}
    B::Basis
    C::Float64
    function NormCacherC1(B, n)
        new(B, 0.0)
    end
end
NormCacher{C1}(B, n) = NormCacherC1(B, n)

function add_column!(Cacher::NormCacherC1, v::AbstractVector, ε::Float64)
    Cacher.C = max(Cacher.C, opnormbound(Cacher.B, C1, v) ⊕₊ ε)
end

"""
Return the norm of the matrix the NormCacher is working on.
"""
function get_norm(Cacher::NormCacherC1)
    n = length(Cacher.B)
    return Cacher.C ⊗₊ Float64(log(n + 1), RoundUp)
end

###############################################################################
# Gram matrix
###############################################################################

@doc raw"""
    gram_matrix(B::Chebyshev; T = Float64)
    inv_gram_matrix(B::Chebyshev; T = Float64)

Gram matrix ``G_{ij} = \langle φ_i, φ_j \rangle`` of the Chebyshev basis
``φ_i(x) = T_{i-1}(2x-1)`` on ``[0,1]``, and its inverse.

The inner product is taken against the **arcsine probability measure**

```math
dμ(x) = \frac{dx}{π\sqrt{x(1-x)}}, \qquad \int_0^1 dμ = 1,
```

the measure that makes the Chebyshev polynomials orthogonal. Substituting
``t = 2x-1`` turns it into the familiar ``dt/(π\sqrt{1-t^2})``, so

```math
G = \operatorname{diag}(1, \tfrac12, \tfrac12, \dots, \tfrac12),
\qquad
G^{-1} = \operatorname{diag}(1, 2, 2, \dots, 2).
```

(Against the unnormalized weight ``1/\sqrt{1-t^2}`` every entry is ``π`` times
these.) Both are returned as `Diagonal` of intervals; every entry is exactly
representable, so the enclosures are thin.

# Transporting an ``\ell^2`` bound to the function space

For ``f = \sum_i c_i φ_i`` we have ``\|f\|_{L^2(μ)}^2 = c^* G c``, so an
operator with coefficient matrix `Q` satisfies

```math
\|Q\|_{L^2(μ) \to L^2(μ)} = \|G^{1/2}\,Q\,G^{-1/2}\|_{\ell^2}.
```

This is the point of using this weight rather than Lebesgue: `G` is diagonal,
so the similarity is a cheap rescaling and `BallArithmetic`'s ``\ell^2``
estimators (`upper_bound_L2_opnorm`, `svd_bound_L2_opnorm_inverse`, …) apply to
`gram_sqrt(B) * Q * inv_gram_sqrt(B)` directly. Under Lebesgue the Gram matrix
is dense and ill-conditioned and no such shortcut exists.

!!! warning
    ``L^2(μ)`` is *not* ``L^2(\mathrm{Leb})``. The rest of the package — the
    Lasota–Yorke constants, `integral_covector`, the invariant-density
    normalization — is set in Lebesgue. Norms computed here refer to the
    weighted space; the spectrum is unchanged by the similarity, but constants
    are not interchangeable.

See also [`gram_sqrt`](@ref), [`inv_gram_sqrt`](@ref).
"""
function gram_matrix(B::Chebyshev; measure::Symbol = :arcsine, T = Float64)
    n = length(B)
    if measure === :arcsine
        d = fill(interval(T, 1) / interval(T, 2), n)
        d[1] = interval(T, 1)
        return LinearAlgebra.Diagonal(d)
    elseif measure === :lebesgue
        return _lebesgue_gram(n, T)
    else
        throw(ArgumentError("measure must be :arcsine or :lebesgue, got $measure"))
    end
end

# G_ij = 1/2 [ 1/(1-(m+n)^2) + 1/(1-(m-n)^2) ] for m+n even, 0 otherwise,
# with m = i-1, n = j-1. The denominators are exact integers, so each entry is
# a single correctly-rounded division. m+n and m-n share parity, so the two
# vanishing denominators (|m±n| = 1) only occur when the entry is 0 anyway.
function _lebesgue_gram(n::Integer, ::Type{T}) where {T}
    G = zeros(Interval{T}, n, n)
    one_T = interval(T, 1)
    half = one_T / interval(T, 2)
    for i = 1:n, j = 1:n
        m, k = i - 1, j - 1
        isodd(m + k) && continue
        G[i, j] =
            half * (one_T / interval(T, 1 - (m + k)^2) + one_T / interval(T, 1 - (m - k)^2))
    end
    return G
end

function inv_gram_matrix(B::Chebyshev; measure::Symbol = :arcsine, T = Float64)
    n = length(B)
    if measure === :arcsine
        d = fill(interval(T, 2), n)
        d[1] = interval(T, 1)
        return LinearAlgebra.Diagonal(d)
    elseif measure === :lebesgue
        return _verified_inverse(_lebesgue_gram(n, T))
    else
        throw(ArgumentError("measure must be :arcsine or :lebesgue, got $measure"))
    end
end

# Verified inverse, one Krawczyk-verified linear solve per column. The Lebesgue
# Gram matrix is only mildly ill-conditioned (cond ~ 1.3n measured), so this
# converges comfortably.
function _verified_inverse(G::Matrix{Interval{T}}) where {T}
    n = size(G, 1)
    GB = BallMatrix(G)
    X = zeros(Interval{T}, n, n)
    for j = 1:n
        c = zeros(T, n)
        c[j] = one(T)
        res = krawczyk_linear_system(GB, BallVector(c, zeros(T, n)))
        res.verified ||
            error("Krawczyk verification failed for column $j of the Gram inverse")
        col = res.solution
        for i = 1:n
            ci, ri = BallArithmetic.mid(col[i]), BallArithmetic.rad(col[i])
            X[i, j] = interval(T, ci - ri, ci + ri)
        end
    end
    return X
end

@doc raw"""
    l2_measure_conversion_bounds(B::Chebyshev; T = Float64) -> (c_leb, C_n)

The two constants relating the Lebesgue and arcsine ``L^2`` norms on the span
of this basis:

```math
\|f\|_{L^2(dx)} \le c_{\mathrm{leb}}\,\|f\|_{L^2(μ)},
\qquad
\|f\|_{L^2(μ)} \le C_n\,\|f\|_{L^2(dx)}.
```

The first is uniform and needs no linear algebra: ``dx/dμ = π\sqrt{x(1-x)}`` is
bounded by ``π/2``, so ``c_{\mathrm{leb}} = \sqrt{π/2}`` on all of ``L^2(μ)``.

The second cannot hold uniformly — ``dμ/dx`` blows up at the endpoints — and is
genuinely finite-dimensional:
``C_n^2 = λ_{\max}(G_μ^{1/2} G_L^{-1} G_μ^{1/2})``. It is evaluated sharply here,
from the verified inverse and the diagonal ``G_μ^{1/2}``. Bounding it instead by
``\|G_L^{-1}\|_2`` — valid since ``\|G_μ^{1/2}\| = 1``, and cheaper in that it
needs no inverse — costs a factor tending to ``\sqrt2``, the extremal direction
lying in the ``\tfrac12``-eigenspace of ``G_μ``. Measured growth of the sharp
constant is ``C_n \approx 0.75\sqrt{n}``.

Use these to carry a resolvent or spectral bound proved in ``L^2(μ)`` — where
Parseval makes it an ``\ell^2`` statement about the coefficient matrix — back to
the ``L^2(dx)`` setting the rest of the package works in.
"""
function l2_measure_conversion_bounds(B::Chebyshev; T = Float64)
    n = length(B)
    c_leb = sqrt(interval(T, π) / interval(T, 2))

    # C_n^2 = ||G_mu^{1/2} G_L^{-1} G_mu^{1/2}||_2, the middle factor verified
    # and the outer ones diagonal, so the product is cheap and the constant sharp.
    Ginv = _verified_inverse(_lebesgue_gram(n, T))
    s = gram_sqrt(B; T = T)
    return (sup(c_leb), sqrt(upper_bound_L2_opnorm(BallMatrix(s * Ginv * s))))
end

@doc raw"""
    gram_sqrt(B::Chebyshev; T = Float64)
    inv_gram_sqrt(B::Chebyshev; T = Float64)

The symmetric factors ``G^{1/2}`` and ``G^{-1/2}`` of the Chebyshev
[`gram_matrix`](@ref) — i.e. ``\operatorname{diag}(1, 1/\sqrt2, \dots)`` and
``\operatorname{diag}(1, \sqrt2, \dots)``.

These, not `G` itself, are what an ``\ell^2`` operator-norm estimator needs:
``\|Q\|_{L^2(μ)} = \|G^{1/2} Q G^{-1/2}\|_{\ell^2}``. Unlike `G`, the entries
are irrational, so the returned intervals are genuine (thin but not exact)
enclosures.
"""
function gram_sqrt(B::Chebyshev; T = Float64)
    n = length(B)
    d = fill(interval(T, 1) / sqrt(interval(T, 2)), n)
    d[1] = interval(T, 1)
    return LinearAlgebra.Diagonal(d)
end

function inv_gram_sqrt(B::Chebyshev; T = Float64)
    n = length(B)
    d = fill(sqrt(interval(T, 2)), n)
    d[1] = interval(T, 1)
    return LinearAlgebra.Diagonal(d)
end

@doc raw"""
    gram_restrict_to_average_zero(B::Chebyshev, BM::BallMatrix; T = Float64)

Restrict `BM` to the average-zero subspace by the Gram change of variables,
returning `(block, chol)` where `block` is an `(n-1) × (n-1)` `BallMatrix`.

Unlike [`restrict_to_average_zero`](@ref), which builds a certified Riesz
projector via a Schur decomposition, this uses the fact that in the
Gram-transformed coordinates the restriction *is* a plain submatrix.

Writing ``G`` for the Lebesgue Gram matrix ([`gram_matrix`](@ref) with
`measure = :lebesgue`) and ``G = U^{*}U`` for its Cholesky factor:

- ``T_0 = 1``, so the constant function is ``e_1`` and the integral covector is
  exactly ``v = G e_1`` — the first column of ``G``;
- ``U`` is upper triangular, hence ``U e_1 \parallel e_1`` and
  ``U^{-*} v = U_{11} e_1``, so average-zero ``\{v \cdot c = 0\}`` becomes
  ``\{y_1 = 0\}`` in the coordinates ``y = Uc``;
- since the transfer operator preserves the integral, ``v^{*}Q = v^{*}``, and
  therefore ``\tilde A = U Q U^{-1}`` has first row exactly ``e_1^{*}``.

The restriction is then ``\tilde A[2:\mathrm{end}, 2:\mathrm{end}]`` — the same
shape as `restrict_to_average_zero(B::Fourier, …)`, of which this is the
special case ``G = I``. Because ``\|M\|_G = \|U M U^{-1}\|_2`` exactly, an
``\ell^2`` bound on the returned block *is* the ``L^2(dx)`` bound on the
average-zero subspace, with no condition-number penalty.

!!! warning "What is and is not certified"
    `chol` is the [`verified_cholesky`](@ref) result for the **midpoint** of
    ``G``, so the enclosure covers the factorization but not the ≤1 ulp
    enclosure radius of the Gram entries themselves. The induced norm is
    therefore that of ``\tilde G = U^{*}U`` rather than of the exact ``G``;
    `chol.residual_norm` and `maximum(radius.(gram_matrix(B; measure = :lebesgue)))`
    quantify the gap. Closing it needs an interval-aware Cholesky, or
    BallArithmetic's `gram_transform`, which reports `gram_residual` directly.
"""
function gram_restrict_to_average_zero(B::Chebyshev, BM::BallMatrix; T = Float64)
    n = length(B)
    size(BM, 1) == n ||
        throw(DimensionMismatch("operator is $(size(BM,1))×$(size(BM,2)), basis has $n"))

    G = _lebesgue_gram(n, T)
    chol = BallArithmetic.verified_cholesky(T.(mid.(G)); use_bigfloat = false)
    chol.success || error("verified Cholesky of the Lebesgue Gram matrix failed")

    # G = U'U with U upper triangular. Move to intervals to invert, then back.
    U = chol.G
    Ui = Interval{T}[
        interval(T, U.c[i, j] - U.r[i, j], U.c[i, j] + U.r[i, j]) for i = 1:n, j = 1:n
    ]
    Ã = BallMatrix(Ui) * BM * BallMatrix(_verified_inverse(Ui))
    return (BallMatrix(Ã.c[2:end, 2:end], Ã.r[2:end, 2:end]), chol)
end


###############################################################################
# Norm interface for the L² weak norms
#
# The T_m are orthogonal against the arcsine measure μ, so with c the coefficient vector
#
#     ||f||_{L²(μ)}² = c₀² + ½ Σ_{m≥1} cₘ²,                                 (*)
#
# i.e. ||c||_{ℓ²}/√2 ≤ ||f||_{L²(μ)} ≤ ||c||_{ℓ²}.  Every bound below follows
# from (*) alone, exactly as the Fourier ones follow from Parseval.
#
# Under `L2` (Lebesgue) each bound is the arcsine one composed with the conversion
# constants of `l2_measure_conversion_bounds`:
#
#     ||f||_{L²(dx)} ≤ √(π/2) ||f||_{L²(μ)},   ||f||_{L²(μ)} ≤ Cₙ ||f||_{L²(dx)},
#
# the first uniform (dx/dμ = π√(x(1-x)) ≤ π/2), the second finite-dimensional.
###############################################################################

# ||v||_{L²} ≤ S₁ ||v||_s + S₂ ||v||_{L¹}: ||v||_{L²} ≤ ||v||_∞ ≤ ||v||_{W^{1,1}}
# ≤ ||v||_{W^{k,1}} by Sobolev embedding on [0,1], for either measure.
weak_by_strong_and_aux_bound(B::Chebyshev{S,L2}) where {S} = (1.0, 0.0)

# ||v||_{L²} ≤ W₁ ||ĉ||_{ℓ¹} + W₂ ||ĉ||_{ℓ∞}
# Both measures are probability measures and |T_m| ≤ 1, so
# ||v||_{L²} ≤ ||v||_∞ ≤ ||ĉ||_{ℓ¹} directly — no conversion factor.
bound_weak_norm_from_linalg_norm(B::Chebyshev{S,L2}) where {S} = (1.0, 0.0)

# ||ĉ||_{ℓ¹} ≤ A ||v||_{L²}:  ||c||_{ℓ¹} ≤ √n ||c||_{ℓ²} ≤ √(2n) ||v||_{L²(μ)}
function bound_linalg_norm_L1_from_weak(B::Chebyshev{S,L2}) where {S}
    # ||c||_{ℓ¹} ≤ √n ||c||_{ℓ²} ≤ √(2n) ||v||_{L²(μ)} ≤ √(2n) Cₙ ||v||_{L²(dx)},
    # the last step by `l2_measure_conversion_bounds`.
    _, C_n = l2_measure_conversion_bounds(B)
    return sqrt_round(2.0 ⊗₊ Float64(length(B), RoundUp), RoundUp) ⊗₊ C_n
end

# ||ĉ||_{ℓ∞} ≤ A ||v||_{L²}:  ||c||_{ℓ∞} ≤ ||c||_{ℓ²} ≤ √2 ||v||_{L²(μ)}
function bound_linalg_norm_L∞_from_weak(B::Chebyshev{S,L2}) where {S}
    # ||c||_{ℓ∞} ≤ ||c||_{ℓ²} ≤ √2 ||v||_{L²(μ)} ≤ √2 Cₙ ||v||_{L²(dx)}.
    _, C_n = l2_measure_conversion_bounds(B)
    return sqrt_round(2.0, RoundUp) ⊗₊ C_n
end

@doc raw"""
    weak_projection_error(B::Chebyshev{S, L2})

``L^2`` projection error for the ``W^{ν,1}`` unit ball, ``ν`` = `B.k`.

This follows from the ``W^{k,1}`` decay of the Chebyshev coefficients, in two
steps, both from Nisoli & Taylor-Crush:

- Theorem 3.12 (their statement of Trefethen, *ATAP*, Thm 7.1): if ``f^{(ν)}``
  has bounded variation ``V`` then ``|\hat b_m| \le 2V/(π\,m(m-1)\cdots(m-ν))``;
- Theorem 3.13: hence ``\|f - π_n f\|_∞ \le 2V/(π ν\, n(n-1)\cdots(n+1-ν))``.

Both ``dx`` and ``dμ`` are *probability* measures on ``[0,1]``, so
``\|\cdot\|_{L^2} \le \|\cdot\|_∞`` with constant 1 and the same bound serves
either measure with no conversion factor. It is the quantity already computed by
[`aux_normalized_projection_error`](@ref).

Theorem 3.13 is stated for the interpolant ``π_n``, which is what this basis
uses (the coefficients come from an FFT at the Chebyshev points), so the
aliasing is already accounted for.

!!! note "A sharper bound is available for the orthogonal projection"
    Applying Parseval to Theorem 3.12 directly gives
    ``\|f - \hat π_n f\|_{L^2(μ)} \le (V\sqrt2/(π\sqrt{2ν+1}))\,(n-ν)^{-(ν+1/2)}``,
    half a power better. That is a bound on the **orthogonal** projection
    ``\hat π_n`` only; transferring it to the interpolant costs an ``\ell^1``
    aliasing estimate which gives the half power straight back, so it is not
    used here.
"""
weak_projection_error(B::Chebyshev{S,L2}) where {S} =
    aux_normalized_projection_error(B)

@doc raw"""
    strong_weak_bound(B::Chebyshev{S, L2})

``\|v\|_{W^{k,1}} \le M \|v\|_{L^2}`` on the span of the basis, obtained by
composing the existing `C1` estimate with ``\|v\|_{C^1} \le M' \|v\|_{L^2}``,
where ``M'`` comes from ``\|v\|_∞ \le \|\hat c\|_{ℓ^1}`` and the Markov bound
``\|v'\|_∞ \le 2(n-1)^2 \|\hat c\|_{ℓ^1}`` already used by
[`bound_weak_norm_from_linalg_norm`](@ref).
"""
function strong_weak_bound(B::Chebyshev{S,L2}) where {S}
    B_c1 = Chebyshev(B.p, B.k, C1)
    W₁, _ = bound_weak_norm_from_linalg_norm(B_c1)      # ||v||_{C1} ≤ W₁ ||ĉ||_{ℓ¹}
    A = bound_linalg_norm_L1_from_weak(B)               # ||ĉ||_{ℓ¹} ≤ A ||v||_{L²}
    return strong_weak_bound(B_c1) ⊗₊ W₁ ⊗₊ A
end

###############################################################################
# Norm interface for the Bernstein-ellipse strong norm
#
# Trefethen, ATAP, Thm 8.1: f analytic in E_ρ with |f| ≤ M has Chebyshev
# coefficients |a_0| ≤ M, |a_k| ≤ 2Mρ^{-k}; Thm 8.2 (8.3): the *interpolant*
# through n+1 Chebyshev points then satisfies ||f - p_n||_∞ ≤ 4Mρ^{-n}/(ρ-1).
# The basis interpolates, so (8.3) is the relevant one.
###############################################################################

# Degree of the interpolant: length(B) points carry degree length(B) - 1.
_cheb_degree(B::Chebyshev) = Float64(length(B) - 1, RoundDown)

@doc raw"""
    weak_projection_error(B::Chebyshev{Eρ, L2})

``\|f - p_n\|_{L^2(dx)} \le 4ρ^{-n}/(ρ-1)`` on the unit ball of the
[`Eρ`](@ref) norm, from Trefethen, *ATAP*, Thm 8.2 (8.3) — the interpolation
form, since this basis interpolates at the Chebyshev points. `dx` is a
probability measure on ``[0,1]``, so the sup-norm bound carries to ``L^2`` with
constant 1.
"""
function _bernstein_projection_error(B::Chebyshev)
    ρ = B.strong.ρ
    ρ > 1 || return Inf
    n = _cheb_degree(B)
    return (4.0 ⊘₊ (ρ ⊖₋ 1.0)) ⊗₊ (ρ^(-n))
end
weak_projection_error(B::Chebyshev{Eρ,L2}) = _bernstein_projection_error(B)

aux_normalized_projection_error(B::Chebyshev{Eρ,L2}) = weak_projection_error(B)

# ||v||_{L²(dx)} ≤ ||v||_∞ ≤ ||v||_{E_ρ}, since [-1,1] ⊂ E_ρ.
weak_by_strong_and_aux_bound(B::Chebyshev{Eρ,L2}) = (1.0, 0.0)

@doc raw"""
    strong_weak_bound(B::Chebyshev{Eρ, L2})

``\|v\|_{E_ρ} \le M\,\|v\|_{L^2(dx)}`` on the span of the basis. On ``E_ρ`` one
has ``|T_k| \le (ρ^k + ρ^{-k})/2 \le ρ^k``, so
``\|v\|_{E_ρ} \le ρ^{n}\|\hat c\|_{ℓ^1}``, and
[`bound_linalg_norm_L1_from_weak`](@ref) supplies ``\|\hat c\|_{ℓ^1}`` in terms
of the weak norm.
"""
function strong_weak_bound(B::Chebyshev{Eρ,L2})
    ρ = B.strong.ρ
    n = _cheb_degree(B)
    return (ρ^n) ⊗₊ bound_linalg_norm_L1_from_weak(B)
end

###############################################################################
# Bernstein ellipses
###############################################################################

@doc raw"""
    bernstein_point(ρ, θ)

The point of the Bernstein ellipse ``∂E_ρ`` at parameter `θ ∈ [0,1]`, namely
``z = (w + w^{-1})/2`` with ``w = ρ\,e^{2πiθ}``, which expands to

```math
z = \frac{ρ + ρ^{-1}}{2}\cos 2πθ \;+\; i\,\frac{ρ - ρ^{-1}}{2}\sin 2πθ .
```

`ρ` and `θ` may be intervals, in which case the result encloses the
corresponding arc.
"""
function bernstein_point(ρ, θ)
    a = (ρ + 1 / ρ) / 2
    b = (ρ - 1 / ρ) / 2
    twoπθ = 2 * interval(π) * θ
    return complex(a * cos(twoπθ), b * sin(twoπθ))
end

@doc raw"""
    bernstein_parameter(z)

The parameter ``ρ ≥ 1`` of the Bernstein ellipse through the point `z`.

``E_ρ`` has foci ``\pm 1`` and major axis ``ρ + ρ^{-1}``, so the focal-distance
characterisation of the ellipse gives

```math
s := |z-1| + |z+1| = ρ + ρ^{-1},
\qquad
ρ = \frac{s + \sqrt{s^2-4}}{2}.
```

Working through `s` keeps everything real: no complex square root, and hence no
branch cut to worry about (the two Joukowski preimages ``w`` and ``w^{-1}`` of
`z` give the same answer by construction). `z` may be a complex interval, and
the result then encloses the parameters of all points it contains.
"""
function bernstein_parameter(z)
    s = _cabs(z - 1) + _cabs(z + 1)
    disc = s * s - 4
    # s ≥ 2 always; clip the rounding of s² - 4 at 0 for z on [-1,1], where ρ = 1.
    disc = intersect_interval(disc, interval(0, Inf))
    isempty_interval(disc) && (disc = interval(0, 0))
    return (s + sqrt(disc)) / 2
end

_cabs(z) = sqrt(real(z) * real(z) + imag(z) * imag(z))

@doc raw"""
    bernstein_expansion(f, ρ; n = 1024) -> ρ_image

Rigorous lower bound on

```math
\min_{z \in ∂E_ρ} \; \texttt{bernstein\_parameter}(f(z)) ,
```

obtained by covering the boundary parameter with `n` intervals, enclosing the
image of each arc under `f`, and taking the smallest enclosure. `f` must accept
a `Complex{Interval}`.

`f` maps ``E_ρ`` strictly outside itself — it *expands* the ellipse — exactly
when the returned value exceeds `ρ`; see [`expands_bernstein_ellipse`](@ref).
This is the Chebyshev analogue of enclosing the image of an annulus and reading
off its inner and outer radii, as done for the Fourier/analytic setting; an
ellipse needs only the one number because ``∂E_ρ`` is a single curve.

# Example

The Chebyshev polynomials are exactly the maps ``T_m(E_ρ) = E_{ρ^m}``, so
`bernstein_expansion(z -> 2z^2 - 1, ρ)` returns ``ρ^2``.
"""
function bernstein_expansion(f, ρ; n::Integer = 1024)
    ρi = ρ isa Interval ? ρ : interval(ρ)
    lo = Inf
    for j = 1:n
        θ = interval((j - 1) / n, j / n)
        val = bernstein_parameter(f(bernstein_point(ρi, θ)))
        lo = min(lo, inf(val))
    end
    return lo
end

@doc raw"""
    expands_bernstein_ellipse(f, ρ; n = 1024) -> (expands, ρ_image)

Whether `f` maps ``∂E_ρ`` strictly outside ``E_ρ``, together with the certified
image parameter from [`bernstein_expansion`](@ref).

For an interval map this is the analytic expansion condition behind an
[`Eρ`](@ref) Lasota–Yorke inequality: if the forward map expands the ellipse
then its inverse branches contract into it, so the transfer operator preserves
analyticity on ``E_ρ``.
"""
function expands_bernstein_ellipse(f, ρ; n::Integer = 1024)
    ρ_image = bernstein_expansion(f, ρ; n = n)
    return (ρ_image > sup(ρ isa Interval ? ρ : interval(ρ)), ρ_image)
end

@doc raw"""
    to_symmetric_interval(f)

Conjugate a map of ``[0,1]`` into a map of ``[-1,1]``, `t ↦ 2f((t+1)/2) - 1`.

Bernstein ellipses live around ``[-1,1]`` while the dynamics in this package
live on ``[0,1]``, so this is what to wrap a branch in before handing it to
[`bernstein_expansion`](@ref).
"""
to_symmetric_interval(f) = t -> 2 * f((t + 1) / 2) - 1

###############################################################################
# Norm interface for weak L²(μ) / auxiliary L¹(μ)
#
# Parseval: ||f||²_{L²(μ)} = b̂₀² + ½ Σ_{k≥1} b̂ₖ², so
# ||ĉ||_{ℓ²}/√2 ≤ ||f||_{L²(μ)} ≤ ||ĉ||_{ℓ²}.  Every constant below follows from
# that alone, exactly as the Fourier ones follow from Parseval on the circle --
# which is the point of measuring against μ rather than dx.
###############################################################################

# ||v||_{L¹(μ)} ≤ ||v||_{L²(μ)}: Cauchy-Schwarz, μ a probability measure.
aux_weak_bound(B::Chebyshev{S,L2μ}) where {S<:NormKind} = 1.0

# ||v||_{L²(μ)} ≤ ||v||_∞ ≤ ||v||_{W^{k,1}} by Sobolev embedding on [0,1].
weak_by_strong_and_aux_bound(B::Chebyshev{S,L2μ}) where {S} = (1.0, 0.0)
# ...and ≤ ||v||_{E_ρ}, since [-1,1] ⊂ E_ρ.
weak_by_strong_and_aux_bound(B::Chebyshev{Eρ,L2μ}) = (1.0, 0.0)

# ||v||_{L²(μ)} ≤ ||ĉ||_{ℓ²} ≤ ||ĉ||_{ℓ¹}
bound_weak_norm_from_linalg_norm(B::Chebyshev{S,L2μ}) where {S} = (1.0, 0.0)

# ||ĉ||_{ℓ¹} ≤ √n ||ĉ||_{ℓ²} ≤ √(2n) ||v||_{L²(μ)} — no conversion factor, this
# is where measuring against μ pays for itself.
bound_linalg_norm_L1_from_weak(B::Chebyshev{S,L2μ}) where {S} =
    sqrt_round(2.0 ⊗₊ Float64(length(B), RoundUp), RoundUp)

# ||ĉ||_{ℓ∞} ≤ ||ĉ||_{ℓ²} ≤ √2 ||v||_{L²(μ)}
bound_linalg_norm_L∞_from_weak(B::Chebyshev{S,L2μ}) where {S} = sqrt_round(2.0, RoundUp)

# μ is a probability measure, so ||·||_{L²(μ)} ≤ ||·||_∞ and the sup-norm
# projection bounds carry over unchanged.
weak_projection_error(B::Chebyshev{S,L2μ}) where {S} =
    aux_normalized_projection_error(B)
weak_projection_error(B::Chebyshev{Eρ,L2μ}) = _bernstein_projection_error(B)
aux_normalized_projection_error(B::Chebyshev{Eρ,L2μ}) = _bernstein_projection_error(B)

function strong_weak_bound(B::Chebyshev{S,L2μ}) where {S}
    B_c1 = Chebyshev(B.p, B.k, C1)
    W₁, _ = bound_weak_norm_from_linalg_norm(B_c1)
    return strong_weak_bound(B_c1) ⊗₊ W₁ ⊗₊ bound_linalg_norm_L1_from_weak(B)
end

function strong_weak_bound(B::Chebyshev{Eρ,L2μ})
    return (B.strong.ρ^_cheb_degree(B)) ⊗₊ bound_linalg_norm_L1_from_weak(B)
end
