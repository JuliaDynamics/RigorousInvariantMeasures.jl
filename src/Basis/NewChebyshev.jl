struct Chebyshev{T<:AbstractVector} <: Basis
    p::T
    k::Integer
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

function Chebyshev(n::Integer, k::Integer; T = Float64)
    return Chebyshev(ChebPoints(n, T), k)
end
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
    a, r = midpoint_radius(x)
    m = length(coeff)
    u = zeros(Interval{T}, m + 1)
    ϵ = zeros(Interval{T}, m + 1)
    u[m] = coeff_a[m]
    for k in reverse(2:m-1)
        u_temp = 2 * a * u[k+1] - u[k+2] + Interval{T}(coeff_a[k])
        u[k], ϵ[k] = midpoint_radius(u_temp)
    end
    u_temp = a * u[2] - u[3] + Interval{T}(coeff_a[1])
    u[1], ϵ[1] = midpoint_radius(u_temp)

    e = zeros(Interval{T}, m + 1)
    e[m] = coeff_r[m]
    for k in reverse(2:m-1)
        e[k] = e[k+1] + 2 * r * abs(u[k+1]) + ϵ[k] + coeff_r[k]
    end
    e[1] = e[2] + r * abs(u[1]) + ϵ[1] + coeff_r[1]
    γ = e[1].hi
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
    a, r = midpoint_radius(x)
    m = length(coeff)
    u = zeros(Interval{T}, m + 1)
    ϵ = zeros(Interval{T}, m + 1)
    u[m] = coeff_a[m]
    for k in reverse(2:m-1)
        u_temp = 2 * a * u[k+1] - u[k+2] + Interval{T}(coeff_a[k])
        u[k], ϵ[k] = midpoint_radius(u_temp)
    end
    u_temp = 2 * a * u[2] - u[3] + Interval{T}(coeff_a[1])
    u[1], ϵ[1] = midpoint_radius(u_temp)

    e = zeros(Interval{T}, m + 1)
    e[m] = coeff_r[m]
    for k in reverse(2:m-1)
        e[k] = e[k+1] + (k + 1) * (2 * r * abs(u[k+1]) + ϵ[k] + coeff_r[k])
    end
    e[1] = e[2] + 2 * r * abs(u[1]) + ϵ[1] + coeff_r[1]
    γ = e[1].hi
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

function weak_projection_error(B::Chebyshev)
    n = Float64(length(B), RoundUp)
    ν = B.k
    νf = Float64(B.k, RoundUp)
    den = n ⊗₋ (νf ⊖₋ 2.0) ⊗₋ π ⊗₋ reduce(⊗₋, [n - i for i = 2:ν-1])
    return (4.0 ⊗₊ (n + 1)) ⊘₊ den
end
function aux_normalized_projection_error(B::Chebyshev)
    n = Float64(length(B), RoundUp)
    ν = B.k
    νf = Float64(B.k, RoundUp)
    den = π ⊗₋ νf ⊗₋ n ⊗₋ reduce(⊗₋, [n - i for i = 1:ν-1])
    return 2.0 ⊘₊ den
end

"""
	strong_weak_bound(B::Chebyshev)

V.A. Markov estimate from GRADIMIR MILOVANOVIC EXTREMAL PROBLEMS AND 
INEQUALITIES OF MARKOV-BERNSTEIN TYPE FOR POLYNOMIALS
"""
# TODO: Check the indexes

function strong_weak_bound(B::Chebyshev)
    n = length(B) - 1
    k = B.k - 1
    # we want to estimate the norm of f^(k) by the C1 norm of f, so we use Markov estimate
    # for derivative k-1
    den = reduce(⊗₋, [Float64(2k - 1 - 2 * i, RoundDown) for i = 0:k-1])
    num = reduce(⊗₊, [Float64(n^2 - i^2, RoundDown) for i = 0:k-1])
    return num ⊘₊ den ⊕₊ 1.0 # the 1.0 is to take into account the L1 norm of f
end
aux_weak_bound(B::Chebyshev) = 1.0

"""
    weak_by_strong_and_aux_bound(B::Chebyshev)

Returns `(S₁, S₂)` such that `||f||_{C1} ≤ S₁·||f||_{W^{k,1}} + S₂·||f||_{L1}`.

For k ≥ 2: Sobolev embedding gives ||f||_∞ ≤ ||f||_{L1} + ||f'||_{L1} and
||f'||_∞ ≤ ||f'||_{L1} + ||f''||_{L1}, so ||f||_{C1} ≤ 2·||f||_{W^{k,1}}.

For k = 1: uses Markov inequality ||p'||_∞ ≤ 2n²·||p||_∞ for polynomials of
degree n on [0,1], giving ||f||_{C1} ≤ (1 + 2n²)·||f||_{W^{1,1}}.
"""
function weak_by_strong_and_aux_bound(B::Chebyshev)
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
function bound_weak_norm_from_linalg_norm(B::Chebyshev)
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
function bound_linalg_norm_L1_from_weak(B::Chebyshev)
    n = Float64(length(B), RoundUp)
    return 2.0 ⊗₊ n ⊖₋ 1.0
end

"""
    bound_linalg_norm_L∞_from_weak(B::Chebyshev)

Returns `A` such that `||ĉ||_{ℓ∞} ≤ A·||f||_{C1}`.

max_j |cⱼ| ≤ 2·||f||_∞ ≤ 2·||f||_{C1}.
"""
function bound_linalg_norm_L∞_from_weak(B::Chebyshev)
    return 2.0
end
weak_norm(B::Chebyshev) = C1
aux_norm(B::Chebyshev) = L1
strong_norm(B::Chebyshev) = W{B.k,1}

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


function chebtransform(w)
    n = length(w) - 1
    z = interval_fft([reverse(w); w[2:end-1]]) / n
    t = real.(z[1:length(w)])
    t[1] /= 2
    t[end] /= 2
    return Interval.(t)
end

using ProgressMeter
function assemble(B::Chebyshev, D::Dynamic; ϵ = 1e-13, max_iter = 100, T = Float64)
    n = length(B.p)
    M = zeros(Interval{T}, (n, n))
    x, labels, x′ = Dual(B, D; ϵ, max_iter)
    #@showprogress enabled=SHOW_PROGRESS_BARS  
    for i = 1:n
        ϕ = B[i]
        w = zeros(Interval{Float64}, n)
        for j = 1:length(x)
            w[labels[j]] += ϕ(x[j]) / abs(x′[j])
        end
        #@info w
        M[:, i] = chebtransform(w)
    end
    return M
end

using IntervalOptimisation

function infnormoffunction(B::Chebyshev, v)
    val = 0
    try
        val = maximize(x -> abs(evalChebyschevCentered(v, x)), interval(0, 1))[1]
    catch
        print("Refining grid")
        f(x) = abs(evalChebyshevCentered(v, x))
        ran = range_estimate(f, interval(0, 1), 5)
        Bval = union(val, ran)
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
        val = union(val, ran)
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

normbound(B::Chebyshev{T}, N::Type{C1}, v) where {T} =
    Float64((infnormoffunction(B, v) + infnormofderivative(B, v)).hi, RoundUp)

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
