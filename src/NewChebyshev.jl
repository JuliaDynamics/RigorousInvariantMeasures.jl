using ..BasisDefinition, ..DynamicDefinition


struct C1 <: NormKind end
struct W{k,l} <: NormKind end
order(::Type{W{k,l}}) where {k,l} = k
regularity(::Type{W{k,l}}) where {k,l} = l

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
    return u[1] + Interval(-γ, γ)
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
    return u[1] + Interval(-γ, γ)
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
    eval_Clenshaw_BackwardFirst(coeff, Interval(mid.(2 * x - 1)))
evalChebyshevDerivative(coeff, x::Interval) = 2 * ChebyshevDerivative(coeff, 2 * x - 1)
function evalChebyschevCentered(coeff, x::Interval)
    m = Interval(mid.(x))
    return evalChebyshev(coeff, m) + evalChebyshevDerivative(coeff, x) * (x - m)
end
###############################################################################
###############################################################################

function BasisDefinition.weak_projection_error(B::Chebyshev)
    n = Float64(length(B), RoundUp)
    ν = B.k
    νf = Float64(B.k, RoundUp)
    den = n ⊗₋ (νf ⊖₋ 2.0) ⊗₋ π ⊗₋ reduce(⊗₋, [n - i for i = 2:ν-1])
    return (4.0 ⊗₊ (n + 1)) ⊘₊ den
end
function BasisDefinition.aux_normalized_projection_error(B::Chebyshev)
    n = Float64(length(B), RoundUp)
    ν = B.k
    νf = Float64(B.k, RoundUp)
    den = π ⊗₋ νf ⊗₋ n ⊗₋ reduce(⊗₋, [n - i for i = 1:ν-1])
    return 2.0 ⊘₊ den
end

"""
	BasisDefinition.strong_weak_bound(B::Chebyshev)

V.A. Markov estimate from GRADIMIR MILOVANOVIC EXTREMAL PROBLEMS AND 
INEQUALITIES OF MARKOV-BERNSTEIN TYPE FOR POLYNOMIALS
"""
# TODO: Check the indexes

function BasisDefinition.strong_weak_bound(B::Chebyshev)
    n = length(B) - 1
    k = B.k - 1
    # we want to estimate the norm of f^(k) by the C1 norm of f, so we use Markov estimate
    # for derivative k-1
    den = reduce(⊗₋, [Float64(2k - 1 - 2 * i, RoundDown) for i = 0:k-1])
    num = reduce(⊗₊, [Float64(n^2 - i^2, RoundDown) for i = 0:k-1])
    return num ⊘₊ den ⊕₊ 1.0 # the 1.0 is to take into account the L1 norm of f
end
BasisDefinition.aux_weak_bound(B::Chebyshev) = 1.0

# Check this!!!
function BasisDefinition.weak_by_strong_and_aux_bound(B::Chebyshev)
    @error "TODO"
    ν = B.k
    return (ν, 1.0)
end
BasisDefinition.bound_weak_norm_from_linalg_norm(B::Chebyshev) = @error "TODO"
BasisDefinition.bound_linalg_norm_L1_from_weak(B::Chebyshev) = @error "TODO"
BasisDefinition.bound_linalg_norm_L∞_from_weak(B::Chebyshev) = @error "TODO"
BasisDefinition.weak_norm(B::Chebyshev) = C1
BasisDefinition.aux_norm(B::Chebyshev) = L1
BasisDefinition.strong_norm(B::Chebyshev) = W{B.k,1}

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

BasisDefinition.is_refinement(Bc::Chebyshev, Bf::Chebyshev) = length(Bc) < length(Bf)
BasisDefinition.integral_covector(B::Chebyshev; T = Float64) =
    [Interval{T}(1); 0; [0.5 * Interval{T}((-1)^n + 1) / (1 - n^2) for n = 2:length(B)-1]]'
BasisDefinition.one_vector(B::Chebyshev) = [1; zeros(length(B) - 1)]

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
    return v / ((i + 1)^2), state + 1
end

struct ChebyshevDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    x′::Vector{Interval}
end

function ChebDualBranch(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter)
    if br.increasing
        endpoint_X = br.X[2]
        der = Contractors.derivative(br.f)(endpoint_X)
        preim_der = preimages_and_derivatives(y, br, ylabel; ϵ, max_iter)
        return [preim_der[1]; endpoint_X],
        [preim_der[2]; length(preim_der[2]) + 1],
        [preim_der[3]; der]
    else
        endpoint_X = br.X[2]
        der = Contractors.derivative(br.f)(endpoint_X)
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
    @showprogress for i = 1:n
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
        val = maximize(x -> abs(evalChebyschevCentered(v, x)), Interval(0, 1))[1]
    catch
        print("Refining grid")
        f(x) = abs(evalChebyshevCentered(v, x))
        ran = range_estimate(f, Interval(0, 1), 5)
        Bval = union(val, ran)
    end
    return val
end

function infnormofderivative(B::Chebyshev, v)
    val = Interval(0)
    try
        val = maximize(x -> abs(evalChebyshevDerivative(v, x)), Interval(0, 1))[1]
    catch
        print("Refining grid")
        f(x) = abs(evalChebyshevDerivative(v, x))
        ran = range_estimate(f, Interval(0, 1), 5)
        val = union(val, ran)
    end
    return val
end


BasisDefinition.is_integral_preserving(B::Chebyshev) = false
function BasisDefinition.opnormbound(B::Chebyshev, N::Type{C1}, v::Vector{S}) where {S}
    return normbound(B, N, v)
end

function BasisDefinition.opnormbound(B::Chebyshev, N::Type{C1}, w::LinearAlgebra.Adjoint)
    return normbound(B, N, w')
end

function BasisDefinition.opnormbound(B::Chebyshev, N::Type{C1}, A::Matrix{S}) where {S}
    n, m = size(A)
    norm = 0.0
    for i = 1:m
        norm = max(
            norm,
            Float64(opnormbound(B, N, A[:, i]), RoundUp) ⊘₊ Float64(1 + i^2, RoundDown),
        )
        # since these are the images of the basis elements, we can use the fact that
        # |Tₖ(x)|= 1, Tₖ' = k*U_{k-1}, |U_{k-1}|=k
    end
    return norm ⊗₊ log(m + 2)
end

BasisDefinition.normbound(B::Chebyshev{T}, N::Type{C1}, v) where {T} =
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
    return Cacher.C ⊗₊ log(n + 2)
end
