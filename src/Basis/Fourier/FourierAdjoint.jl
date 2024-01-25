
using IntervalArithmetic
using ..RigorousInvariantMeasures: MonotonicBranch, PwMap, Dual, C1, NormCacher
import ..RigorousInvariantMeasures: NormKind, assemble
using LinearAlgebra

export FourierAdjoint

struct Cω <: NormKind end

# we will store the Fourier expansion in the following way,
# to make it as coherent possible as the output of the FFT
# [0:N; -N:-1]

struct FourierAdjoint{T<:AbstractVector} <: Fourier
    p::T
    k::Integer
end

function FourierAdjoint(k::Integer, n::Integer; T = Float64)
    @warn "This basis breaks the usual interface of the package, i.e., 
    the dynamic is input as a function instead than a PwMap"
    return FourierAdjoint(FourierPoints(n, T), k)
end
Base.show(io::IO, B::FourierAdjoint) =
    print(io, "FFT on $(length(B.p)) points restricted to highest frequency $(B.k)")

"""
Return the size of the Fourier basis
"""
Base.length(B::FourierAdjoint) = 2 * B.k + 1
#Base.lastindex(B::FourierAdjoint) = length(B)

###############################################################################
###############################################################################

# function weak_projection_error(B::Fourier)
#     n = Float64(length(B), RoundUp)
#     ν = B.k
#     νf = Float64(B.k, RoundUp)
#     den = n ⊗₋ (νf ⊖₋ 2.0) ⊗₋ π ⊗₋ reduce(⊗₋, [n - i for i in 2:ν-1])
#     return (4.0 ⊗₊ (n + 1)) ⊘₊ den
# end

# function aux_normalized_projection_error(B::Fourier)
#     n = Float64(length(B), RoundUp)
#     ν = B.k
#     νf = Float64(B.k, RoundUp)
#     den = π ⊗₋ νf ⊗₋ n ⊗₋ reduce(⊗₋, [n - i for i in 1:ν-1])
#     return 2.0 ⊘₊ den
# end

# """
# 	strong_weak_bound(B::Chebyshev)

# V.A. Markov estimate from GRADIMIR MILOVANOVIC EXTREMAL PROBLEMS AND 
# INEQUALITIES OF MARKOV-BERNSTEIN TYPE FOR POLYNOMIALS
# """
# TODO: Check the indexes

# function strong_weak_bound(B::Fourier)
#     n = length(B) - 1
#     k = B.k - 1
#     # we want to estimate the norm of f^(k) by the C1 norm of f, so we use Markov estimate
#     # for derivative k-1
#     den = reduce(⊗₋, [Float64(2k - 1 - 2 * i, RoundDown) for i in 0:k-1])
#     num = reduce(⊗₊, [Float64(n^2 - i^2, RoundDown) for i in 0:k-1])
#     return num ⊘₊ den ⊕₊ 1.0 # the 1.0 is to take into account the L1 norm of f
# end
# aux_weak_bound(B::Fourier) = 1.0

# Check this!!!
# function weak_by_strong_and_aux_bound(B::Fourier)
#     @error "TODO"
#     ν = B.k
#     return (ν, 1.0)
# end
# bound_weak_norm_from_linalg_norm(B::Fourier) = @error "TODO"
# bound_linalg_norm_L1_from_weak(B::Fourier) = @error "TODO"
# bound_linalg_norm_L∞_from_weak(B::Fourier) = @error "TODO"
# weak_norm(B::Fourier) = C1
# aux_norm(B::Fourier) = L1
# strong_norm(B::Fourier) = W{B.k,1}

# """
# Make so that B[j] returns the basis function of coordinate j
# """
# function Base.getindex(B::FourierAdjoint, i::Int)

#     K = length(B)
#     @boundscheck 1 <= i <= K || throw(BoundsError(B, i))

#     N = (K - 1) ÷ 2

#     if i == 1
#         return x -> Interval(1.0)
#     end

#     if 1 < i <= N + 1
#         return x -> exp(2 * pi * im * (i - 1) * x)
#     end

#     if N + 2 <= i <= K
#         L = i + (-2 * N - 2)
#         return x -> exp(2 * pi * im * L * x)
#     end
# end


# is_refinement(Bc::FourierAdjoint, Bf::FourierAdjoint) = length(Bc) < length(Bf)
# integral_covector(B::FourierAdjoint; T=Float64) = [Interval{T}(1); zeros(length(B) - 1)]'
# one_vector(B::FourierAdjoint) = [1; zeros(length(B) - 1)]

# Base.length(S::AverageZero{T}) where {T<:FourierAdjoint} = length(S.basis) - 1

# function Base.iterate(S::AverageZero{T}, state=1) where {T<:FourierAdjoint}
#     B = S.basis
#     i = state
#     if i == length(B)
#         return nothing
#     end
#     v = zeros(length(B))
#     v[i+1] = 1
#     return v, state + 1
# end

struct FourierAdjointDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
end

function FourierAdjointDualBranch(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter)

    mask = [br.X[1] <= x < br.X[2] for x in y]

    return [br.f(p) for p in y[mask]]
end

function Dual(B::FourierAdjoint, D::PwMap; ϵ, max_iter)
    #@assert is_full_branch(D)

    #T = eltype(B.p)
    #x = []

    results = collect(
        FourierAdjointDualBranch(B.p, b, 1:length(B.p); ϵ, max_iter) for b in branches(D)
    )

    @info typeof(results)

    #@info results

    x = vcat((result for result in results)...)

    #for b in branches(D)
    #    x = vcat(x, )
    #end

    #@info typeof(x)
    return FourierAdjointDual(x)
end

function Dual(B::FourierAdjoint, T::Function; ϵ, max_iter)
    #@assert is_full_branch(D)

    #T = eltype(B.p)
    #x = []

    #results = collect(FourierAdjointDualBranch(B.p, b, 1:length(B.p); ϵ, max_iter) for b in branches(D))

    #@info typeof(results)

    #@info results

    x = T.(B.p)

    #for b in branches(D)
    #    x = vcat(x, )
    #end

    #@info typeof(x)
    return FourierAdjointDual(x)
end


# Base.length(dual::FourierDual) = length(dual.x)
# Base.eltype(dual::FourierDual) = Tuple{eltype(dual.xlabel),Tuple{eltype(dual.x),eltype(dual.x')}}
# function Base.iterate(dual::FourierDual, state=1)
#     if state <= length(dual.x)
#         return ((dual.xlabel[state], (dual.x[state], abs(dual.x'[state]))), state + 1)
#     else
#         return nothing
#     end
# end


# function FourierTransform(w)
#     #@info sum(w)
#     n = length(w)
#     z = fft(w) / (n)
#     #t = real.(z[1:length(w)])
#     return z#Interval.(t)
# end

function eval_on_dual(B::FourierAdjoint, computed_dual::FourierAdjointDual, ϕ)

    x = computed_dual.x

    #@info x, maximum(diam.(x))

    return ϕ.(Interval.(x))
end

function assemble(B::FourierAdjoint, D; ϵ = 0.0, max_iter = 100, T = Float64)
    return assemble_common(B, D; ϵ = 0.0, max_iter = 100, T = Float64)'
end

# using ProgressMeter
# function assemble_standard(B::FourierAdjoint, D::Dynamic; ϵ=0.0, max_iter=100, T=Float64)
#     n = length(B)

#     @info n

#     k = (n - 1) ÷ 2

#     @info k

#     M = zeros(Complex{Interval{Float64}}, (n, n))
#     computed_dual = Dual(B, D; ϵ, max_iter)
#     @showprogress for i in 1:n
#         ϕ = B[i]
#         w = eval_on_dual(B, computed_dual, ϕ)
#         #@info w

#         FFTw = interval_fft(w)

#         M[:, i] = [FFTw[1:k+1]; FFTw[end-k+1:end]]
#     end
#     return M'
# end

# function assemble(B::Fourier, D::Dynamic, ϵ=0.0; T=Float64)
#     n = length(B)
#     #n_odd = (n%2 == 1)
#     #n_less = Int64(floor(n/2))
#     #B_expanded = Fourier(B.k,Int64(2*B.k))
#     #M_big = assemble_standard(B_expanded,D)
#     #M_upper = hcat(M_big[1:n_odd+n_less,1:n_odd+n_less],M_big[1:n_odd+n_less,end-n_less+1:end])
#     #M_lower = hcat(M_big[end-n_less+1:end,1:n_odd+n_less],M_big[end-n_less+1:end,end-n_less+1:end])
#     #return vcat(M_upper,M_lower)
#     return assemble_standard(B, D)
# end

# function ConvolutionMatrix(C, ϵ=0.0; T=Float64)
#     c1 = fft(C)
#     return diagm(c1 / c1[1])
# end

# function FourierCoefficientsError(B, f::Vector{Complex{Interval{T}}}) where {T}
#     u = Float64(eps(T), RoundUp)
#     gamma4 = 4 * u / (1 - 4 * u)
#     mu = u
#     eta = mu + gamma4 * (sqrt(2) + mu)
#     N = length(B)
#     fm2 = sum((mid.(abs.(f))) .^ 2)^0.5
#     fr2 = sum((radius.(abs.(f))) .^ 2)^0.5
#     e = log2(N) / sqrt(2 * N + 1) * (eta / (1 - eta) * fm2 + fr2)
#     return e
# end

# function ProjectorError(B)
#     V = 1#νth variation
#     ν = 6#Number of derivatives
#     n = length(B)
#     println(n)
#     n_prod1 = 1
#     n_prod2 = 1
#     for i = 0:ν-1
#         n_prod1 *= n - i
#         n_prod2 *= n - i
#         println(n_prod1)
#         println(n_prod2)
#     end

#     n_prod2 /= (n - 1)
#     #||πLπ-L||≤||πLπ-Lπ||+||Lπ-L||≤ϵ||Lπ||+a||π-1||≤ϵap+aϵ=ϵa(1+p)
#     return 4 * V / (π * ν * n_prod1) + 4 * (n + 1) * V / (π * (ν - 2) * n_prod2)
# end

# function OperatorError(B)
#     p = 1#Projection operator norm
#     a = 1.5#Operator operator norm
#     ϵ = ProjectorError(B)
#     return ϵ * a * (1 + p)
# end

# using IntervalOptimisation

# function infnormoffunction(B::Fourier, v)
#     val = 0
#     try
#         val = maximize(x -> abs(evalFourierCentered(v, x)), Interval(0, 1))[1]
#     catch
#         print("Refining grid")
#         f(x) = abs(evalFourierCentered(v, x))
#         ran = range_estimate(f, Interval(0, 1), 5)
#         Bval = union(val, ran)
#     end
#     return val
# end

# function infnormofderivative(B::Fourier, v)
#     val = Interval(0)
#     try
#         val = maximize(x -> abs(evalFourierDerivative(v, interval(x))), Interval(0, 1))[1]
#     catch
#         print("Refining grid")
#         f(x) = abs(evalFourierDerivative(v, x))
#         ran = range_estimate(f, Interval(0, 1), 5)
#         Bval = union(val, ran)
#     end
#     return val
# end


# is_integral_preserving(B::Fourier) = false
# function opnormbound(B::Fourier, N::Type{C1}, v::Vector{S}) where {T,S}
#     return normbound(B, N, v)
# end

# function opnormbound(B::Fourier, N::Type{C1}, w::LinearAlgebra.Adjoint) where {T,S}
#     return normbound(B, N, w')
# end

# function opnormbound(B::Fourier, N::Type{C1}, A::Matrix{S}) where {T,S}
#     n, m = size(A)
#     norm = 0.0
#     for i in 1:m
#         norm = max(norm, Float64(opnormbound(B, N, A[:, i]), RoundUp) ⊘₊ Float64(1 + i^2, RoundDown))
#         # since these are the images of the basis elements, we can use the fact that
#         # |Tₖ(x)|= 1, Tₖ' = k*U_{k-1}, |U_{k-1}|=k
#     end
#     return norm ⊗₊ log(m + 2)
# end

# normbound(B::Fourier{T}, N::Type{C1}, v) where {T} = Float64((infnormoffunction(B, v) + infnormofderivative(B, v)).hi, RoundUp)

# # mutable struct NormCacherC1 <: NormCacher{C1}
# # 	B::Basis
# #     C::Float64
# #     function NormCacherC1(B, n)
# #         new(B, 0.0)
# #     end
# # end
# # NormCacher{C1}(B, n) = NormCacherC1(B, n)

# # function add_column!(Cacher::NormCacherC1, v::AbstractVector, ε::Float64)
# #     Cacher.C = max(Cacher.C, opnormbound(Cacher.B, C1, v) ⊕₊ ε)
# # end

# # """
# # Return the norm of the matrix the NormCacher is working on.
# # """
# # function get_norm(Cacher::NormCacherC1)
# #     n = length(Cacher.B)
# # 	return Cacher.C ⊗₊ log(n+2)
# # end
