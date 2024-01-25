using IntervalArithmetic
using ..RigorousInvariantMeasures: MonotonicBranch, PwMap, Dual, C1, NormCacher
import ..RigorousInvariantMeasures: NormKind, derivative, interval_fft, assemble
using LinearAlgebra

export FourierAnalytic

struct Aη <: NormKind
    η::Any
end
#struct L2 <: NormKind end

# we will store the Fourier expansion in the following way,
# to make it as coherent possible as the output of the FFT
# [0:N; -N:-1]

struct FourierAnalytic{T<:AbstractVector} <: Fourier
    p::T
    k::Integer
end

Base.length(B::FourierAnalytic) = 2 * B.k + 1


function FourierAnalytic(k::Integer, n::Integer; T=Float64)
    return FourierAnalytic(FourierPoints(n, T), k)
end
Base.show(io::IO, B::FourierAnalytic) =
    print(io, "FFT on $(length(B.p)) points restricted to highest frequency $(B.k)")



###############################################################################
###############################################################################

function weak_projection_error(B::FourierAnalytic)
    N = B.k
    η = Interval(strong_norm(B).η)
    λ = -2 * Interval(pi) * N * η

    return exp(-λ).hi
end

function aux_normalized_projection_error(B::FourierAnalytic)
    N = B.k
    η = Interval(strong_norm(B).η)
    λ = -2 * Interval(pi) * N * η

    return exp(-λ).hi
end

@doc raw"""
Return the weak-strong norm bound when restricted on the finite dimensional subspace
We use here
`` \sum_{i=-K}^K (|g_i|\exp^{2\pi \eta |i|})^2 \leq \sum_{-K}^K |g_i|^2 \sum_{-K}^K \exp^{2(2\pi \eta |i|)}
\leq ||g||_2 2\frac{1-\exp^{2K+2(2\pi \eta |i|)}}{1-\exp^{2(2\pi \eta |i|)}}`` 
"""
function strong_weak_bound(B::FourierAnalytic)
    k = B.k
    η = Interval(strong_norm(B).η)
    λ = 2 * 2 * Interval(pi) * η

    return 2 * (1 - exp(λ * (k + 1)) / (1 - exp(λ)))
end
aux_weak_bound(B::FourierAnalytic) = 1.0

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

struct FourierAnalyticDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    xp::Vector{Interval}
end

using ..RigorousInvariantMeasures: preimages_and_derivatives

function FourierAnalyticDualBranch(
    y,
    br::MonotonicBranch,
    ylabel=1:length(y);
    ϵ,
    max_iter,
)
    if br.increasing
        endpoint_X = br.X[2]
        der = derivative(br.f, endpoint_X)
        preim_der = preimages_and_derivatives(y, br, ylabel; ϵ, max_iter)
        return [preim_der[1];],#endpoint_X is calculated in preimages
        [preim_der[2];],
        [preim_der[3];]#same as previous comment
    else
        endpoint_X = br.X[2]
        der = derivative(br.f, endpoint_X)
        preim_der = preimages_and_derivatives(y, br, 1:length(y)-1; ϵ, max_iter)
        return [preim_der[1]; endpoint_X],
        [preim_der[2]; length(preim_der[2]) + 1],
        [preim_der[3]; der]
    end
end

function Dual(B::FourierAnalytic, D::PwMap; ϵ, max_iter)
    #@assert is_full_branch(D)
    results = collect(
        FourierAnalyticDualBranch(B.p, b, 1:length(B.p); ϵ, max_iter) for b in branches(D)
    )
    x = vcat((result[1] for result in results)...)
    xlabel = vcat((result[2] for result in results)...)
    xp = vcat((result[3] for result in results)...)
    return FourierAnalyticDual(x, xlabel, xp)
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

function eval_on_dual(B::FourierAnalytic, computed_dual::FourierAnalyticDual, ϕ)

    x, labels, xp = computed_dual.x, computed_dual.xlabel, computed_dual.xp

    #@info x, maximum(diam.(x))
    w = zeros(Complex{Interval{Float64}}, length(B.p))
    for j = 1:length(x)
        C = ϕ(x[j]) / abs(xp[j])
        if real(C) != ∅ && imag(C) != ∅
            w[labels[j]] += C
        end
    end

    return w

end



using ProgressMeter
function assemble(B::FourierAnalytic, D::Dynamic; ϵ=0.0, max_iter=100, T=Float64)
    return assemble_common(B, D; ϵ, max_iter, T)
    # n = length(B)

    # k = (n - 1) ÷ 2

    # M = zeros(Complex{Interval{Float64}}, (n, n))
    # computed_dual = Dual(B, D; ϵ, max_iter)
    # @showprogress for i = 1:n
    #     ϕ = B[i]
    #     w = eval_on_dual(B, computed_dual, ϕ)

    #     # zeros(Complex{Interval{Float64}}, length(B.p))
    #     # for j in 1:length(x)
    #     #     C = ϕ(x[j]) / abs(xp[j])
    #     #     if real(C) != ∅ && imag(C) != ∅
    #     #         w[labels[j]] += C
    #     #     end
    #     # end
    #     # #@info w
    #     FFTw = interval_fft(w)

    #     M[:, i] = [FFTw[1:k+1]; FFTw[end-k+1:end]]
    # end
    # return M
end

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
