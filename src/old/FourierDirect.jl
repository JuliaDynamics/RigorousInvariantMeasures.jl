# COV_EXCL_START

using IntervalArithmetic
using ..RigorousInvariantMeasures: MonotonicBranch, PwMap, Dual, C1, NormCacher






using LinearAlgebra


struct Cω <: NormKind end
struct L2 <: NormKind end
#struct W{k, l} <:NormKind end
#order(::Type{W{k, l}}) where {k, l} = k
#regularity(::Type{W{k, l}}) where {k, l} = l

struct FourierAnalytic{T<:AbstractVector} <: Basis
    p::T
    k::Integer
end

FourierPoints(n, T) = [Interval{T}(i) / (n) for i = 0:n-1]

function FourierAnalytic(n::Integer, k::Integer; T = Float64)
    return Fourier(FourierPoints(n, T), k)
end
Base.show(io::IO, B::FourierAnalytic) = print(
    io,
    "Fourier basis on $(length(B)) points, highest frequency $(Int64(floor((length(B)-1)/2)))",
)

"""
Return the size of the Fourier basis
"""
Base.length(B::FourierAnalytic) = length(B.p)

##########################################################
# These are the methods to rigorously enclose the image of a trigonometric polynomial
# they need refactoring, but not now
#########################################################
"""
	_eval_T

Eval the Fourier basis up to degree n on an array of 
points in [0, 1].

"""
function _eval_Fourier_T(n, x::Array{T}) where {T}
    k = length(x)
    M = zeros(Complex, k, n + 1)
    M[:, 1] = ones(k)
    for i = 1:Int64(floor(n / 2))
        M[:, i+1] = exp.(i * 2im * pi * x)
        M[:, end-i+1] = exp.(-i * 2im * pi * x)
    end
    return M
end

function _eval_Fourier_T_derivative(n, x::Array{T}) where {T}
    k = length(x)
    M = zeros(Complex, k, n + 1)

    for i = 1:Int64(floor(n / 2))
        M[:, i+1] = i * 2im * pi * exp.(i * 2im * pi * x)
        M[:, end-i+1] = -i * 2im * pi * exp.(-i * 2im * pi * x)
    end
    return M
end

# ModFreq(i, j, x) = 2im * pi * (i * x[1] + j * x[2])
# """
# We Create an array M of size (k,m,n,l) M[:,:,i,j] is a matrix of basis vectors evaluated at (i,j), M[k,m,:,:] is the basis vector k,m
# """
# function _eval_Fourier_T_2D(n, l, x::Matrix{T}) where {T}
#     (k, m) = size(x)
#     M = zeros(Complex, k, m, n + 1, l + 1)
#     M[:, :, 1, 1] = ones(k, m)
#     for i in 1:Int64(floor(n / 2))
#         for j in 1:Int64(floor(m / 2))
#             M[:, :, i+1, j+1] = exp.(ModFreq.(i, j, x))
#             M[:, :, i+1, end-j+1] = exp.(ModFreq.(i, -j, x))
#             M[:, :, end-i+1, j+1] = exp.(ModFreq.(-i, j, x))
#             M[:, :, end-i+1, end-j+1] = exp.(ModFreq.(-i, -j, x))
#         end
#     end
#     return M
# end

evalFourierVector(coeff, B::FourierDirect) =
    real.(_eval_Fourier_T(length(coeff) - 1, B.p) * coeff)
evalFourierDerivativeVector(coeff, B::FourierDirect) =
    real.(_eval_Fourier_T_derivative(length(coeff) - 1, B.p) * coeff)

evalFourier(coeff, x) = mid(real.(_eval_Fourier_T(length(coeff) - 1, [x]) * coeff)[1])
evalFourierImag(coeff, x) = mid(imag.(_eval_Fourier_T(length(coeff) - 1, [x]) * coeff)[1])
evalFourierDerivative(coeff, x) =
    mid(real.(_eval_Fourier_T_derivative(length(coeff) - 1, [x]) * coeff)[1])

"""
Maybe it is better to state that you want the mid of your interval in coding rather than here
"""
evalFourier(coeff::Interval, x) = evalFourier(mid.(coeff), x)
evalFourierImag(coeff::Interval, x) = evalFourierImag(mid.(coeff), x)
evalFourierDerivative(coeff::Interval, x) = evalFourier(mid.(coeff), x)

function evalFourierCentered(coeff, x::Interval)
    m = Interval(mid.(x))
    return evalFourier(coeff, m) + evalFourierDerivative(coeff, x) * (x - m)
end

###############################################################################
###############################################################################

function weak_projection_error(B::Fourier)
    n = Float64(length(B), RoundUp)
    ν = B.k
    νf = Float64(B.k, RoundUp)
    den = n ⊗₋ (νf ⊖₋ 2.0) ⊗₋ π ⊗₋ reduce(⊗₋, [n - i for i = 2:ν-1])
    return (4.0 ⊗₊ (n + 1)) ⊘₊ den
end

function aux_normalized_projection_error(B::Fourier)
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

function strong_weak_bound(B::Fourier)
    n = length(B) - 1
    k = B.k - 1
    # we want to estimate the norm of f^(k) by the C1 norm of f, so we use Markov estimate
    # for derivative k-1
    den = reduce(⊗₋, [Float64(2k - 1 - 2 * i, RoundDown) for i = 0:k-1])
    num = reduce(⊗₊, [Float64(n^2 - i^2, RoundDown) for i = 0:k-1])
    return num ⊘₊ den ⊕₊ 1.0 # the 1.0 is to take into account the L1 norm of f
end
aux_weak_bound(B::Fourier) = 1.0

# Check this!!!
function weak_by_strong_and_aux_bound(B::Fourier)
    @error "TODO"
    ν = B.k
    return (ν, 1.0)
end
bound_weak_norm_from_linalg_norm(B::Fourier) = @error "TODO"
bound_linalg_norm_L1_from_weak(B::Fourier) = @error "TODO"
bound_linalg_norm_L∞_from_weak(B::Fourier) = @error "TODO"
weak_norm(B::Fourier) = C1
aux_norm(B::Fourier) = L1
strong_norm(B::Fourier) = W{B.k,1}

"""
Make so that B[j] returns the basis function of coordinate j
"""
function Base.getindex(B::Fourier, i::Int)
    n = length(B)
    v = zeros(n)
    v[i] = 1
    @boundscheck 1 <= i <= n || throw(BoundsError(B, i))
    return x -> evalFourier(v, x)
end

is_refinement(Bc::Fourier, Bf::Fourier) = length(Bc) < length(Bf)
integral_covector(B::Fourier; T = Float64) = [Interval{T}(1); zeros(length(B) - 1)]'
one_vector(B::Fourier) = [1; zeros(length(B) - 1)]

Base.length(S::AverageZero{T}) where {T<:Fourier} = length(S.basis) - 1

function Base.iterate(S::AverageZero{T}, state = 1) where {T<:Fourier}
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

struct FourierDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    xp::Vector{Interval}
end

using ..RigorousInvariantMeasures: preimages_and_derivatives

function FourierDualBranch(y, br::MonotonicBranch, ylabel = 1:length(y), ϵ = 0.0)
    if br.increasing
        endpoint_X = br.X[2]
        der = derivative(br.f)(endpoint_X)
        preim_der = preimages_and_derivatives(y, br, ylabel, ϵ)
        return [preim_der[1];],#endpoint_X is calculated in preimages
        [preim_der[2];],
        [preim_der[3];]#same as previous comment
    else
        endpoint_X = br.X[2]
        der = derivative(br.f)(endpoint_X)
        preim_der = preimages_and_derivatives(y, br, 1:length(y)-1, ϵ)
        return [preim_der[1]; endpoint_X],
        [preim_der[2]; length(preim_der[2]) + 1],
        [preim_der[3]; der]
    end
end

function Dual(B::Fourier, D::PwMap, ϵ)
    #@assert is_full_branch(D)
    results = collect(FourierDualBranch(B.p, b, 1:length(B.p), ϵ) for b in branches(D))
    x = vcat((result[1] for result in results)...)
    xlabel = vcat((result[2] for result in results)...)
    xp = vcat((result[3] for result in results)...)
    return x, xlabel, xp
end

Base.length(dual::FourierDual) = length(dual.x)
Base.eltype(dual::FourierDual) =
    Tuple{eltype(dual.xlabel),Tuple{eltype(dual.x),eltype(dual.x')}}
function Base.iterate(dual::FourierDual, state = 1)
    if state <= length(dual.x)
        return ((dual.xlabel[state], (dual.x[state], abs(dual.x'[state]))), state + 1)
    else
        return nothing
    end
end


function FourierTransform(w)
    #@info sum(w)
    n = length(w)
    z = fft(w) / (n)
    #t = real.(z[1:length(w)])
    return z#Interval.(t)
end

using ProgressMeter
function assemble_standard(B::Fourier, D::Dynamic, ϵ = 0.0; T = Float64)
    n = length(B.p)
    M = zeros(Complex{Interval{Float64}}, (n, n))
    x, labels, xp = Dual(B, D, ϵ)
    @showprogress for i = 1:n
        ϕ = B[i]
        w = zeros(Complex{Interval{Float64}}, n)
        for j = 1:length(x)
            C = ϕ(x[j]) / abs(xp[j])
            if real(C) != ∅ && imag(C) != ∅
                w[labels[j]] += C
            end
        end
        #@info w
        M[:, i] = FourierTransform(mid.(w))
    end
    return M
end

function assemble(B::Fourier, D::Dynamic, ϵ = 0.0; T = Float64)
    n = length(B)
    #n_odd = (n%2 == 1)
    #n_less = Int64(floor(n/2))
    #B_expanded = Fourier(B.k,Int64(2*B.k))
    #M_big = assemble_standard(B_expanded,D)
    #M_upper = hcat(M_big[1:n_odd+n_less,1:n_odd+n_less],M_big[1:n_odd+n_less,end-n_less+1:end])
    #M_lower = hcat(M_big[end-n_less+1:end,1:n_odd+n_less],M_big[end-n_less+1:end,end-n_less+1:end])
    #return vcat(M_upper,M_lower)
    return assemble_standard(B, D)
end

function ConvolutionMatrix(C, ϵ = 0.0; T = Float64)
    c1 = fft(C)
    return diagm(c1 / c1[1])
end

function FourierCoefficientsError(B, f::Vector{Complex{Interval{T}}}) where {T}
    u = Float64(eps(T), RoundUp)
    gamma4 = 4 * u / (1 - 4 * u)
    mu = u
    eta = mu + gamma4 * (sqrt(2) + mu)
    N = length(B)
    fm2 = sum((mid.(abs.(f))) .^ 2)^0.5
    fr2 = sum((radius.(abs.(f))) .^ 2)^0.5
    e = log2(N) / sqrt(2 * N + 1) * (eta / (1 - eta) * fm2 + fr2)
    return e
end

function ProjectorError(B)
    V = 1#νth variation
    ν = 6#Number of derivatives
    n = length(B)
    println(n)
    n_prod1 = 1
    n_prod2 = 1
    for i = 0:ν-1
        n_prod1 *= n - i
        n_prod2 *= n - i
        println(n_prod1)
        println(n_prod2)
    end

    n_prod2 /= (n - 1)
    #||πLπ-L||≤||πLπ-Lπ||+||Lπ-L||≤ϵ||Lπ||+a||π-1||≤ϵap+aϵ=ϵa(1+p)
    return 4 * V / (π * ν * n_prod1) + 4 * (n + 1) * V / (π * (ν - 2) * n_prod2)
end

function OperatorError(B)
    p = 1#Projection operator norm
    a = 1.5#Operator operator norm
    ϵ = ProjectorError(B)
    return ϵ * a * (1 + p)
end

using IntervalOptimisation

function infnormoffunction(B::Fourier, v)
    val = 0
    try
        val = maximize(x -> abs(evalFourierCentered(v, x)), Interval(0, 1))[1]
    catch
        print("Refining grid")
        f(x) = abs(evalFourierCentered(v, x))
        ran = range_estimate(f, Interval(0, 1), 5)
        Bval = union(val, ran)
    end
    return val
end

function infnormofderivative(B::Fourier, v)
    val = Interval(0)
    try
        val = maximize(x -> abs(evalFourierDerivative(v, interval(x))), Interval(0, 1))[1]
    catch
        print("Refining grid")
        f(x) = abs(evalFourierDerivative(v, x))
        ran = range_estimate(f, Interval(0, 1), 5)
        Bval = union(val, ran)
    end
    return val
end


is_integral_preserving(B::Fourier) = false
function opnormbound(B::Fourier, N::Type{C1}, v::Vector{S}) where {T,S}
    return normbound(B, N, v)
end

function opnormbound(B::Fourier, N::Type{C1}, w::LinearAlgebra.Adjoint) where {T,S}
    return normbound(B, N, w')
end

function opnormbound(B::Fourier, N::Type{C1}, A::Matrix{S}) where {T,S}
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

normbound(B::Fourier{T}, N::Type{C1}, v) where {T} =
    Float64((infnormoffunction(B, v) + infnormofderivative(B, v)).hi, RoundUp)

# mutable struct NormCacherC1 <: NormCacher{C1}
# 	B::Basis
#     C::Float64
#     function NormCacherC1(B, n)
#         new(B, 0.0)
#     end
# end
# NormCacher{C1}(B, n) = NormCacherC1(B, n)

# function add_column!(Cacher::NormCacherC1, v::AbstractVector, ε::Float64)
#     Cacher.C = max(Cacher.C, opnormbound(Cacher.B, C1, v) ⊕₊ ε)
# end

# """
# Return the norm of the matrix the NormCacher is working on.
# """
# function get_norm(Cacher::NormCacherC1)
#     n = length(Cacher.B)
# 	return Cacher.C ⊗₊ log(n+2)
# end

# COV_EXCL_STOP
