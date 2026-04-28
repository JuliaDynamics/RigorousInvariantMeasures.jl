using IntervalArithmetic
using ..RigorousInvariantMeasures: MonotonicBranch, PwMap, Dual, C1, NormCacher, Aη, W, L2, NormKind
import ..RigorousInvariantMeasures: derivative
using LinearAlgebra

export FourierAnalytic

# we will store the Fourier expansion in the following way,
# to make it as coherent possible as the output of the FFT
# [0:N; -N:-1]

struct FourierAnalytic{S<:NormKind,W<:NormKind,T<:AbstractVector} <: Fourier
    p::T
    k::Integer
    strong::S
    weak::W
end

Base.length(B::FourierAnalytic) = 2 * B.k + 1

# Default constructors (backward-compatible): Aη strong, L2 weak
function FourierAnalytic(k::Integer, n::Integer; η = 0.1, T = Float64)
    return FourierAnalytic(FourierPoints(n, T), k, Aη(η), L2())
end
function FourierAnalytic(p::AbstractVector, k::Integer)
    return FourierAnalytic(p, k, Aη(0.1), L2())
end

# Sobolev constructor
function FourierAnalytic(k_freq::Integer, n::Integer, ::Type{W{j,l}}; T = Float64) where {j,l}
    return FourierAnalytic(FourierPoints(n, T), k_freq, W{j,l}(), L2())
end

Base.show(io::IO, B::FourierAnalytic) =
    print(io, "FourierAnalytic on $(length(B.p)) points, freq $(B.k), strong=$(B.strong), weak=$(B.weak)")

###############################################################################
# Norm interface
###############################################################################

strong_norm(B::FourierAnalytic) = B.strong
weak_norm(B::FourierAnalytic) = typeof(B.weak)
aux_norm(B::FourierAnalytic) = L1

# --- Projection errors ---

# Aη strong norm: exponential decay of Fourier tail
function weak_projection_error(B::FourierAnalytic{Aη})
    N = B.k
    η = interval(B.strong.η)
    return sup(exp(2 * interval(pi) * N * η))
end

function aux_normalized_projection_error(B::FourierAnalytic{Aη})
    N = B.k
    η = interval(B.strong.η)
    return sup(exp(2 * interval(pi) * N * η))
end

# W{k,1} strong norm: polynomial decay O(N^{-(k-1)}) for Fourier coefficient decay
# ||P_h f - f||_{L²} ≤ (2πN)^{-k} ||f^{(k)}||_{L²} ≤ (2πN)^{-k} ||f||_{W^{k,1}}
function weak_projection_error(B::FourierAnalytic{W{k,l}}) where {k,l}
    N = B.k
    if k == 0
        return 1.0
    end
    denom = (2 * interval(pi) * N)^k
    return sup(1 / denom)
end

function aux_normalized_projection_error(B::FourierAnalytic{W{k,l}}) where {k,l}
    N = B.k
    if k == 0
        return 1.0
    end
    denom = (2 * interval(pi) * N)^k
    return sup(1 / denom)
end

# --- Strong-weak bound: ||v||_s ≤ M₁n · ||v||_{L²} ---

@doc raw"""
Return the weak-strong norm bound when restricted on the finite dimensional subspace.
For Aη, by Cauchy-Schwarz:
``||v||_{Aη} = \sum |ĉ_k| e^{2πη|k|} \leq (\sum e^{4πη|k|})^{1/2} \cdot ||ĉ||_{ℓ²}``
Geometric series: ``\sum_{|k|\leq N} e^{4πη|k|} = 1 + 2(e^{4πη} - e^{4πη(N+1)})/(1 - e^{4πη})``
"""
function strong_weak_bound(B::FourierAnalytic{Aη})
    k = B.k
    η = interval(B.strong.η)
    λ = 4 * interval(pi) * η

    # Sum of geometric series: 1 + 2*(e^λ - e^{λ(N+1)})/(1 - e^λ)
    # = 1 + 2*e^λ*(1 - e^{λN})/(1 - e^λ)
    geo_sum = 1 + 2 * exp(λ) * (1 - exp(λ * k)) / (1 - exp(λ))
    return sup(sqrt(geo_sum))
end

@doc raw"""
For W^{k,1}: Bernstein + Cauchy-Schwarz on [0,1]:
``||v^{(j)}||_{L¹} ≤ ||v^{(j)}||_{L²} ≤ (2πN)^j ||v||_{L²}``
Sum: ``M₁n = \sum_{j=0}^{k} (2πN)^j``
"""
function strong_weak_bound(B::FourierAnalytic{W{k,l}}) where {k,l}
    N = B.k
    M = interval(1.0)
    for j = 1:k
        M = M + (2 * interval(pi) * N)^j
    end
    return sup(M)
end

###############################################################################
# Dual and assembly (unchanged)
###############################################################################

struct FourierAnalyticDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    xp::Vector{Interval}
end

using ..RigorousInvariantMeasures: preimages_and_derivatives

function FourierAnalyticDualBranch(
    y,
    br::MonotonicBranch,
    ylabel = 1:length(y);
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

function eval_on_dual(B::FourierAnalytic, computed_dual::FourierAnalyticDual, ϕ)

    x, labels, xp = computed_dual.x, computed_dual.xlabel, computed_dual.xp

    w = zeros(Complex{Interval{Float64}}, length(B.p))
    for j = 1:length(x)
        C = ϕ(x[j]) / abs(xp[j])
        if !isempty_interval(real(C)) && !isempty_interval(imag(C))
            w[labels[j]] += C
        end
    end

    return w

end

# assemble(::FourierAnalytic, D; …) is provided by the FFTWExt extension.
