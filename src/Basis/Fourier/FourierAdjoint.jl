
using IntervalArithmetic
using ..RigorousInvariantMeasures: MonotonicBranch, PwMap, Dual, C1, NormCacher, Cω, W, L2, NormKind
import ..RigorousInvariantMeasures: assemble
using LinearAlgebra

export FourierAdjoint

# we will store the Fourier expansion in the following way,
# to make it as coherent possible as the output of the FFT
# [0:N; -N:-1]

struct FourierAdjoint{S<:NormKind,WN<:NormKind,T<:AbstractVector} <: Fourier
    p::T
    k::Integer
    strong::S
    weak::WN
end

# Default constructor (backward-compatible): Cω strong, L2 weak
function FourierAdjoint(k::Integer, n::Integer; T = Float64)
    @warn "This basis breaks the usual interface of the package, i.e.,
    the dynamic is input as a function instead than a PwMap"
    return FourierAdjoint(FourierPoints(n, T), k, Cω(), L2())
end
function FourierAdjoint(p::AbstractVector, k::Integer)
    return FourierAdjoint(p, k, Cω(), L2())
end

Base.show(io::IO, B::FourierAdjoint) =
    print(io, "FourierAdjoint on $(length(B.p)) points, freq $(B.k), strong=$(B.strong), weak=$(B.weak)")

"""
Return the size of the Fourier basis
"""
Base.length(B::FourierAdjoint) = 2 * B.k + 1

###############################################################################
# Norm interface
###############################################################################

strong_norm(B::FourierAdjoint) = B.strong
weak_norm(B::FourierAdjoint) = typeof(B.weak)
aux_norm(B::FourierAdjoint) = L1

# Projection errors for Cω (same structure as Aη with a default strip width)
function weak_projection_error(B::FourierAdjoint{Cω})
    N = B.k
    # Conservative bound: use exponential decay with a small default η
    η = Interval(0.1)
    return exp(2 * Interval(pi) * N * η).hi
end

function aux_normalized_projection_error(B::FourierAdjoint{Cω})
    N = B.k
    η = Interval(0.1)
    return exp(2 * Interval(pi) * N * η).hi
end

function strong_weak_bound(B::FourierAdjoint{Cω})
    k = B.k
    η = Interval(0.1)
    λ = 4 * Interval(pi) * η
    geo_sum = 1 + 2 * exp(λ) * (1 - exp(λ * k)) / (1 - exp(λ))
    return sqrt(geo_sum).hi
end

###############################################################################
# Dual and assembly (unchanged)
###############################################################################

struct FourierAdjointDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
end

function FourierAdjointDualBranch(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter)

    mask = [br.X[1] <= x < br.X[2] for x in y]

    return [br.f(p) for p in y[mask]]
end

function Dual(B::FourierAdjoint, D::PwMap; ϵ, max_iter)
    results = collect(
        FourierAdjointDualBranch(B.p, b, 1:length(B.p); ϵ, max_iter) for b in branches(D)
    )

    @info typeof(results)

    x = vcat((result for result in results)...)

    return FourierAdjointDual(x)
end

function Dual(B::FourierAdjoint, T::Function; ϵ, max_iter)
    x = T.(B.p)
    return FourierAdjointDual(x)
end

function eval_on_dual(B::FourierAdjoint, computed_dual::FourierAdjointDual, ϕ)

    x = computed_dual.x

    return ϕ.(Interval.(x))
end

function assemble(B::FourierAdjoint, D; ϵ = 0.0, max_iter = 100, T = Float64)
    return assemble_common(B, D; ϵ = 0.0, max_iter = 100, T = Float64)'
end
