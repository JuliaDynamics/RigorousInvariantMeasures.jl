
using IntervalArithmetic
using ..RigorousInvariantMeasures: MonotonicBranch, PwMap, Dual, C1, NormCacher, Cω, W, L2, NormKind, WeightedSobolev
import ..RigorousInvariantMeasures: assemble
using LinearAlgebra

export FourierAdjoint

# we will store the Fourier expansion in the following way,
# to make it as coherent possible as the output of the FFT
# [0:N; -N:-1]

@doc raw"""
    FourierAdjoint(k, n; T = Float64)
    FourierAdjoint(p, k)

The adjoint realisation of the Fourier basis of `2k + 1` frequencies on `n` points, with
strong norm `Cω` by default.  With a Sobolev strong norm (`W{s,l}`, or the weighted norm of
`FourierAdjoint(k, n, D::CircleMap, s)`) the order must be
``s \ge 2``: W{1,1} Fourier bases are not supported, since no projection error is derived
for order 1. Use `Ulam` for ``W^{1,1}``/BV, or a Fourier basis with Sobolev order ``\ge 2``
and `L2` weak norm.
"""
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

# FourierAdjoint(k_freq, n, a::WeightedSobolev) and FourierAdjoint(k_freq, n, D::CircleMap, s):
# see the docstring of the FourierAnalytic methods of the same signatures
function FourierAdjoint(k_freq::Integer, n::Integer, a::WeightedSobolev; T = Float64)
    _fourier_weighted_order(a)
    return FourierAdjoint(FourierPoints(n, T), k_freq, a, L2())
end
function FourierAdjoint(
    k_freq::Integer,
    n::Integer,
    D::CircleMap,
    s::Integer;
    A_target = nothing,
    T = Float64,
)
    _require_symbolics("FourierAdjoint")
    return FourierAdjoint(k_freq, n, _fourier_weighted_norm(D, s, A_target); T = T)
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

# Projection errors for Cω: the tail of the unit ball of Σ|ĉ_k| e^{2πη|k|} with the
# default η = 0.1, Σ_{|k|>N} e^{-2πη|k|} = 2q^{N+1}/(1-q), q = e^{-2πη}, as for Aη in
# `_aeta_tail_bound` (FourierAnalytic.jl). Until 2026-11-25 this returned exp(+2πNη),
# which grows with N: the sign error the comment above `_aeta_tail_bound` records as
# fixed for Aη.
function _comega_tail_bound(B::FourierAdjoint{Cω})
    q = exp(-2 * interval(pi) * interval(0.1))
    return sup(2 * q^(B.k + 1) / (1 - q))
end
weak_projection_error(B::FourierAdjoint{Cω}) = _comega_tail_bound(B)
aux_normalized_projection_error(B::FourierAdjoint{Cω}) = _comega_tail_bound(B)

function strong_weak_bound(B::FourierAdjoint{Cω})
    k = B.k
    η = interval(0.1)
    λ = 4 * interval(pi) * η
    geo_sum = 1 + 2 * exp(λ) * (1 - exp(λ * k)) / (1 - exp(λ))
    return sup(sqrt(geo_sum))
end

# W^{k,1} strong norm: same polynomial Fourier-tail estimates as the
# `FourierAnalytic` realization — the basis is the same, only the matrix
# representation differs (the adjoint is its transpose), and the basis-level
# constants depend only on the truncation index `B.k` and the norm choice.
weak_projection_error(B::FourierAdjoint{<:Union{W,WeightedSobolev}}) =
    _fourier_sobolev_projection_error(B)
aux_normalized_projection_error(B::FourierAdjoint{<:Union{W,WeightedSobolev}}) =
    _fourier_sobolev_projection_error(B)

function strong_weak_bound(B::FourierAdjoint{W{k,l}}) where {k,l}
    N = B.k
    M = interval(1.0)
    for j = 1:k
        M = M + (2 * interval(pi) * N)^j
    end
    return sup(M)
end


# WeightedSobolev: ||v||_a = Σ_{j>=1} a_j ||v^(j)||_1 <= Σ_{j>=1} a_j (2πK)^j ||v||_{L²}, from
# ||v^(j)||_{L¹} <= ||v^(j)||_{L²} <= (2πK)^j ||v||_{L²} on trigonometric polynomials of degree K.
function strong_weak_bound(B::FourierAdjoint{WeightedSobolev})
    M = interval(0.0)
    for (j, a) in enumerate(B.strong.a)
        M = M + interval(a) * (2 * interval(pi) * B.k)^j
    end
    return sup(M)
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

    @debug typeof(results)

    x = vcat((result for result in results)...)

    return FourierAdjointDual(x)
end

function Dual(B::FourierAdjoint, T::Function; ϵ, max_iter)
    x = T.(B.p)
    return FourierAdjointDual(x)
end

# The adjoint dual carries no weights and no relabelling: node `j` contributes
# ϕ(x_j) to entry `j`. Cf. `eval_on_dual` below.
dual_nodes(B::FourierAdjoint, computed_dual::FourierAdjointDual) =
    (computed_dual.x, 1:length(computed_dual.x), ones(length(computed_dual.x)))

function eval_on_dual(B::FourierAdjoint, computed_dual::FourierAdjointDual, ϕ)

    x = computed_dual.x

    return ϕ.(Interval.(x))
end

# assemble(::FourierAdjoint, D; …) is provided by the FFTWExt extension.
