
#export DiscretizedNoiseKernel, UniformNoise

using LinearAlgebra

abstract type NoiseKernel end
function Base.:*(M::NoiseKernel, v)
    @error "Not implemented"
end

opnormbound(N::NormKind, M::NoiseKernel) = @error "Not Implemented"
opradius(N::NormKind, M::NoiseKernel) = @error "Not Implemented"
nonzero_per_row(M::NoiseKernel) = @error "Not Implemented"

# DiscretizedNoiseKernelFFT, UniformNoiseFFT and Mfft are defined in the
# FFTWExt extension (loaded via `using FFTW`). They were never wired into
# any test or example, but kept in case someone wants the FFT-based
# uniform-noise kernel later.

"""
    DiscretizedNoiseKernelUlam{T,S}

Ulam discretization of a noise kernel whose bounds are of type `T`, one of the
types `Interval{T}` is built on.

`T` is a parameter rather than fixed to `Float64` because `w` and `z` are the
scratch buffers the multiplication writes through: they are touched once per
output entry, so at the abstract `Vector` they cost a boxed load every
iteration, 2.2 MB per application at `k = 16384`; and fixed to `Vector{Float64}`
they would silently narrow a `Vector{Interval{BigFloat}}` argument to double
precision, which the error term added afterwards does not account for.
"""
struct DiscretizedNoiseKernelUlam{T<:Real,S<:AbstractVector{T}} <: NoiseKernel
    B::Ulam
    ξ::Any
    v::S
    rad::T
    boundarycondition::Symbol
    w::Vector{T}
    z::Vector{T}
end

"""
    UniformNoiseUlam([T,] ξ, B::Ulam, boundarycondition = :periodic)

Ulam discretization of the uniform noise kernel of half-width `ξ`, with bounds
of type `T` (`Float64` by default). Apply it to a `Vector{Interval{T}}` of the
matching type: a vector of a wider type would be narrowed to `T` on the way into
the scratch buffers, and the error term is not written to cover that.
"""
UniformNoiseUlam(ξ, B::Ulam, boundarycondition = :periodic) =
    UniformNoiseUlam(Float64, ξ, B, boundarycondition)

function UniformNoiseUlam(
    ::Type{T},
    ξ,
    B::Ulam,
    boundarycondition = :periodic,
) where {T<:Real}
    k = length(B)
    n = 2 * Int64(ceil(ξ * k))
    v = zeros(Interval{T}, n)
    a = 1 / (2 * interval(T, ξ))
    v[2:n-1] = a * ones(Interval{T}, n - 2)
    v[1] = (k - sum(v)) / 2
    v[n] = v[1]
    nw = boundarycondition == :reflecting ? k + n + 2 : k + n
    boundarycondition ∈ (:periodic, :reflecting) || return throw(
        ArgumentError("boundary condition must be :periodic or :reflecting, got $boundarycondition"),
    )
    return DiscretizedNoiseKernelUlam(
        B,
        interval(T, ξ),
        mid.(v),
        T(opnormbound(L1, v) - k, RoundUp),
        boundarycondition,
        zeros(T, nw),
        zeros(T, k),
    )
end

#TODO, but at the moment this is fine, it is a Markov operator
opnormbound(B::Ulam, ::Type{L1}, M::DiscretizedNoiseKernelUlam) = 1.0
opradius(::Type{L1}, M::DiscretizedNoiseKernelUlam) = M.rad
nonzero_per_row(M::DiscretizedNoiseKernelUlam) = length(M.v)
# Lemma 47 of Galatolo-Monge-Nisoli gives ‖N_ξ‖_{L¹→Var} ≤ Var(ρ_ξ), and for the
# uniform kernel of half-width ξ the density is 1/(2ξ) on a support of length
# 2ξ, so it rises once and falls once and Var(ρ_ξ) = 1/ξ. This returned 1/(2ξ),
# which is half of that and so not an upper bound: a spike of unit L¹ mass is
# spread to a box of height 1/(2ξ), whose variation is exactly 1/ξ.
dfly(::Type{TotalVariation}, ::Type{L1}, N::DiscretizedNoiseKernelUlam) =
    (0.0, sup(1 / N.ξ))


function Base.:*(M::DiscretizedNoiseKernelUlam, v)
    mult(M, v, Val(M.boundarycondition))
end

"""
Wrap `src` periodically into the workspace `w`, so that `w[j] = src[j-l mod k]`.
"""
function _extend_periodic!(w, src, l, k)
    w .= 0
    @views w[l+1:l+k] .= src
    @views w[1:l] .= src[end-l+1:end]
    @views w[end-l+1:end] .= src[1:l]
    return w
end

# `dot` rather than `sum(M.v .* h)`: the latter allocated a fresh length-n
# vector on each of the k iterations, 222 MB per application at k = 16384.
# The error bound below is unaffected, since γ_n bounds the error of any
# summation order.
function mult(M::DiscretizedNoiseKernelUlam, v, ::Val{:periodic})
    n = length(M.v)
    k = length(v)
    l = n ÷ 2

    _extend_periodic!(M.w, v, l, k)

    @inbounds for i = 1:k
        h = @view M.w[i:i+n-1]
        v[i] = dot(M.v, h) / k
    end

    return v
end

function mult(
    M::DiscretizedNoiseKernelUlam,
    v::Vector{Interval{T}},
    ::Val{:periodic},
) where {T}
    n = length(M.v)
    k = length(v)
    l = n ÷ 2

    nrmv = opnormbound(M.B, L1, v)
    midv = mid.(v)
    radv = radius.(v)
    nrmrad = opnormbound(M.B, L1, radv)

    _extend_periodic!(M.w, midv, l, k)

    @inbounds for i = 1:k
        h = @view M.w[i:i+n-1]
        v[i] = dot(M.v, h) / k
    end
    δₖ = opradius(L1, M)
    γₖ = gamma(T, nonzero_per_row(M))
    nrm_MK = opnormbound(M.B, L1, M)
    normMK = nrm_MK ⊕₊ δₖ

    ϵ = (γₖ ⊗₊ normMK) ⊗₊ nrmv ⊕₊ normMK ⊗₊ nrmrad

    return v + fill(interval(-ϵ, ϵ), length(v))
end



"""
Accumulate the reflecting-boundary convolution of the float vector `src` into
`dest`, which must have length `k` and is overwritten.

This is definition 8 of Galatolo-Monge-Nisoli read literally: extend `src` by
zero outside `1:k`, convolve on ℤ, then push the mass that fell outside back in
with the mirror map `π` of `reflect_outward_idx`. Every step preserves mass, so
the result is Markov whatever the parity of the window.

Reflecting the input instead and reading off `1:k` is the adjoint of this, not
the same operator; with the even window this kernel builds it is not even
mass-preserving, losing O(1/k) of the mass on a density that is not symmetric.
"""
function _convolve_reflecting!(dest, M, src, l, k, n)
    fill!(dest, zero(eltype(dest)))
    @inbounds for j = (1-l):(k+l)
        s = zero(eltype(dest))
        for t = 1:n
            idx = j + t - 1 - l
            if 1 <= idx <= k
                s += M.v[t] * src[idx]
            end
        end
        dest[reflect_outward_idx(j, k)] += s / k
    end
    return dest
end

function mult(M::DiscretizedNoiseKernelUlam, v, ::Val{:reflecting})
    n = length(M.v)
    k = length(v)
    l = n ÷ 2

    _convolve_reflecting!(M.z, M, v, l, k, n)
    v .= M.z

    return v
end

function mult(
    M::DiscretizedNoiseKernelUlam,
    v::Vector{Interval{T}},
    ::Val{:reflecting},
) where {T}
    n = length(M.v)
    k = length(v)
    l = n ÷ 2

    nrmv = opnormbound(M.B, L1, v)
    midv = mid.(v)
    radv = radius.(v)
    nrmrad = opnormbound(M.B, L1, radv)

    _convolve_reflecting!(M.z, M, midv, l, k, n)

    δₖ = opradius(L1, M)
    γₖ = gamma(T, nonzero_per_row(M))
    nrm_MK = opnormbound(M.B, L1, M)
    normMK = nrm_MK ⊕₊ δₖ

    ϵ = (γₖ ⊗₊ normMK) ⊗₊ nrmv ⊕₊ normMK ⊗₊ nrmrad

    for i = 1:k
        v[i] = interval(M.z[i])
    end

    return v + fill(interval(-ϵ, ϵ), length(v))
end
