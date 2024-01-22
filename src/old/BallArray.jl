# COV_EXCL_START
"""
Implement some ball matrix arithmetic, i.e., verified matrix operations
for (abstract) interval matrices stored in midpoint-radius format.

The 'radius' for a matrix can be either a matrix of the same dimension, or
a scalar which represents a norm.

A type parameter specifies the "norm type" for the radius. It can be `:Abs` for
entrywise radii, or `OpnormL1`, `OpnormLinf` for normwise radii (specifying which norm to use)
"""

using FastRounding

abstract type NormKind end
struct OpnormL1 <: NormKind end
struct OpnormLinf <: NormKind end
struct Abs <: NormKind end

"""
Certified upper bound to ||A||, abs(A), or whatever is used in the specified radius
"""
function normbound(A::AbstractVecOrMat{T}, ::OpnormL1) where {T}
    # partly taken from JuliaLang's LinearAlgebra/src/generic.jl
    Tnorm = typeof(float(real(zero(T))))
    Tsum = promote_type(Float64, Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for j = 1:size(A, 2)
            nrmj::Tsum = 0
            for i = 1:size(A, 1)
                nrmj = nrmj ⊕₊ abs(A[i, j])
            end
            nrm = max(nrm, nrmj)
        end
    end
    return convert(Tnorm, nrm)
end

function normbound(A::AbstractVecOrMat{T}, ::OpnormLinf) where {T}
    # partly taken from JuliaLang's LinearAlgebra/src/generic.jl
    Tnorm = typeof(float(real(zero(T))))
    Tsum = promote_type(Float64, Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for i = 1:size(A, 1)
            nrmi::Tsum = 0
            for j = 1:size(A, 2)
                nrmi = nrmi ⊕₊ abs(A[i, j])
            end
            nrm = max(nrm, nrmi)
        end
    end
    return convert(Tnorm, nrm)
end

function normbound(A::AbstractArray{T}, ::Abs) where {T}
    return abs.(A)
end

radiustype(MatrixType, ::OpnormLinf) = Float64
radiustype(MatrixType, ::OpnormL1) = Float64
radiustype(MatrixType, ::Abs) = MatrixType

"""
Type for ball arrays. Comes in two variants, one with cached normbound(midpoint)
(handy if you have to compute many products with it) and one without.
"""
abstract type BallArray{
    N<:NormKind,
    MatrixType<:Union{AbstractArray,Number},
    RadiusType<:Union{AbstractArray,Number},
} end

struct BallArrayWithoutCachedNorm{N,MatrixType,RadiusType} <:
       BallArray{N,MatrixType,RadiusType}
    midpoint::MatrixType
    radius::RadiusType
    function BallArrayWithoutCachedNorm{N,MatrixType,RadiusType}(
        midpoint,
        radius = zero(RadiusType),
    ) where {N,MatrixType,RadiusType}
        if !isa(radius, radiustype(MatrixType, N()))
            error("Wrong radius type specified")
        end
        return new(midpoint, radius)
    end
end

struct BallArrayWithCachedNorm{N,MatrixType,RadiusType} <:
       BallArray{N,MatrixType,RadiusType}
    midpoint::MatrixType
    radius::RadiusType
    norm::RadiusType
    function BallArrayWithCachedNorm{N,MatrixType,RadiusType}(
        midpoint,
        radius = zero(RadiusType),
    ) where {N,MatrixType,RadiusType}
        if !isa(radius, radiustype(MatrixType, N()))
            error("Wrong radius type specified")
        end
        return new{N,MatrixType,RadiusType}(midpoint, radius, normbound(midpoint, N()))
    end
end

function BallArray{N}(midpoint, radius) where {N}
    return BallArrayWithoutCachedNorm{N,typeof(midpoint),typeof(radius)}(midpoint, radius)
end

function BallArrayWithCachedNorm{N}(midpoint, radius) where {N}
    return BallArrayWithCachedNorm{N,typeof(midpoint),typeof(radius)}(midpoint, radius)
end

function midpoint_norm(A::BallArrayWithCachedNorm)
    return A.norm
end

function midpoint_norm(A::BallArrayWithoutCachedNorm{N,MT,RT}) where {N,MT,RT}
    return normbound(A, N())
end

"""
Product between ball matrices / vectors.
"""
function Base.:*(A::BallArray, B::BallArray)
    # A*maximum(B, dims=2) TODO
end

# TODO: using this file for tests

function matvec(A, v)
    w = zero(v)
    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)
    for j = 1:n
        for k in nzrange(A, j)
            @inbounds i = rows[k]
            @inbounds Aij = vals[k]
            @inbounds w[i] += Aij * v[j]
        end
    end
end

function matvec_roundup(A, v)
    w = zero(v)
    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)
    for j = 1:n
        for k in nzrange(A, j)
            @inbounds i = rows[k]
            @inbounds Aij = vals[k]
            @inbounds w[i] = w[i] ⊕₊ Aij ⊗₊ v[j]
        end
    end
end
# COV_EXCL_STOP
