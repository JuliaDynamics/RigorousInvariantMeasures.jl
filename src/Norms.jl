"""
Functions to deal with various types of norms and seminorms
"""

using FastRounding
using IntervalOptimisation
using SparseArrays: getcolptr
using .DynamicDefinition

using .DynamicDefinition: derivative

"""
'Absolute value' definition that returns mag(I) for an interval and abs(x) for a real
"""
abs_or_mag(x::Number) = abs(x)
abs_or_mag(x::Interval) = mag(x)

"""
Certified upper bound to ||A|| (of specified NormKind)
"""
function opnormbound(::Type{L1}, A::AbstractVecOrMat{T}) where {T}
    # partly taken from JuliaLang's LinearAlgebra/src/generic.jl
    Tnorm = typeof(abs_or_mag(float(real(zero(T)))))
    Tsum = promote_type(Float64, Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for j = 1:size(A, 2)
            nrmj::Tsum = 0
            for i = 1:size(A, 1)
                nrmj = nrmj ⊕₊ abs_or_mag(A[i,j])
            end
            nrm = max(nrm,nrmj)
        end
    end
    return convert(Tnorm, nrm)
end

function opnormbound(::Type{Linf}, A::AbstractVecOrMat{T}) where {T}
    # partly taken from JuliaLang's LinearAlgebra/src/generic.jl
    Tnorm = typeof(abs_or_mag(float(real(zero(T)))))
    Tsum = promote_type(Float64, Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for i = 1:size(A, 1)
            nrmi::Tsum = 0
            for j = 1:size(A, 2)
                nrmi = nrmi ⊕₊ abs_or_mag(A[i,j])
            end
            nrm = max(nrm,nrmi)
        end
    end
    return convert(Tnorm, nrm)
end

function opnormbound(::Type{L1}, A::SparseMatrixCSC) where {T}
    # partly taken from JuliaLang's Sparsearray/src/linalg.jl
    m, n = size(A)
    Tnorm = typeof(abs_or_mag(float(real(zero(eltype(A))))))
    Tsum = promote_type(Float64,Tnorm)
    nA::Tsum = 0
    @inbounds begin
        for j=1:n
            colSum::Tsum = 0
            for i = getcolptr(A)[j]:getcolptr(A)[j+1]-1
                colSum = colSum ⊕₊ abs_or_mag(nonzeros(A)[i])
            end
            nA = max(nA, colSum)
        end
    end
    return convert(Tnorm, nA)
end

function opnormbound(::Type{Linf}, A::SparseMatrixCSC) where {T}
    # partly taken from JuliaLang's Sparsearray/src/linalg.jl
    m, n = size(A)
    Tnorm = typeof(abs_or_mag(float(real(zero(eltype(A))))))
    Tsum = promote_type(Float64,Tnorm)
    rowSum = zeros(Tsum,m)
    @inbounds begin
        for i=1:length(nonzeros(A))
            rowSum[rowvals(A)[i]] = rowSum[rowvals(A)[i]] ⊕₊ abs_or_mag(nonzeros(A)[i])
        end
    end
    return convert(Tnorm, maximum(rowSum))
end

"""
Rigorous upper bound on a vector norm. Note that Linf, L1 are the "analyst's" norms
"""
normbound(N::Type{L1}, v::AbstractVector) = opnormbound(L1, v) ⊘₊ Float64(length(v), RoundDown)
normbound(N::Type{Linf}, v::AbstractVector) = opnormbound(Linf, v)

"""
Types to compute norms iteratively by "adding a column at a time".
"""
abstract type NormCacher{T} end

mutable struct NormCacherL1 <: NormCacher{L1}
    C::Float64
    function NormCacherL1(n)
        new(0)
    end
end
"""
Create a new NormCacher to compute the normbound of the empty matrix with n rows
"""
NormCacher{L1}(n) = NormCacherL1(n)

mutable struct NormCacherLinf <: NormCacher{Linf}
    C::Vector{Float64}
    function NormCacherLinf(n)
        new(zeros(n))
    end
end
NormCacher{Linf}(n) = NormCacherLinf(n)

"""
Update a NormCacher to add one column to the matrix it is computing a norm of.
This column may be affected by an error ε (in the same norm).
"""
function add_column!(Cacher::NormCacherL1, v::AbstractVector, ε::Float64)
    Cacher.C = max(Cacher.C, opnormbound(L1, v) ⊕₊ ε)
end

function add_column!(Cacher::NormCacherLinf, v::AbstractVector, ε::Float64)
    Cacher.C = Cacher.C .⊕₊ abs.(v) .⊕₊ ε
end

"""
Return the norm of the matrix the NormCacher is working on.
"""
function get_norm(Cacher::NormCacherL1)
    return Cacher.C
end

function get_norm(Cacher::NormCacherLinf)
    return maximum(Cacher.C)
end

"""
(A, B) = dfly(strongnorm, auxnorm, dynamic)

Constants (A, B) such that ||Lf||_s ≦ A||f||_s + B||f||_aux
"""
dfly(::Type{<:NormKind}, ::Type{<:NormKind}, ::Dynamic) = @error "Not implemented"

function dfly(::Type{TotalVariation}, ::Type{L1}, D::Dynamic)
    dist = maximise(x -> distorsion(D, x), D.domain)[1]
    lam = maximise(x-> abs(1/derivative(D, x)), D.domain)[1]

    if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
        @error "The function is not expanding"
    end

    if is_full_branch(D)
        return lam.hi, dist.hi
    else
        if !(abs(lam) < 0.5)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        # We need a way to estimate the branch widths
        @error "Not implemented"
    end
end

function dfly(::Type{TotalVariation}, ::Type{L1}, D::PwMap)
    dist = 0
    lam = 0
    disc = 0
    for i in 1:nbranches(D)
        f(x) = D.Ts[i](x)
        domain = hull(D.endpoints[i], D.endpoints[i+1])
        fprime(x) = f(TaylorSeries.Taylor1([x, 1], 1))[1]
        fsecond(x) = f(TaylorSeries.Taylor1([x, 1], 2))[2]/2
        distorsion(x)=abs(fsecond(x)/(fprime(x)^2))
        lambda(x) = abs(1/fprime(x))
        dist = max(dist, maximise(distorsion, domain)[1].hi)
        lam = max(lam, maximise(lambda, domain)[1].hi)
        low_rad = (abs(D.endpoints[i]-D.endpoints[i+1])/2).lo
        disc = max(disc, ((1/Interval(low_rad)).hi))
    end

    if is_full_branch(D)
        if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
            @error "The function is not expanding"
        end
        return lam, dist
    else
        if !(abs(2⊗₊lam) < 1)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        return 2⊗₊lam, dist ⊕₊ disc
    end
end

function dfly(::Type{Lipschitz}, ::Type{L1}, D::Dynamic)
    # TODO: should assert that D is globally C2 instead, but we don't have that kind of infrastructure yet.
    @assert is_full_branch(D)

    dist = maximise(x -> distorsion(D, x), D.domain)[1]
    lam = maximise(x-> abs(1/derivative(D, x)), D.domain)[1]

    lam = lam.hi
    dist = dist.hi

    return lam ⊗₊ (2. ⊗₊  dist ⊕₊ 1.), dist ⊗₊ (2. ⊗₊ dist ⊕₊ 1.)
end
