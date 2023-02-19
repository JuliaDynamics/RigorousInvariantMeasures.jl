"""
Functions to deal with various types of norms and seminorms
"""

using FastRounding
using IntervalOptimisation
using TaylorSeries: Taylor1
using SparseArrays: getcolptr
using .DynamicDefinition


"""
'Absolute value' definition that returns mag(I) for an interval and abs(x) for a real
"""
abs_or_mag(x::Number) = Float64(abs(x), RoundUp)
abs_or_mag(x::Interval) = Float64(mag(x), RoundUp)

"""
Computes a rigorous upper bound for z*z'
"""
z_times_conjz(z::Complex) = square_round(abs_or_mag(real(z)), RoundUp) ⊕₊ square_round(abs_or_mag(imag(z)), RoundUp)
abs_or_mag(z::Complex) = sqrt_round(z_times_conjz(z), RoundUp)

"""
Certified upper bound to ||A|| (of specified NormKind)
"""
function BasisDefinition.opnormbound(::Type{L1}, A::AbstractVecOrMat{T}) where {T}
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

function BasisDefinition.opnormbound(::Type{Linf}, A::AbstractVecOrMat{T}) where {T}
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

"""
These functions compute a rigorous upper bound for the 2-norm of a vector;
we have a specialized version for complex numbers to avoid taking
the sqrt root and squaring again 
"""
function BasisDefinition.opnormbound(::Type{L2}, v::Vector{T}) where {T<:Real}
    # partly taken from JuliaLang's LinearAlgebra/src/generic.jl
    Tnorm = typeof(abs_or_mag(float(real(zero(T)))))
    Tsum = promote_type(Float64, Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for j = 1:length(v)
            nrm = nrm ⊕₊ square_round(abs_or_mag(v[j]), RoundUp)
        end
    end
    return convert(Tnorm, sqrt_round(nrm, RoundUp))
end

function BasisDefinition.opnormbound(::Type{L2}, v::Vector{T}) where {T<:Complex}
    # partly taken from JuliaLang's LinearAlgebra/src/generic.jl
    Tnorm = typeof(abs_or_mag(float(real(zero(T)))))
    Tsum = promote_type(Float64, Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for j = 1:length(v)
            nrm = nrm ⊕₊ z_times_conjz(v[j])
        end
    end
    return convert(Tnorm, sqrt_round(nrm, RoundUp))
end

function BasisDefinition.opnormbound(::Type{L1}, A::SparseMatrixCSC)
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

function BasisDefinition.opnormbound(::Type{Linf}, A::SparseMatrixCSC)
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
BasisDefinition.normbound(N::Type{L1}, v::AbstractVector) = opnormbound(L1, v) ⊘₊ Float64(length(v), RoundDown)
BasisDefinition.normbound(N::Type{Linf}, v::AbstractVector) = opnormbound(Linf, v)

"""
Types to compute norms iteratively by "adding a column at a time".
"""
abstract type NormCacher{T} end

mutable struct NormCacherL1 <: NormCacher{L1}
    B::Basis
    C::Float64
    function NormCacherL1(B, n)
        new(B, 0)
    end
end
"""
Create a new NormCacher to compute the normbound of the empty matrix with n rows
"""
NormCacher{L1}(B, n) = NormCacherL1(B, n)

mutable struct NormCacherLinf <: NormCacher{Linf}
    B::Basis
    C::Vector{Float64}
    function NormCacherLinf(B, n)
        new(B, zeros(n))
    end
end
NormCacher{Linf}(B, n) = NormCacherLinf(B, n)

"""
Update a NormCacher to add one column to the matrix it is computing a norm of.
This column may be affected by an error ε (in the same norm).
"""
function add_column!(Cacher::NormCacherL1, v::AbstractVector, ε::Float64)
    Cacher.C = max(Cacher.C, opnormbound(Cacher.B, L1, v) ⊕₊ ε)
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

# I don't think this is used in production anymore
function dfly(N1::Type{TotalVariation}, N2::Type{L1}, D::Dynamic)
    dist = max_distortion(D)
    lam = max_inverse_derivative(D)

    if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
        @error "The function is not expanding"
    end

    if is_full_branch(D)
        return lam.hi, dist.hi
    else
        if !(abs(lam) < 0.5)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        endpts = endpoints(D)
        min_width = minimum([endpts[i+1]-endpts[i] for i in 1:length(endpts)-1])
        return lam.hi, dist.hi⊕₊(2/min_width).hi
    end
end

function dfly(N1::Type{TotalVariation}, N2::Type{L1}, D::PwMap)
    if has_infinite_derivative_at_endpoints(D)
        return dfly_inf_der(N1, N2, D, 10^-3)
    end
    
    dist = max_distortion(D)
    lam = max_inverse_derivative(D)
    vec = endpoints(D)
    disc = maximum(2/abs(vec[i]-vec[i+1]) for i in 1:nbranches(D))

    if is_full_branch(D)
        if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
            @error "The function is not expanding"
        end
        return lam.hi, dist.hi
    else
        if !(abs(2*lam) < 1)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        return (2*lam).hi, (dist + disc).hi
    end
end

function dfly(::Type{Lipschitz}, ::Type{L1}, D::Dynamic)
    # TODO: should assert that D is globally C2 instead, but we don't have that kind of infrastructure yet.
    @assert is_full_branch(D)

    dist = max_distortion(D)
    #@info dist
    lam = max_inverse_derivative(D)
    #@info lam

    return ((lam*(2*dist+1)).hi, (dist*(dist+1)).hi)
end
