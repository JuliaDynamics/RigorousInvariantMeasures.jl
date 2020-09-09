"""
Functions to deal with various types of norms and seminorms
"""

using FastRounding
using ValidatedNumerics

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

# TODO: specialize for sparse matrices

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

	distorsion(x)=der_der(D, x)/(der(D, x)^2)
	lambda(x) = 1/der(D, x)
    dist = range_estimate(distorsion, D.domain)
    lam = range_estimate(lambda, D.domain)

    if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) <= 1.
        @error "The function is not expanding"
    end

    if is_full_branch(D)
        return abs(lam).hi, abs(dist).hi
    else
        if !(abs(lam) < 0.5)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        # We need a way to estimate the branch widths
        @error "Not implemented"
    end
end

function dfly(::Type{Lipschitz}, ::Type{L1}, D::Dynamic)
    distorsion(x)=der_der(D, x)/(der(D, x)^2)
	lambda(x) = 1/der(D, x)
    dist = range_estimate(distorsion, D.domain)
    lam = range_estimate(lambda, D.domain)

    lam = abs(lam).hi
    dist = abs(dist).hi

    return round_expr(lam*(2*dist+1), RoundUp), round_expr(dist*(2*dist+1), RoundUp)
end
