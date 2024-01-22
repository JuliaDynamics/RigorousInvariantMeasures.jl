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
                nrmj = nrmj ⊕₊ abs_or_mag(A[i, j])
            end
            nrm = max(nrm, nrmj)
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
                nrmi = nrmi ⊕₊ abs_or_mag(A[i, j])
            end
            nrm = max(nrm, nrmi)
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

import SparseArrays

function BasisDefinition.opnormbound(::Type{L1}, A::SparseArrays.SparseMatrixCSC)
    # partly taken from JuliaLang's Sparsearray/src/linalg.jl
    m, n = size(A)
    Tnorm = typeof(abs_or_mag(float(real(zero(eltype(A))))))
    Tsum = promote_type(Float64, Tnorm)
    nA::Tsum = 0
    @inbounds begin
        for j = 1:n
            colSum::Tsum = 0
            for i = getcolptr(A)[j]:getcolptr(A)[j+1]-1
                colSum = colSum ⊕₊ abs_or_mag(nonzeros(A)[i])
            end
            nA = max(nA, colSum)
        end
    end
    return convert(Tnorm, nA)
end

function BasisDefinition.opnormbound(::Type{Linf}, A::SparseArrays.SparseMatrixCSC)
    # partly taken from JuliaLang's Sparsearray/src/linalg.jl
    m, n = size(A)
    Tnorm = typeof(abs_or_mag(float(real(zero(eltype(A))))))
    Tsum = promote_type(Float64, Tnorm)
    rowSum = zeros(Tsum, m)
    @inbounds begin
        for i = 1:length(nonzeros(A))
            rowSum[rowvals(A)[i]] = rowSum[rowvals(A)[i]] ⊕₊ abs_or_mag(nonzeros(A)[i])
        end
    end
    return convert(Tnorm, maximum(rowSum))
end

"""
Rigorous upper bound on a vector norm. Note that Linf, L1 are the "analyst's" norms
"""
BasisDefinition.normbound(N::Type{L1}, v::AbstractVector) =
    opnormbound(L1, v) ⊘₊ Float64(length(v), RoundDown)
BasisDefinition.normbound(N::Type{Linf}, v::AbstractVector) = opnormbound(Linf, v)
