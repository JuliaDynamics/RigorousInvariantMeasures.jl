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
                nrmj = nrmj ⊕₊ abs_or_mag(A[i, j])
            end
            nrm = max(nrm, nrmj)
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
function opnormbound(::Type{L2}, v::Vector{T}) where {T<:Real}
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

function opnormbound(::Type{L2}, v::Vector{T}) where {T<:Complex}
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

"""
Certified upper bound to the L2 operator norm of a matrix via BallArithmetic.

`upper_bound_L2_opnorm` is `min(Collatz, √(‖·‖₁‖·‖_∞))`: cheap, but on the
matrices that arise here it overestimates the spectral norm by 1.4–1.7×. The
verified-SVD enclosure is sharp to ~1e-8 relative and costs a few ms at the
sizes we use, and the weak norm feeds straight into γ_N, so we prefer it and
keep the `min` with the cheap bound — the result can only improve.

Midpoints may be real or complex, and the element type may be `BigFloat` —
`svdbox` handles it through BallArithmetic's `GenericSchurExt`, so load
`GenericSchur` alongside this package to get the sharp bound in extended
precision. Should the SVD be unavailable or fail, we fall back to the cheap
bound; since we return the `min` of the two, the result is a valid upper bound
either way.
"""
function opnormbound(::Type{L2}, A::AbstractMatrix{T}) where {T}
    BM = BallMatrix(Matrix(A))  # materialize to handle Adjoint/Transpose types
    return _l2_opnorm_ball(BM)
end

function _l2_opnorm_ball(BM::BallMatrix)
    cheap = upper_bound_L2_opnorm(BM)
    sharp = try
        BallArithmetic.svd_bound_L2_opnorm(BM)
    catch
        return cheap
    end
    return (isfinite(sharp) && sharp > 0) ? min(cheap, sharp) : cheap
end

"""
Rigorous upper bound on the L2 norm of a vector, using Parseval identity.
"""
function normbound(::Type{L2}, v::AbstractVector)
    return opnormbound(L2, collect(v))
end

import SparseArrays

function opnormbound(::Type{L1}, A::SparseArrays.SparseMatrixCSC)
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

function opnormbound(::Type{Linf}, A::SparseArrays.SparseMatrixCSC)
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
normbound(N::Type{L1}, v::AbstractVector) =
    opnormbound(L1, v) ⊘₊ Float64(length(v), RoundDown)
normbound(N::Type{Linf}, v::AbstractVector) = opnormbound(Linf, v)
