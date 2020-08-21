"""
Functions to estimate Q|_{U^0}. See our paper for details.
"""

using LinearAlgebra
using SparseArrays
using FastRounding
using ValidatedNumerics
using ValidatedNumerics.IntervalArithmetic: round_expr

"""
Returns the maximum number of (structural) nonzeros in a row of A
"""
function max_nonzeros_per_row(A::SparseMatrixCSC)
    rows = rowvals(A)
    m, n = size(A)
    nonzeros_in_each_row = zeros(eltype(rows), m)
    for i in rows
        nonzeros_in_each_row[i] += 1
    end
    return maximum(nonzeros_in_each_row)
end

"""
γₙ constants for floating point error estimation, as in [Higham, Accuracy and Stability of Numerical Algorithms]
"""
function gamma(T, n::Integer)
    u = eps(T)
    nu = u ⊗₊ T(n)
    return nu ⊘₊ (one(T) ⊖₋ nu)
end

"""
Estimates the norms ||Q||, ||Q^2||, ... ||Q^m|| on U^0.

Q is the matrix M if is_integral_preserving==true, or
M + e*(f-f*M) otherwise. Here M is a matrix ∈ LL.

U is the matrix [ones(1,n-1); I_(n-1,n-1)]. It is currently assumed that
f*U==0, because the initial vectors v are supposed to have f*v==0.
This assumption would better be relaxed in future.

The following constants may be specified as keyword arguments:

normQ, normE, normv0, normEF, normIEF, normN

otherwise they are computed (which may be slower).

e and f must be specified in case is_integral_preserving==false
In case is_integral_preserving is true, they may be specified but they are then ignored.
(TODO: maybe this should be better integrated in the syntax).
"""
function norm_of_powers(N::NormKind, m::Integer, LL::SparseMatrixCSC{Interval{RealType}, IndexType}, is_integral_preserving::Bool ;
        e::Vector=[0.],
        f::Adjoint=adjoint([0.]),
        normv0::Real=-1., #used as "missing" value
        normQ::Real=-1.,
        normE::Real=-1.,
        normEF::Real=-1.,
        normIEF::Real=-1.,
        normN::Real=-1.) where {RealType, IndexType}

    n = size(LL, 1)
    M = mid.(LL)
    R = radius.(LL)
    δ = opnormbound(R, N)
    γz = gamma(RealType, max_nonzeros_per_row(LL))
    γn = gamma(RealType, n+3) # not n+2 like in the paper, because we wish to allow for f to be the result of rounding
    ϵ = zero(RealType)

    nrmM = opnormbound(M, N)
    if normQ == -1
        normQ = is_integral_preserving ? nrmM ⊕₊ δ : nrmM ⊕₊ δ ⊕₊ normE * opnormbound(f - f*LL, N)
    end

    # precompute norms
    if !is_integral_preserving
        if normE == -1.
            normE = opnormbound(e, N)
        end
        if normQ == -1.
            if is_integral_preserving
                normQ = nrmM ⊕₊ δ
            else
                defect = opnormbound(f - f*LL, N)
                normQ = nrmM ⊕₊ δ ⊕₊ normE * defect
            end
        end
        if normEF == -1.
            normEF = opnormbound(e*f, N)
        end
        if normIEF == -1.
            normIEF =  opnormbound([Matrix(UniformScaling{Float64}(1),n,n) e*f], N)
        end
        if normN == -1.
            normN = opnormbound(Matrix(UniformScaling{Float64}(1),n,n) - e*f, N)
        end
    end

    # initialize normcachers
    normcachers = [NormCacher{typeof(N)}(n) for j in 1:m]

    # main loop

    for j = 1:n-1
        v = zeros(n) # TODO: check for type stability in cases with unusual types
        v[1] = 1. # TODO: in full generality, this should contain entries of f rather than ±1
        v[j+1] = -1.
        if normv0 == -1.
            nrmv = opnormbound(v, N)
        else
            nrmv = normv0
        end
        ϵ = 0.
        nrmw = nrmv # we assume that initial vectors are already integral-preserving
        for k = 1:m
            w = M * v
            if is_integral_preserving
                v = w
                ϵ = round_expr((γz * nrmM + δ)*nrmv + normQ*ϵ, RoundUp)
            else
                v = w - e * (f*w)
                new_nrmw = opnormbound(w, N)
                ϵ = round_expr(γn*normIEF*(new_nrmw + normEF*nrmw) + normN*(γz*nrmM + δ)*nrmv + normQ*ϵ, RoundUp)
                nrmw = new_nrmw
            end
            nrmv = opnormbound(v, N)
            add_column!(normcachers[k], v, ϵ) #TODO: Could pass and reuse nrmv in the case of norm-1
        end
    end
    return map(get_norm, normcachers)
end
