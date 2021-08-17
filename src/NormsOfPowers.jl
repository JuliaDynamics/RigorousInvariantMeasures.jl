"""
Functions to estimate Q|_{U^0}. See our paper for details.
"""

using LinearAlgebra
using SparseArrays
using FastRounding
using ValidatedNumerics

export norms_of_powers, refine_norms_of_powers, norms_of_powers_dfly, norms_of_powers_trivial, norms_of_powers_from_coarser_grid

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
    nu = u ⊗₊ T(n) #TODO: in theory, this should be rounded up/down. In practice, all integers up to 2^53 or so fit in Float64, so it won't be needed.
    return nu ⊘₊ (one(T) ⊖₋ nu)
end

using ProgressMeter
"""
Estimates the norms ||Q||, ||Q^2||, ... ||Q^m|| on U^0.

U is the matrix [ones(1,n-1); -I_(n-1,n-1)]. It is currently assumed that
f*U==0 (i.e., all elements of f are equal).

f must be an interval vector.

The following constants may be specified as keyword arguments:

normQ, normE, normv0, normEF, normIEF, normN

otherwise they are computed (which may be slower).

e and f must be specified in case is_integral_preserving==false
In case is_integral_preserving is true, they may be specified but they are then ignored.

Implementation note: currently we perform this computation one column at a time,
to be able to scale (slowly) to cases with large size; for moderate sizes, it would
indeed be better to do the computation all columns at the same time, in BLAS level 3.
"""
function norms_of_powers(N::Type{<:NormKind}, m::Integer, Q::DiscretizedOperator, f::AbstractArray;
        normv0::Real=-1., #used as "missing" value
        normQ::Real=-1.,
        normE::Real=-1.,
        normEF::Real=-1.,
        normIEF::Real=-1.,
        normN::Real=-1.)

    @assert eltype(f) <: Interval
    T = typeof(zero(eltype(Q.L)).hi) # gets "Float64" from Q.L
    n = size(Q.L, 1)
    M = mid.(Q.L)
    R = radius.(Q.L)
    δ = opnormbound(N, R)
    γz = gamma(T, max_nonzeros_per_row(Q.L))
    γn = gamma(T, n+3) # not n+2 like in the paper, because we wish to allow for f to be the result of rounding
    ϵ = zero(T)

    nrmM = opnormbound(N, M)

    # precompute norms
    if !is_integral_preserving(Q)
        if normE == -1.
            normE = opnormbound(N, Q.e)
        end
        if normEF == -1.
            normEF = opnormbound(N, Q.e*f)
        end
        if normIEF == -1.
            normIEF =  opnormbound(N, [Matrix(UniformScaling{Float64}(1),n,n) Q.e*f])
        end
        if normN == -1.
            normN = opnormbound(N, Matrix(UniformScaling{Float64}(1),n,n) - Q.e*f)
        end
    end

    if normQ == -1.
        if is_integral_preserving(Q)
            normQ = nrmM ⊕₊ δ
        else
            defect = opnormbound(N, Q.w)
            normQ = nrmM ⊕₊ δ ⊕₊ normE ⊗₊ defect
        end
    end

    # initialize normcachers
    normcachers = [NormCacher{N}(n) for j in 1:m]

    midf = map(mid, f)

    # main loop

    v = zeros(T, n)
    @showprogress for j = 1:n-1
        v .= zero(T) # TODO: check for type stability in cases with unusual types
        v[1] = one(T) # TODO: in full generality, this should contain entries of f rather than ±1
        v[j+1] = -one(T)
        if normv0 == -1.
            nrmv = opnormbound(N, v)
        else
            nrmv = normv0
        end
        ϵ = 0.
        nrmw = nrmv # we assume that initial vectors are already integral-preserving
        for k = 1:m
            w = M * v
            if is_integral_preserving(Q)
                v = w
                ϵ = (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
            else
                v = w - Q.e * (midf*w)  # TODO: we are currently assuming that f is not too large, to estimate the error (the result of only one floating point operation)
                new_nrmw = opnormbound(N, w)
                ϵ = γn ⊗₊ normIEF ⊗₊ (new_nrmw ⊕₊ normEF ⊗₊ nrmw) ⊕₊ normN ⊗₊ (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
                nrmw = new_nrmw
            end
            nrmv = opnormbound(N, v)
            add_column!(normcachers[k], v, ϵ) #TODO: Could pass and reuse nrmv in the case of norm-1
        end
    end
    return map(get_norm, normcachers)
end

"""
Array of "trivial" bounds for the powers of a DiscretizedOperator (on the whole space)
coming from from ||Q^k|| ≤ ||Q||^k
"""
function norms_of_powers_trivial(N::Type{<:NormKind}, Q::DiscretizedOperator, m::Integer)
    norms = fill(NaN, m)
    norms[1] = opnormbound(N, Q)
    for i = 2:m
        norms[i] = norms[i-1] ⊗₀ norms[1]
    end
    return norms
end

"""
Arrays of bounds to ||Q^k||_{w → s} = sup_{||f||_w=1} ||Q^k f||_s
and to ||Q^k||_{w}
coming theoretically from iterated DFLY inequalities (the "small matrix method").

Returns two arrays (strongs, norms) of length m:
strongs[k] bounds ||Q^k f||_s, norms[k] bounds ||Q^k f||)
"""
function norms_of_powers_dfly(Bas::Basis, D::Dynamic, m)
    A, B = dfly(strong_norm(Bas), aux_norm(Bas), D)
    Eh = BasisDefinition.aux_normalized_projection_error(Bas)
    M₁n = BasisDefinition.strong_weak_bound(Bas)
    M₂ = BasisDefinition.aux_weak_bound(Bas)
    S₁, S₂ = BasisDefinition.weak_by_strong_and_aux_bound(Bas)

    norms = fill(NaN, m)
    strongs = fill(NaN, m)

    v = Array{Float64}([M₁n; M₂])
    # We evaluate [S₁ S₂] * ([1 0; Eh 1]*[A B; 0 1])^k * [M₁n; M₂] (with correct rounding)
    for k = 1:m
        # invariant: v[1] bounds ||Q^kf||_s for ||f||_w=1
        # v[2] bounds |||Q^kf||| for ||f||_w=1
        v[1] = A ⊗₊ v[1] ⊕₊ B ⊗₊ v[2]
        v[2] = Eh ⊗₊ v[1] ⊕₊ v[2]
        strongs[k] = v[1]
        norms[k] = S₁ ⊗₊ v[1] ⊕₊ S₂ ⊗₊ v[2]
    end
    return strongs, norms
end


"""
Compute better and/or more estimates of power norms using
the fact that ||Q^{k+h}|| ≤ ||Q^k|| * ||Q^h||.
This uses multiplicativity, so it will not work for mixed norms,
e.g., ||Q^k||_{s → w}, or ||M^k|_{U^0}||
(unless M preserves U^0, which is the case for Q|_{U^0}).
"""
function refine_norms_of_powers(norms::Vector, m)
    better_norms = fill(NaN, m)
    better_norms[1] = norms[1]
    for k in 2:m
        better_norms[k] = minimum(better_norms[i] ⊗₊ better_norms[k-i] for i in 1:k-1)
        if k <= length(norms)
            better_norms[k] = min(better_norms[k], norms[k])
        end
    end
    return better_norms
end
refine_norms_of_powers(norms::Vector) = refine_norms_of_powers(norms, length(norms))

"""
Estimate norms of powers from those on a coarser grid (see paper for details)
"""
function norms_of_powers_from_coarser_grid(fine_basis::Basis, coarse_basis::Basis, D::Dynamic, coarse_norms::Vector, normQ::Real)
    if !BasisDefinition.is_refinement(fine_basis, coarse_basis)
        @error "The fine basis is not a refinement of the coarse basis"
    end
    m = length(coarse_norms)
    fine_norms = fill(NaN, m)
    (strongs, norms) = norms_of_powers_dfly(fine_basis, D, m)

    # adds a 0th element to strongs
    strongs0(k::Integer) = k==0 ? BasisDefinition.strong_weak_bound(fine_basis) : strongs[k]
    coarse_norms0(k::Integer) = k==0 ? 1. : coarse_norms[k]

    Kh =  BasisDefinition.weak_projection_error(coarse_basis)
    for k in 1:m
		temp = 0.
		for k in 0:m-1
			temp = temp ⊕₊ coarse_norms0(m-1-k) ⊗₊ (normQ ⊗₊ strongs0(k) ⊕₊ strongs0(k+1))
		end
		fine_norms[k] = coarse_norms[k] ⊕₊ 2. ⊗₊ Kh ⊗₊ temp
	end
    return fine_norms
end

"""
Estimate ||I + Q + Q^2 + … || (infinite sum) using a list of computed norm bounds norms[k] ≥ ||Q^k||.
"""
function infinite_sum_norms(norms::Vector)
    m = length(norms)
    if norms[m] >= 1
        @error "The last norm is >= 1: the bounds are not sufficient to show that the sum converges"
    end
    S = 1.
    for k = 1:m-1
        S = S ⊕₊ norms[k]
    end
    S = S ⊘₊ (1. ⊖₋ norms[m])
    return S
end
