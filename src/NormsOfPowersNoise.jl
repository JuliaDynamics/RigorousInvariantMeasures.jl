"""
Functions to estimate Q|_{U^0}. See our paper for details.
"""

using LinearAlgebra
using SparseArrays
using FastRounding
using ValidatedNumerics
using FFTW

export norms_of_powers_noise

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
"""
function norms_of_powers_noise(N::Type{L1}, 
        m::Integer, 
        Q::DiscretizedOperator,
        MK::NoiseKernel, 
        f::AbstractArray;
        normv0::Real=-1., #used as "missing" value
        normQ::Real=-1.,
        normE::Real=-1.,
        normEF::Real=-1.,
        normIEF::Real=-1.,
        normN::Real=-1.,
        normρ::Real=-1.)

    @assert eltype(f) <: Interval
    T = typeof(zero(eltype(Q.L)).hi) # gets "Float64" from Q.L
    n = size(Q.L, 1)
    M = mid.(Q.L)
        
    R = radius.(Q.L)
    δ = opnormbound(N, R)
    
    # this is the operator radius of the noise operator,
    # we use this function since the operator may be defined implictly
    
    γz = gamma(T, max_nonzeros_per_row(Q.L))
    γn = gamma(T, n+3) # not n+2 like in the paper, because we wish to allow for f to be the result of rounding
    
    
    ϵ = zero(T)

    nrmM = opnormbound(N, M)

    # these need to be implemented for each kernel
    δₖ = opradius(N, MK) #this is the opnorm of the radius matrix
    
    # essentially nonzero_per_row(MK) is
    # the  number of mul_add needed to obtain one
    # entry of MK*v
    γₖ = gamma(T, nonzero_per_row(MK)) 
    
    nrm_MK = opnormbound(N, MK)
    normMK = nrm_MK ⊕₊ δₖ

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
            
            # the noise step
            nrmv = opnormbound(N, v)
            
            v = w
            w = MK * v
            v = w
            
            ϵ = (γₖ ⊗₊nrm_MK ⊕₊δₖ) ⊗₊ nrmv ⊕₊ normMK ⊗₊ ϵ
            
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
function norms_of_powers_trivial_noise(N::Type{<:NormKind}, 
                                       Q::DiscretizedOperator, 
                                       MK::NoiseKernel,
                                       m::Integer)
    norms = fill(NaN, m)
    δₖ = opradius(N, MK)
    nrm_MK = opnormbound(N, MK)
    norms[1] = opnormbound(N, Q)⊗₊ nrm_MK⊕₊δₖ
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
function norms_of_powers_abstract_noise(Bas::Basis, N::NoiseKernel, m)
    A, B = dfly(strong_norm(Bas), aux_norm(Bas), N)
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
Estimate norms of powers from those on a coarser grid (see paper for details)
"""
function norms_of_powers_from_coarser_grid_noise(fine_basis::Basis, coarse_basis::Basis, var_ρ::Real, coarse_norms::Vector, normQ::Real)
    if !BasisDefinition.is_refinement(fine_basis, coarse_basis)
        @error "The fine basis is not a refinement of the coarse basis"
    end
    m = length(coarse_norms)
    fine_norms = fill(NaN, m+1)
    fine_norms[1] = 1.0

    Kh =  BasisDefinition.weak_projection_error(coarse_basis)
    
    coarse_norms0(k::Integer) = k==0 ? 1. : coarse_norms[k]

    for k in 1:m
        temp = 0.0
        for j in 0:k-1
            temp = temp ⊕₊ coarse_norms0(j) ⊗₊ var_ρ
        end
        temp = 2*(temp ⊕₊ 1.0)
        fine_norms[k+1] = coarse_norms[k] ⊕₊ Kh ⊗₊ temp
    end

    return fine_norms
end

function powernormboundsnoise(B; Q=DiscretizedOperator(B, D), NK = NK::NoiseKernel)
	m = 8
	computed_norms = []
	while true
		computed_norms = norms_of_powers_noise(weak_norm(B), m, Q, NK, integral_covector(B))
		if any(computed_norms .< 0.1)
			break
		end
		m = 2*m
	end
	trivial_norms = norms_of_powers_trivial_noise(weak_norm(B), Q, NK, m)
	(dfly_strongs, dfly_norms) = norms_of_powers_abstract_noise(B, NK, m)
	# in the current version, dfly_norms seem to be always larger and could be omitted
	# however they do not cost much to compute
	norms = min.(trivial_norms, computed_norms, dfly_norms)

	m_extend = 2*m
	better_norms = []
	while true
		better_norms = refine_norms_of_powers(norms, m_extend)
		if better_norms[end] < 1e-8
			break
		end
		m_extend = 2*m_extend
	end

	return better_norms

end
