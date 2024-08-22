"""
Functions to estimate Q|_{U^0}. See our paper for details.
"""

using LinearAlgebra
using SparseArrays
using FastRounding

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
function norms_of_powers_noise(
        B::Ulam,
        N::Type{L1}, 
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
    δ = opnormbound(B, N, R)
    
    # this is the operator radius of the noise operator,
    # we use this function since the operator may be defined implictly
    
    γz = gamma(T, max_nonzeros_per_row(Q.L))
    γn = gamma(T, n+3) # not n+2 like in the paper, because we wish to allow for f to be the result of rounding
    
    
    ϵ = zero(T)

    nrmM = opnormbound(B, N, M)

    # these need to be implemented for each kernel
    δₖ = opradius(N, MK) #this is the opnorm of the radius matrix
    
    # essentially nonzero_per_row(MK) is
    # the  number of mul_add needed to obtain one
    # entry of MK*v
    γₖ = gamma(T, nonzero_per_row(MK)) 
    
    nrm_MK = opnormbound(B, N, MK)
    normMK = nrm_MK ⊕₊ δₖ

    # precompute norms
    if !is_integral_preserving(Q)
        if normE == -1.
            normE = opnormbound(B, N, Q.e)
        end
        if normEF == -1.
            normEF = opnormbound(B, N, Q.e*f)
        end
        if normIEF == -1.
            normIEF =  opnormbound(B, N, [Matrix(UniformScaling{Float64}(1),n,n) Q.e*f])
        end
        if normN == -1.
            normN = opnormbound(B, N, Matrix(UniformScaling{Float64}(1),n,n) - Q.e*f)
        end
    end

    if normQ == -1.
        if is_integral_preserving(Q)
            normQ = nrmM ⊕₊ δ
        else
            defect = opnormbound(B, N, Q.w)
            normQ = nrmM ⊕₊ δ ⊕₊ normE ⊗₊ defect
        end
    end

    # initialize normcachers
    normcachers = [NormCacher{N}(B, n) for j in 1:m]

    midf = map(mid, f)

    # main loop

    v = zeros(T, n)

    @showprogress for j = 1:n-1
        v .= zero(T) # TODO: check for type stability in cases with unusual types
        v[1] = one(T) # TODO: in full generality, this should contain entries of f rather than ±1
        v[j+1] = -one(T)
        if normv0 == -1.
            nrmv = opnormbound(B, N, v)
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
                new_nrmw = opnormbound(B, N, w)
                ϵ = γn ⊗₊ normIEF ⊗₊ (new_nrmw ⊕₊ normEF ⊗₊ nrmw) ⊕₊ normN ⊗₊ (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
                nrmw = new_nrmw
            end
            
            # the noise step
            nrmv = opnormbound(B, N, v)
            
            v = w
            w = MK * v
            v = w
            
            ϵ = (γₖ ⊗₊nrm_MK ⊕₊δₖ) ⊗₊ nrmv ⊕₊ normMK ⊗₊ ϵ
            
            nrmv = opnormbound(B, N, v)
            add_column!(normcachers[k], v, ϵ) #TODO: Could pass and reuse nrmv in the case of norm-1
        end
    end
    return map(get_norm, normcachers)
end

"""
Array of "trivial" bounds for the powers of a DiscretizedOperator (on the whole space)
coming from from ||Q^k|| ≤ ||Q||^k
"""
function norms_of_powers_trivial_noise(B::Basis,
                                       N::Type{<:NormKind}, 
                                       Q::DiscretizedOperator, 
                                       MK::NoiseKernel,
                                       m::Integer)
    norms = fill(NaN, m)
    δₖ = opradius(N, MK)
    nrm_MK = opnormbound(B, N, MK)
    norms[1] = opnormbound(B, N, Q) ⊗₊ nrm_MK⊕₊δₖ
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
TODO: Check if it works for other basis types
"""
function norms_of_powers_from_coarser_grid_noise(fine_basis::Ulam, 
                                            coarse_basis::Ulam, 
                                            Q::DiscretizedOperator,
                                            NK::NoiseKernel, 
                                            coarse_norms::Vector)
    if !BasisDefinition.is_refinement(fine_basis, coarse_basis)
        @error "The fine basis is not a refinement of the coarse basis"
    end
    m = length(coarse_norms)
    fine_norms = fill(NaN, m+1)
    trivial_norms = norms_of_powers_trivial_noise(fine_basis, weak_norm(fine_basis), Q, NK, m+1)

    A, B = dfly(strong_norm(fine_basis), aux_norm(fine_basis), NK)

    # adds a 0th element to strongs
    trivial_norms0(k::Integer) = k==0 ? 1. : trivial_norms[k]
    coarse_norms0(k::Integer) = k==0 ? 1. : coarse_norms[k]

    Kh =  BasisDefinition.weak_projection_error(coarse_basis)
    
    fine_norms[1] = trivial_norms0(1)

    for k in 1:m
		temp = 0.
		for i in 0:k-1
			temp = temp ⊕₊ coarse_norms0(i) ⊗₊ (trivial_norms0(k-i) ⊕₊ trivial_norms0(k-i-1))
		end
        fine_norms[k+1] = coarse_norms0(k) ⊕₊ B ⊗₊ Kh ⊗₊ (temp ⊕₊ 1.0)  
	end
    return fine_norms
end

function norms_of_powers_from_coarser_grid_noise_abstract(coarse_basis::Ulam, 
                                                        NK::NoiseKernel,
                                                        coarse_norms::Vector)
    
    m = length(coarse_norms)
    abstract_norms = fill(NaN, m+1)
    trivial_norms = ones(m+1)

    A, B = dfly(strong_norm(coarse_basis), aux_norm(coarse_basis), NK)

    # adds a 0th element to strongs
    trivial_norms0(k::Integer) = k==0 ? 1. : trivial_norms[k]
    coarse_norms0(k::Integer) = k==0 ? 1. : coarse_norms[k]

    Kh =  BasisDefinition.weak_projection_error(coarse_basis)

    fine_norms = fill(NaN, m+1)
    fine_norms[1] = trivial_norms0(1)

    for k in 1:m
        temp = 0.
        for i in 0:k-1
            temp = temp ⊕₊ coarse_norms0(i) ⊗₊ (trivial_norms0(k-i) ⊕₊ trivial_norms0(k-i-1))
        end
    fine_norms[k+1] = coarse_norms0(k) ⊕₊ B ⊗₊ Kh ⊗₊ (temp ⊕₊ 1.0)  
    end
    return fine_norms
end


function powernormboundsnoise(B; Q=DiscretizedOperator(B, D), NK = NK::NoiseKernel, threshold = 0.1)
	m = 8
	computed_norms = []
	while true
		computed_norms = norms_of_powers_noise(B, weak_norm(B), m, Q, NK, integral_covector(B))
		if any(computed_norms .< threshold)
			break
		end
		m = 2*m
	end
	trivial_norms = norms_of_powers_trivial_noise(B, weak_norm(B), Q, NK, m)
	norms = min.(trivial_norms, computed_norms)

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

"""
Uses power norm bounds already computed for a coarse operator to estimate
the same norms for a finer operator
"""
function finepowernormboundsnoise(B, 
                                  Bfine,  
                                  coarse_norms; 
                                  Qfine::DiscretizedOperator,
                                  NKfine::NoiseKernel)
	                              m = length(coarse_norms)

	trivial_norms_fine = norms_of_powers_trivial_noise(Bfine, weak_norm(B), Qfine, NKfine, m+1)
	twogrid_norms_fine = norms_of_powers_from_coarser_grid_noise(Bfine, 
                                                                 B, 
                                                                 Qfine,
                                                                 NKfine, 
                                                                 coarse_norms)
    
	norms_fine = min.(trivial_norms_fine, twogrid_norms_fine)

	better_norms_fine = refine_norms_of_powers(norms_fine, m+1)
	return better_norms_fine
end

function abstractpowernormboundsnoise(B, NK, coarse_norms; m = length(coarse_norms))
    trivial_norms_fine = ones(m+1)
    abstract_norms_from_coarse = norms_of_powers_from_coarser_grid_noise_abstract(B, NK, coarse_norms)

    norms_fine = min.(trivial_norms_fine, abstract_norms_from_coarse)
    better_norms_abstract = refine_norms_of_powers(norms_fine, m+1)
    return better_norms_abstract
end

"""
Return a numerical approximation to the (hopefully unique) invariant vector
of the dynamic with discretized operator Q.

The vector is normalized so that integral_covector(B)*w ≈ 1
"""
function invariant_vector_noise(B::Basis, Q::DiscretizedOperator, NK::NoiseKernel; tol = 0.0, iter =10)
	mQ = mid(Q)
	v = one_vector(B)
    for _ in 1:iter
        v = mQ*v
        v = NK*v
    end
    v = v ./ (mid.(integral_covector(B))*v)
    return v
end

"""
Return an upper bound to Q_h*w - w in the given norm
"""
function residualboundnoise(B::Basis, N::Type{<:NormKind}, Q::DiscretizedOperator, NK::NoiseKernel, w::AbstractVector)
	v = Q*w
    v = NK*v # actually, a priori estimate... TODO modify mult NK*v
    return normbound(B, N, v - w)
end

"""
Bounds rigorously the distance of w from the fixed point of Q (normalized with integral = 1),
using a vector of bounds norms[k] ≥ ||Q_h^k|_{U_h^0}||.
"""
function distance_from_invariant_noise(B::Basis, 
                                 Q::DiscretizedOperator, 
                                 NK::NoiseKernel, 
                                 w::AbstractVector, 
                                 norms::Vector; 
                                 ε₁::Float64 = residualboundnoise(B, weak_norm(B), Q, NK, w), 
                                 ε₂::Float64 = mag(integral_covector(B) * w - 1), 
                                 normQ::Float64 = opnormbound(B, weak_norm(B), Q))
	if ε₂ > 1e-8
		@error "w does not seem normalized correctly"
	end
    _ , us = dfly(strong_norm(B), aux_norm(B), NK)
	Cs = infinite_sum_norms(norms)
	Kh =  BasisDefinition.weak_projection_error(B)
	normw = normbound(B, weak_norm(B), w)

	return Cs ⊗₊ (2. ⊗₊ Kh ⊗₊ (1. ⊕₊ normQ) ⊗₊ us ⊕₊ ε₁ ⊘₊ (1. ⊖₋ ε₂)) ⊕₊ ε₂ ⊘₊ (1. ⊖₋ ε₂) ⊗₊ normw
end