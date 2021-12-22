using LinearAlgebra, Arpack, FastRounding, ValidatedNumerics
using ..DynamicDefinition, ..BasisDefinition

export invariant_vector, finepowernormbounds, powernormbounds, distance_from_invariant

"""
Return a numerical approximation to the (hopefully unique) invariant vector
of the dynamic with discretized operator Q.

The vector is normalized so that integral_covector(B)*w ≈ 1
"""
function invariant_vector(B::Basis, Q::DiscretizedOperator; tol = 0.0)
	mQ = mid(Q)
	n = size(Q)[1]
	# setting a larger nev seems to slow things down
	F = eigs(mQ; tol=tol, nev=1, ritzvec=true, v0=ones((n,)))
	w = F[2][:, 1]
	@assert imag(w) ≈ zeros(n)
	w = real(w) # this seems a pretty safe assumption.
	# In the Ulam case, in principle we could enforce w >= 0, but in practice
	# it will hardly ever be relevant.
	w = w ./ (mid.(integral_covector(B))*w) #normalization
	return w
end

"""
Return an upper bound to Q_h*w - w in the given norm
"""
function residualbound(B::Basis, N::Type{<:NormKind}, Q::DiscretizedOperator, w::AbstractVector)
	return normbound(B, N, Q*w - w)
end

"""
Bounds rigorously the distance of w from the fixed point of Q (normalized with integral = 1),
using a vector of bounds norms[k] ≥ ||Q_h^k|_{U_h^0}||.
"""
function distance_from_invariant(B::Basis, D::Dynamic, Q::DiscretizedOperator, w::AbstractVector, norms::Vector; ε₁::Float64 = residualbound(B, weak_norm(B), Q, w), ε₂::Float64 = mag(integral_covector(B) * w - 1), normQ::Float64 = opnormbound(B, weak_norm(B), Q))
	if ε₂ > 1e-8
		@error "w does not seem normalized correctly"
	end
	us = BasisDefinition.invariant_measure_strong_norm_bound(B, D)
	Cs = infinite_sum_norms(norms)
	Kh =  BasisDefinition.weak_projection_error(B)
	normw = normbound(B, weak_norm(B), w)

	return Cs ⊗₊ (2. ⊗₊ Kh ⊗₊ (1. ⊕₊ normQ) ⊗₊ us ⊕₊ ε₁ ⊘₊ (1. ⊖₋ ε₂)) ⊕₊ ε₂ ⊘₊ (1. ⊖₋ ε₂) ⊗₊ normw
end

# """
# This function returns a sequence of Cᵢ, \\tilde{C}ᵢ for a matrix P
# on a subspace V such that ||P^i|_V||_1\\leq C_i and
# ||P^i|_V||_{\\infty}\\leq \\tilde{C}_i, with respect to the
# """
# function contractmatrix(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
# 	# vector of the Cᵢ for P
# 	C = zeros(m)
# 	S = zeros((length(B), m))
# 	tilde_C = zeros(m)
#
# 	PP = mid.(P)
#
# 	for v in BasisDefinition.AverageZero(B)
# 		λ₁, λ₂ = norm(v, 1), norm(v, Inf)
# 		for i in 1:m
# 			v = PP*v
# 			C[i] = max(C[i], norm(v, 1)/λ₁)
# 			S[:, i]+=abs.(v)/λ₂
# 		end
# 	end
#
#
# 	for k in 1:m
# 		tilde_C[k] = maximum(S[:,k])
# 	end
#
# 	# we keep track of the error due to the basis we have chosen
# 	η₁ = BasisDefinition.spaceconstant(B, Val(:L1))
# 	η₂ = BasisDefinition.spaceconstant(B, Val(:L∞))
#
#  	return η₁*C, η₂*tilde_C
# end
#
# @deprecate contractmatrix norms_of_powers

"""
This function returns the bound on the weak norm of the discretized operator
"""
function boundnorm(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
	W₁, W₂ = BasisDefinition.bound_weak_norm_from_linalg_norm(B)
	α₁ = BasisDefinition.bound_linalg_norm_L1_from_weak(B)
	α₂ = BasisDefinition.bound_linalg_norm_L∞_from_weak(B)
	C, tilde_C = contractmatrix(B, P, m)
	return (W₁/α₁)*C+(W₂/α₂)*tilde_C
end

"""
Uses different strategies to compute power norm bounds.

If specified, `m` norms of powers are estimated computationally, and then
`m_extend` norms are obtained with a cheaper refinement process. Otherwise
these numbers are selected automatically.

A vector of length m_extend is returned, such that norms[k] ≥ ||Q_h^k|_{U_h^0}||
"""
function powernormbounds(B, D, m, m_extend; Q=DiscretizedOperator(B, D))
	normQ = opnormbound(B, weak_norm(B), Q)
	trivial_norms = norms_of_powers_trivial(normQ, m)
	computed_norms = norms_of_powers(weak_norm(B), m, Q, integral_covector(B))

	# not interesting at the moment
	#(dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, m)
	# in the current version, dfly_norms seem to be always larger and could be omitted
	# however they do not cost much to compute
	#norms = min.(trivial_norms, computed_norms, dfly_norms)
	norms = min.(trivial_norms, computed_norms)

	better_norms = refine_norms_of_powers(norms, m_extend)

	return better_norms
end

function powernormbounds(B, D; Q=DiscretizedOperator(B, D))
	m = 8
	computed_norms = []
	while true
		computed_norms = norms_of_powers(B, weak_norm(B), m, Q, integral_covector(B))
		if any(computed_norms .< 0.1)
			break
		end
		m = 2*m
	end
	normQ = opnormbound(B, weak_norm(B), Q)
	trivial_norms = norms_of_powers_trivial(normQ, m)
	# (dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, m)
	# in the current version, dfly_norms seem to be always larger and could be omitted
	# however they do not cost much to compute
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
function finepowernormbounds(B, B_fine, D, coarse_norms; Q_fine=DiscretizedOperator(B_fine, D))
	m = length(coarse_norms)

	norm_Q_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
	
	trivial_norms_fine = norms_of_powers_trivial(norm_Q_fine, m)
	twogrid_norms_fine = norms_of_powers_from_coarser_grid(B_fine, B, D, coarse_norms, norm_Q_fine)

	(dfly_strongs_fine, dfly_norms_fine) = norms_of_powers_dfly(B_fine, D, m)

	norms_fine = min.(trivial_norms_fine, twogrid_norms_fine, dfly_norms_fine)

	better_norms_fine = refine_norms_of_powers(norms_fine, m)
	return better_norms_fine
end
