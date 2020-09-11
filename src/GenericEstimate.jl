using LinearAlgebra, Arpack, FastRounding, ValidatedNumerics
using ..DynamicDefinition, ..BasisDefinition

export invariant_vector

"""
Return the (hopefully unique) invariant vector of the dynamic with discretized operator P.

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
function residualbound(N::Type{<:NormKind}, Q::DiscretizedOperator, w::AbstractVector)
	return normbound(N, Q*w - w)
end

"""
Bounds rigorously the distance of w from the fixed point of Q (normalized with integral = 1),
using a vector of bounds norms[k] ≥ ||Q_h^k|_{U_h^0}||.
"""
function distance_from_invariant(B::Basis, D::Dynamic, Q::DiscretizedOperator, w::AbstractVector, norms::Vector; ε₁::Float64 = residualbound(weak_norm(B), Q, w), ε₂::Float64 = mag(integral_covector(B) * w - 1), normQ::Float64 = opnormbound(weak_norm(B), Q))
	if ε₂ > 1e-8
		@error "w does not seem normalized correctly"
	end
	us = BasisDefinition.invariant_measure_strong_norm_bound(B, D)
	Cs = infinite_sum_norms(norms)
	Kh =  BasisDefinition.weak_projection_error(B)
	normw = normbound(weak_norm(B), w)

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
