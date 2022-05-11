using LinearAlgebra, Arpack, FastRounding, ValidatedNumerics
using ..DynamicDefinition, ..BasisDefinition

export invariant_vector, finepowernormbounds, powernormbounds, distance_from_invariant, 
	compute_coarse_grid_quantities, compute_fine_grid_quantities, one_grid_estimate, two_grid_estimate

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
If ε₁ and normQ are given, then Q can be omitted
"""
function distance_from_invariant(B::Basis, D::Dynamic, Q::Union{DiscretizedOperator,Nothing}, w::AbstractVector, norms::Vector;
	ε₁::Float64 = residualbound(B, weak_norm(B), Q, w), ε₂::Float64 = mag(integral_covector(B) * w - 1), dfly_coefficients=dfly(strong_norm(B), aux_norm(B), D))
	if ε₂ > 1e-8
		@error "w does not seem normalized correctly"
	end
	us = BasisDefinition.invariant_measure_strong_norm_bound(B, D; dfly_coefficients=dfly_coefficients)
	Cs = infinite_sum_norms(norms)
	Kh =  BasisDefinition.weak_projection_error(B)
	normw = normbound(B, weak_norm(B), w)
	normL = BasisDefinition.bound_weak_norm_abstract(B, D; dfly_coefficients=dfly_coefficients)

	return Cs ⊗₊ (2. ⊗₊ Kh ⊗₊ (1. ⊕₊ normL) ⊗₊ us ⊕₊ ε₁ ⊘₊ (1. ⊖₋ ε₂)) ⊕₊ ε₂ ⊘₊ (1. ⊖₋ ε₂) ⊗₊ normw
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

"""
Computes bounds for norms of powers, taking (optionally) minimum values for the number of norms to compute
"""
function powernormbounds(B, D; Q=DiscretizedOperator(B, D), m=8)
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
function finepowernormbounds(B, B_fine, D, coarse_norms; normQ_fine=opnormbound(B_fine, weak_norm(B_fine)DiscretizedOperator(B_fine, D)), dfly_coefficients=dfly(strong_norm(B_fine), aux_norm(B_fine), D))
	m = length(coarse_norms)
	
	trivial_norms_fine = norms_of_powers_trivial(normQ_fine, m)
	twogrid_norms_fine = norms_of_powers_from_coarser_grid(B_fine, B, D, coarse_norms, normQ_fine; dfly_coefficients=dfly_coefficients)

	(dfly_strongs_fine, dfly_norms_fine) = norms_of_powers_dfly(B_fine, D, m; dfly_coefficients=dfly_coefficients)

	norms_fine = min.(trivial_norms_fine, twogrid_norms_fine, dfly_norms_fine)

	better_norms_fine = refine_norms_of_powers(norms_fine, m)
	return better_norms_fine
end

"""
Struct that encapsulates all the quantities computed from the fine basis that are needed in the two-grid estimate.
It is meant as an intermediate quantity that can be saved on the disk to avoid recomputing Q all the times
"""
struct FineGridQuantities
	B::Basis
	D::Dynamic
	normQ
	w
	ε₁
	ε₂
	time_assembling # time to compute B, D, Q, normQ
	time_eigen      # time to compute w, ε₁, ε₂
end
"""
Struct that encapsulates the additional quantities needed on the coarse basis for a two-grid estimate,
or on the (only) basis for a one-grid estimate. 
It is meant as an intermediate quantity that can be saved on the disk to avoid recomputing Q all the times.
"""
struct CoarseGridQuantities
	B::Basis
	D::Dynamic
	norms
	dfly_coefficients
	time_assembling # time to compute B, D, Q (without normQ, needed in the two-grid estimate)
	time_norms      # time to compute norms
	time_dfly       # time to compute dfly_coefficients
end

"""
Compute FineGridQuantities, given a function f(n) that computes B, D, Q = f(n)
"""
function compute_fine_grid_quantities(f, n)
	time_assembling1 = @elapsed B, D, Q = f(n)
	time_assembling2 = @elapsed normQ = opnormbound(B, weak_norm(B), Q)
	time_eigen1 = @elapsed w = invariant_vector(B, Q)
	time_eigen2 = @elapsed ε₁, ε₂ = residualbound(B, weak_norm(B), Q, w), mag(integral_covector(B) * w - 1)
	return FineGridQuantities(B, D, normQ, w, ε₁, ε₂, time_assembling1 + time_assembling2, time_eigen1 + time_eigen2)
end

"""
Compute FineGridQuantities _and_ CoarseGridQuantities, given a function f(n) that computes B, D, Q = f(n)
"""
function compute_coarse_grid_quantities(f, n; m = 8)
	time_assembling1 = @elapsed B, D, Q = f(n)
	time_assembling2 = @elapsed normQ = opnormbound(B, weak_norm(B), Q)
	time_eigen1 = @elapsed w = invariant_vector(B, Q)
	time_eigen2 = @elapsed ε₁, ε₂ = residualbound(B, weak_norm(B), Q, w), mag(integral_covector(B) * w - 1)
	time_norms = @elapsed norms = powernormbounds(B, D, Q=Q, m=m)
    time_dfly = @elapsed dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D)
	return FineGridQuantities(B, D, normQ, w, ε₁, ε₂, time_assembling1 + time_assembling2, time_eigen1 + time_eigen2),
		CoarseGridQuantities(B, D, norms, dfly_coefficients, time_assembling1, time_norms, time_dfly)
end

"""
Compute a one-grid error estimate.

The first return argument is the error, the second is the time breakdown according to ["dfly", "assembling", "eigen", "norms", "error"]. 
(The sum of that vector is the total time taken)
"""
function one_grid_estimate(C::CoarseGridQuantities, F::FineGridQuantities)
	@assert(C.B==F.B)
	# Dynamics don't compare unfortunately
	time_error = @elapsed error = distance_from_invariant(F.B, F.D, nothing, F.w, C.norms; dfly_coefficients=C.dfly_coefficients, ε₁=F.ε₁, ε₂=F.ε₂)
	return error, [C.time_dfly, F.time_assembling, F.time_eigen, C.time_norms, time_error]
end

"""
Compute a two-grid error estimate.

The first return argument is the error, the second is the time breakdown according to ["dfly", "coarse", "assembling", "eigen", "norms", "error"]. 
(The sum of that vector is the total time taken)
"""
function two_grid_estimate(C::CoarseGridQuantities, F::FineGridQuantities; m_extend=400)
	@assert(is_refinement(F.B, C.B))
	# Dynamics don't compare unfortunately
	time_error_fine = @elapsed error_fine = distance_from_invariant(F.B, F.D, nothing, F.w, 
		finepowernormbounds(C.B, F.B, F.D, refine_norms_of_powers(C.norms, m_extend); normQ_fine=F.normQ, dfly_coefficients=C.dfly_coefficients);
		dfly_coefficients=C.dfly_coefficients, ε₁=F.ε₁, ε₂=F.ε₂)
	return error_fine, [C.time_dfly, C.time_assembling+C.time_norms, F.time_assembling, F.time_eigen, time_error_fine]
end
