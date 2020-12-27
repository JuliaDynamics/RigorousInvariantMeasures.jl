"""
Generic function to estimate decay time on a single vector

Args:
	PP (scipy matrix): approximation of the discretized matrix to use
	v (numpy vector):
	my_norm (function):
	M, K (reals with RNDU): constants such that ||P-PP|| ≤ K` and ||P^n|| ≤ M`,
	where P is the exact discretized matrix.

Returns:
	n such that ||PP^n v|| ≤ target`.

Raises:
	ValueError if insufficient precision is detected
"""
function vector_decay_time(B::Basis, PP, v, M=0, K=0, target = 0.5)
	current_norm = weak_norm(B, v)
	error_on_computed_norm = zero(M)
	n = 0
	MK = M * K
	error_propagation_constant = K #this constant is K at the first step and MK afterwards; see notes
	while current_norm + error_on_computed_norm >= target
		n += 1
		v = PP * v
		current_norm = weak_norm(B, v)
		error_on_computed_norm += error_propagation_constant * current_norm
		error_propagation_constant = MK
		if error_on_computed_norm > target
			raise() # ValueError, 'Insufficient precision'
		end
	end
	return n
end

"""
Number of iterations needed to contract all vectors in `basis.contracting_pairs()` to a given target alpha
"""
function decay_time(D::Dynamic, B::Basis, P::AbstractMatrix{Interval{T}}, alpha = 0.5, n_jobs = 1) where {T}

	PP = mid.(P)

	M = bound_on_norms_of_powers(B, D, project_left=True, project_right=True)
	K = numerical_error(B, D, P, PP)
	my_norm(v) = weak_norm(B, v)

	#alpha = Rdown(alpha) # we need a lower bound to alpha*s to make sure we contract below it
	#decay_times = Parallel(n_jobs=n_jobs, verbose=1)(delayed(vector_decay_time)(PP, v, my_norm, M, K, alpha*s) for v, s in basis.contracting_pairs())

	return max(decay_times)
end
