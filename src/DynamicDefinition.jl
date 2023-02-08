"""
defines generic Dynamic type
"""

module DynamicDefinition
export Dynamic, MarkovDynamic, preim, branches, nbranches, plottable, is_full_branch, domain, endpoints, branch, max_distortion, max_expansivity

abstract type Dynamic end
abstract type MarkovDynamic <: Dynamic end
using IntervalArithmetic, IntervalOptimisation

domain(S::Dynamic) = @error "Not implemented"
nbranches(S::Dynamic) = @error "Not implemented"
branch(S::Dynamic, k) = @error "Not implemented"
branches(S::Dynamic, k) = @error "Not implemented"

"""
Return a non-interval version of the map as a function. This can be used, for instance, for plot(plottable(D)).
"""
plottable(S::Dynamic) = @error "Not implemented"

"""
	preim(S::Dynamic, k, y, ϵ)

Computes the preim of y in branch k of a dynamic, with accuracy ϵ
"""
function preim end

is_full_branch(S::Dynamic) = @error "Not implemented"

"""
	endpoints(S::Dynamic)

Endpoints of the branches, in increasing order (returned as a vector of intervals)
"""
endpoints(S::Dynamic) = @error "Not implemented"

"""
Maximum of |1/T'|
"""
function max_expansivity(D::Dynamic, tol=1e-3)
	@error "Not implemented"
	# v = endpoints(D)
	# # due to the fact that D(x::Taylor1 ) is defined, this calls all the right methods 
	# # the call to the branch method was superfluous
	# return maximum(maximise(x -> abs(1/derivative(D, x)), hull(v[k], v[k+1]), tol=tol)[1] for k in 1:nbranches(D))
end

"""
Maximum of distortion(D, x) = |T''| / (T')^2, over all branches
"""
function max_distortion(D::Dynamic, tol=1e-3)
	@error "Not implemented"
	# v = endpoints(D)
	# # due to the fact that D(x::Taylor1 ) is defined, this calls all the right methods 
	# # the call to the branch method was superfluous
	# return maximum(maximise(x -> distortion(D, x), hull(v[k], v[k+1]), tol=tol)[1] for k in 1:nbranches(D))
end

end
