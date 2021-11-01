"""
defines generic Dynamic type
"""

module DynamicDefinition
export Dynamic, MarkovDynamic, preim, nbranches, plottable, is_full_branch, domain, derivative, distorsion, endpoints, branch, expansivity, max_distorsion, orientation

abstract type Dynamic end
abstract type MarkovDynamic <: Dynamic end
using IntervalArithmetic, IntervalOptimisation

domain(S::Dynamic) = @error "Not implemented"
nbranches(S::Dynamic) = @error "Not implemented"
branch(S::Dynamic, k) = @error "Not implemented"
plottable(S::Dynamic) = @error "Not implemented"

"""
preim(S::Dynamic, k, y, ϵ)

Computes the preim of y in branch k of a dynamic, with accuracy ϵ
"""
function preim end

is_full_branch(S::Dynamic) = @error "Not implemented"

"""
Endpoints of the branches, in increasing order (returned as a vector of intervals)
"""
endpoints(S::Dynamic) = @error "Not implemented"

# Derivative and distorsion of a generic function (*not* a dynamic). Here for convenience,
# the isempty check is required because otherwise derivative(x -> 4*x, ∅) == 4.

#
import TaylorSeries
"""
Nth derivative of a function (or a dynamic)
"""
derivative(f, x) = derivative(1, f, x)
derivative(n, f, x) = isempty(x) ?  ∅ : f(TaylorSeries.Taylor1([x, 1], n))[n] * factorial(n)

"""
Distorsion of a function (or a dynamic), i.e., |f′′ / f′^2|
"""
function distorsion(f, x)
	if isempty(x)
		return ∅
	end
	series = f(TaylorSeries.Taylor1([x, 1], 2))
	f′ = series[1]
	f′′ = 2*series[2]
	return abs(f′′ / f′^2)
end

"""
Maximum of |1/T'|
"""
function expansivity(D::Dynamic, tol=1e-3)
	v = endpoints(D)
	# due to the fact that D(x::Taylor1 ) is defined, this calls all the right methods 
	# the call to the branch method was superfluous
	return maximum(maximise(x -> abs(1/derivative(D, x)), hull(v[k], v[k+1]), tol=tol)[1] for k in 1:nbranches(D))
end

"""
Maximum of distorsion(D, x) = |T''| / (T')^2, over all branches
"""
function max_distorsion(D::Dynamic, tol=1e-3)
	v = endpoints(D)
	# due to the fact that D(x::Taylor1 ) is defined, this calls all the right methods 
	# the call to the branch method was superfluous
	return maximum(maximise(x -> distorsion(D, x), hull(v[k], v[k+1]), tol=tol)[1] for k in 1:nbranches(D))
end

"""
orientation(D, k)

Orientation of branch k: 1. for increasing, -1. for decreasing
"""
function orientation end

end
