"""
defines generic Dynamic type
"""

module DynamicDefinition

using ValidatedNumerics
using TaylorSeries:Taylor1
export Dynamic, MarkovDynamic, preim, nbranches, plottable, is_full_branch, domain, derivative, distorsion, endpoints, branch, expansivity, max_distorsion

abstract type Dynamic end
abstract type MarkovDynamic <: Dynamic end

domain(S::Dynamic) = @error "Not implemented"
nbranches(S::Dynamic) = @error "Not implemented"
branch(S::Dynamic, k) = @error "Not implemented"
plottable(S::Dynamic) = @error "Not implemented"
preim(S::Dynamic, k, y, ϵ) = @error "Not implemented"
is_full_branch(S::Dynamic) = @error "Not implemented"

"""
Endpoints of the branches, in increasing order (returned as a vector of intervals)
"""
endpoints(S::Dynamic) = @error "Not implemented"

# Derivative and distorsion of a generic function (*not* a dynamic). Here for convenience,

# the isempty check is required because otherwise derivative(x -> 4*x, ∅) == 4.
"""
Nth derivative of a function (or a dynamic)
"""
derivative(f, x) = derivative(1, f, x)
derivative(n, f, x) = isempty(x) ?  ∅ : f(Taylor1([x, 1], n))[n] * factorial(n)

"""
Distorsion of a function (or a dynamic), i.e., |f′′ / f′^2|
"""
function distorsion(f, x)
	if isempty(x)
		return ∅
	end
	series = f(Taylor1([x, 1], 2))
	f′ = series[1]
	f′′ = 2*series[2]
	return abs(f′′ / f′^2)
end

"""
Maximum of |1/T'|
"""
function expansivity(D::Dynamic, tol=1e-8)
	v = endpoints(D)
	return maximum(maximise(x -> abs(1/derivative(branch(D,k), x)), hull(v[k], v[k+1]), tol=tol)[1] for k in 1:nbranches(D))
end

"""
Maximum of distorsion(D, x) = |T''| / (T')^2, over all branches
"""
function max_distorsion(D::Dynamic, tol=1e-8)
	v = endpoints(D)
	return maximum(maximise(x -> distorsion(branch(D, k) , x), hull(v[k], v[k+1]), tol=tol)[1] for k in 1:nbranches(D))
end

# these do not work properly, since T() for a Mod1Dynamic contains the unquotiented map. Better not to use them at all.

# function iterate(T, x, ::Val{n}) where {n}
# 	for i in 1:n
# 		x=T(x)
# 	end
# 	return x
# end
#
# iterate(T, x, n)=iterate(T, x, Val(n))

end
