export Dynamic,
    MarkovDynamic,
    preim,
    branches,
    nbranches,
    plottable,
    is_full_branch,
    domain,
    endpoints,
    branch,
    max_distortion,
    max_inverse_derivative,
    is_increasing

abstract type Dynamic end
abstract type MarkovDynamic <: Dynamic end

function domain end
function nbranches end
function branch end
function branches end
function is_increasing end

"""
Return a non-interval version of the map as a function. This can be used, for instance, for plot(plottable(D)).
"""
function plottable end

@doc raw"""
	preim(S::Dynamic, k, y, ϵ)

Computes the preim of y in branch k of a dynamic, with accuracy ϵ
"""
function preim end

function is_full_branch end

"""
	endpoints(S::Dynamic)

Endpoints of the branches, in increasing order (returned as a vector of intervals)
"""
function endpoints end

@doc raw"""
Maximum of ``|1/T'|``
"""
function max_inverse_derivative end

@doc raw"""
Maximum of distortion(D, x) = ``|T''| / (T')^2``, over all branches
"""
function max_distortion end
