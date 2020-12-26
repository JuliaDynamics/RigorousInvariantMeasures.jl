module PwDynamicDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors
using TaylorSeries: Taylor1

export PwMap, preim, nbranches, plottable

"""
Dynamic based on a piecewise defined map.

The map is defined as T(x) = Ts[k](x) if x ∈ [endpoints[k], endpoints[k+1]).

We set is_full[k]==true if we can prove/assume that the kth branch is full, i.e., Ts[k](endpoints[k])==0 and Ts[k](endpoints[k+1])==1 or vice versa.
"""
struct PwMap <: Dynamic
	Ts::Array{Function, 1}
	endpoints::Array{Interval, 1}
	is_full
	orientations # these will be filled in automatically, usually
end


# the isempty check is required because otherwise derivative(x -> 4*x, ∅) == 4.
derivative(f, x) = isempty(x) ?  ∅ : f(Taylor1([x, 1], 1))[1]

PwMap(Ts, endpoints) = PwMap(Ts, endpoints, fill(false, length(endpoints)-1))
PwMap(Ts, endpoints, is_full) = PwMap(Ts, map(Interval, endpoints), is_full, [derivative(Ts[i], (endpoints[i]+endpoints[i+1])/2) for i in 1:length(endpoints)-1])

DynamicDefinition.nbranches(D::PwMap)=length(D.endpoints)-1

DynamicDefinition.is_full_branch(D::PwMap) = all(D.is_full)


function DynamicDefinition.preim(D::PwMap, k, y, ϵ)
	@assert 1 <= k <= nbranches(D)
	domain = hull(D.endpoints[k], D.endpoints[k+1])
	root(x->D.Ts[k](x)-y, domain, ϵ)
end

function DynamicDefinition.plottable(D::PwMap, x)
	for k in nbranches(D)
		domain = hull(D.endpoints[k], D.endpoints[k+1])
		if x in domain
			return D.Ts[k](x)
		end
	end
end

"""
Computes the hull of an iterable of intervals
"""
common_hull(S) = interval(minimum(x.lo for x in S), maximum(x.hi for x in S))

function DynamicDefinition.der(D::PwMap, x)
	common_hull(derivative(D.Ts[i], x ∩ hull(D.endpoints[i],D.endpoints[i+1])) for i in 1:length(D.endpoints)-1)
end

end
