module PwDynamicDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors

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

DynamicDefinition.nbranches(D::PwMap)=length(D.endpoints)-1

DynamicDefinition.is_full_branch(D::PwMap) = all(D.is_full)

import DualNumbers: Dual

PwMap(Ts, endpoints) = PwMap(Ts, endpoints, fill(false, length(endpoints)-1))

PwMap(Ts, endpoints, is_full) = PwMap(Ts, endpoints, is_full, [sign(Ts[i](Dual((endpoints[i]+endpoints[i+1])/2, 1)).epsilon) for i in 1:length(endpoints)-1])

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

end
