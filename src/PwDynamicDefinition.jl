module PwDynamicDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors

export PwMap, preim, nbranches, plottable

struct PwMap <: Dynamic
	Ts::Array{Function, 1}
	endpoints::Array{Interval, 1}
	orientations
end

DynamicDefinition.nbranches(D::PwMap)=length(D.endpoints)-1

DynamicDefinition.is_full_branch(D::PwMap) = false

import DualNumbers: Dual
PwMap(Ts, endpoints) = PwMap(Ts, endpoints, [sign(Ts[i](Dual((endpoints[i]+endpoints[i+1])/2, 1)).epsilon) for i in 1:length(endpoints)-1])



function DynamicDefinition.preim(D::PwMap, k, y, ϵ)
	# we need to treat the case with the other orientation, 0 not fixed point...
	@assert 1 <= k <= nbranches(D)
	domain = hull(D.endpoints[k], D.endpoints[k+1])
	root(x->D.Ts[k](x)-y, domain, ϵ)
end

import Base: iterate
import ..BasisDefinition: DualComposedWithDynamic
import InvariantMeasures: Ulam

function Base.iterate(S::DualComposedWithDynamic{Ulam, PwMap}, state = (1, 1))
	i, k = state

	if i == length(S.basis)+1
			return nothing
	end

	# remark that this version supposes that for each i there exists a preimage
	# another more specific version should be implemented for maps with
	# incomplete branches

	x₁ = preim(S.dynamic, k, getindex(S.basis, i-1), S.ϵ)
	#@info "x₁" x₁
	x₂ = preim(S.dynamic, k, getindex(S.basis, i), S.ϵ)
	#@info "x₂" x₂
	

	if S.dynamic.orientations[k]>0
			if isempty(x₁) && !isempty(x₁)
				x₁ = S.dynamic.endpoints[k] 
			end
			if isempty(x₂) && !isempty(x₁)
				x₂ = S.dynamic.endpoints[k+1] 
			end
			lower, upper = x₁, x₂
	elseif	S.dynamic.orientations[k]<0
			if isempty(x₂) && !isempty(x₁)
				x₂ = S.dynamic.endpoints[k] 
			end
			if isempty(x₁) && !isempty(x₂)
				x₁ = S.dynamic.endpoints[k+1] 
			end
			lower, upper = x₂, x₁
	end

	if k == nbranches(S.dynamic)
		if isempty(lower) && isempty(upper)
			return ((i, :∅), (i+1, 1))
		end
		return ((i, (lower, upper)), (i+1, 1))
	else
		if isempty(lower) && isempty(upper)
			return ((i, :∅), (i, k+1))
		end
		return ((i, (lower, upper)), (i, k+1))
	end
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

import TaylorSeries
import InvariantMeasures: dfly

function dfly(::Type{TotalVariation}, ::Type{L1}, D::InvariantMeasures.PwDynamicDefinition.PwMap)
	dist = 0
	lam = 0
	disc = 0
	for i in 1:nbranches(D)
		f(x) = D.Ts[i](x)
		domain = hull(D.endpoints[i], D.endpoints[i+1])
		fprime(x) = f(TaylorSeries.Taylor1([x, 1], 1))[1]
		fsecond(x) = f(TaylorSeries.Taylor1([x, 1], 2))[2]/2
		distorsion(x)=abs(fsecond(x)/(fprime(x)^2))
		lambda(x) = abs(1/fprime(x))
    	dist = max(dist, maximise(distorsion, domain)[1].hi)
    	lam = max(lam, maximise(lambda, domain)[1].hi)
    	low_rad = (abs(D.endpoints[i]-D.endpoints[i+1])/2).lo
		disc = max(disc, ((1/Interval(low_rad)).hi))
	end
	return 2*lam, dist+disc
end
