module PwDynamicDefinitionIsaia
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors

export PwIsa, preim, nbranches, plottable

struct PwIsa <: Dynamic
	Ts::Array{Function, 1}
	domains::Array{Interval, 1}
	orientations
end

DynamicDefinition.nbranches(D::PwIsa)=length(D.domains)

DynamicDefinition.is_full_branch(D::PwIsa) = D.is_full_branch

DynamicDefinition.plottable(D::PwIsa, x) = @error "Not Implemented"

import DualNumbers: Dual
PwIsa(Ts, domains) = PwIsa(Ts, domains, [sign(Ts[i](Dual(mid(domains[i]), 1)).epsilon) for i in 1:length(domains)])



function DynamicDefinition.preim(D::PwIsa, k, y, ϵ)
	# we need to treat the case with the other orientation, 0 not fixed point...
	@assert 1 <= k <= length(D.domains)
	root(x->D.Ts[k](x)-y, D.domains[k], ϵ)
end

import Base: iterate
import ..BasisDefinition: DualComposedWithDynamic
import InvariantMeasures: Ulam

function Base.iterate(S::DualComposedWithDynamic{Ulam, PwIsa}, state = (1, 1))
	i, k = state

	if i == length(S.basis)+1
			return nothing
	end

	# remark that this version supposes that for each i there exists a preimage
	# another more specific version should be implemented for maps with
	# incomplete branches

	x₁ = preim(S.dynamic, k, getindex(S.basis, i-1), S.ϵ)
	@info "x₁" x₁
	x₂ = preim(S.dynamic, k, getindex(S.basis, i), S.ϵ)
	@info "x₂" x₂
	

	if S.dynamic.orientations[k]>0
			if isempty(x₁) && !isempty(x₁)
				x₁ = typeof(x₂)(S.dynamic.domains[k].lo) 
			end
			if isempty(x₂) && !isempty(x₁)
				x₂ = typeof(x₁)(S.dynamic.domains[k].hi) 
			end
			lower, upper = x₁, x₂
	elseif	S.dynamic.orientations[k]<0
			if isempty(x₂) && !isempty(x₁)
				x₂ = typeof(x₁)(S.dynamic.domains[k].lo) 
			end
			if isempty(x₁) && !isempty(x₂)
				x₁ = typeof(x₂)(S.dynamic.domains[k].hi) 
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

end

import TaylorSeries
import InvariantMeasures: dfly

function dfly(::Type{TotalVariation}, ::Type{L1}, D::InvariantMeasures.PwDynamicDefinitionIsaia.PwIsa)
	dist = 0
	lam = 0
	disc = 0
	for i in 1:nbranches(D)
		f(x) = D.Ts[i](x)
		domain = D.domains[i]
		fprime(x) = f(TaylorSeries.Taylor1([x, 1], 1))[1]
		fsecond(x) = f(TaylorSeries.Taylor1([x, 1], 2))[2]/2
		distorsion(x)=abs(fsecond(x)/(fprime(x)^2))
		lambda(x) = abs(1/fprime(x))
    	dist = max(dist, maximise(distorsion, domain)[1].hi)
    	lam = max(lam, maximise(lambda, domain)[1].hi)
		disc = max(disc, (1/Interval(radius(domain)).hi))
	end
	return 2*lam, dist+disc
end
