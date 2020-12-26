module Mod1DynamicDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors

using ..DynamicDefinition: derivative

export Mod1Dynamic, preim, nbranches, plottable

struct Mod1Dynamic <: MarkovDynamic
	T::Function
	nbranches::Integer
	orientation::Integer
	domain::Interval
	is_full_branch::Bool
end


function Mod1Dynamic(T::Function, nbranches = undef, domain::Interval{S} = Interval{Float64}(0,1)) where {S}
	@assert domain == 0..1 # TODO: this only works for domain == 0..1, for now
	range_diff = T(1..1)-T(0..0)
	@assert 0 ∉ range_diff # TODO: check that the function is always increasing OR decreasing with the derivative instead?
	orientation = range_diff > 0 ? 1 : -1
	nbranches = ceil(orientation * range_diff).hi
	is_full_branch = isinteger(range_diff)
	return Mod1Dynamic(T, nbranches, orientation, domain, is_full_branch)
end

DynamicDefinition.domain(S::Mod1Dynamic) = S.domain

DynamicDefinition.nbranches(S::Mod1Dynamic)=S.nbranches

DynamicDefinition.is_full_branch(S::Mod1Dynamic) = S.is_full_branch

DynamicDefinition.plottable(S::Mod1Dynamic, x) = mod1(S.T(x), 1.0)

function DynamicDefinition.preim(D::Mod1Dynamic, k, y, ϵ)
	# we need to treat the case with the other orientation, 0 not fixed point...
	@assert 1 <= k <= D.nbranches
	f(x) = D.T(x)-D.T(0)-(y-D.T(0)+(k-1)*D.orientation)
	root(f, D.domain, ϵ)
end

DynamicDefinition.derivative(n, D::Mod1Dynamic, x) = derivative(n, D.T, x)
DynamicDefinition.distorsion(D::Mod1Dynamic, x) = distorsion(D.T, x)

end
