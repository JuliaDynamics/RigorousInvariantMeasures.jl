module Mod1DynamicDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors

export Mod1Dynamic, preim, nbranches, plottable

struct Mod1Dynamic <: Dynamic
	T::Function
	nbranches::Integer
	orientation::Integer
	domain::Interval
end

function Mod1Dynamic(T::Function, nbranches = undef, domain::Interval{S} = Interval{Float64}(0,1)) where {S}
	bound_values = T(1)-T(0)
	nbranches = floor(Int64, abs(bound_values)) ### buggy
	orientation = bound_values > 0 ? 1 : -1
	
	return Mod1Dynamic(T, nbranches, orientation, domain)
end

DynamicDefinition.nbranches(S::Mod1Dynamic)=S.nbranches
### for the moment we suppose T(0)=0

DynamicDefinition.plottable(S::Mod1Dynamic, x) = mod1(S.T(x), 1.0)

function DynamicDefinition.preim(D::Mod1Dynamic, k, y, ϵ) 
	# we need to treat the case with the other orientation, 0 not fixed point...
	@assert k<= D.nbranches
	f(x) = D.T(x)-D.T(0)-(y-D.T(0)+(k-1)*D.orientation)
	root(f, D.domain, ϵ)
end

end