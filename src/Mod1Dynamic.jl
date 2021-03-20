module Mod1DynamicDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors

using ..DynamicDefinition: derivative

export Mod1Dynamic, preim, nbranches, plottable

"""
Defines a Dynamic on [0,1] as the Mod-1 quotient of a given map.

An alternative newer implementation relying on piecewise-defined functions is in `mod1_dynamic`
"""
struct Mod1Dynamic{FT} <: MarkovDynamic
	T::FT
	nbranches::Int
	orientation::Int
	domain::Interval{Float64}
	is_full_branch::Bool
end

Mod1Dynamic(T::FT, nbranches = undef, domain = Interval{Float64}(0,1)) where {FT} = Mod1Dynamic{FT}(T, nbranches, domain)

# TODO: serious doubts that this works if T(0) is not an integer...
function Mod1Dynamic{FT}(T, nbranches = undef, domain = Interval{Float64}(0,1)) where {FT}
	@assert domain == 0..1 # TODO: this only works for domain == 0..1, for now
	range_diff = T(1..1)-T(0..0)
	@assert 0 ∉ range_diff # TODO: check that the function is always increasing OR decreasing with the derivative instead?
	orientation = range_diff > 0 ? 1 : -1
	nbranches = ceil(orientation * range_diff).hi
	is_full_branch = isinteger(range_diff)
	return Mod1Dynamic{FT}(T, nbranches, orientation, domain, is_full_branch)
end

DynamicDefinition.domain(S::Mod1Dynamic{FT}) where {FT} = S.domain

DynamicDefinition.nbranches(S::Mod1Dynamic{FT}) where {FT} =S.nbranches

DynamicDefinition.is_full_branch(S::Mod1Dynamic{FT}) where {FT} = S.is_full_branch

function DynamicDefinition.preim(D::Mod1Dynamic{FT}, k, y, ϵ) where {FT}
	# we need to treat the case with the other orientation, 0 not fixed point...
	@assert 1 <= k <= D.nbranches
	f(x) = D.T(x)-D.T(0)-(y-D.T(0)+(k-1)*D.orientation)
	root(f, D.domain, ϵ)
end

DynamicDefinition.derivative(n, D::Mod1Dynamic{FT}, x) where {FT} = derivative(n, D.T, x)
DynamicDefinition.distorsion(D::Mod1Dynamic{FT}, x) where {FT} = distorsion(D.T, x)

DynamicDefinition.max_distorsion(D::Mod1Dynamic{FT}, tol=1e-3) where {FT} = maximise(x -> distorsion(D.T, x), domain(D), tol=tol)[1]
DynamicDefinition.expansivity(D::Mod1Dynamic{FT}, tol=1e-3) where {FT} = maximise(x -> abs(1/derivative(D, x)), domain(D), tol=tol)[1]

function DynamicDefinition.plottable(D::Mod1Dynamic{FT}, x) where {FT}
	@assert 0 <= x <= 1
	return mod(D.T(x), 1.)
end

using RecipesBase
@recipe f(::Type{Mod1Dynamic{FT}}, D::Mod1Dynamic{FT}) where {FT} = x -> plottable(D, x)

end
