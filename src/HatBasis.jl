module HatBasis
"""
Hat basis on the Torus [0,1]
"""

using ..BasisDefinition, ..Mod1DynamicDefinition
using ValidatedNumerics

export Hat, HatFunction

"""
Equispaced partition of size n of [0,1]
"""
struct EquispacedPartition{T} <: AbstractVector{T}
	n::Integer
end
function Base.getindex(p::EquispacedPartition{T}, i::Int) where {T}
	@boundscheck 1 <= i <= p.n || throw(BoundsError(p, i))
	return convert(T, i-1) / p.n
end
Base.size(p::EquispacedPartition) = (p.n,)
Base.IndexStyle(::Type{<:EquispacedPartition}) = IndexLinear()

struct Hat{T<:AbstractVector} <:Basis
	p::T
end

"""
Return the size of the Hat basis
"""
Base.length(b::Hat) = length(b.p)

"""
Hat function (on the reals)

This is a piecewise linear function such that:
	f(x) = 0 if x <= lo
	f(mi) = 1
	f(hi)
"""
struct HatFunction{T<: Number}
	lo::T
	mi::T
	hi::T
	function HatFunction{T}(lo, mi, hi) where T <: Number
		@assert lo <= mi <= hi
		new{T}(lo, mi, hi)
	end
end
HatFunction(lo::T, mi::T, hi::T) where {T<: Number} = HatFunction{T}(lo,mi,hi);
HatFunction(lo::Number, mi::Number, hi::Number) = HatFunction(promote(lo,mi,hi)...)

"""
Evaluate a HatFunction.

Must work correctly when S is an interval.
"""
function (f::HatFunction{T})(x::S) where T where S <: Number
	lo = convert(S, f.lo)
	mi = convert(S, f.mi)
	hi = convert(S, f.hi)

	left_branch = (x-lo)/(mi-lo)
	right_branch = (hi-x)/(hi-mi)
	return max(min(left_branch, right_branch, 1), 0)
end

"""
A separate type for intervals on the torus (mod 1) to "remind" us of the quotient

The interval is normalized in the constructor: the caller may assume that
* i.lo < 1
* i.hi < i.lo + 1 OR i==Interval(0,1)
"""
struct IntervalOnTorus{T <: Real}
	I::Interval{T}
	function IntervalOnTorus{T}(I) where T<:Real
		# Note that this test may "expand" intervals such as 1e-30..1, to 0..1, but it is not a big issue anyway
		if diam(I) >= 1.
			new{T}(0..1)
		else
			# Note that I - floor(I.lo) may return something smaller than zero in some rounding modes
			new{T}(max(I - floor(I.lo), 0))
		end
	end
end
IntervalOnTorus(I::Interval{T}) where {T} = IntervalOnTorus{T}(I)

"""
Evaluate a HatFunction correctly on an interval on the Torus.

Note that IntervalOnTorus is not a subtype of number, so this does not overlap with the previous definition
"""
function (f::HatFunction{T})(x::S) where T where S <: IntervalOnTorus
	return max(f(x.I), f(x.I-1))
end

# Older parts

Base.iterate(b::Hat, state = 1) = state < length(b)+1 ? (b[state-1], state+1) : nothing

"""
Returns the point on which the i-th element of the Hat basis has value 1
TODO: should probably return a custom type
"""
function Base.getindex(b::Hat, i::Int)
	@boundscheck 0 <= i < length(b) || throw(BoundsError(b, i))
	return Float64(i)/b.n
end

"""
This iterator returns the preimages of the endpoints
of the intervals defining the Ulam basis through the dynamic
"""
function Base.iterate(S::DualComposedWithDynamic{Hat, Mod1Dynamic}, state = (1, 1))
	i, k = state

	if i == length(S.basis)+1
			return nothing
	end

	# remark that this version supposes that for each i there exists a preimage
	# another more specific version should be implemented for maps with
	# incomplete branches


	x = preim(S.dynamic, k, getindex(S.basis, i-1), S.ϵ)
	der = der(S.dynamic, x)

	if k == nbranches(S.dynamic)
		return ((i, (x, der)), (i+1, 1))
	else
		return ((i, (x, der)), (i, k+1))
	end
end

"""
Returns the indices of the elements of the Ulam basis that intersect with the interval y
"""
nonzero_on(B::Hat, y) = max(floor(Int64, y[1].lo*length(B)), 1), min(ceil(Int64, y[2].hi*length(B)), length(B))

"""
Constructor for the ProjectDualElement iterator
"""
function ProjectDualElement(B::Hat, y)
	j_min, j_max = nonzero_on(B, y)
	return ProjectDualElement(B, j_min, j_max, y)
end


ϕ(B::Hat, j, y) = @error Not Implemented


# """
# Given a preimage ```y``` of a point ```x```, this iterator returns
# ```\\phi_j(y)/T'(y) ```
# """
# function Base.iterate(S::ProjectDualElement{Hat}, state = S.j_min)
# 	if state == S.j_max+1
# 		return nothing
# 	end
# 	y, der = S.dual_element
#
# 	return ((state, ϕ(S.basis, state, y)/der),
# 		    state+1)
# end

evaluate_integral(B::Hat, i, T) = T(i)/length(B)

function Base.iterate(S::AverageZero{Hat}, state = 1)
	n = length(S.basis)
	if state == n
		return nothing
	end
	v = zeros(Float64, n)
	v[1] = 1
	v[state+1]=-1
	return (v, state+1)
end

using RecipesBase

@userplot PlotHat
@recipe function f(h::PlotHat)
	if length(h.args)!= 2 || (typeof(h.args[1])!= Ulam) || !(typeof(h.args[2])<:AbstractVector)
		error("Plot Ulam needs as an input a Ulam Basis and a vector")
	end

	B = h.args[1]
	w = h.args[2]

	seriestype := :path
	collect(B), mid.(w)
end

end
