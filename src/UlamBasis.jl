using ..BasisDefinition, ..DynamicDefinition, ..Contractors, ..Mod1DynamicDefinition, ..PwDynamicDefinition
using ValidatedNumerics, LinearAlgebra
import Base: iterate
import ..BasisDefinition: one_vector, integral_covector, is_integral_preserving, strong_norm, weak_norm, aux_norm

"""
Equispaced Ulam basis on [0,1] of size n
"""
struct Ulam{T<:AbstractVector} <:Basis
	p::T
	# TODO: check in constructor that p is sorted, starts with 0 and ends with 1
end

Ulam(n::Integer) = Ulam(LinRange(0., 1., n+1))
Base.length(B::Ulam) = length(B.p) - 1

function BasisDefinition.is_dual_element_empty(::Ulam, d)
	return isempty(d[1])
end

Base.length(S::DualComposedWithDynamic{<:Ulam, <:Dynamic}) = length(S.basis) * nbranches(S.dynamic)

"""
Returns dual elements which are pairs (i, (a,b))
i is an interval index, and (a,b) are the endpoints of its preimage
"""
function Base.iterate(S::DualComposedWithDynamic{<:Ulam, <:Dynamic}, state = (1, 1))
	i, k = state

	if i == length(S.basis)+1
			return nothing
	end

	a = preim(S.dynamic, k, S.basis.p[i], S.ϵ)
	b = preim(S.dynamic, k, S.basis.p[i+1], S.ϵ)

	if isempty(a) && !isempty(b)
		ep = endpoints(S.dynamic)
		if orientation(S.dynamic, k) > 0
			a = convert(typeof(b), ep[k])
		else
			a = convert(typeof(b), ep[k+1])
		end
	elseif isempty(b) && !isempty(a)
		ep = endpoints(S.dynamic)
		if orientation(S.dynamic, k) > 0
			b = convert(typeof(a), ep[k+1])
		else
			b = convert(typeof(a), ep[k])
		end
	end

	# this is a bit lazy because we could compute orientations and tell which is which
	# but it won't matter unless a,b are computed with very low precision
	a, b = min(a,b), max(a,b)

	if k == nbranches(S.dynamic)
		return ((i, (a, b)), (i+1, 1))
	else
		return ((i, (a, b)), (i, k+1))
	end
end

Base.eltype(f::DualComposedWithDynamic{<:Ulam, <:Dynamic}) = Tuple{Int64,Tuple{Interval{Float64},Interval{Float64}}}

"""
Returns the indices of the elements of the Ulam basis that intersect with the interval y
We do not assume an order of a and b; this should not matter unless
the preimages are computed with very low accuracy
We assume, though, that y comes from the (possibly inexact) numerical approximation
of an interval in [0,1], i.e., we restrict to y ∩ [0,1]
"""
function BasisDefinition.nonzero_on(B::Ulam, (a, b))
	y = hull(a, b)

	# finds in which semi-open interval [p[k], p[k+1]) y.lo and y.hi fall
	lo = searchsortedlast(B.p, y.lo)
	hi = searchsortedlast(B.p, y.hi)

	# they may be n+1 if y.hi==1
	lo = clamp(lo, 1, length(B))
	hi = clamp(hi, 1, length(B))

	return (lo, hi)
end

"""
Relative measure of the intersection of (a,b) wrt the whole interval (c,d)
Assumes that a,b and c,d are sorted correctly
"""
function relative_measure((a,b)::Tuple{<:Interval,<:Interval}, (c,d)::Tuple{<:Interval,<:Interval})
	lower = max(a, c)
	upper = min(b, d)
	intersection = max(upper - lower, 0) / (d-c)
	return intersection
end

"""
Given a preimage of an interval ```I_i```, this iterator returns
its relative intersection with all the elements of the Ulam basis that
have nonzero intersection with it
"""
function Base.iterate(S::ProjectDualElement{BT,DT}, state = S.j_min) where {BT<:Ulam,DT}
	if state == S.j_max+1
		return nothing
	end
	j = state
	x = relative_measure(S.dual_element,
			(@interval(S.basis.p[j]),
			@interval(S.basis.p[j+1])))
	return (j, x), state+1
end
Base.eltype(f::ProjectDualElement{<:Ulam, DT}) where{DT} = Tuple{Int64,Interval{Float64}}

BasisDefinition.evaluate(B::Ulam{T}, i, x) where {T} = (x>(i-1)/n) && (x<i/n) ? 1 : 0

BasisDefinition.evaluate_integral(B::Ulam{S}, i, T::Type) where{S} = T(i)/length(B)

function Base.iterate(S::AverageZero{Ulam{T}}, state = 1) where{T}
	n = length(S.basis)
	if state == n
		return nothing
	end
	v = zeros(Float64, n)
	v[1] = 1
	v[state+1]=-1
	return (v, state+1)
end

BasisDefinition.is_refinement(Bf::Ulam, Bc::Ulam) = Bc.p ⊆ Bf.p

function integral_covector(B::Ulam{T}) where{T}
	n = length(B)
	return 1/n * ones(Interval{Float64}, n)'
end

is_integral_preserving(B::Ulam{T}) where {T} = true

function one_vector(B::Ulam)
	return ones(length(B))
end

BasisDefinition.strong_norm(B::Ulam) = TotalVariation
BasisDefinition.weak_norm(B::Ulam) = L1
BasisDefinition.aux_norm(B::Ulam) = L1

# See BasisDefinition for docs on these constants
# These must be rounded up correctly!

BasisDefinition.weak_projection_error(B::Ulam) = 0.5 ⊘₊ Float64(length(B), RoundDown)
BasisDefinition.aux_normalized_projection_error(B::Ulam) = 0.
BasisDefinition.strong_weak_bound(B::Ulam) = Float64(length(B), RoundUp)
BasisDefinition.aux_weak_bound(B::Ulam) = 1.
BasisDefinition.weak_by_strong_and_aux_bound(B::Ulam) = (0., 1.)
BasisDefinition.bound_weak_norm_from_linalg_norm(B::Ulam) = (1., 0.)
BasisDefinition.bound_linalg_norm_L1_from_weak(B::Ulam) = 1.
BasisDefinition.bound_linalg_norm_L∞_from_weak(B::Ulam) = Float64(length(B), RoundUp)

function BasisDefinition.invariant_measure_strong_norm_bound(B::Ulam, D::Dynamic)
	A, B = dfly(strong_norm(B), aux_norm(B), D)
	@assert A < 1.
	return B ⊘₊ (1. ⊖₋ A)
end

using RecipesBase
using LaTeXStrings

"""
Plots a function in the Ulam basis
"""
@recipe function f(B::Ulam, w::AbstractVector)

	legend --> :bottomright

	if eltype(w) <: Interval
		w = mid.(w)
	end

	@series begin
		seriestype --> :steppost
		label --> L"f_{\delta}"
		ylims --> (0, NaN)
		B.p, vcat(w, w[end])
	end
end

"""
Displays error on a function in the Ulam basis

The w argument is unused, but kept for compatibility with other functions
for different bases
"""
@recipe function f(B::Ulam, error::Number, w=nothing)

	if isfinite(error)
		@series begin
			seriestype --> :path
			seriesalpha --> 0.5
			fillrange --> 0
			label --> "L1 Error"
			[0; sqrt(error)], [sqrt(error); sqrt(error)]
		end
	end
end
