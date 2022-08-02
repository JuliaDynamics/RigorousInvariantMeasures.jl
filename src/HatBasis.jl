"""
Hat basis on the Torus [0,1]
"""

using ..BasisDefinition, ..DynamicDefinition

import ..BasisDefinition: one_vector, integral_covector, is_integral_preserving

@doc """
	Hat

This type represents a Hat basis on ``S^1``. It contains a vector with the midpoints of the hats.
""" 
struct Hat{T<:AbstractVector} <:Basis
	p::T
	# TODO: check in constructor that p is sorted, starts with 0 and ends with 1
end

@doc """
	Hat(n::Integer)

This constructs a Hat basis on ``S^1`` on equispaced points
"""
Hat(n::Integer) = Hat(LinRange(0., 1., n+1))

"""
	Base.length(B::Hat{T})

Return the size of the Hat basis
"""
Base.length(B::Hat{T}) where {T} = Base.length(B.p) - 1

"""
	IntervalOnTorus

A separate type for intervals on the torus (mod 1) to "remind" us of the quotient

The interval is normalized in the constructor: the caller may assume that
* 0 <= i.lo < 1
* i.hi < i.lo + 1 OR i==Interval(0,1)
"""
struct IntervalOnTorus{T <: Real}
	I::Interval{T}
	function IntervalOnTorus{T}(I::Interval) where {T<:Real}
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

@doc raw"""
	HatFunctionOnTorus

Hat function (on the torus)

This is a piecewise linear function such that:
	f(x) = 0 if x <= lo
	f(mi) = 1
	f(x) = 0 if x >= ho
"""
struct HatFunctionOnTorus{T<: Number}
	lo::T
	mi::T
	hi::T
	function HatFunctionOnTorus{T}(lo, mi, hi) where {T <: Number}
		@assert 0 <= lo < 1
		@assert 0 <= mi < 1
		@assert 0 <= hi < 1
		new{T}(lo, mi, hi)
	end
end
HatFunctionOnTorus(lo::T, mi::T, hi::T) where {T<: Number} = HatFunctionOnTorus{T}(lo,mi,hi);
HatFunctionOnTorus(lo::Number, mi::Number, hi::Number) = HatFunctionOnTorus(promote(lo,mi,hi)...)

"""
Evaluate a HatFunctionOnTorus correctly on an IntervalOnTorus

Assumption: this is only called on functions defined on our partition,
so either mi==0, hi==0, or the three values are in increasing order
"""
function (f::HatFunctionOnTorus{T})(x::IntervalOnTorus) where {T}
	lo = Interval{T}(f.lo)
	mi = Interval{T}(f.mi)
	hi = Interval{T}(f.hi)
	# Since I is normalized, we only need two "hats": the one centered in mi,
	# and the one centered in mi+1. It is not complicated to check this
	# (one needs to treat mi==0 separately)

	# constructs correct endpoints for the hat centered in mi, also in the edge cases
	if mi == 0
		lo = lo - 1
	end
	if hi == 0
		hi = 1
	end

	I = x.I
	left_branch = (I-lo)/(mi-lo)
	right_branch = (hi-I)/(hi-mi)
	# 1 is not necessary in practice, but it avoids spurious results from rounding
	first_hat = max(min(left_branch, right_branch, 1), 0)

	I = x.I - 1
	left_branch = (I-lo)/(mi-lo)
	right_branch = (hi-I)/(hi-mi)
	# 1 is not necessary in practice, but it avoids spurious results from rounding
	second_hat = max(min(left_branch, right_branch, 1), 0)

	return max(first_hat, second_hat)
end

"""
	Base.getindex(B::Hat, i::Int)

Make so that B[j] returns a HatFunctionOnTorus with the j-th basis element
"""
function Base.getindex(B::Hat, i::Int)
	n = length(B)
	@boundscheck 1 <= i <= n || throw(BoundsError(B, i))
	return HatFunctionOnTorus(B.p[mod(i-1, 1:n)], B.p[mod(i, 1:n)], B.p[mod(i+1, 1:n)])
end

# """
# Return (in an iterator) the pairs (i, (x, |T'(x)|)) where x is a preimage of p[i], which
# describe the "dual" L* evaluation(p[i])
# """
# function Base.iterate(S::DualComposedWithDynamic{T, D}, state = (1, 1)) where T<:Hat where D<:Dynamic
# 	@assert is_full_branch(S.dynamic)

# 	i, k = state

# 	if i == length(S.basis)+1
# 			return nothing
# 	end

# 	x = preim(S.dynamic, k, S.basis.p[i], S.ϵ)
# 	absT′ = abs(derivative(S.dynamic, x))

# 	if k == nbranches(S.dynamic)
# 		return ((i, (x, absT′)), (i+1, 1))
# 	else
# 		return ((i, (x, absT′)), (i, k+1))
# 	end
# end

function BasisDefinition.is_dual_element_empty(::Hat, d)
	# TODO: the preim() may indeed be empty, so there could be an additional check here
	return false
end

BasisDefinition.is_refinement(Bf::Hat, Bc::Hat) = Bc.p ⊆ Bf.p

function BasisDefinition.integral_covector(B::Hat)
	n = length(B)
	return 1/n * ones(Interval{Float64}, n)'
end

function BasisDefinition.one_vector(B::Hat)
	return ones(length(B))
end


"""
	BasisDefinition.nonzero_on(B::Hat, dual_element)

Return the range of indices of the elements of the basis whose support intersects
with the given dual element (i.e., a pair (y, absT')).
The range may end with length(B)+1; this must be interpreted "mod length(B)":
it means that it intersects with the hat function peaked in 0 as well
(think for instance y = 0.9999).
"""
function BasisDefinition.nonzero_on(B::Hat, dual_element)
	y, absT′ = dual_element
	# Note that this cannot rely on arithmetic unless it is verified

	y = y ∩ Interval(0.,1.) # we assume it's bona-fide interval in [0,1]
	# this should work for preims(), since they are supposed to return
	# a number in [0,1]

	# finds in which semi-open interval [p[k], p[k+1]) y.lo and y.hi fall
	lo = searchsortedlast(B.p, y.lo)
	hi = searchsortedlast(B.p, y.hi)
	lo = min(lo, length(B)) # lo may be n+1 if y.lo==1
	hi = min(hi, length(B)) # hi may be n+1 if y.hi==1
	hi = hi + 1 # because the hat centered in p[k] is also nonzero in the interval before

	if lo == 1 # 1:N+1 does not make sense and would mean that the first interval is counted twice
		hi = min(hi, length(B))
	end
	return (lo, hi)
end

"""
Given a preimage ```y``` of a point ```x```, this iterator returns
```\\phi_j(y)/T'(y) ```
"""
function Base.iterate(S::ProjectDualElement{T,DT}, state = S.j_min) where {T <: Hat,DT}
	if state == S.j_max+1
		return nothing
	end
	y, absT′ = S.dual_element
	j = state
	y_normalized = IntervalOnTorus(y)
	n = length(S.basis)

	return ((j, S.basis[mod(j, 1:n)](y_normalized) / absT′),
		    state+1)
end

BasisDefinition.strong_norm(B::Hat) = Lipschitz
BasisDefinition.weak_norm(B::Hat) = Linf
BasisDefinition.aux_norm(B::Hat) = L1

BasisDefinition.evaluate_integral(B::Hat, i, T) = T(i)/length(B)

function Base.iterate(S::AverageZero{Hat{T}}, state = 1) where {T}
	n = length(S.basis)
	if state == n
		return nothing
	end
	v = zeros(Float64, n)
	v[1] = 1
	v[state+1]=-1
	return (v, state+1)
end

Base.length(S::AverageZero{Hat{T}}) where {T}= length(S.basis)-1

BasisDefinition.weak_projection_error(B::Hat) = 0.5 ⊘₊ Float64(length(B), RoundDown)
BasisDefinition.aux_normalized_projection_error(B::Hat) = 0.5 ⊘₊ Float64(length(B), RoundDown)
BasisDefinition.strong_weak_bound(B::Hat) = 2. ⊗₊ Float64(length(B), RoundDown)
BasisDefinition.aux_weak_bound(B::Hat) = 1.
BasisDefinition.weak_by_strong_and_aux_bound(B::Hat) = (1., 1.)
BasisDefinition.bound_weak_norm_from_linalg_norm(B::Hat) = @error "TODO"
BasisDefinition.bound_linalg_norm_L1_from_weak(B::Hat) = @error "TODO"
BasisDefinition.bound_linalg_norm_L∞_from_weak(B::Hat) = @error "TODO"
BasisDefinition.bound_weak_norm_abstract(B::Hat, D=nothing; dfly_coefficients=dfly(strong_norm(B), aux_norm(B), D)) = dfly_coefficients[2] ⊕₊ 1.

BasisDefinition.opnormbound(B::Hat{T}, N::Type{Linf}, A::AbstractVecOrMat{S}) where {S,T} = opnormbound(N, A)
BasisDefinition.normbound(B::Hat{T}, N::Type{Linf}, v) where {T} = normbound(N, v)

function BasisDefinition.invariant_measure_strong_norm_bound(B::Hat, D::Dynamic; dfly_coefficients=dfly(strong_norm(B), aux_norm(B), D))
	A, B = dfly_coefficients
	@assert A < 1.
	return B ⊘₊ (1. ⊖₋ A)
end


using RecipesBase

"""
Plots a function in the Hat basis
"""
@recipe function f(B::Hat, w::AbstractVector)

	legend --> :bottomright

	if eltype(w) <: Interval
		w = mid.(w)
	end

	@series begin
		seriestype --> :path
		label --> L"f_{\delta}"
		ylims --> (0, NaN)
		B.p, vcat(w, w[end])
	end
end

"""
Displays error on a function in the Hat basis
"""
@recipe function f(B::Hat, error::Number, w)

	if eltype(w) <: Interval
		w = mid.(w)
	end

	if isfinite(error)
		@series begin
			seriestype --> :path
			seriesalpha --> 0.5
			fillrange --> vcat(w, w[end]) .- error
			label --> "Error area"
			B.p, vcat(w, w[end]) .+ error
		end
	end
end

struct HatDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    x′::Vector{Interval}
end

Dual(B::Hat, D, ϵ) = HatDual(preimages_and_derivatives(B.p, D, 1:length(B.p)-1, ϵ)...)
Base.length(dual::HatDual) = length(dual.x)
Base.eltype(dual::HatDual) = Tuple{eltype(dual.xlabel), Tuple{eltype(dual.x), eltype(dual.x′)}}
function Base.iterate(dual::HatDual, state=1)
    if state <= length(dual.x)
        return ((dual.xlabel[state], (dual.x[state], abs(dual.x′[state]))), state+1)
    else
        return nothing
    end
end
