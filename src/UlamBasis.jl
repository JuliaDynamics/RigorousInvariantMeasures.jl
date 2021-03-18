using ..BasisDefinition, ..DynamicDefinition, ..Contractors, ..Mod1DynamicDefinition, ..PwDynamicDefinition
using ValidatedNumerics, LinearAlgebra
import Base: iterate
import ..BasisDefinition: one_vector, integral_covector, is_integral_preserving, strong_norm, weak_norm, aux_norm

"""
Equispaced Ulam basis on [0,1] of size n
"""
struct Ulam <:Basis
	n::Integer #TODO, change to partition
end

Base.length(B::Ulam) = B.n

Base.iterate(B::Ulam, state = 1) = state < length(B)+1 ? (B[state-1], state+1) : nothing

"""
Returns the left endpoint of the i-th element of the Ulam basis
"""
Base.getindex(B::Ulam, i) = Float64(i)/B.n


"""
This iterator returns the preimages of the endpoints
of the intervals defining the Ulam basis through the dynamic
"""
function Base.iterate(S::DualComposedWithDynamic{Ulam, <:Dynamic}, state = (1, 1))
	i, k = state

	if i == length(S.basis)+1
			return nothing
	end

	# remark that this version supposes that for each i there exists a preimage
	# another more specific version should be implemented for maps with
	# incomplete branches

	x₁ = preim(S.dynamic, k, getindex(S.basis, i-1), S.ϵ)
	x₂ = preim(S.dynamic, k, getindex(S.basis, i), S.ϵ)

	lower, upper = x₁, x₂

	if k == nbranches(S.dynamic)
		return ((i, (lower, upper)), (i+1, 1))
	else
		return ((i, (lower, upper)), (i, k+1))
	end
end

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


"""
Returns the indices of the elements of the Ulam basis that intersect with the interval y
"""
BasisDefinition.nonzero_on(B::Ulam, y) = max(floor(Int64, y[1].lo*length(B)), 1), min(ceil(Int64, y[2].hi*length(B)), length(B))



function relative_measure(S::Ulam, y, a, b)
	lower = max(y[1], a)
	upper = min(y[2], b)
	return max(upper-lower, 0)/(b-a)
end

"""
Given a preimage of an interval ```I_i```, this iterator returns
its relative intersection with all the elements of the Ulam basis that
have nonzero intersection with him
"""
function Base.iterate(S::ProjectDualElement{Ulam}, state = S.j_min)
	if state == S.j_max+1
		return nothing
	end

	return ((state, relative_measure(S.basis, S.dual_element,
			getindex(S.basis, state-1),
			getindex(S.basis, state))),
		    state+1)
end

BasisDefinition.evaluate(B::Ulam, i, x) = (x>(i-1)/n) && (x<i/n) ? 1 : 0

BasisDefinition.evaluate_integral(B::Ulam, i, T::Type)  = T(i)/length(B)

function Base.iterate(S::AverageZero{Ulam}, state = 1)
	n = length(S.basis)
	if state == n
		return nothing
	end
	v = zeros(Float64, n)
	v[1] = 1
	v[state+1]=-1
	return (v, state+1)
end

BasisDefinition.is_refinement(Bf::Ulam, Bc::Ulam) = rem(Bf.n, Bc.n) == 0

function integral_covector(B::Ulam)
	n = length(B)
	return 1/n * ones(Interval{Float64}, n)'
end

is_integral_preserving(B::Ulam) = true

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
		vcat(collect(Float64, B), 1.), vcat(w, w[end])
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
