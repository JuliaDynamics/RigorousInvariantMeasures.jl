module UlamBasis

using ..BasisDefinition, ..DynamicDefinition, ..Contractors, ..Mod1DynamicDefinition
using ValidatedNumerics
import Base: iterate

export Ulam

"""
Equispaced Ulam basis on [0,1] of size n 
"""
struct Ulam <:Basis
	n::Integer #TODO, change to partition
end

"""
Returns the size of the Ulam basis
"""
Base.length(b::Ulam) = b.n


Base.iterate(b::Ulam, state = 1) = state < length(b)+1 ? (b[state-1], state+1) : nothing

"""
Returns the left endpoint of the i-th element of the Ulam basis
"""
Base.getindex(b::Ulam, i) = Float64(i)/b.n

"""
This iterator returns the preimages of the endpoints 
of the intervals defining the Ulam basis through the dynamic
"""
function Base.iterate(S::DualComposedWithDynamic{Ulam, Dyn}, state = (1, 1)) where {Dyn<:Dynamic}
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

SpaceConstant(B::Ulam, ::Val{:L1}) = 1
SpaceConstant(B::Ulam, ::Val{:L∞}) = 1


bound_weak_norm_from_linalg_norm(B::Ulam) = (1/length(B), 0.0)
bound_linalg_norm_L1_from_weak(B::Ulam) = 1/length(B)
bound_linalg_norm_L∞_from_weak(B::Ulam) = length(B) #this is defined for coherence

function dfly(B::Ulam, D::Dynamic) 
	distorsion(x) = der_der(D, x)/(der(D, x))^2
	der_range = range_estimate(x->der(D, x), D.domain) 
	distorsion_range = range_estimate(distorsion, D.domain)
 	return (1/abs(der_range)).hi, abs(distorsion_range).hi
end

using RecipesBase
using LaTeXStrings

@userplot PlotUlam
@recipe function f(h::PlotUlam)
	if length(h.args)<2 || (typeof(h.args[1])!= Ulam) || !(typeof(h.args[2])<:AbstractVector)
		error("Plot Ulam needs as an input a Ulam Basis and a vector")
	end

	B = h.args[1]
	w = h.args[2]
	D = h.args[3]

	layout := (3, 1)
	link := :both
	
	if eltype(w) <: Interval
		w = mid.(w)
	end


	# bar plot
	@series begin
		linecolor := :blue
		linealpha := 0.8
		fillalpha := 0.8
		seriestype := :bar
		label := L"f_{\delta}"
		collect(B), w
	end

	#dynamic plot
	@series begin
		linecolor := :green
	    seriestype := :path
	    label := L"T(x)"
	    G = x-> plottable(D, x)
	    collect(B), G.(collect(B))
	end

	@series begin
		linecolor := :red
	    seriestype := :path
	    label := L"T'(x)"
	    G = x-> der(D, x)
	    collect(B), G.(collect(B))
	end
end

end