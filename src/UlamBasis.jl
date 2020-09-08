module UlamBasis

using ..BasisDefinition, ..DynamicDefinition, ..Contractors, ..Mod1DynamicDefinition
using ValidatedNumerics, LinearAlgebra
import Base: iterate

export Ulam

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

"""
Returns constant η such that ``\\min ||Uv||\\geq η ||v||``, where
U = collect(AverageZero) (the matrix whose columns are the vectors in AverageZero)
TODO: introduce a default method that uses rigorous SVD to bound these constants
from below.
This is used in the following estimate
```math:
	\\max_x ||Ax|_V||/||x|| = \\max ||AUz||/||Uz|| \\leq \\max_z \\eta||AUz||/||z||
```
"""
BasisDefinition.spaceconstant(B::Ulam, ::Val{:L1}) = 1
BasisDefinition.spaceconstant(B::Ulam, ::Val{:L∞}) = 1

"""
 	Rigorous estimate (from above) of ||v||_w

	Args:
		B basis
	 	v (numpy vector):

	Returns: x such that ``||v||_w \\leq x``
"""
BasisDefinition.norm_estimate(B::Ulam, v) = norm(v, 1)

"""
	Rigorous norm of a vector.

	Args:
		B basis
	 	v

	Returns:
	 		its (weak) norm.
"""
BasisDefinition.rigorous_weak_norm(B::Ulam, v) = norm(v,1)

"""
	Rigorous estimate (from above) of the matrix norm

	Args:
		B Basis
		PP Matrix

	 	Returns: x such that ``||PP||_w``
"""
BasisDefinition.matrix_norm_estimate(B::Ulam, P) = opnorm(P, 1)

"""
	Diameter (in the matrix norm) of an interval matrix.

	Must be rigorous.

	Returns:
	 		M such that :math:`\\|P_1-P_2\\|_w \\leq M` for all :math:`P_1,P_2 \\in P`.
"""
BasisDefinition.matrix_norm_diameter(B::Ulam, P) = opnorm(diam.(P), 1)

"""
	Computes the residual (in norm) of the computed Perron vector

	Args:
	 		P (interval matrix):
	 		v (numpy vector):

	Returns:
	 		res (real RNDU): an upper bound to :math:`\\|Pv-v\\|`
"""
BasisDefinition.residual_estimate(B::Ulam, P, v) = rigorous_weak_norm(P*v-v)

"""
	Returns a constant K such that `||P_h f-f||\\leq K h ||f||_s`

	Arg:
		B::Basis
"""
BasisDefinition.normapprox(B::Ulam) = 1/2

"""
	Returns a constant E such that `|||P_h f|||\\leq |||f|||+E h ||f||_s`
	Arg:
		B::Basis
"""
BasisDefinition.boundweak(B::Ulam) = 0

"""
	Returns a constant M₁ such that for a vector v in Uₕ `||v||_s\\leq \\frac{M_1}{h}||v||`
"""
BasisDefinition.boundstrongbyweak(B::Ulam) = 1

"""
	Returns a constant M₂ such that for a vector v in Uₕ `|||v|||\\leq M_2||v||`
"""
BasisDefinition.boundauxiliarybyweak(B::Ulam) = 1

"""
	Returns constants S₁, S₂ such that for a vector v in Uₕ `||v||\\leq S_1||v||_s+S_2|||v|||`
"""
BasisDefinition.boundweakbystrongauxiliary(B::Ulam) = (0, 1)

"""
	Returns constants W₁, W₂ such that for a vector v in Uₕ `||v||\\leq W_1||v||_1+W_2||v||_{\\infty}`
"""
BasisDefinition.bound_weak_norm_from_linalg_norm(B::Ulam) = (1, 0)

"""
	Returns constant A such that for a vector v in Uₕ `||v||_1\\leq A||v||`
"""
BasisDefinition.bound_linalg_norm_L1_from_weak(B::Ulam) = 1

"""
	Returns constant A such that for a vector v in Uₕ `||v||_\\infty \\leq A||v||`
"""
BasisDefinition.bound_linalg_norm_L∞_from_weak(B::Ulam) = length(B)


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
