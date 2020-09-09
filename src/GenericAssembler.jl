using ValidatedNumerics, SparseArrays
using LinearAlgebra: Adjoint
using ..DynamicDefinition, ..BasisDefinition

"""
Very generic assembler function
"""
function assemble(B::Basis, D::Dynamic, ϵ=2^(-40); T = Float64)
	n = length(B)
	# TODO: this should be assembled using `sparse`, instead, for speed
 	P = spzeros(Interval{T}, n, n)

	for (i, dual_element) in DualComposedWithDynamic(B, D, ϵ)
		for (j, x) in ProjectDualElement(B, dual_element)
			P[i,mod(j,1:n)] += x
		end
	end

	return P
end

abstract type DiscretizedOperator end

struct IntegralPreservingDiscretizedOperator{T<:AbstractMatrix} <: DiscretizedOperator
	L:: T
end
IntegralPreservingDiscretizedOperator(L) = IntegralPreservingDiscretizedOperator{typeof(L)}(L)

"""
An operator of the form Q = L + e*w (sparse plus rank-1)
"""
struct NonIntegralPreservingDiscretizedOperator{T<:AbstractMatrix, S<:AbstractVector, U<:Adjoint} <: DiscretizedOperator
	L:: T
	e:: S
	w:: U
end
NonIntegralPreservingDiscretizedOperator(L, e, w) = NonIntegralPreservingDiscretizedOperator{typeof(L), typeof(e), typeof(w)}(L, e, w)

opnormbound(N::Type{<:NormKind}, Q::IntegralPreservingDiscretizedOperator) = opnormbound(N, Q.L)

function DiscretizedOperator(B::Basis, D::Dynamic, ϵ=2^(-40); T = Float64)
	L = assemble(B, D, ϵ; T)
	if is_integral_preserving(B)
		return IntegralPreservingDiscretizedOperator(L)
	else
		f = integral_covector(B)
		e = one_vector(B)
		w = f - f*L #will use interval arithmetic when L is an interval matrix
		return NonIntegralPreservingDiscretizedOperator(L, e, w)
	end
end

function opnormbound(N::Type{<:NormKind}, Q::NonIntegralPreservingDiscretizedOperator)
	normL = opnormbound(N, Q.L)
	norme = opnormbound(N, Q.e)
	normw = opnormbound(N, Q.w)
	return round_expr(normL + norme * normw, RoundUp)
end
