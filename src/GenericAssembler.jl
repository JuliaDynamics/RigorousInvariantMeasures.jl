using ValidatedNumerics, SparseArrays
using ..DynamicDefinition, ..BasisDefinition

using LinearAlgebra

import ValidatedNumerics.IntervalArithmetic: mid
import Base: size, eltype
import LinearAlgebra: mul!

"""
Very generic assembler function
"""
function assemble(B::Basis, D::Dynamic, ϵ=2^(-40); T = Float64)
	I = Int64[]
	J = Int64[]
	nzvals = Interval{T}[]
	n = length(B)

	# TODO: reasonable size hint?

	for (i, dual_element) in DualComposedWithDynamic(B, D, ϵ)
		if dual_element != :∅
			for (j, x) in ProjectDualElement(B, dual_element)
				push!(I, i)
				push!(J, mod(j,1:n))
				push!(nzvals, x)
			end
		end
	end

	return sparse(I, J, nzvals, n, n)
end

abstract type DiscretizedOperator end

struct IntegralPreservingDiscretizedOperator{T<:AbstractMatrix} <: DiscretizedOperator
	L:: T
end
IntegralPreservingDiscretizedOperator(L) = IntegralPreservingDiscretizedOperator{typeof(L)}(L)

"""
An operator of the form Q = L + e*w (sparse plus rank-1).
"""
struct NonIntegralPreservingDiscretizedOperator{T<:AbstractMatrix, S<:AbstractVector, U<:AbstractMatrix} <: DiscretizedOperator
	L:: T
	e:: S
	w:: U
end
NonIntegralPreservingDiscretizedOperator(L, e, w) = NonIntegralPreservingDiscretizedOperator{typeof(L), typeof(e), typeof(w)}(L, e, w)

# Some cruft needed for eigs
Base.size(Q::DiscretizedOperator) = size(Q.L)
Base.eltype(Q::DiscretizedOperator) = eltype(Q.L)
LinearAlgebra.issymmetric(Q::IntegralPreservingDiscretizedOperator) = issymmetric(Q.L)
LinearAlgebra.issymmetric(Q::NonIntegralPreservingDiscretizedOperator) = issymmetric(Q.L) && Q.e' == Q.w

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
	return normL ⊕₊ norme ⊗₊ normw
end

function IntervalArithmetic.mid(Q::IntegralPreservingDiscretizedOperator)
	return IntegralPreservingDiscretizedOperator(map(mid, Q.L))
end
function IntervalArithmetic.mid(Q::NonIntegralPreservingDiscretizedOperator)
	# we are assuming that e is *not* an interval, for now. Types will break otherwise;
	# this may need to be fixed in an API.
	return NonIntegralPreservingDiscretizedOperator(map(mid, Q.L), Q.e, map(mid, Q.w))
end

function LinearAlgebra.mul!(Y, Q::IntegralPreservingDiscretizedOperator, v::AbstractArray, α, β)
	mul!(Y, Q.L, v, α, β)
end

function LinearAlgebra.mul!(Y, Q::NonIntegralPreservingDiscretizedOperator, v::AbstractArray, α, β)
	mul!(Y, Q.L, v, α, β)
	T = Base.promote_eltype(Q.e, Q.w)
	Y .+= convert(Array{T},Q.e) * (Q.w * v) * α
end

# these should be defined in terms of mul!, but it's simpler for now
function Base.:*(Q::IntegralPreservingDiscretizedOperator, v::Array)
	return Q.L * v
end

function Base.:*(Q::NonIntegralPreservingDiscretizedOperator, v::Array)
	T = Base.promote_eltype(Q.e, Q.w)
	return Q.L * v + convert(Array{T},Q.e) * (Q.w * v)
end

BasisDefinition.is_integral_preserving(Q::NonIntegralPreservingDiscretizedOperator) = false
BasisDefinition.is_integral_preserving(Q::IntegralPreservingDiscretizedOperator) = true
