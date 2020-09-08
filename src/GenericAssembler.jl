module GenericAssembler
using ValidatedNumerics, SparseArrays
using ..DynamicDefinition, ..BasisDefinition

export assemble

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

# Kept here for possible future use
#
# abstract type DiscretizedOperator end
#
# struct IntegralPreservingDiscretizedOperator{T<:AbstractMatrix} <: DiscretizedOperator
# 	L:: T
# end
#
# """
# An operator of the form Q = L + e*w (sparse plus rank-1)
# """
# struct NonIntegralPreservingDiscretizedOperator{T<:AbstractMatrix, S<:AbstractVector, U<:LinearAlgebra.Adjoint} <: DiscretizedOperator
# 	L:: T
# 	e:: S
# 	w:: U

end
