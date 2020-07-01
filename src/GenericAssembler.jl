module GenericAssembler
using ValidatedNumerics, SparseArrays
using ..DynamicDefinition, ..BasisDefinition 

export assemble

function assemble(B::Basis, D::Dynamic, ϵ=2^(-20); T = Float64, prec = 53)
"""
Very generic assembler function
"""
	n = length(B)
 	P = spzeros(Interval{T}, n, n)


	for (i, dual_element) in DualComposedWithDynamic(B, D, ϵ)
		for (j, x) in ProjectDualElement(B, dual_element)
			P[i,j] += x
		end
	end

	return P
end

end