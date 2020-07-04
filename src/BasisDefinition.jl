module BasisDefinition
using ..DynamicDefinition

import Base

export Basis, DualComposedWithDynamic, ProjectDualElement, AverageZero, assemble


abstract type Basis end

length(B::Basis) = @error "Not Implemented"

struct DualComposedWithDynamic{B<:Basis, D<:Dynamic} 
	basis::B
	dynamic::D
	ϵ::Float64
end

Base.iterate(S::DualComposedWithDynamic{B, D}, state) where {B<:Basis, D<:Dynamic} = @error "Not implemented"

struct ProjectDualElement{B<:Basis}
	basis::B
	j_min
	j_max
	dual_element
end
ProjectDualElement(basis::B, j_min, j_max, y) where {B} = ProjectDualElement{B}(basis, j_min, j_max, y)

nonzero_on(B::Basis, I) = @error "Not Implemented"

function ProjectDualElement(B::Basis, y)
 	j_min, j_max = nonzero_on(B, y)
 	return ProjectDualElement(B, j_min, j_max, y)
end

Base.iterate(S::ProjectDualElement{B}, state) where {B} = @error "Not Implemented"

evaluate(B::Basis, i, x) = @error "Not Implemented"
evaluate_integral(B::Basis, i; T = Float64) = @error "Not Implemented"

function integral(B::Basis, v; T = Float64)
	"""
	 	Integral of a function in U_h

	 	Args:
	 		v (any type of vector):

	 	Returns:
	 		the integral, computed with the arithmetic of v.
	"""
	return sum([T(v[i])*evaluate_integral(B, i, T) for i in 1:length(B)])
end

"""
	Yields a basis of the space of average zero vectors
"""
struct AverageZero{B<:Basis}
	basis::B
end

Base.iterate(S::AverageZero{B}, state) where {B} = @error "Not Implemented"
Base.length(S::AverageZero) = length(S.basis)-1

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
SpaceConstant(B::Basis, ::Val{:L1}) = @error "Not Implemented"
SpaceConstant(B::Basis, ::Val{:L∞}) = @error "Not Implemented"


"""
 	Rigorous estimate (from above) of ||v||_w
		
	Args:
		B basis
	 	v (numpy vector):
			
	Returns: x such that ``||v||_w \\leq x``
"""
norm_estimate(B::Basis, v) = @error "Not Implemented"

"""
	Rigorous norm of a vector.
				
	Args:
		B basis
	 	v 
			
	Returns:
	 		its (weak) norm. 
"""
rigorous_weak_norm(B::Basis, v) = @error "Not Implemented" 	

"""
	Rigorous estimate (from above) of the matrix norm
		
	Args:
		B Basis
		PP Matrix
			
	 	Returns: x such that ``||PP||_w``
""" 	
matrix_norm_estimate(B::Basis, P) = @error "Not Implemented"

"""
	Diameter (in the matrix norm) of an interval matrix.
		
	Must be rigorous.
		
	Returns:
	 		M such that :math:`\\|P_1-P_2\\|_w \\leq M` for all :math:`P_1,P_2 \\in P`.
"""
matrix_norm_diameter(B::Basis, P) = @error "Not Implemented"

"""
	Computes the residual (in norm) of the computed Perron vector
		
	Args:
	 		P (interval matrix):
	 		v (numpy vector):
		
	Returns:
	 		res (real RNDU): an upper bound to :math:`\\|Pv-v\\|`
"""
residual_estimate(B::Basis, P, v) = @error "Not Implemented" 

"""
	Returns a constant K such that `||P_h f-f||\\leq K h ||f||_s`

	Arg:
		B::Basis
"""
normapprox(B::Basis) = @error "Not implemented"

"""
	Returns a constant E such that `|||P_h f|||\\leq |||f|||+E h ||f||_s`
	Arg:
		B::Basis
"""
boundweak(B::Basis) = @error "Not implemented"

"""
	Returns a constant M₁ such that for a vector v in Uₕ `||v||_s\\leq \\frac{M_1}{h}||v||` 
"""
boundstrongbyweak(B::Basis) = @error "Not implemented"

"""
	Returns a constant M₂ such that for a vector v in Uₕ `|||v|||\\leq M_2||v||` 
"""
boundauxiliarybyweak(B::Basis) = @error "Not implemented"

"""
	Returns constants S₁, S₂ such that for a vector v in Uₕ `||v||\\leq S_1||v||_s+S_2|||v|||` 
"""
boundweakbystrongauxiliary(B::Basis) = @error "Not implemented"

"""
	Returns constants W₁, W₂ such that for a vector v in Uₕ `||v||\\leq W_1||v||_1+W_2||v||_{\\infty}`
"""
bound_weak_norm_from_linalg_norm(B::Basis) = @error "Not implemented"

"""
	Returns constant A such that for a vector v in Uₕ `||v||_1\\leq A||v||`
"""
bound_linalg_norm_L1_from_weak(B::Basis) = @error "Not implemented"

"""
	Returns constant A such that for a vector v in Uₕ `||v||_\\infty \\leq A||v||`
"""
bound_linalg_norm_L∞_from_weak(B::Basis) = @error "Not implemented"

"""
	Returns constants A, B such that `||Lf||_s\\leq A||f||_s+B|||f|||`
"""
dfly(B::Basis, D::Dynamic) = @error "Not implemented"

end