module BasisDefinition
using ..DynamicDefinition

import Base

export Basis, DualComposedWithDynamic, ProjectDualElement, AverageZero, assemble, integral_covector, one_vector, is_integral_preserving, strong_norm, weak_norm, aux_norm

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

strong_norm(B::Basis) = @error "Must be specialized"
weak_norm(B::Basis) = @error "Must be specialized"
aux_norm(B::Basis) = @error "Must be specialized"

"""
Check if Bfine is a refinement of Bcoarse
"""
is_refinement(Bfine::Basis, Bcoarse::Basis) = @error "Not implemented"

"""
	Covector that represents the integral
"""
integral_covector(B::Basis) = @error "Must be specialized"

"""
	Vector that represents the function 1
"""
one_vector(B::Basis) = @error "Must be specialized"

"""
	Integral-preserving discretizations may specialize this to "true"
"""
is_integral_preserving(B::Basis) = false

"""
Integral of a function in U_h

Args:
	v (any type of vector):

Returns:
	the integral, computed with the arithmetic of v.
"""
function integral(B::Basis, v; T = Float64)
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
	Return a constant Kh (typically scales as h ~ 1/n) such that `||P_h f-f||\\leq Kh ||f||_s`
	Must be rounded up correctly!

	Arg:
		B::Basis
"""
weak_projection_error(B::Basis) = @error "Not implemented"

"""
	Return a constant Eh (typically scales as h ~ 1/n) such that `|||P_h f|||\\leq |||f|||+ Eh * ||f||_s`
	Must be rounded up correctly!

	Arg:
		B::Basis
"""
aux_normalized_projection_error(B::Basis) = @error "Not implemented"

"""
	Return a constant M₁n such that for a vector v in Uₕ `||v||_s\\leq M1n*||v||`
	Must be rounded up correctly!
"""
strong_weak_bound(B::Basis) = @error "Not implemented"

"""
	Return a constant M₂ such that for a vector v in Uₕ `|||v|||\\leq M_2||v||`
	Must be rounded up correctly!
"""
aux_weak_bound(B::Basis) = @error "Not implemented"

"""
	Return constants S₁, S₂ such that for a vector v in Uₕ `||v||\\leq S_1||v||_s+S_2|||v|||`
	Must be rounded up correctly!
"""
weak_by_strong_and_aux_bound(B::Basis) = @error "Not implemented"

"""
	Return constants W₁, W₂ such that for a vector v in Uₕ `||v||\\leq W_1||v||_1+W_2||v||_{\\infty}`
	Must be rounded up correctly!
"""
bound_weak_norm_from_linalg_norm(B::Basis) = @error "Not implemented"

"""
	Return a constant A such that for a vector v in Uₕ `||v||_1\\leq A||v||`
	Must be rounded up correctly!
"""
bound_linalg_norm_L1_from_weak(B::Basis) = @error "Not implemented"

"""
	Return a constant A such that for a vector v in Uₕ `||v||_\\infty \\leq A||v||`
	Must be rounded up correctly!
"""
bound_linalg_norm_L∞_from_weak(B::Basis) = @error "Not implemented"

end
