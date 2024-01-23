export Basis,
    DualComposedWithDynamic,
    ProjectDualElement,
    AverageZero,
    assemble,
    integral_covector,
    one_vector,
    is_integral_preserving,
    strong_norm,
    weak_norm,
    aux_norm,
    is_dual_element_empty,
    nonzero_on,
    is_refinement,
    opnormbound,
    normbound,
    bound_weak_norm_abstract,
    Dual

abstract type Basis end

Base.length(B::Basis) = @error "Not Implemented"

struct DualComposedWithDynamic{B<:Basis,D<:Dynamic}
    basis::B
    dynamic::D
    ϵ::Float64
end

#Base.iterate(S::DualComposedWithDynamic{B, D}, state) where {B<:Basis, D<:Dynamic} = @error "Not implemented"

struct ProjectDualElement{B<:Basis,DT}
    basis::B
    j_min::Int64
    j_max::Int64
    dual_element::DT
end
ProjectDualElement(basis::B, j_min, j_max, y::DT) where {B,DT} =
    ProjectDualElement{B,DT}(basis, j_min, j_max, y)
Base.length(S::ProjectDualElement{B,DT}) where {B,DT} = S.j_max - S.j_min + 1

is_dual_element_empty(B::Basis, I) = @error "Not Implemented"
nonzero_on(B::Basis, I) = @error "Not Implemented"

function ProjectDualElement(B::Basis, y)
    j_min, j_max = nonzero_on(B, y)
    return ProjectDualElement(B, j_min, j_max, y)
end

# Base.iterate(S::ProjectDualElement, state) = @error "Not Implemented"

"""
	evaluate(B::Basis, i, x)

Evaluate the i-th basis element at x
"""
evaluate(B::Basis, i, x) = @error "Not Implemented"

"""
	evaluate_integral(B::Basis, i; T = Float64)

Value of the integral on [0,1] of the i-th basis element
"""
evaluate_integral(B::Basis, i; T = Float64) = @error "Not Implemented"

"""
	strong_norm(B::Basis)

Return the type of the strong norm of the basis

# Example

```jldoctest
julia> using RigorousInvariantMeasures

julia> B = Ulam(1024)
Ulam{LinRange{Float64, Int64}}(LinRange{Float64}(0.0, 1.0, 1025))

julia> strong_norm(B)
TotalVariation
```
"""
strong_norm(B::Basis) = @error "Must be specialized"
"""
	weak_norm(B::Basis)

Return the type of the weak norm of the basis

# Example

```jldoctest
julia> using RigorousInvariantMeasures

julia> B = Ulam(1024)
Ulam{LinRange{Float64, Int64}}(LinRange{Float64}(0.0, 1.0, 1025))

julia> weak_norm(B)
L1
```
"""
weak_norm(B::Basis) = @error "Must be specialized"
aux_norm(B::Basis) = @error "Must be specialized"

"""
	is_refinement(Bfine::Basis, Bcoarse::Basis)
Check if Bfine is a refinement of Bcoarse

# Example

```jldoctest
julia> using RigorousInvariantMeasures

julia> B = Ulam(1024)
Ulam{LinRange{Float64, Int64}}(LinRange{Float64}(0.0, 1.0, 1025))

julia> Bfine = Ulam(2048)
Ulam{LinRange{Float64, Int64}}(LinRange{Float64}(0.0, 1.0, 2049))

julia> is_refinement(Bfine, B)
true

julia> Bfine = Ulam(2049)
Ulam{LinRange{Float64, Int64}}(LinRange{Float64}(0.0, 1.0, 2050))

julia> is_refinement(Bfine, B)
false
```

"""
is_refinement(Bfine::Basis, Bcoarse::Basis) = @error "Not Implemented"

"""
	integral_covector(B::Basis)

Return a covector that represents the integral in the basis B
"""
integral_covector(B::Basis) = @error "Must be specialized"

"""
	one_vector(B::Basis)

Vector that represents the function 1 in the basis B
"""
one_vector(B::Basis) = @error "Must be specialized"

"""
	is_integral_preserving(B::Basis)

Integral-preserving discretizations may specialize this to "true"
"""
is_integral_preserving(B::Basis) = false

"""
	integral(B::Basis, v; T = Float64)

Return the integral of the function with coefficients v in the basis B 
"""
function integral(B::Basis, v; T = Float64)
    return sum([T(v[i]) * evaluate_integral(B, i, T) for i = 1:length(B)])
end

"""
	AverageZero{B<:Basis}

Yield a basis of the space of average zero vectors
"""
struct AverageZero{B<:Basis}
    basis::B
end

Base.iterate(S::AverageZero{B}, state) where {B} = @error "Not Implemented"
Base.length(S::AverageZero{T}) where {T} = length(S.basis) - 1

"""
	weak_projection_error(B::Basis)

Return a constant Kh (typically scales as h ~ 1/n) such that 
```||P_h f-f||\\leq Kh ||f||_s```
Must be rounded up correctly!
This function is not exported explictly but is used in all the estimates.

# Example

```jldoctest
julia> using RigorousInvariantMeasures;

julia> B = Ulam(1024)
Ulam{LinRange{Float64, Int64}}(range(0.0, stop=1.0, length=1025))

julia> RigorousInvariantMeasures.weak_projection_error(B)

0.00048828125
````
"""
weak_projection_error(B::Basis) = @error "Not Implemented"

"""
	aux_normalized_projection_error(B::Basis)

Return a constant Eh (typically scales as h ~ 1/n) such that 
```|||P_h f|||\\leq |||f|||+ Eh * ||f||_s```
Must be rounded up correctly!
"""
aux_normalized_projection_error(B::Basis) = @error "Not Implemented"

"""
	strong_weak_bound(Bstruct HatDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    x′::Vector{Interval}
end

Dual(B::Hat, D; ϵ, max_iter) =
    HatDual(preimages_and_derivatives(B.p, D, 1:length(B.p)-1; ϵ, max_iter)...)
Base.length(dual::HatDual) = length(dual.x)
Base.eltype(dual::HatDual) =
    Tuple{eltype(dual.xlabel),Tuple{eltype(dual.x),eltype(dual.x′)}}
function Base.iterate(dual::HatDual, state = 1)
    if state <= length(dual.x)
        return ((dual.xlabel[state], (dual.x[state], abs(dual.x′[state]))), state + 1)
    else
        return nothing
    end
end::Basis)
Return a constant ``M₁n`` such that for a vector ``v ∈ Uₕ`` 
```||v||_s\\leq M1n*||v||```
Must be rounded up correctly!
"""
strong_weak_bound(B::Basis) = @error "Not Implemented"

"""
	aux_weak_bound(B::Basis)
Return a constant ``M₂`` such that for a vector ``v ∈ Uₕ`` 
```|||v|||\\leq M_2||v||```
Must be rounded up correctly!
"""
aux_weak_bound(B::Basis) = @error "Not Implemented"

"""
	Return constants ``S₁, S₂`` such that for a vector ``v ∈ Uₕ`` 
	```||v||\\leq S_1||v||_s+S_2|||v|||```
Must be rounded up correctly!
"""
weak_by_strong_and_aux_bound(B::Basis) = @error "Not Implemented"

"""
	bound_weak_norm_from_linalg_norm(B::Basis)
Return constants W₁, W₂ such that for a vector ``v ∈ Uₕ`` 
```||v||\\leq W_1||v||_1+W_2||v||_{\\infty}```
Must be rounded up correctly!
"""
bound_weak_norm_from_linalg_norm(B::Basis) = @error "Not Implemented"

"""
	bound_linalg_norm_L1_from_weak(B::Basis)

Return a constant ``A`` such that for a vector ``v ∈ Uₕ`` 
```||v||_1\\leq A||v||```
Must be rounded up correctly!
"""
bound_linalg_norm_L1_from_weak(B::Basis) = @error "Not Implemented"

"""
	bound_linalg_norm_L∞_from_weak(B::Basis)

Return a constant ``A`` such that for a vector ``v ∈ Uₕ`` 
```||v||_\\infty \\leq A||v||```
Must be rounded up correctly!
"""
bound_linalg_norm_L∞_from_weak(B::Basis) = @error "Not Implemented"

"""
	invariant_measure_strong_norm_bound(B::Basis, D::Dynamic)

Bounds ``||u||_s``, where ``u`` is the invariant measure normalized with 
``i(u)=1``.
"""
invariant_measure_strong_norm_bound(B::Basis, D::Dynamic) = @error "Must be specialized"


"""
	bound_weak_norm_abstract(B::Basis, D=nothing; dfly_coefficients=nothing)

Returns an a priori bound on the weak norm of the abstract operator ``L``
"""
bound_weak_norm_abstract(B::Basis, D = nothing; dfly_coefficients = nothing) =
    @error "Must be specialized"

using ..RigorousInvariantMeasures: NormKind
opnormbound(B::Basis, N::Type{<:NormKind}, M::AbstractVecOrMat{S}) where {S} =
    @error "Must be specialized"
normbound(B::Basis, N::Type{<:NormKind}, v) = @error "Must be specialized"

# careful, this is defined outside the module!!!

"""
Replacement of DualComposedWithDynamic.
"""
abstract type Dual end
