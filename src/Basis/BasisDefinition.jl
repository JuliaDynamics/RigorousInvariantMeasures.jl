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
    projection,
    Dual

abstract type Basis end

# `Base.length(::Basis)` is left to the concrete basis types — Julia's
# native `MethodError` if a subtype forgets to provide one is the right
# signal (was an `@error "Not Implemented"` stub that silently returned
# nothing, see also the migration note for the other stubs below).

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

function is_dual_element_empty end
function nonzero_on end

"""
    projection(B::Basis, f::Function; kwargs...)

Project a function `f` onto the basis `B`. The result type and available
keyword arguments are basis-dependent. Currently implemented via the
`TaylorModels` extension for the `Ulam` basis.
"""
function projection end

function ProjectDualElement(B::Basis, y)
    j_min, j_max = nonzero_on(B, y)
    return ProjectDualElement(B, j_min, j_max, y)
end

# Base.iterate(S::ProjectDualElement, state) = @error "Not Implemented"

"""
	evaluate(B::Basis, i, x)

Evaluate the i-th basis element at x
"""
function evaluate end

"""
	evaluate_integral(B::Basis, i; T = Float64)

Value of the integral on [0,1] of the i-th basis element
"""
function evaluate_integral end

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
function strong_norm end

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
function weak_norm end

function aux_norm end

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
function is_refinement end

"""
	integral_covector(B::Basis)

Return a covector that represents the integral in the basis B
"""
function integral_covector end

"""
	one_vector(B::Basis)

Vector that represents the function 1 in the basis B
"""
function one_vector end

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

# `Base.iterate(::AverageZero, state)` is left to the concrete bases.
Base.length(S::AverageZero{T}) where {T} = length(S.basis) - 1

"""
	weak_projection_error(B::Basis)

Return a constant Kh (typically scales as h ~ 1/n) such that 

``||P_h f-f||\\leq Kh ||f||_s``

Must be rounded up correctly!
This function is not exported explictly but is used in all the estimates.

# Example

```jldoctest
julia> using RigorousInvariantMeasures;

julia> B = Ulam(1024)
Ulam{LinRange{Float64, Int64}}(LinRange{Float64}(0.0, 1.0, 1025))

julia> RigorousInvariantMeasures.weak_projection_error(B)
0.00048828125
```
"""
function weak_projection_error end

"""
	aux_normalized_projection_error(B::Basis)

Return a constant Eh (typically scales as h ~ 1/n) such that

``|||P_h f|||\\leq |||f|||+ Eh * ||f||_s``

Must be rounded up correctly!
"""
function aux_normalized_projection_error end

"""
    strong_weak_bound(B::Basis)
Return a constant ``M₁n`` such that for a vector ``v ∈ Uₕ``

``||v||_s\\leq M1n*||v||``

Must be rounded up correctly!
"""
function strong_weak_bound end

"""
	aux_weak_bound(B::Basis)
Return a constant ``M₂`` such that for a vector ``v ∈ Uₕ``

``|||v|||\\leq M_2||v||``

Must be rounded up correctly!
"""
function aux_weak_bound end

"""
Return constants ``S₁, S₂`` such that for a vector ``v ∈ Uₕ``

``||v||\\leq S_1||v||_s+S_2|||v|||``

Must be rounded up correctly!
"""
function weak_by_strong_and_aux_bound end

"""
	bound_weak_norm_from_linalg_norm(B::Basis)
Return constants W₁, W₂ such that for a vector ``v ∈ Uₕ``

``||v||\\leq W_1||v||_1+W_2||v||_{\\infty}``

Must be rounded up correctly!
"""
function bound_weak_norm_from_linalg_norm end

@doc raw"""
	bound_linalg_norm_L1_from_weak(B::Basis)

Return a constant ``A`` such that for a vector ``v ∈ Uₕ``

``||v||_1\leq A||v||``

Must be rounded up correctly!
"""
function bound_linalg_norm_L1_from_weak end

"""
	bound_linalg_norm_L∞_from_weak(B::Basis)

Return a constant ``A`` such that for a vector ``v ∈ Uₕ``
```||v||_\\infty \\leq A||v||```
Must be rounded up correctly!
"""
function bound_linalg_norm_L∞_from_weak end

"""
	invariant_measure_strong_norm_bound(B::Basis, D::Dynamic)

Bounds ``||u||_s``, where ``u`` is the invariant measure normalized with
``i(u)=1``.
"""
function invariant_measure_strong_norm_bound end


"""
	bound_weak_norm_abstract(B::Basis, D=nothing; dfly_coefficients=nothing)

Returns an a priori bound on the weak norm of the abstract operator ``L``
"""
function bound_weak_norm_abstract end

# `opnormbound` and `normbound` are first defined in `src/NormBounds.jl` (which
# is `included` before this file). We don't need a per-`Basis` fallback method
# here — concrete bases provide their own three-arg specializations, and any
# call without a matching method now raises a `MethodError` instead of the
# previous silent `@error "Must be specialized"`.

using ..RigorousInvariantMeasures: NormKind

# careful, this is defined outside the module!!!

"""
Replacement of DualComposedWithDynamic.
"""
abstract type Dual end
