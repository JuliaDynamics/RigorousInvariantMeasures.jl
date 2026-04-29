export ProjectedFunction,
    Observable,
    integrateobservable,
    integral_pairing,
    weak_dual_norm_bound

@doc raw"""
    ProjectedFunction{TB,TV,TWDB,TPE}

A function `f` projected onto a basis `B`, packaged together with the
two bounds needed for rigorous integration against another projected
function:

- `B` — the basis.
- `v` — the basis-specific discrete coefficient vector.
- `weak_dual_bound` — upper bound on ``\|f\|_{w^*}`` (dual of the weak
  norm of `B`). For weak `L²` (`Fourier`) this is ``\|f\|_{L²}``; for
  weak `L¹` (`Ulam`) this is ``\|f\|_{L^∞}``.
- `proj_error` — upper bound on the weak-norm projection error
  ``\|f - φ_v\|_w``.

The integration formula

```math
\Bigl|\int f\,g\,dx - \langle f_N, g_N\rangle\Bigr|
  \;\leq\; p_f.\text{proj\_error} \cdot \|p_g.v\|_{w^*}
        + p_f.\text{weak\_dual\_bound} \cdot p_g.\text{proj\_error}
```

is computed by [`integral_pairing`](@ref). For `Ulam` (cell-wise exactness)
the first term is identically zero — see the Ulam method docstring.

`Observable` is provided as a `const` alias to `ProjectedFunction`: prior
code using `Observable(B, ϕ; …)` keeps working. The two roles (observable
vs density) are now distinguished only by which arguments the caller
emphasizes, not by the type system.

Constructors are provided per basis by extensions:

- `ProjectedFunction(B::Ulam, f; tol, var_bound, weak_dual_bound)` —
  `TaylorModelsExt` (load with `using TaylorModels`). `weak_dual_bound`
  defaults to a Taylor-model bound on ``\|f\|_{L^∞}``; `var_bound`
  defaults to a `VariationBound`-computed total variation; `proj_error
  = var_bound / length(B)`.
- `ProjectedFunction(B::FourierAnalytic{W{k,l},…}, f;
   weak_dual_bound, Wk1_seminorm, M = 4*B.k)` and
  `ProjectedFunction(B::FourierAdjoint{W{k,l},…}, …)` —
  `FFTWExt` (load with `using FFTW`).

The struct itself has no inherent dependency on either extension; concrete
field types are inferred from the constructor.
"""
struct ProjectedFunction{TB<:Basis,TV<:AbstractVector,TWDB,TPE}
    B::TB
    v::TV
    weak_dual_bound::TWDB
    proj_error::TPE
end

"""
    Observable

Alias for [`ProjectedFunction`](@ref). The two were separate structs in
earlier versions of the package; the unified type carries both an
observable-side (`weak_dual_bound`) and a projection-side (`proj_error`)
bound, so a single struct serves both roles.
"""
const Observable = ProjectedFunction

# Legacy 3-arg form (`proj_error` left as `nothing`). Used by callers that
# build a discrete object without a known projection-error bound.
ProjectedFunction(B::Basis, v::AbstractVector, weak_dual_bound) =
    ProjectedFunction(B, v, weak_dual_bound, nothing)

"""
    integrateobservable(B, ϕ::ProjectedFunction, f, error)

Integrate the discretized observable `ϕ` against a density coefficient vector
`f`, with an upper bound `error` on the basis-side density projection error.
Legacy two-term bound; prefer [`integral_pairing`](@ref). Methods provided by
extensions.
"""
function integrateobservable end

"""
    weak_dual_norm_bound(B::Basis, v::AbstractVector)

Upper bound on ``\\|φ_v\\|_{w^*}`` where ``φ_v`` is the basis-`B` reconstruction
of the coefficient vector `v` and ``w`` is `weak_norm(B)`.
"""
function weak_dual_norm_bound end

@doc raw"""
    integral_pairing(ϕ::ProjectedFunction, ρ::ProjectedFunction)
    integral_pairing(ϕ::ProjectedFunction, ρ::AbstractVector, ρ_w_error;
                     ρ_dual_weak_bound = weak_dual_norm_bound(ϕ.B, ρ))

Compute a rigorous enclosure of ``\int_0^1 ϕ(x)\,ρ(x)\,dx``.

The error decomposes as

```math
\Bigl|\int ϕ\,ρ - \langle ϕ_N,\,ρ_N\rangle\Bigr|
\;\leq\; \|ϕ - ϕ_N\|_w\,\|ρ_N\|_{w^*}
       + \|ϕ\|_{w^*}\,\|ρ - ρ_N\|_w
```

so the result is

    ⟨ϕ.v, ρ.v⟩_basis
        + [-(ϕ.proj_error · ‖ρ.v‖_{w*}
             + ϕ.weak_dual_bound · ρ.proj_error),
           +(…)]

The two-vector form lets callers pass a density that wasn't constructed
through this package (e.g. a result of `invariant_vector`); they supply
`ρ_w_error` and optionally a tighter `ρ_dual_weak_bound` than what
`weak_dual_norm_bound(ϕ.B, ρ)` would compute.

Methods are provided per basis. For `Ulam` the cell-wise integration of
`ϕ` is exact against any piecewise-constant `ρ_N`, so the first term
collapses and only `ϕ.weak_dual_bound · ρ.proj_error` contributes.
"""
function integral_pairing end

# Two-PF form: unpack the second one into the vector form.
integral_pairing(ϕ::ProjectedFunction, ρ::ProjectedFunction; kwargs...) =
    integral_pairing(ϕ, ρ.v, ρ.proj_error; kwargs...)

# Lift a real- or complex-typed coefficient vector to interval-typed so dot
# products in `integral_pairing` run in interval arithmetic and the result
# is a rigorous enclosure of the floating-point rounding error.
_lift_to_interval(v::AbstractVector{<:Real}) = interval.(v)
_lift_to_interval(v::AbstractVector{<:Complex{<:Real}}) =
    [interval(real(z)) + im * interval(imag(z)) for z in v]
_lift_to_interval(v::AbstractVector) = v  # Already interval-typed.

# ---------------------------------------------------------------------------
# Algebraic operations on `ProjectedFunction`
# ---------------------------------------------------------------------------

# Combine two bounds rigorously. If either operand is `nothing` (the legacy
# 3-arg constructor leaves `proj_error` unspecified), the result is also
# `nothing` — the caller hasn't supplied enough information.
_add_bound(::Nothing, _) = nothing
_add_bound(_, ::Nothing) = nothing
_add_bound(::Nothing, ::Nothing) = nothing
_add_bound(a, b) = interval(a) + interval(b)

_mul_bound(::Nothing, _) = nothing
_mul_bound(_, ::Nothing) = nothing
_mul_bound(::Nothing, ::Nothing) = nothing
_mul_bound(a, b) = interval(a) * interval(b)

@doc raw"""
    p1 + p2
    p1 - p2

Add or subtract two `ProjectedFunction`s on the same basis. The
coefficient vector is combined componentwise; both bounds combine via
the triangle inequality:

    weak_dual_bound = p1.weak_dual_bound + p2.weak_dual_bound
    proj_error      = p1.proj_error      + p2.proj_error

Subtraction uses the same combine rule for bounds (still triangle
inequality on the dual norm and the weak norm). Bounds combine in
interval arithmetic; if either input has `proj_error == nothing` the
result inherits `nothing` for that field.
"""
function Base.:+(p1::ProjectedFunction, p2::ProjectedFunction)
    @assert length(p1.v) == length(p2.v) "Sum requires same coefficient length"
    return ProjectedFunction(
        p1.B,
        p1.v .+ p2.v,
        _add_bound(p1.weak_dual_bound, p2.weak_dual_bound),
        _add_bound(p1.proj_error, p2.proj_error),
    )
end

function Base.:-(p1::ProjectedFunction, p2::ProjectedFunction)
    @assert length(p1.v) == length(p2.v) "Subtraction requires same coefficient length"
    return ProjectedFunction(
        p1.B,
        p1.v .- p2.v,
        _add_bound(p1.weak_dual_bound, p2.weak_dual_bound),
        _add_bound(p1.proj_error, p2.proj_error),
    )
end
