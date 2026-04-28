export Observable,
    ProjectedFunction,
    integrateobservable,
    integral_pairing,
    weak_dual_norm_bound

@doc raw"""
    Observable{TB,TV,TIB,TPE}

Discretization of an observable ``ϕ`` on a basis `B`.

Fields:

- `B` — the basis.
- `v` — the basis-specific coefficient vector.
- `inf_bound` — an upper bound on ``\|ϕ\|_{w^*}`` (the dual of the basis's
  weak norm). For weak `L²` (the default for `FourierAnalytic` /
  `FourierAdjoint`) this is ``\|ϕ\|_{L²}``; for weak `L¹` (`Ulam`)
  this is ``\|ϕ\|_{L^∞}``.
- `proj_error` — an upper bound on the **weak-norm** projection error
  ``\|ϕ - \tilde P_N ϕ\|_w`` of the observable itself (analog of
  `ProjectedFunction.err_bound` but for the observable side). May be
  `nothing` if the constructor wasn't given enough information to bound
  it — e.g. the legacy 3-arg `Observable(B, v, inf_bound)` call.

The intended use is the integration estimate

```math
\Big| \int ϕ\,ρ\,dx - \int ϕ_N\,ρ_N\,dx \Big|
\;\leq\; \|ϕ - ϕ_N\|_w\,\|ρ\|_{w^*}
       + \|ϕ_N\|_{w^*}\,\|ρ - ρ_N\|_w
```

so a caller bounds the integration error using `proj_error · ‖ρ‖_{w*}`
plus `inf_bound · ProjectedFunction(ρ).err_bound`.

Constructors are provided per basis by extensions:

- `Observable(B::Ulam, ϕ; tol)` — `TaylorModelsExt` (load with `using TaylorModels`).
- `Observable(B::FourierAnalytic{W{k,l},…}, ϕ; inf_bound, Wk1_seminorm, oversample)`
  and `Observable(B::FourierAdjoint{W{k,l},…}, …)` — `FFTWExt`
  (load with `using FFTW`).

The struct itself has no inherent dependency on either extension; concrete
field types are inferred from the constructor.
"""
struct Observable{TB<:Basis,TV<:AbstractVector,TIB,TPE}
    B::TB
    v::TV
    inf_bound::TIB
    proj_error::TPE
end

# Legacy 3-arg form (`proj_error` left as `nothing`). Used by extensions that
# don't yet compute the observable-side projection error (e.g. the Ulam
# constructor in `TaylorModelsExt`).
Observable(B::Basis, v::AbstractVector, inf_bound) =
    Observable(B, v, inf_bound, nothing)

@doc raw"""
    ProjectedFunction{TB,TV,TEB}

A function `f` projected onto a basis. `v` is the basis-specific coefficient
vector and `err_bound` is an upper bound on the L¹ projection error
``\|f - P_N f\|_{L^1}``.

Constructors are provided per basis by extensions, with the same scheme as
`Observable`. The Fourier extensions use the `Wk1_seminorm` kwarg for
``\|f^{(k)}\|_{L^1}`` (W^{k,1} seminorm); the Ulam extension uses
`VarBound` (a total-variation bound on `f`).
"""
struct ProjectedFunction{TB<:Basis,TV<:AbstractVector,TEB}
    B::TB
    v::TV
    err_bound::TEB
end

"""
    integrateobservable(B, ϕ::Observable, f, error)

Integrate the discretized observable ``ϕ`` against a density coefficient vector
`f`, with an upper bound on the basis-side projection error `error`. Legacy
two-term bound; prefer [`integral_pairing`](@ref). Methods provided by
extensions.
"""
function integrateobservable end

"""
    weak_dual_norm_bound(B::Basis, v::AbstractVector)

Upper bound on ``\\|φ_v\\|_{w^*}`` where ``φ_v`` is the basis-`B` reconstruction
of the coefficient vector `v` and ``w`` is `weak_norm(B)`. Methods are
provided per basis (`Ulam`, `Fourier`, …); extensions add more.
"""
function weak_dual_norm_bound end

@doc raw"""
    integral_pairing(ϕ::Observable, ρ::AbstractVector, ρ_w_error;
                     ρ_dual_weak_bound = weak_dual_norm_bound(ϕ.B, ρ))
    integral_pairing(ϕ::Observable, ρ::ProjectedFunction)

Compute a rigorous enclosure of ``\int_0^1 ϕ(x)\,ρ(x)\,dx``, where `ϕ` is an
`Observable` discretized on basis `ϕ.B` and `ρ` is a density coefficient
vector in the same basis (with `ρ_w_error` an upper bound on
``\|ρ - ρ_N\|_w``).

The error is decomposed as

```math
\Bigl|\int (ϕ\,ρ - ϕ_N\,ρ_N)\,dx\Bigr|
\;\leq\; \|ϕ - ϕ_N\|_w\,\|ρ_N\|_{w^*}
       + \|ϕ\|_{w^*}\,\|ρ - ρ_N\|_w
```

so the result is

    ⟨ϕ.v, ρ⟩_basis
        + [-(ϕ.proj_error · ρ_dual_weak_bound + ϕ.inf_bound · ρ_w_error),
           +(…)]

`ρ_dual_weak_bound` defaults to `weak_dual_norm_bound(ϕ.B, ρ)` (computed
from the finite-dimensional coefficient vector); pass a tighter
user-provided bound to override.

Requires `ϕ.proj_error !== nothing` — i.e. an `Observable` constructed with
enough information to bound ``\|ϕ - ϕ_N\|_w``. Methods are provided per
basis.
"""
function integral_pairing end

integral_pairing(ϕ::Observable, ρ::ProjectedFunction; kwargs...) =
    integral_pairing(ϕ, ρ.v, ρ.err_bound; kwargs...)

# Lift a real- or complex-typed coefficient vector to interval-typed so dot
# products in `integral_pairing` run in interval arithmetic and the result
# is a rigorous enclosure of the floating-point rounding error.
_lift_to_interval(v::AbstractVector{<:Real}) = interval.(v)
_lift_to_interval(v::AbstractVector{<:Complex{<:Real}}) =
    [interval(real(z)) + im * interval(imag(z)) for z in v]
_lift_to_interval(v::AbstractVector) = v  # Already interval-typed.
