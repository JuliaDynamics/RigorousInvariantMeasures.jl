@doc raw"""
Constructors for [`Observable`](@ref) and [`ProjectedFunction`](@ref) on
Fourier bases. They use the FFT to compute the coefficient vector and a
user-supplied seminorm to bound the discretization error in the basis's
weak norm.

Observables and projected functions may live in a more regular space than
the one preserved by the transfer operator; the user supplies the bound on
the appropriate seminorm of the input. For a `Fourier{W{k,l},…}` basis
with ``k ≥ 2``, the user passes `Wk1_seminorm` — an upper bound on
``\|f^{(k)}\|_{L^1}``.

For both `ProjectedFunction.err_bound` and `Observable.proj_error`, the
quantity returned is an upper bound on the **weak-norm** projection error
``\|f - \tilde P_N f\|_w`` (with weak norm `weak_norm(B)`); for the `L²`
weak norm used by both `FourierAnalytic` and `FourierAdjoint`, this is
the L² norm.

`Observable.inf_bound` is an upper bound on the observable in the **dual
of the weak norm** — the quantity that controls the integration error
``|\int (\phi - \tilde P_N \phi)\,d\mu|`` against an invariant density
``\mu`` expressed in the weak norm. For weak L², this is ``\|ϕ\|_{L²}``.

# Aliasing control

The DFT computes ``\tilde f_n = \hat f_n + \sum_{m\neq 0} \hat f_{n+mM}``
where ``M`` is the FFT grid size. To control the aliasing error, callers
choose `M` via a kwarg (default `M = 4 * B.k`, i.e. four times the basis
cutoff frequency). With grid size `M`, the FFT samples at `M` equispaced
points, the central `length(B) = 2*B.k + 1` Fourier modes are kept, and
the aliasing bound uses `M` in place of `length(B)`. Larger `M` divides
the per-coefficient aliasing bound by `(M / length(B))^k`; for the default
`M = 4 * B.k`, that is roughly `2^k`.

The aliasing bound used (loose but explicit): for ``k ≥ 2``,

```math
α \;\leq\; S\,\frac{2π^2}{3}\,\Bigl(\frac{1}{π M}\Bigr)^k
```

(``S = `` `Wk1_seminorm`). The L² truncation tail

```math
\|f - P_N f\|_{L^2}
\;\leq\; S\,\sqrt{\tfrac{2}{2k-1}}\,\Bigl(\tfrac{1}{2π}\Bigr)^k
\,\Bigl(\tfrac{1}{k_\text{freq}}\Bigr)^{k-½}
```

does not depend on the FFT grid size — only on the kept basis modes
``k_\text{freq} = `` `B.k`.
"""

using RigorousInvariantMeasures: FourierPoints

# -- Per-coefficient aliasing inflation for f ∈ W^{k,1} on an M-point FFT grid.
function _fourier_aliasing_bound_Wk1(M::Integer, k::Integer, seminorm)
    k ≥ 2 ||
        throw(ArgumentError("W^{k,1} aliasing bound requires k ≥ 2; got k = $k"))
    Si = interval(Float64(seminorm))
    Mi = interval(Float64(M))
    return sup(Si * (2 * interval(pi)^2 / 3) / (interval(pi) * Mi)^k)
end

# -- L² truncation-tail bound for f ∈ W^{k,1} truncated at frequency k_freq.
function _fourier_L2_tail_bound_Wk1(k_freq::Integer, k::Integer, seminorm)
    k ≥ 2 ||
        throw(ArgumentError("W^{k,1} truncation bound requires k ≥ 2; got k = $k"))
    Si = interval(Float64(seminorm))
    kfi = interval(Float64(k_freq))
    inner = 2 * Si^2 / (2 * interval(pi))^(2k) /
            (interval(Float64(2k - 1)) * kfi^(2k - 1))
    return sup(sqrt(inner))
end

# -- Common path: oversampled FFT, truncate to basis modes, inflate by aliasing.
function _fourier_coeffs_with_aliasing(
    B::Fourier,
    f::Function,
    k::Integer,
    seminorm,
    M::Integer,
)
    M ≥ length(B) ||
        throw(ArgumentError("FFT grid size M = $M must be ≥ length(B) = $(length(B))"))
    samples = [f(p) for p in FourierPoints(M, Float64)]
    raw = interval_fft(samples)               # length M
    k_freq = B.k
    coeffs =
        M == length(B) ? raw :
        [raw[1:k_freq+1]; raw[M-k_freq+1:M]]
    α = _fourier_aliasing_bound_Wk1(M, k, seminorm)
    α_box = interval(-α, α) + im * interval(-α, α)
    return coeffs .+ Ref(α_box)
end

# ---------------------------------------------------------------------------
# Observable / ProjectedFunction for Fourier{W{k,l}, …}
# ---------------------------------------------------------------------------

@doc raw"""
    Observable(B::FourierAnalytic{W{k,l},…}, ϕ::Function;
               inf_bound, Wk1_seminorm, M = 4 * B.k)
    Observable(B::FourierAdjoint{W{k,l},…}, ϕ::Function;
               inf_bound, Wk1_seminorm, M = 4 * B.k)

Discretize an observable `ϕ` on a Fourier basis whose strong norm is
``W^{k,1}`` (``k ≥ 2``).

`Wk1_seminorm` is an upper bound on ``\|ϕ^{(k)}\|_{L^1}``. `inf_bound`
is an upper bound on `ϕ` in the dual of the basis's weak norm. `M` is
the FFT grid size; default `4 * B.k` gives a roughly `2^k` reduction in
aliasing relative to the un-oversampled `M = length(B)`. The
`proj_error` field of the result is the L² truncation tail bound
``\|ϕ - P_N ϕ\|_{L²}`` derived from `Wk1_seminorm`.
"""
function RigorousInvariantMeasures.Observable(
    B::FourierAnalytic{W{k,l}},
    ϕ::Function;
    inf_bound,
    Wk1_seminorm,
    M::Integer = 4 * B.k,
) where {k,l}
    coeffs = _fourier_coeffs_with_aliasing(B, ϕ, k, Wk1_seminorm, M)
    proj_err = _fourier_L2_tail_bound_Wk1(B.k, k, Wk1_seminorm)
    return RigorousInvariantMeasures.Observable(B, coeffs, inf_bound, proj_err)
end

function RigorousInvariantMeasures.Observable(
    B::FourierAdjoint{W{k,l}},
    ϕ::Function;
    inf_bound,
    Wk1_seminorm,
    M::Integer = 4 * B.k,
) where {k,l}
    coeffs = _fourier_coeffs_with_aliasing(B, ϕ, k, Wk1_seminorm, M)
    proj_err = _fourier_L2_tail_bound_Wk1(B.k, k, Wk1_seminorm)
    return RigorousInvariantMeasures.Observable(B, coeffs, inf_bound, proj_err)
end

@doc raw"""
    ProjectedFunction(B::FourierAnalytic{W{k,l},…}, f::Function;
                      Wk1_seminorm, M = 4 * B.k)
    ProjectedFunction(B::FourierAdjoint{W{k,l},…}, f::Function;
                      Wk1_seminorm, M = 4 * B.k)

Project `f` onto a Fourier basis with strong norm ``W^{k,1}``
(``k ≥ 2``). `Wk1_seminorm` is an upper bound on ``\|f^{(k)}\|_{L^1}``.
`M` is the FFT grid size (see the module docstring for the aliasing
formula).

`err_bound` is an upper bound on the **weak-norm** projection error
``\|f - \tilde P_N f\|_{w}``; for `L²` weak, that is the L² truncation
tail derived above.
"""
function RigorousInvariantMeasures.ProjectedFunction(
    B::FourierAnalytic{W{k,l}},
    f::Function;
    Wk1_seminorm,
    M::Integer = 4 * B.k,
) where {k,l}
    coeffs = _fourier_coeffs_with_aliasing(B, f, k, Wk1_seminorm, M)
    err = _fourier_L2_tail_bound_Wk1(B.k, k, Wk1_seminorm)
    return RigorousInvariantMeasures.ProjectedFunction(B, coeffs, err)
end

function RigorousInvariantMeasures.ProjectedFunction(
    B::FourierAdjoint{W{k,l}},
    f::Function;
    Wk1_seminorm,
    M::Integer = 4 * B.k,
) where {k,l}
    coeffs = _fourier_coeffs_with_aliasing(B, f, k, Wk1_seminorm, M)
    err = _fourier_L2_tail_bound_Wk1(B.k, k, Wk1_seminorm)
    return RigorousInvariantMeasures.ProjectedFunction(B, coeffs, err)
end

# `projection` dispatches to `ProjectedFunction` for any Fourier basis;
# the concrete subtype's constructor (W{k,l} above) handles the math
# and the kwargs. Bases without a defined `ProjectedFunction` constructor
# (e.g. `Aη`/`Cω` strong norm — not yet implemented) raise `MethodError`.
RigorousInvariantMeasures.projection(B::Fourier, f::Function; kwargs...) =
    RigorousInvariantMeasures.ProjectedFunction(B, f; kwargs...)
