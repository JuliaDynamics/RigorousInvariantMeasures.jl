@doc raw"""
Constructors for [`Observable`](@ref) and [`ProjectedFunction`](@ref) on
Fourier bases. They use the FFT to compute the coefficient vector and a
user-supplied seminorm to bound the discretization error in the basis's
weak norm.

For a `Fourier{W{k,l},…}` basis with ``k ≥ 1`` (so this also covers BV /
``W^{1,1}``), the user passes `Wk1_seminorm` — an upper bound on
``\|f^{(k)}\|_{L^1}`` (for ``k=1`` this is the BV seminorm, the total
variation). The coefficient bound ``|\hat f_n| \leq S/(2π|n|)^k`` (where
``S`` is the seminorm) drives every error formula on this page.

# Error model

Let ``\varphi_N(x) = \sum_{|n|\leq k_\text{freq}} \tilde f_n e^{2π i n x}``
be the discrete Fourier reconstruction returned by the constructor (the
``\tilde f_n`` are the FFT-computed coefficients, encoded as complex
intervals carrying only the FFT roundoff bound from `interval_fft`).
Decompose

```math
f - \varphi_N \;=\; \underbrace{(f - P_N f)}_{\text{truncation}}
\;+\; \underbrace{(P_N f - \varphi_N)}_{\text{aliasing}}.
```

Both pieces have closed-form L² bounds:

```math
\|f - P_N f\|_{L^2}^2
  \;=\; \sum_{|n|>k_\text{freq}} |\hat f_n|^2
  \;\leq\; \frac{2 S^2}{(2π)^{2k}\,(2k-1)\,k_\text{freq}^{2k-1}}
  \;=:\; T_1^2,
```

```math
\|P_N f - \varphi_N\|_{L^2}^2
  \;\leq\; \|f - P_M f\|_{L^2}^2
  \;\leq\; \frac{2 S^2}{(2π)^{2k}\,(2k-1)\,(M/2)^{2k-1}}
  \;=:\; T_2^2,
```

(the aliasing energy on the kept ``|n| \leq k_\text{freq}`` modes is
bounded by Parseval-aliasing's total energy, which is the M-point
truncation tail). Combining,

```math
\|f - \varphi_N\|_{L^2} \;\leq\; \sqrt{T_1^2 + T_2^2}.
```

This works for any ``k ≥ 1``. For oversampling factor ``s = M/N \geq 1``
the aliasing tail ``T_2`` shrinks like ``s^{-(k-½)}`` relative to
``T_1``; the default `M = 4 * B.k` halves the relative aliasing
contribution.

# Coefficient interpretation

`v[n]` encloses the **FFT-computed** coefficient ``\tilde f_n`` (with
FFT roundoff but no aliasing inflation). It does **not** enclose the
true continuous Fourier coefficient ``\hat f_n``. Aliasing is folded
into `err_bound` / `proj_error` at the L² level — appropriate for
integration via [`integral_pairing`](@ref), but if you need an enclosure
of ``\hat f_n`` itself you have to inflate `v[n]` separately.
"""

using RigorousInvariantMeasures: FourierPoints

# T₁: truncation-tail L² bound for f ∈ W^{k,1} truncated at frequency k_freq.
# Works for any k ≥ 1.
function _fourier_L2_tail_bound_Wk1(k_freq::Integer, k::Integer, seminorm)
    k ≥ 1 ||
        throw(ArgumentError("W^{k,1} truncation bound requires k ≥ 1; got k = $k"))
    Si = interval(Float64(seminorm))
    kfi = interval(Float64(k_freq))
    inner = 2 * Si^2 / (2 * interval(pi))^(2k) /
            (interval(Float64(2k - 1)) * kfi^(2k - 1))
    return sup(sqrt(inner))
end

# Combined L² error: √(T₁² + T₂²), where T₁ is the truncation tail at k_freq
# and T₂ is the M-grid truncation tail bounding the L² aliasing energy.
function _fourier_L2_total_error_Wk1(
    k_freq::Integer,
    k::Integer,
    seminorm,
    M_grid::Integer,
)
    k ≥ 1 ||
        throw(ArgumentError("W^{k,1} total-error bound requires k ≥ 1; got k = $k"))
    Si = interval(Float64(seminorm))
    kfi = interval(Float64(k_freq))
    half_grid = interval(Float64(M_grid ÷ 2))
    T1_sq = 2 * Si^2 / (2 * interval(pi))^(2k) /
            (interval(Float64(2k - 1)) * kfi^(2k - 1))
    T2_sq = 2 * Si^2 / (2 * interval(pi))^(2k) /
            (interval(Float64(2k - 1)) * half_grid^(2k - 1))
    return sup(sqrt(T1_sq + T2_sq))
end

# Sample, FFT, take the central `length(B)` coefficients (in standard FFT
# layout: DC at index 1, positives, then negatives). No per-coefficient
# aliasing inflation — aliasing is folded into the L² total-error bound
# returned alongside.
function _fourier_coeffs_no_aliasing(B::Fourier, f::Function, M_grid::Integer)
    M_grid ≥ length(B) || throw(
        ArgumentError(
            "FFT grid size M = $M_grid must be ≥ length(B) = $(length(B))",
        ),
    )
    # The transform is taken on a larger grid, of size a power of two, for which the
    # error bound of the transform is proved; only the frequencies
    # 0, 1, …, K, -K, …, -1 are kept. The caller's error bound is stated for the
    # grid size it asked for, which is not larger.
    M_fft = nextpow(2, M_grid)
    samples = [f(p) for p in FourierPoints(M_fft, Float64)]
    raw = interval_fft(samples)
    k_freq = B.k
    return [raw[1:k_freq+1]; raw[M_fft-k_freq+1:M_fft]]
end

# ---------------------------------------------------------------------------
# ProjectedFunction for Fourier{W{k,l}, …}
# ---------------------------------------------------------------------------

@doc raw"""
    ProjectedFunction(B::FourierAnalytic{W{k,l},…}, f::Function;
                      weak_dual_bound, Wk1_seminorm, M = 4 * B.k)
    ProjectedFunction(B::FourierAdjoint{W{k,l},…}, f::Function;
                      weak_dual_bound, Wk1_seminorm, M = 4 * B.k)

Discretize `f` on a Fourier basis with strong norm ``W^{k,1}`` (``k ≥ 1``,
including BV / ``W^{1,1}``).

- `Wk1_seminorm` — upper bound on ``\|f^{(k)}\|_{L^1}`` (total variation
  for ``k = 1``). Drives the L² truncation/aliasing tail formulas.
- `weak_dual_bound` — upper bound on ``\|f\|_{w^*}`` (for the default
  weak `L²` Fourier basis, this is ``\|f\|_{L^2}``). Required because
  computing it from a callable `f` isn't generally cheap.
- `M::Integer` — FFT grid size. Default `4 * B.k` controls the
  aliasing tail (see module docstring).

`proj_error` of the result is the combined L² error
``\|f - φ_N\|_{L²} \leq \sqrt{T_1^2 + T_2^2}``.

(`Observable` is a `const` alias for `ProjectedFunction`, so
`Observable(B, ϕ; …)` calls the same constructor.)
"""
function RigorousInvariantMeasures.ProjectedFunction(
    B::FourierAnalytic{W{k,l}},
    f::Function;
    weak_dual_bound,
    Wk1_seminorm,
    M::Integer = 4 * B.k,
) where {k,l}
    coeffs = _fourier_coeffs_no_aliasing(B, f, M)
    proj_err = _fourier_L2_total_error_Wk1(B.k, k, Wk1_seminorm, M)
    return RigorousInvariantMeasures.ProjectedFunction(
        B,
        coeffs,
        weak_dual_bound,
        proj_err,
    )
end

function RigorousInvariantMeasures.ProjectedFunction(
    B::FourierAdjoint{W{k,l}},
    f::Function;
    weak_dual_bound,
    Wk1_seminorm,
    M::Integer = 4 * B.k,
) where {k,l}
    coeffs = _fourier_coeffs_no_aliasing(B, f, M)
    proj_err = _fourier_L2_total_error_Wk1(B.k, k, Wk1_seminorm, M)
    return RigorousInvariantMeasures.ProjectedFunction(
        B,
        coeffs,
        weak_dual_bound,
        proj_err,
    )
end

# `projection` dispatches to `ProjectedFunction` for any Fourier basis;
# the concrete subtype's constructor handles the math and the kwargs.
# Bases without a defined `ProjectedFunction` constructor (e.g. `Aη`/`Cω`
# strong norm — not yet implemented) raise `MethodError`.
RigorousInvariantMeasures.projection(B::Fourier, f::Function; kwargs...) =
    RigorousInvariantMeasures.ProjectedFunction(B, f; kwargs...)
