
export coarse_fine_weak_resolvent,
    strong_resolvent_lift,
    coarse_fine_weak_resolvent_auto_N,
    bN_constant,
    projector_distance_bound

# Float64 upper / lower bounds for x^n via repeated FastRounding multiplication.
function _pow_round_up(x::Real, n::Integer)
    n ≥ 0 || throw(ArgumentError("n must be ≥ 0; got n = $n"))
    r = 1.0
    for _ = 1:n
        r = r ⊗₊ x
    end
    return r
end

function _pow_round_down(x::Real, n::Integer)
    n ≥ 0 || throw(ArgumentError("n must be ≥ 0; got n = $n"))
    r = 1.0
    for _ = 1:n
        r = r ⊗₋ x
    end
    return r
end

@doc raw"""
    bN_constant(a, b, M, N)

Float64 upper bound on the Lasota–Yorke iteration constant

```math
b_N := b\sum_{j=0}^{N-1} a^{N-1-j} M^j,
```

from Lemma A.4 of [Nisoli, "Certified spectral approximation of transfer
operators and the Gauss map"]. With `a, b, M ≥ 0`, this is the same
``b_N`` appearing in

```math
\|L^N u\|_s \;\le\; a^N \|u\|_s + b_N \|u\|_w,
```

and inside the perturbation factor of Proposition A.14. `N ≥ 1`; for
`N = 1` the sum has the single term ``j = 0``, giving ``b_1 = b``.

# Caller note (package convention)

This helper is paper-faithful, with a single ``M``. In the package the
DFLY is on `(strong, aux)`, so the *correct* value of `M` here is the
**aux-norm contractivity** of ``L`` — typically `1` for L¹-preserving
transfer operators with `aux_norm = L¹`. The `b` passed in should be
``b^{\text{pkg}} \cdot M_2`` where ``M_2 = \mathrm{aux\_weak\_bound}(B)``.

In particular, do *not* pass `bound_weak_norm_abstract(B, D)` as `M`:
that quantity (which can exceed 1) leads to `b_N = O(M^N)` blowing up
and to Prop A.14 never closing. See [`coarse_fine_weak_resolvent`](@ref)
for the basis-aware wiring that handles this conversion automatically.
"""
function bN_constant(a::Real, b::Real, M::Real, N::Integer)
    N ≥ 1 || throw(ArgumentError("N must be ≥ 1; got N = $N"))
    s = 0.0
    for j = 0:N-1
        s = s ⊕₊ (_pow_round_up(a, N - 1 - j) ⊗₊ _pow_round_up(M, j))
    end
    return b ⊗₊ s
end

# S_N^{(k)}(z) upper bound (eq. (32)). Internal — callers reach this via the
# basis-aware entry points below.
function _SN_bound(
    abs_z::Real,
    N::Integer;
    norm_powers::Union{AbstractVector,Nothing} = nothing,
    M::Union{Real,Nothing} = nothing,
)
    N ≥ 1 || throw(ArgumentError("N must be ≥ 1; got N = $N"))
    abs_z > 0 || throw(ArgumentError("abs_z must be > 0"))
    (norm_powers === nothing) == (M === nothing) &&
        throw(ArgumentError("supply exactly one of `norm_powers` or `M`"))

    if norm_powers !== nothing
        length(norm_powers) ≥ N || throw(
            ArgumentError(
                "norm_powers needs ≥ N = $N entries; got $(length(norm_powers))",
            ),
        )
        s = 0.0
        for ℓ = 0:N-1
            abszℓ_lower = _pow_round_down(abs_z, ℓ)
            abszℓ_lower > 0 || return Inf
            s = s ⊕₊ (Float64(norm_powers[ℓ+1]) ⊘₊ abszℓ_lower)
        end
        return s ⊘₊ abs_z
    else
        Mf = Float64(M)
        s = 0.0
        rℓ = 1.0
        r = Mf ⊘₊ abs_z
        for _ = 0:N-1
            s = s ⊕₊ rℓ
            rℓ = rℓ ⊗₊ r
        end
        return s ⊘₊ abs_z
    end
end

@doc raw"""
    abstract_weak_norm_bound(B::Basis; dfly_coefficients)

Float64 upper bound on ``\|L_k u\|_w / \|u\|_w`` for ``u \in U_h``,
derived entirely from the basis interface via the vector inequality
``\|v\|_w \le S_1\|v\|_s + S_2 |||v|||`` (`weak_by_strong_and_aux_bound`)
combined with the DFLY inequality and the norm-equivalence constants on
``U_h``:

```math
\|L_k u\|_w
   \;\le\; S_1\bigl(a\,M_{1n} + b\,M_2\bigr)\|u\|_w + S_2\,M_2\,\|u\|_w,
\qquad u \in U_h.
```

Assumes the aux norm is contracted by ``L`` (``|||L u||| \le |||u|||``),
which holds for transfer operators of integral-preserving maps with the
typical aux = L¹ convention. For bases where the standard
`bound_weak_norm_abstract` gives a sharper estimate, prefer that or pass
a computed bound from `powernormbounds` directly.
"""
function abstract_weak_norm_bound(
    B::Basis;
    dfly_coefficients,
)
    a, b = dfly_coefficients
    S1, S2 = weak_by_strong_and_aux_bound(B)
    M1n = strong_weak_bound(B)
    M2 = aux_weak_bound(B)
    return (S1 ⊗₊ ((a ⊗₊ M1n) ⊕₊ (b ⊗₊ M2))) ⊕₊ (S2 ⊗₊ M2)
end

# Core Prop A.14 scalar formula — internal helper isolating the arithmetic
# so the tests can pin down each component without constructing a basis.
#
# Note on the two `M`s in Prop A.14:
# - `M_aux` enters `b_N` and is the aux-norm contractivity of L (= 1 for
#   L¹-preserving transfer operators with aux = L¹). The caller is
#   responsible for passing `b` already scaled by `M_2 = aux_weak_bound(B)`
#   so the paper-formula `b_N = b ∑ a^{N-1-j} M^j` matches the iterated
#   (s, w) DFLY in the package convention.
# - `M_k_weak` enters `S_N^{(k)}(z)` and is the *weak* operator norm of
#   the fine discretized operator L_k (e.g., the `upper_bound_L2_opnorm`
#   for L²-weak Fourier).
function _coarse_fine_weak_resolvent_scalar(
    abs_z::Real,
    K_z::Real,
    M_k_weak;
    a::Real,
    b::Real,
    M_aux::Real,
    Δ_k::Real,
    E_sw::Real,
    E_k_ws::Real,
    N::Integer,
)
    N ≥ 1 || throw(ArgumentError("N must be ≥ 1; got N = $N"))
    abs_z > 0 || throw(ArgumentError("abs_z must be > 0"))

    bN = bN_constant(a, b, M_aux, N)
    SN = if M_k_weak isa AbstractVector
        _SN_bound(abs_z, N; norm_powers = M_k_weak)
    else
        _SN_bound(abs_z, N; M = M_k_weak)
    end

    aN = _pow_round_up(a, N)
    factor = (aN ⊗₊ E_k_ws) ⊕₊ bN

    abszN_lower = _pow_round_down(abs_z, N)
    abszN_lower > 0 || return Inf
    inv_abszN = 1.0 ⊘₊ abszN_lower

    β̃ = inv_abszN ⊗₊ Δ_k ⊗₊ K_z ⊗₊ factor
    β̃ < 1.0 || return Inf

    num_upper = SN ⊕₊ (inv_abszN ⊗₊ E_sw ⊗₊ K_z ⊗₊ factor)
    denom_lower = 1.0 ⊖₋ β̃
    denom_lower > 0 || return Inf

    return num_upper ⊘₊ denom_lower
end

@doc raw"""
    strong_resolvent_lift(B::Basis, D::Dynamic, abs_z, R_w_coarse; kwargs...)

Float64 upper bound on the *strong* resolvent of the infinite-dimensional
operator,

```math
\mathcal{R}_s(z, L) \;\le\; \frac{b\,\mathcal{R}_w(z, L_{k_0})\,E_{s\to w} + 1}
        {|z| - a - b\,\mathcal{R}_w(z, L_{k_0})\,\delta_{k_0}},
```

derived from Proposition A.7 of the paper. The numerator bounds
``\|u\|_s`` for the unique solution of ``(zI - L)u = f`` with
``\|f\|_s = 1``; the denominator is the perturbation gap

```math
|z| > a + b\,\mathcal{R}_w(z, L_{k_0})\,\delta_{k_0}.
```

Use this to lift a *coarse-level* weak resolvent certificate
``R_{w,\text{coarse}} \ge \mathcal{R}_w(z, L_{k_0})`` to a strong
resolvent bound for the infinite-dimensional ``L``. Returns `Inf` when
the perturbation gap fails to close.

Standing constants are pulled from the basis interface; override via
keyword for sharper specific-case bounds:

| paper            | default                                       |
|------------------|-----------------------------------------------|
| ``a, b``         | `dfly(strong_norm(B), aux_norm(B), D)`        |
| ``\delta_{k_0}`` | `weak_projection_error(B)`                    |
| ``E_{s\to w}``   | `1.0`                                         |
"""
function strong_resolvent_lift(
    B::Basis,
    D::Dynamic,
    abs_z::Real,
    R_w_coarse::Real;
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
    δ::Real = weak_projection_error(B),
    E_sw::Real = 1.0,
)
    abs_z > 0 || throw(ArgumentError("abs_z must be > 0"))
    a, b = dfly_coefficients

    bRδ = b ⊗₊ R_w_coarse ⊗₊ δ
    denom_lower = (abs_z ⊖₋ a) ⊖₋ bRδ
    denom_lower > 0 || return Inf

    num_upper = (b ⊗₊ R_w_coarse ⊗₊ E_sw) ⊕₊ 1.0
    return num_upper ⊘₊ denom_lower
end

@doc raw"""
    coarse_fine_weak_resolvent(B::Basis, D::Dynamic, abs_z, K_z;
                               N, M_k_weak=..., kwargs...)

Float64 upper bound on the *fine* weak resolvent
``\mathcal{R}_w(z, L_k)`` of the finite-rank discretization, from
Proposition A.14 (DFLY appendix) of [Nisoli, "Certified spectral
approximation of transfer operators and the Gauss map"]:

```math
\mathcal{R}_w(z, L_k) \;\le\;
   \frac{S_N^{(k)}(z) + |z|^{-N}\,E_{s\to w}\,K(z)\,\bigl(a^N E_{k,w\to s} + b_N\bigr)}
        {1 - \widetilde{\beta}_N(z)},
```

```math
\widetilde{\beta}_N(z) \;:=\; |z|^{-N}\,\Delta_k\,K(z)\,
        \bigl(a^N E_{k,w\to s} + b_N\bigr),
\qquad
S_N^{(k)}(z) := \frac{1}{|z|}\sum_{\ell=0}^{N-1}\frac{\|L_k^{\ell}\|_w}{|z|^{\ell}}.
```

# Caller-supplied inputs
- `abs_z` — magnitude ``|z|``; must satisfy ``|z| > a`` for a useful bound.
- `K_z`   — upper bound on the *strong* resolvent
            ``\mathcal{R}_s(z, L)`` of the infinite-dimensional ``L`` at
            this ``z`` (typically from [`strong_resolvent_lift`](@ref) at
            a coarse level).
- `N`     — Laurent-truncation parameter, ≥ 1.

# Basis-derived defaults (override via keyword)
| paper / role       | default                                            |
|--------------------|----------------------------------------------------|
| ``a, b``           | `dfly(strong_norm(B), aux_norm(B), D)`             |
| ``M_{\rm aux}``    | `1.0` (aux-norm contractivity of L; correct for    |
|                    | any L¹-preserving transfer operator with `aux = L¹`).|
| ``\Delta_k``       | `weak_projection_error(B)`                         |
| ``E_{s\to w}``     | `1.0`                                              |
| ``E_{k,w\to s}``   | `strong_weak_bound(B)` (the ``M_{1n}`` constant)   |
| ``M_k`` (for ``S_N^{(k)}``) | [`abstract_weak_norm_bound(B; …)`](@ref) |
|                    | (from `weak_by_strong_and_aux_bound` + DFLY). Pass |
|                    | a computed scalar or an `AbstractVector` of upper  |
|                    | bounds on ``[\|L_k^\ell\|_w]_{\ell=0}^{N-1}`` (e.g.|
|                    | `max.(1.0, powernormbounds(B, D))`) for a sharper sum. |

Internally `b` is multiplied by `M_2 = aux_weak_bound(B)` before being
fed into `b_N`, so the paper formula

```math
b_N = b \sum_{j=0}^{N-1} a^{N-1-j} M^j
```

is applied with `b ← b · M_2` and `M ← M_aux`. This is the iterated
(strong, weak) DFLY constant derived from the package's (strong, aux)
DFLY plus the aux-to-weak conversion — and stays bounded by
`b · M_2 / (1 - a)` when `M_aux = 1`.

Returns `Inf` when the perturbation condition
``\widetilde{\beta}_N(z) < 1`` fails (typically because ``N`` is too
small or the discretization is too coarse).

# Example
```julia
B = Ulam(1024)
D = mod1_dynamic(x -> 2 * x)
norms = powernormbounds(B, D)              # ‖Q^ℓ|_U‖_w sequence
M_k_seq = max.(1.0, norms)                 # whole-space upper bounds

# K_z certified at a coarse level via strong_resolvent_lift / a Schur
# resolvent bound on a contour.
R_w_fine = coarse_fine_weak_resolvent(B, D, abs_z, K_z;
                                       N = 8, M_k_weak = M_k_seq)
```
"""
function coarse_fine_weak_resolvent(
    B::Basis,
    D::Dynamic,
    abs_z::Real,
    K_z::Real;
    N::Integer,
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
    Δ_k::Real = weak_projection_error(B),
    E_sw::Real = 1.0,
    E_k_ws::Real = strong_weak_bound(B),
    M_aux::Real = 1.0,
    M_2::Real = aux_weak_bound(B),
    M_k_weak = abstract_weak_norm_bound(B; dfly_coefficients = dfly_coefficients),
)
    a, b = dfly_coefficients
    b_eff = b ⊗₊ M_2
    return _coarse_fine_weak_resolvent_scalar(
        abs_z,
        K_z,
        M_k_weak;
        a = a,
        b = b_eff,
        M_aux = M_aux,
        Δ_k = Δ_k,
        E_sw = E_sw,
        E_k_ws = E_k_ws,
        N = N,
    )
end

@doc raw"""
    coarse_fine_weak_resolvent_auto_N(B::Basis, D::Dynamic, abs_z, K_z;
                                       μ=nothing, kwargs...)

Convenience wrapper applying Corollary A.16: choose

```math
N_k := \left\lceil\frac{\log E_{k,w\to s}}{|\log(a/\mu)|}\right\rceil,
```

where ``\mu \in (a, M)`` is a target spectral radius (defaulting to the
geometric mean ``\sqrt{aM}`` when not supplied). This is the truncation
depth that balances ``a^{N_k}`` against ``(\mu/M)^{N_k}`` in the bound.

Returns `(R_w_fine, N_k)`. Remaining keywords are forwarded to
[`coarse_fine_weak_resolvent`](@ref) and share its basis-interface
defaults.
"""
function coarse_fine_weak_resolvent_auto_N(
    B::Basis,
    D::Dynamic,
    abs_z::Real,
    K_z::Real;
    μ::Union{Real,Nothing} = nothing,
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
    Δ_k::Real = weak_projection_error(B),
    E_sw::Real = 1.0,
    E_k_ws::Real = strong_weak_bound(B),
    M_aux::Real = 1.0,
    M_2::Real = aux_weak_bound(B),
    M_k_weak = abstract_weak_norm_bound(B; dfly_coefficients = dfly_coefficients),
)
    a, _ = dfly_coefficients
    # Cor. A.16's μ-range condition is `μ ∈ (a, M_aux)` in the iterated
    # (s, aux) DFLY — since M_aux = 1 by default and a < 1 by the LY
    # contraction, this is `μ ∈ (a, 1)`.
    0 < a < M_aux ||
        throw(ArgumentError("Cor. A.16 requires 0 < a < M_aux; got a = $a, M_aux = $M_aux"))
    μ_val = μ === nothing ? sqrt(a * M_aux) : Float64(μ)
    a < μ_val < M_aux ||
        throw(ArgumentError("μ must lie in (a, M_aux); got μ = $μ_val"))
    abs_z ≥ μ_val ||
        @warn("Cor. A.16 assumes |z| ≥ μ; got |z| = $abs_z < μ = $μ_val")

    denom = abs(log(a / μ_val))
    denom > 0 || throw(ArgumentError("log(a/μ) is zero — pick a ≠ μ"))
    N_k = max(1, Int(ceil(log(max(E_k_ws, 1.0)) / denom)))

    R = coarse_fine_weak_resolvent(
        B,
        D,
        abs_z,
        K_z;
        N = N_k,
        dfly_coefficients = dfly_coefficients,
        Δ_k = Δ_k,
        E_sw = E_sw,
        E_k_ws = E_k_ws,
        M_aux = M_aux,
        M_2 = M_2,
        M_k_weak = M_k_weak,
    )
    return (R, N_k)
end

@doc raw"""
    projector_distance_bound(B_coarse, B_fine, D, Q_coarse, Q_fine, radius;
                              eigenvalue = 1.0, samples = 256, N, kwargs...)

Float64 upper bound on the *Riesz projector distance*
``\|P_L - P_{L_{\rm fine}}\|_{s \to w}`` for the spectral projector
onto the eigenspace of ``L`` and ``L_{\rm fine}`` enclosed by the small
circle

```math
\Gamma = \{ z \in \mathbb{C} : |z - \mathrm{eigenvalue}| = r \}.
```

Combines Lemma A.11 (strong-to-weak resolvent perturbation) with the
contour-integral representation of the Riesz projector:

```math
\| P_L - P_{L_{\rm fine}} \|_{s \to w}
   \;\le\; \frac{\ell(\Gamma)}{2\pi}
   \sup_{z \in \Gamma} \mathcal{R}_w(z, L_{\rm fine}) \,\delta_{\rm fine}\, \mathcal{R}_s(z, L)
   \;=\; r \cdot
   \sup_{z \in \Gamma} \mathcal{R}_w(z, L_{\rm fine}) \,\delta_{\rm fine}\, K(z).
```

# How the pieces are obtained
The strategy is the coarse-fine one: validated certification is done on
the **coarse** matrix only (small enough for cheap Schur), and both
fine-level quantities are then propagated.

- `CertifScripts.run_certification` on the coarse `Q_coarse.L` with a
  `CertifScripts.CertificationCircle` of radius `r` centred
  at `eigenvalue` gives ``\sup_{z \in \Gamma} \mathcal{R}_w(z, L_{\rm coarse})``.
- ``K(z) := \mathcal{R}_s(z, L)`` is obtained from the coarse weak
  resolvent via [`strong_resolvent_lift`](@ref) (Prop A.7), evaluated
  at the worst contour magnitude ``|z|_{\min} = |\text{eigenvalue}| - r``.
- ``\sup_{z \in \Gamma} \mathcal{R}_w(z, L_{\rm fine})`` is propagated
  to the fine level via [`coarse_fine_weak_resolvent`](@ref) (Prop A.14),
  also at the worst contour magnitude.
- ``\delta_{\rm fine}`` is `weak_projection_error(B_fine)`.

This way `CertifScripts` is run only on the coarse matrix; the fine
level is reached entirely by the propagation pipeline.

# Inputs
- `B_coarse, B_fine` — coarse and fine `Basis` (same family).
- `D::Dynamic`        — the dynamic.
- `Q_coarse, Q_fine`  — `DiscretizedOperator(B, D)` for each level.
- `radius`            — circle radius `r`; must satisfy
                        `0 < r < |eigenvalue|` and in practice
                        `r < 1 − λ₂` (else Γ encloses extra spectrum).
- `eigenvalue::Number` — circle centre (default `1.0`).
- `samples::Integer`   — number of contour samples for `CertifScripts`
                         (default 256).
- `N::Integer`         — Laurent-truncation depth for Prop A.14.
- `M_k_weak`           — bound on ``\|L_{\rm fine}\|_w`` (default
                         `upper_bound_L2_opnorm(BallMatrix(Q_fine.L))`).
- Remaining kwargs (`M_aux`, `M_2`, `Δ_k`, `E_sw`, `E_k_ws`,
  `dfly_coefficients`, `schur_data`) are forwarded — see
  [`coarse_fine_weak_resolvent`](@ref).

# Connection to the eigenvector distance
For an L¹-preserving transfer operator with simple eigenvalue 1, both
Riesz projectors are rank-1 of the form ``P f = v\,\langle\ell, f\rangle``
where ``\ell`` is the integral functional. Hence

```math
\|P_L - P_{L_{\rm fine}}\|_{s \to w}
   \;=\; \|v - v_{\rm fine}\|_w \cdot \|\ell\|_{s \to \mathbb{R}},
```

so the projector-distance bound translates directly into a bound on the
invariant-measure distance up to the factor ``\|\ell\|_{s\to\mathbb R}``
(typically ≤ 1 for the integral functional on regular spaces).

Returns a named tuple with `projector_distance`, `R_w_coarse`, `K_z`,
`R_w_fine`, `δ_k`, `radius`, and the underlying `cert`/`schur_data`.
If Prop A.7 or A.14 fails to close, `projector_distance` is `Inf`.
"""
function projector_distance_bound(
    B_coarse::Basis,
    B_fine::Basis,
    D::Dynamic,
    Q_coarse::DiscretizedOperator,
    Q_fine::DiscretizedOperator,
    radius::Real;
    eigenvalue::Number = 1.0,
    samples::Integer = 256,
    N::Integer,
    schur_data = nothing,
    M_k_weak::Union{Real,AbstractVector,Nothing} = nothing,
    M_aux::Real = 1.0,
    M_2::Real = aux_weak_bound(B_fine),
    Δ_k::Real = weak_projection_error(B_fine),
    E_sw::Real = 1.0,
    E_k_ws::Real = strong_weak_bound(B_fine),
    dfly_coefficients = dfly(strong_norm(B_fine), aux_norm(B_fine), D),
)
    radius > 0 ||
        throw(ArgumentError("radius must be > 0; got radius = $radius"))
    abs_eig = abs(eigenvalue)
    abs_z_min = abs_eig - radius
    abs_z_min > 0 ||
        throw(ArgumentError(
            "radius must be < |eigenvalue| = $abs_eig; got radius = $radius",
        ))

    BM_coarse = BallMatrix(Q_coarse.L)
    sd = schur_data === nothing ?
         CertifScripts.compute_schur_and_error(BM_coarse) : schur_data
    circle = CertifScripts.CertificationCircle(
        Complex{Float64}(eigenvalue), radius; samples = samples,
    )
    cert = CertifScripts.run_certification(BM_coarse, circle; schur_data = sd)
    R_w_coarse = cert.resolvent_original

    K_z = strong_resolvent_lift(
        B_coarse, D, abs_z_min, R_w_coarse;
        dfly_coefficients = dfly_coefficients,
        δ = weak_projection_error(B_coarse),
        E_sw = E_sw,
    )
    isfinite(K_z) || return (
        projector_distance = Inf, R_w_coarse = R_w_coarse, K_z = K_z,
        R_w_fine = Inf, δ_k = Δ_k, radius = radius,
        cert = cert, schur_data = sd,
    )

    M_k_used = M_k_weak === nothing ?
               upper_bound_L2_opnorm(BallMatrix(Q_fine.L)) : M_k_weak
    R_w_fine = coarse_fine_weak_resolvent(
        B_fine, D, abs_z_min, K_z;
        N = N,
        dfly_coefficients = dfly_coefficients,
        Δ_k = Δ_k,
        E_sw = E_sw,
        E_k_ws = E_k_ws,
        M_aux = M_aux,
        M_2 = M_2,
        M_k_weak = M_k_used,
    )
    isfinite(R_w_fine) || return (
        projector_distance = Inf, R_w_coarse = R_w_coarse, K_z = K_z,
        R_w_fine = Inf, δ_k = Δ_k, radius = radius,
        cert = cert, schur_data = sd,
    )

    proj_dist = radius ⊗₊ R_w_fine ⊗₊ Δ_k ⊗₊ K_z

    return (
        projector_distance = proj_dist,
        R_w_coarse = R_w_coarse,
        K_z = K_z,
        R_w_fine = R_w_fine,
        δ_k = Δ_k,
        radius = radius,
        cert = cert,
        schur_data = sd,
    )
end
