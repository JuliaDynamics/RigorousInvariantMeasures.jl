"""
Specialized a posteriori noise error estimator: the bootstrap of Galatolo, Monge, Nisoli,
Nonlinearity 33 (2020) 4237, §3.3 (Lemmas 24, 25, 27 and §3.3.5), with the bounds of their
Proposition 23 and Lemma 56.

Ported from the Python `noise_specialized_measure_error.py` in the
`compinvmeas-python-release-noise` package.  These routines produce
significantly tighter error bounds than the generic DFLY-based approach
in `distance_from_invariant_noise` by exploiting the *actual computed
density vector* to obtain data-dependent bounds.

The lemmas are written for the operator `π_δ N_ξ π_δ L π_δ`, so the kernel is the
[`UlamProjectedUniformKernel`](@ref); `ξ` is read from it. Two misprints of the paper are
corrected here: in §3.3.5 the constant and the coefficient are swapped (with the bounds
`A‖f-f̃‖ + B`, `‖f_ξ - f_{ξ,δ}‖ ≤ B/(1-α) + A/(1-α)‖f_ξ - f̃‖`), and in (30) the first
argument of the minimum is also multiplied by `Var_I(N_ξ L f̃)`, as in the proof.
"""

using IntervalArithmetic

# ──────────────────────────────────────────────────────────────────
# 1.  Total variation of a step-function vector
# ──────────────────────────────────────────────────────────────────

"""
    total_variation(v::AbstractVector)

Compute a rigorous upper bound for the total variation of a
step-function represented by the vector `v`:

    Var(v) = Σᵢ |v[i] − v[i+1]|

Uses directed rounding (round-up) to keep the result rigorous.
"""
function total_variation(v::AbstractVector)
    var = 0.0
    @inbounds for i in 1:length(v)-1
        var = var ⊕₊ Float64(abs(v[i] - v[i+1]), RoundUp)
    end
    return var
end

# ──────────────────────────────────────────────────────────────────
# 2.  Preimage data for the W-norm estimate  (preparation_for_estimate1)
# ──────────────────────────────────────────────────────────────────

"""
    _mag_or_inf(g, x)

An upper bound on `|g(x)|` for an interval `x`, or `Inf` where `g` cannot be bounded on
`x`: a `DomainError` or an empty result, as for the derivative of `(1/8 - x)^(1/3)` at
`x = 1/8`. An infinite bound makes the estimates below fall back to the mass, which is
valid whatever the true value.
"""
function _mag_or_inf(g, x)
    v = try
        g(x)
    catch e
        e isa DomainError || rethrow()
        return Inf
    end
    isempty_interval(v) && return Inf
    return Float64(mag(abs(v)), RoundUp)
end

"""
Information about a single preimage of an interval I under a branch Tₖ.
"""
struct PreimageInfo
    A::Float64            # left endpoint of preimage (mid, rounded down)
    B::Float64            # right endpoint of preimage (mid, rounded up)
    invtp::Float64        # upper bound for |1/T'| on preimage
    distortion::Float64   # upper bound for |T''/(T')²| on preimage
    boundary_terms::Vector{Tuple{Float64,Float64}}  # (endpoint y of the branch domain with T(y) ∈ I, upper bound on |1/T'(y)|)
end

"""
    prepare_Wnorm_estimate(D::PwMap, B_basis::Ulam, B_est::Ulam)

For every interval Iᵢ in the estimation partition `B_est`, compute
preimage data under each branch of `D`.

Returns `Vector{Vector{PreimageInfo}}` of length `length(B_est)`:
  - `prep[i]` contains one `PreimageInfo` per branch of `D` that has
    a preimage overlapping with Iᵢ.

This is the Julia equivalent of the Python `preparation_for_estimate1`.
"""
function prepare_Wnorm_estimate(D::PwMap, B_est::Ulam)
    Kest = length(B_est)
    prep = Vector{Vector{PreimageInfo}}(undef, Kest)

    # Precompute branch limits  T(endpoint) for each branch
    branch_limits = Vector{Vector{Interval{Float64}}}(undef, nbranches(D))
    for k in 1:nbranches(D)
        br = D.branches[k]
        branch_limits[k] = [br.f(br.X[1]), br.f(br.X[2])]
    end

    # Build the Ulam dual (preimages of the estimation partition under D)
    dual = Dual(B_est, D; ϵ=1e-14, max_iter=100)

    # We accumulate preimage infos per interval
    for i in 1:Kest
        prep[i] = PreimageInfo[]
    end

    # The UlamDual iterator yields (ylabel, (a, b)) pairs
    # ylabel is the interval index in B_est, (a,b) are preimage endpoints
    for (ylabel, (a, b)) in dual
        # Determine which branch this preimage comes from
        # The UlamDual concatenates preimages across branches.
        # Each branch produces exactly Kest preimage pairs (since full-branch maps),
        # but non-full-branch maps produce fewer.
        # We track branches by counting: the preimage points are ordered by branch.

        # Figure out which branch produced this preimage
        # We do this by checking which branch domain contains the preimage
        A_val = min(a, b)
        B_val = max(a, b)
        TiI = hull(A_val, B_val)

        bnum = 0
        for k in 1:nbranches(D)
            br_domain = hull(D.branches[k].X[1], D.branches[k].X[2])
            if !isempty_interval(intersect_interval(TiI, br_domain))
                bnum = k
                break
            end
        end
        if bnum == 0
            continue  # skip empty preimages
        end

        br = D.branches[bnum]

        # Compute |1/T'| upper bound on preimage
        invtp_val = _mag_or_inf(inverse_derivative(br.f), TiI)

        # Compute |T''/(T')²| (distortion) upper bound on preimage
        distortion_val = _mag_or_inf(x -> distortion(br.f, x), TiI)

        # Boundary terms of Lemma 56: the endpoints y of the branch domain with T(y) ∈ I,
        # each contributing |g(y)|/|T'(y)|, so 1/T' is taken at y and not at T(y)
        I_interval = hull(interval(B_est.p[ylabel]), interval(B_est.p[ylabel + 1]))
        bterms = Tuple{Float64,Float64}[]
        for (y, Ty) in zip(br.X, branch_limits[bnum])
            if !isempty_interval(intersect_interval(Ty, I_interval))
                y_invtp = _mag_or_inf(inverse_derivative(br.f), y)
                push!(bterms, (mid(y), y_invtp))
            end
        end

        pi = PreimageInfo(
            Float64(inf(interval(A_val)), RoundDown),
            Float64(sup(interval(B_val)), RoundUp),
            invtp_val,
            distortion_val,
            bterms,
        )
        push!(prep[ylabel], pi)
    end

    return prep
end

# ──────────────────────────────────────────────────────────────────
# 3.  Estimate ||(1-π)Lf||_W  (estimate_W_norm_of_1mpi_Lf)
# ──────────────────────────────────────────────────────────────────

"""
    estimate_Wnorm_of_discretization_error(prep, density_vec, B_density::Ulam)

Estimate `||(1-π)Lf||_W` where `π` is the Ulam projection and `f` is
the density represented by `density_vec` on the basis `B_density`.

For each interval `I` and each preimage, computes two alternative bounds
and takes the minimum:
  - Variation-based:  `(δ/4) * estimate_var`
  - Measure-based:    `meas_in_ab`

Returns `(est, meas_Lf, var_Lf)`:
  - `est`: the W-norm estimate (scalar, rounded up)
  - `meas_Lf`: vector of measure of Lf per estimation interval
  - `var_Lf`: vector of variation of Lf per estimation interval
"""
function estimate_Wnorm_of_discretization_error(
    prep::Vector{Vector{PreimageInfo}},
    density_vec::AbstractVector,
    B_density::Ulam,
)
    Kest = length(prep)
    K_density = length(density_vec)
    δ = Float64(1.0, RoundUp) / Float64(K_density, RoundDown)

    estimate = 0.0
    meas_Lf = zeros(Kest)
    var_Lf = zeros(Kest)

    partition = B_density.p

    for i in 1:Kest
        meas_in_Ii = 0.0
        var_in_Ii = 0.0

        for pinfo in prep[i]
            A, B = pinfo.A, pinfo.B
            invtp_ab = pinfo.invtp
            distortion_ab = pinfo.distortion

            meas_in_ab = 0.0
            var_in_ab = 0.0

            # Binary search: find which density partition intervals overlap [A, B]
            jmin = searchsortedlast(partition, A)
            jmax = searchsortedlast(partition, B)
            jmin = clamp(jmin, 1, K_density)
            jmax = clamp(jmax, 1, K_density)

            @inbounds for j in jmin:jmax
                x = Float64(partition[j])
                y = Float64(partition[j + 1])

                # Intersection of [A,B] with [x,y]
                lower = Float64(max(A, x), RoundUp)
                upper = Float64(min(B, y), RoundUp)
                prop_num = Float64(max(upper - lower, 0.0), RoundUp)
                prop_den = Float64(y - x, RoundDown)
                proportion = prop_num / prop_den  # rounded up since num↑ den↓

                # Lf integral over this piece
                meas_in_ab = meas_in_ab ⊕₊ (proportion ⊗₊ Float64(abs(density_vec[j]), RoundUp) ⊗₊ δ)

                # Variation contribution
                if j < jmax
                    var_in_ab = var_in_ab ⊕₊ Float64(abs(Float64(density_vec[j], RoundUp) - Float64(density_vec[j + 1], RoundDown)), RoundUp)
                end
            end

            # Estimate of the variation of L_i f in I, Lemma 56; where |1/T'| or the
            # distortion is unbounded the lemma gives nothing, the estimate is infinite and
            # the minimum below falls back to the mass of L_i f in I
            estimate_var_Ii_from_ab = if isfinite(invtp_ab) && isfinite(distortion_ab)
                invtp_ab ⊗₊ var_in_ab ⊕₊ distortion_ab ⊗₊ meas_in_ab
            else
                Inf
            end

            # Boundary terms |f(y)|/|T'(y)|; f is a step function, so |f(y)| is bounded by
            # the cells on both sides of y when y is a partition point
            for (y, y_invtp) in pinfo.boundary_terms
                j = clamp(searchsortedlast(partition, y), 1, K_density)
                fy = abs(density_vec[j])
                if j > 1 && partition[j] == y
                    fy = max(fy, abs(density_vec[j-1]))
                end
                if fy != 0             # 0·∞ would be NaN; a zero value contributes nothing
                    estimate_var_Ii_from_ab = estimate_var_Ii_from_ab ⊕₊ Float64(fy, RoundUp) ⊗₊ y_invtp
                end
            end

            meas_in_Ii = meas_in_Ii ⊕₊ meas_in_ab
            var_in_Ii = var_in_Ii ⊕₊ estimate_var_Ii_from_ab

            # Key: take minimum of two alternative bounds
            alt_est = (δ / 4) ⊗₊ estimate_var_Ii_from_ab
            estimate = estimate ⊕₊ min(alt_est, meas_in_ab)
        end

        meas_Lf[i] = meas_in_Ii
        var_Lf[i] = var_in_Ii
    end

    # Multiply the final estimate by δ/2
    return (estimate ⊗₊ δ / 2), meas_Lf, var_Lf
end

# ──────────────────────────────────────────────────────────────────
# 4.  Reflecting index helper
# ──────────────────────────────────────────────────────────────────

# The reflection of cell indices is `reflect_outward_idx` (UniformNoiseUlam.jl), the mirror
# at 0 and 1 of the reflecting boundary condition. The `_reflect` that was here mirrored
# about the centre of the first cell on the left (cell 0 to cell 2) and about the edge on
# the right (cell k+1 to cell k).

"The noise support `ξ` of a kernel, as an exact rational."
noise_support(K::UlamProjectedUniformKernel) = Rational{BigInt}(K.s) // length(K.B)

# ──────────────────────────────────────────────────────────────────
# 5.  Estimate Var(NLf)  (estimate_var_of_NLf)
# ──────────────────────────────────────────────────────────────────

"""
    estimate_variation_after_noise(meas_Lf, var_Lf, K::UlamProjectedUniformKernel)

Upper bounds on `Var_I(N_ξ g)` for the cells `I` of a partition of `Kest` cells of width
`h`, from the masses `meas_Lf` and variations `var_Lf` of `g` on those cells. Since
`(N_ξ g)' = (g(· + ξ/2) - g(· - ξ/2))/ξ`,

    Var_I(N_ξ g) ≤ (‖g‖_{L¹(I+ξ/2)} + ‖g‖_{L¹(I-ξ/2)}) / ξ,

and since `N_ξ g` is an average of translates of `g`,

    Var_I(N_ξ g) ≤ (1/ξ) ∫_{|t|≤ξ/2} Var_{I+t}(g) dt ≤ (h/ξ) Var_{I ± ξ/2}(g);

the minimum of the two is returned, cell by cell, with the reflecting boundary.
"""
function estimate_variation_after_noise(
    meas_Lf::AbstractVector,
    var_Lf::AbstractVector,
    K::UlamProjectedUniformKernel,
)
    Kest = length(meas_Lf)
    ξKest = noise_support(K) * Kest              # ξ/h, exact
    ξ_dn = Float64(noise_support(K), RoundDown)
    ξKest_dn = Float64(ξKest, RoundDown)
    f = Int(floor(ξKest / 2))
    c = Int(ceil(ξKest / 2))

    var_NLf = zeros(Kest)
    for i in 1:Kest
        # I ± ξ/2 meet the cells i ± f and i ± c
        est1 = (
            meas_Lf[reflect_outward_idx(i - c, Kest)] ⊕₊
            meas_Lf[reflect_outward_idx(i + c, Kest)] ⊕₊
            meas_Lf[reflect_outward_idx(i - f, Kest)] ⊕₊
            meas_Lf[reflect_outward_idx(i + f, Kest)]
        ) ⊘₊ ξ_dn

        # I + [-ξ/2, ξ/2] is covered by the cells i-c … i+c
        est2 = 0.0
        for j in (i - c):(i + c)
            est2 = est2 ⊕₊ var_Lf[reflect_outward_idx(j, Kest)]
        end
        est2 = est2 ⊘₊ ξKest_dn

        var_NLf[i] = min(est1, est2)
    end

    return var_NLf
end

# ──────────────────────────────────────────────────────────────────
# 6.  Estimate measure of NLf  (estimate_meas_of_NLf)
# ──────────────────────────────────────────────────────────────────

"""
    estimate_measure_after_noise(meas_Lf, K::UlamProjectedUniformKernel)

Upper bounds on `‖N_ξ g‖_{L¹(I)}` for the cells `I` of width `h` of a partition of `Kest`
cells, from the masses `meas_Lf` of `g`: a point `p` contributes to `∫_I N_ξ g` for the
`t` with `p - t ∈ I`, a set of measure at most `h`, so

    ‖N_ξ g‖_{L¹(I)} ≤ (h/ξ) ‖g‖_{L¹(I ± ξ/2)},

and `I ± ξ/2` is covered by the cells `i - c … i + c`, `c = ⌈ξ/(2h)⌉`.
"""
function estimate_measure_after_noise(
    meas_Lf::AbstractVector,
    K::UlamProjectedUniformKernel,
)
    Kest = length(meas_Lf)
    ξKest = noise_support(K) * Kest
    ξKest_dn = Float64(ξKest, RoundDown)
    c = Int(ceil(ξKest / 2))

    meas_NLf = zeros(Kest)
    for i in 1:Kest
        s = 0.0
        for j in (i - c):(i + c)
            s = s ⊕₊ meas_Lf[reflect_outward_idx(j, Kest)]
        end
        meas_NLf[i] = s ⊘₊ ξKest_dn
    end
    return meas_NLf
end

# ──────────────────────────────────────────────────────────────────
# 7.  Derivative bounds  (preparation_for_estimate2)
# ──────────────────────────────────────────────────────────────────

"""
    prepare_derivative_bounds(D::PwMap, B::Ulam)

For each interval in the partition of `B`, compute an upper bound
for `|T'|`.  This is used to estimate the `W → W` norm.

Returns a `Vector{Float64}` of length `length(B)`.
"""
function prepare_derivative_bounds(D::PwMap, B::Ulam)
    K = length(B)
    derivs = zeros(K)
    for i in 1:K
        I = hull(interval(B.p[i]), interval(B.p[i + 1]))
        # Evaluate |T'| on I by trying each branch that overlaps
        deriv_bound = 0.0
        for br in D.branches
            br_domain = hull(br.X[1], br.X[2])
            J = intersect_interval(I, br_domain)
            if !isempty_interval(J)
                d = _mag_or_inf(x -> derivative(br.f, x), J)
                deriv_bound = max(deriv_bound, d)
            end
        end
        derivs[i] = deriv_bound
    end
    return derivs
end

# ──────────────────────────────────────────────────────────────────
# 8.  L1 norm of NL(1-π)NLf  (estimate_L1_norm_of_NL_1mpi_NLf)
# ──────────────────────────────────────────────────────────────────

"""
    estimate_L1_NL_discretization_error(var_NLf, derivs, δ, varRho)

Estimate `||NL(1-π)NLf||₁` — the L1 norm of the second-order
discretization error.

For each interval: `min(δ*varRho/4*deriv[i], 1) * var_NLf[i]`,
summed and scaled by `δ/2`.
"""
function estimate_L1_NL_discretization_error(
    var_NLf::AbstractVector,
    derivs::AbstractVector,
    δ::Float64,
    varRho::Float64,
)
    retv = 0.0
    K = length(var_NLf)
    for i in 1:K
        factor = min(δ ⊗₊ varRho / 4 ⊗₊ derivs[i], 1.0)
        retv = retv ⊕₊ factor ⊗₊ var_NLf[i]
    end
    return retv ⊗₊ δ / 2
end

# ──────────────────────────────────────────────────────────────────
# 9.  A priori error bound  (estimate_L1_error_apriori)
# ──────────────────────────────────────────────────────────────────

"""
    noise_error_apriori(K_size, noise_size_rel, α, sumCi)

Simple a priori error bound for the noisy system.

    error ≤ (1/(1-α)) * (2*sumCi + 1) * (δ/2) * varRho

where `δ = 1/K_size` and `varRho = 2/noise_size_rel`.
"""
function noise_error_apriori(K_size::Integer, noise_size_rel::Float64, α::Float64, sumCi::Float64)
    δ = Float64(1.0, RoundUp) / Float64(K_size, RoundDown)
    varRho = 2.0 ⊘₊ Float64(noise_size_rel, RoundDown)
    return (1.0 ⊘₊ (1.0 ⊖₋ α)) ⊗₊ (2.0 ⊗₊ sumCi ⊕₊ 1.0) ⊗₊ (δ / 2) ⊗₊ varRho
end

# ──────────────────────────────────────────────────────────────────
# 10.  Main a posteriori error estimate  (estimate_L1_error_aposteriori)
# ──────────────────────────────────────────────────────────────────

"""
    noise_error_aposteriori(prep1, prep2, meas_vec, numeric_error,
                            B_density::Ulam, K::UlamProjectedUniformKernel,
                            α, sumCi)

Compute the a posteriori L1 error bound for a noisy dynamical system
using the actual computed density `meas_vec`.

# Arguments
- `prep1`: output of `prepare_Wnorm_estimate`
- `prep2`: output of `prepare_derivative_bounds`
- `meas_vec`: the computed density vector (approximation of the invariant density)
- `numeric_error`: the numerical error from the eigenvector computation
  (e.g., `residualboundnoise`)
- `B_density`: the `Ulam` basis on which `meas_vec` lives
- `K`: the noise kernel, `π_δ N_ξ π_δ` on the partition of `meas_vec`
- `α`, `sumCi`: for one `n`, `α ≥ ‖L_δⁿ|_V‖ < 1` and `sumCi ≥ Σ_{i=0}^{n-1} ‖L_δⁱ|_V‖`
  (with `‖L_δ⁰|_V‖ = 1`), as Lemma 15 and (14) require; a per-step rate is not an `α`

The estimation partition of `prep1` must have a number of cells dividing that of
`B_density`, its cells being unions of cells of the density partition (Lemmas 25, 27).

# Returns
`(L1Error, apriori_error, Linf_estimate)`:
- `L1Error`: the a posteriori L1 error bound (should be tighter than generic)
- `apriori_error`: the a priori error bound for comparison
- `Linf_estimate`: per-interval L∞ error estimates (vector)
"""
function noise_error_aposteriori(
    prep1::Vector{Vector{PreimageInfo}},
    prep2::AbstractVector,
    meas_vec::AbstractVector,
    numeric_error::Float64,
    B_density::Ulam,
    K::UlamProjectedUniformKernel,
    α::Float64,
    sumCi::Float64,
)
    K_size = length(meas_vec)
    length(K.B) == K_size ||
        throw(ArgumentError("the kernel is on $(length(K.B)) cells, the density on $K_size"))
    mod(K_size, length(prep1)) == 0 ||
        throw(ArgumentError("the estimation partition ($(length(prep1)) cells) does not refine into the density partition ($K_size cells)"))
    0 <= α < 1 || throw(ArgumentError("α = $α must be in [0, 1)"))
    noise_size_rel = Float64(noise_support(K), RoundDown)
    δ = Float64(1.0, RoundUp) / Float64(K_size, RoundDown)
    varRho = 2.0 ⊘₊ Float64(noise_size_rel, RoundDown)

    # A priori bound for comparison
    apriori_error = noise_error_apriori(K_size, noise_size_rel, α, sumCi) ⊕₊ numeric_error

    # Step 1: estimate ||(1-π)Lf||_W
    est1, meas_Lf, var_Lf = estimate_Wnorm_of_discretization_error(prep1, meas_vec, B_density)

    # Step 2: estimate Var(NLf)
    var_NLf = estimate_variation_after_noise(meas_Lf, var_Lf, K)
    tot_var_NLf = 0.0
    for v in var_NLf
        tot_var_NLf = tot_var_NLf ⊕₊ v
    end

    # Step 3: decompose error into three types

    # Type 1: (1-π)f, estimated as (1-π)NLf_δ + ...
    A1 = (δ / 2) ⊗₊ tot_var_NLf
    B1 = (δ / 2) ⊗₊ varRho

    # Type 2: W-norm projection error
    A2 = varRho ⊗₊ est1
    B2 = (δ / 2) ⊗₊ varRho

    # Type 3: second-order term
    est2 = estimate_L1_NL_discretization_error(var_NLf, prep2, δ, varRho)
    A3 = est2 ⊕₊ ((δ^2 / 4) ⊗₊ varRho ⊗₊ tot_var_NLf)
    B3 = (δ / 2) ⊗₊ varRho

    # Combine
    A = A1 ⊕₊ sumCi ⊗₊ (A2 ⊕₊ A3)
    Bval = B1 ⊕₊ sumCi ⊗₊ (B2 ⊕₊ B3)
    C = A ⊘₊ (1.0 ⊖₋ α)
    Dval = Bval ⊘₊ (1.0 ⊖₋ α)

    if Dval >= 1.0
        @warn "Denominator 1-D is non-positive (D=$Dval); bound does not converge. Returning Inf."
        return Inf, apriori_error, fill(Inf, length(prep1))
    end

    L1Error = (C ⊕₊ numeric_error) ⊘₊ (1.0 ⊖₋ Dval)

    # L∞ estimates per interval
    Kest = length(prep1)
    meas_NLf = estimate_measure_after_noise(meas_Lf, K)
    Linf_to_real_f = L1Error ⊘₊ Float64(noise_size_rel, RoundDown)
    Linf_estimate = zeros(Kest)
    for i in 1:Kest
        Linf_estimate[i] = var_NLf[i] ⊕₊ meas_NLf[i] ⊗₊ Float64(Kest, RoundUp) ⊕₊ Linf_to_real_f
    end

    return L1Error, apriori_error, Linf_estimate
end
