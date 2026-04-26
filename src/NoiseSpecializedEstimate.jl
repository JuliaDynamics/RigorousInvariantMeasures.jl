"""
Specialized a posteriori noise error estimator.

Ported from the Python `noise_specialized_measure_error.py` in the
`compinvmeas-python-release-noise` package.  These routines produce
significantly tighter error bounds than the generic DFLY-based approach
in `distance_from_invariant_noise` by exploiting the *actual computed
density vector* to obtain data-dependent bounds.
"""

using FastRounding
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
Information about a single preimage of an interval I under a branch Tₖ.
"""
struct PreimageInfo
    A::Float64            # left endpoint of preimage (mid, rounded down)
    B::Float64            # right endpoint of preimage (mid, rounded up)
    invtp::Float64        # upper bound for |1/T'| on preimage
    distortion::Float64   # upper bound for |T''/(T')²| on preimage
    boundary_terms::Vector{Tuple{Float64,Float64}}  # (branch-limit value inside I, |1/T'| there)
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
        invtp_val = Float64(mag(abs(inverse_derivative(br.f, TiI))), RoundUp)

        # Compute |T''/(T')²| (distortion) upper bound on preimage
        distortion_val = Float64(mag(abs(distortion(br.f, TiI))), RoundUp)

        # Check boundary terms: branch limits that overlap with interval I
        I_interval = hull(interval(B_est.p[ylabel]), interval(B_est.p[ylabel + 1]))
        bterms = Tuple{Float64,Float64}[]
        for bl in branch_limits[bnum]
            if !isempty_interval(intersect_interval(bl, I_interval))
                bl_invtp = Float64(mag(abs(inverse_derivative(br.f, bl))), RoundUp)
                push!(bterms, (mid(bl), bl_invtp))
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

            # Estimate of variation of Lf in I from this branch
            if isfinite(invtp_ab) && isfinite(distortion_ab)
                estimate_var_Ii_from_ab = invtp_ab ⊗₊ var_in_ab ⊕₊ distortion_ab ⊗₊ meas_in_ab
            else
                # Fallback for infinite derivatives (e.g. singularities):
                # use the trivial bound Var(Lf|_preimage) ≤ 2 * meas / δ
                estimate_var_Ii_from_ab = 2.0 ⊗₊ meas_in_ab ⊘₊ δ
            end

            # Boundary term contributions
            for (bl, b_invtp) in pinfo.boundary_terms
                bl_jmin = searchsortedlast(partition, bl)
                bl_jmin = clamp(bl_jmin, 1, K_density)
                b_contrib = Float64(abs(density_vec[bl_jmin]), RoundUp) ⊗₊ b_invtp
                if isfinite(b_contrib)
                    estimate_var_Ii_from_ab = estimate_var_Ii_from_ab ⊕₊ b_contrib
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

"""
    _reflect(i, size)

Reflecting boundary condition index mapping (0-based internally).
Maps index `i` (1-based) into `1:size` using reflection at boundaries.
"""
function _reflect(i::Int, sz::Int)
    # Convert to 0-based for the Python-style reflect
    i0 = i - 1
    if i0 < 0
        i0 = -i0
    end
    if i0 >= sz
        i0 = 2 * sz - i0 - 1
    end
    # Handle periodic wrapping for very large offsets
    i0 = mod(i0, 2 * sz)
    if i0 >= sz
        i0 = 2 * sz - i0 - 1
    end
    return i0 + 1  # back to 1-based
end

# ──────────────────────────────────────────────────────────────────
# 5.  Estimate Var(NLf)  (estimate_var_of_NLf)
# ──────────────────────────────────────────────────────────────────

"""
    estimate_variation_after_noise(meas_Lf, var_Lf, K::UniformKernelUlam)

Estimate `Var(NLf)` — the variation of Lf after applying the noise
operator N.

For each interval `i`, takes the minimum of:
  - Measure-based: boundary values divided by noise width
  - Variation-based: averaged variation over noise window

Returns `var_NLf` vector of length `Kest`.
"""
function estimate_variation_after_noise(
    meas_Lf::AbstractVector,
    var_Lf::AbstractVector,
    K::UniformKernelUlam,
)
    Kest = length(meas_Lf)
    k_basis = length(K.B)
    l = K.l
    rel_noise_size = Float64(2 * l + 1) / Float64(k_basis)

    var_NLf = zeros(Kest)

    abs_xiq2f = Int(floor(rel_noise_size * Kest / 2))
    abs_xiq2c = Int(ceil(rel_noise_size * Kest / 2))

    for i in 1:Kest
        # Measure-based estimate: boundary terms / noise width
        est1 = (
            meas_Lf[_reflect(i - abs_xiq2c, Kest)] ⊕₊
            meas_Lf[_reflect(i + abs_xiq2c, Kest)] ⊕₊
            meas_Lf[_reflect(i - abs_xiq2f, Kest)] ⊕₊
            meas_Lf[_reflect(i + abs_xiq2f, Kest)]
        ) ⊘₊ Float64(rel_noise_size, RoundDown)

        # Variation-based estimate: averaged variation over window
        est2 = 0.0
        for j in (i - abs_xiq2c):(i + abs_xiq2c)
            est2 = est2 ⊕₊ var_Lf[_reflect(j, Kest)]
        end
        if abs_xiq2c > 0
            est2 = est2 ⊘₊ Float64(2 * abs_xiq2c, RoundDown)
        end

        var_NLf[i] = min(est1, est2)
    end

    return var_NLf
end

# ──────────────────────────────────────────────────────────────────
# 6.  Estimate measure of NLf  (estimate_meas_of_NLf)
# ──────────────────────────────────────────────────────────────────

"""
    estimate_measure_after_noise(meas_Lf, K::UniformKernelUlam)

Estimate the measure of `NLf` per interval by averaging `meas_Lf`
over the noise window.

Returns `meas_NLf` vector of length `Kest`.
"""
function estimate_measure_after_noise(
    meas_Lf::AbstractVector,
    K::UniformKernelUlam,
)
    Kest = length(meas_Lf)
    k_basis = length(K.B)
    l = K.l
    rel_noise_size = Float64(2 * l + 1) / Float64(k_basis)

    abs_noise_size = Int(round(rel_noise_size * Kest))
    abs_xiq2 = abs_noise_size ÷ 2

    meas_NLf = zeros(Kest)

    for i in 1:Kest
        s = 0.0
        for j in (i - abs_xiq2):(i + abs_xiq2)
            s = s ⊕₊ meas_Lf[_reflect(j, Kest)]
        end
        meas_NLf[i] = s ⊘₊ Float64(abs_noise_size, RoundDown)
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
                d = Float64(mag(abs(derivative(br.f, J))), RoundUp)
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
                            B_density::Ulam, K::UniformKernelUlam,
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
- `K`: the `UniformKernelUlam` noise kernel
- `α`: spectral contraction rate (should be < 1)
- `sumCi`: `Σ Cᵢ` from the norms of powers

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
    K::UniformKernelUlam,
    α::Float64,
    sumCi::Float64,
)
    K_size = length(meas_vec)
    k_basis = length(K.B)
    l = K.l
    noise_size_rel = Float64(2 * l + 1) / Float64(k_basis)
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
