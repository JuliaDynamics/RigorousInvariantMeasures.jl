# Example: Comparing generic vs finer noise error estimates
#
# This example demonstrates the full workflow for computing both the
# standard DFLY-based error bound (distance_from_invariant_noise) and
# the finer a posteriori bound (noise_error_aposteriori) for a noisy
# dynamical system, then compares them.

using RigorousInvariantMeasures
using IntervalArithmetic

# ── 1. Define the dynamics ──────────────────────────────────────────

# A perturbed 4x map on [0,1]
D = mod1_dynamic(x -> 4x + RigorousInvariantMeasures.sinpi(8x) / 100)
println("Dynamic: 4x + sin(8πx)/100  mod 1")
println("  $(nbranches(D)) branches")

# ── 2. Assemble operator + noise kernel ─────────────────────────────

K_size = 512
B = Ulam(K_size)
noise_half_width = 16   # noise window = 2*16+1 = 33 bins
NK = UniformKernelUlamReflecting(B, noise_half_width)

noise_size_rel = (2 * noise_half_width + 1) / K_size
println("Partition size: $K_size")
println("Noise relative size: $noise_size_rel")

Q = DiscretizedOperator(B, D)
println("Operator assembled.")

# ── 3. Compute invariant measure via power iteration ────────────────

w = invariant_vector_noise(B, Q, NK; iter=30)
println("Invariant vector computed (power iteration, 30 steps).")
println("  sum(w)/K = $(sum(w)/K_size)  (should be ≈ 1)")

# ── 4. Compute norms of powers ──────────────────────────────────────

println("\nComputing norms of powers...")
norms = powernormboundsnoise(B; Q=Q, NK=NK)
println("  $(length(norms)) norms computed")
println("  last norm: $(norms[end])")

# ── 5. Generic DFLY error bound ─────────────────────────────────────

generic_error = distance_from_invariant_noise(B, Q, NK, w, norms)
println("\n── Generic DFLY error bound ──")
println("  ||f - f*||₁ ≤ $generic_error")

# ── 6. Finer a posteriori error bound ──────────────────────────────

println("\nPreparing finer estimate...")

# Spectral contraction rate
m = length(norms)
α = norms[end]^(1.0 / m)
if α >= 1.0
    α = 0.99
end

# Sum of norms Σ ||Q^i||
using RigorousInvariantMeasures: infinite_sum_norms
sumCi = infinite_sum_norms(norms)

# Numerical residual error: floating-point precision of eigenvector computation
# (NOT the interval-arithmetic residual, which includes discretization error
# that the specialized estimator already accounts for)
mQ = mid(Q)
v_fp = mQ.L * w
v_fp = NK * v_fp
numeric_error = sum(abs.(v_fp .- w)) / K_size

# Preimage data for the finer estimate
println("  Computing preimage data...")
prep1 = prepare_Wnorm_estimate(D, B)
println("  Computing derivative bounds...")
prep2 = prepare_derivative_bounds(D, B)

# The a posteriori estimate
println("  Computing a posteriori bound...")
L1err, apriori_err, Linf_est = noise_error_aposteriori(
    prep1, prep2, w, numeric_error, B, NK, α, sumCi,
)

println("\n── A posteriori finer error bound ──")
println("  ||f - f*||₁ ≤ $L1err  (finer)")
println("  ||f - f*||₁ ≤ $apriori_err  (a priori)")

# ── 7. Compare the two bounds ──────────────────────────────────────

println("\n── Comparison ──")
println("  Generic DFLY bound:       $generic_error")
println("  A posteriori finer bound: $L1err")
if L1err < generic_error
    ratio = generic_error / L1err
    println("  Improvement factor:       $(round(ratio, digits=2))x tighter")
else
    println("  (Finer bound is not tighter in this case)")
end

println("\n  Max L∞ estimate: $(maximum(Linf_est))")
println("  Numeric error:   $numeric_error")
println("  α (contraction): $α")
println("  Σ Cᵢ:            $sumCi")
