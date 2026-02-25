# Example: Visualization of a dynamical system with uniform noise
#
# This example demonstrates the plot_noisy_system function which creates
# a 3-panel visualization showing:
# 1. The dynamic with noise point cloud: x → T(x) + ξₙ
# 2. The approximated stationary measure
# 3. The spectrum of the discretized operator with unit circle

using RigorousInvariantMeasures
using IntervalArithmetic
using Plots
using LaTeXStrings

# Create the dynamic T(x) = 4x + 0.1*sin(2πx) mod 1
# This is an expanding map with small nonlinear perturbation
D = mod1_dynamic(x -> 4*x + 0.1*RigorousInvariantMeasures.sinpi(2*x), full_branch=true)

# Create Ulam basis with 256 elements
B = Ulam(256)

# Create noise kernel with half-width l=64
# This gives noise amplitude ξ = (2*64+1)/256 ≈ 0.504
l = 64
K = UniformKernelUlamPeriodic(B, l)

# Assemble the discretized transfer operator
Q = DiscretizedOperator(B, D)

# Compute the invariant vector with noise using power iteration
w = invariant_vector_noise(B, Q, K; iter=100)

# Get the midpoint of the operator matrix for spectrum computation
L = mid(Q).L

# Create the 3-panel visualization
p = plot_noisy_system(D, K, w, L; n_samples=100, n_points=300)

# Save to PNG
savefig(p, joinpath(@__DIR__, "noisy_system_plot.png"))

println("Plot saved to examples/noisy_system_plot.png")
