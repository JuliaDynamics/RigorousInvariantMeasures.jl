# Example: Visualization of the Lorenz map with uniform noise

using RigorousInvariantMeasures
using IntervalArithmetic
using Plots
using LaTeXStrings

# Create the Lorenz map with default parameters (θ=109/64, α=51/64)
# This is an intermittent map with a neutral fixed point at x=0.5
D = Lorenz()

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
savefig(p, joinpath(@__DIR__, "Lorenz_noisy_system_plot.png"))

println("Plot saved to examples/Lorenz_noisy_system_plot.png")
