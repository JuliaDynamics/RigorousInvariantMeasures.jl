# Example: Visualization of the BZ map with uniform noise

using RigorousInvariantMeasures
using IntervalArithmetic
using Plots
using LaTeXStrings
using TaylorSeries

# Define BZ map with proper handling of the singularity at x=1/8
# Uses exp(α*log(x)) for fractional powers with special handling for x=0
function BZ_fixed()
    A_big = Interval{BigFloat}(BigFloat("0.5060735690368223513195993710530479569801417368282037493809901142182256388277772"))
    BZA = Interval{Float64}(A_big)
    B_big = Interval{BigFloat}(
        BigFloat("0.02328852830307032296813220750095076307514120085284790788241725176646123060202265"),
        BigFloat("0.02328852830307032296813220750095076307514120085284790788241849636951680283042265"),
    )
    BZB = Interval{Float64}(B_big)
    C_big = Interval{BigFloat}(BigFloat("0.121205692738975111744666848150620569782497212127938371936404761693002104361654"))
    BZC = Interval{Float64}(C_big)

    α = Interval(1.0) / 3

    # Safe power function for Interval type
    function safe_power(x::Interval, α)
        if inf(x) <= 0 && sup(x) >= 0
            if sup(x) > 0
                return hull(Interval(0), exp(α * log(Interval(max(1e-300, inf(x)), sup(x)))))
            else
                return Interval(0)
            end
        else
            return exp(α * log(x))
        end
    end

    # For Taylor series, just use exp(α*log(x)) directly
    # The constant term check handles the singularity
    function safe_power(x::Taylor1, α)
        x0 = x[0]  # constant term (an Interval)
        if inf(x0) <= 0 && sup(x0) >= 0
            # Near singularity - use a regularized version
            # Add small positive value to avoid log(0)
            x_reg = x + Interval(1e-15)
            return exp(α * log(x_reg))
        else
            return exp(α * log(x))
        end
    end

    T_left_leq_1_8(x) = (BZA - safe_power(Interval(1.0) / 8 - x, α)) * exp(-x) + BZB
    T_left_geq_1_8(x) = (BZA + safe_power(x - Interval(1.0) / 8, α)) * exp(-x) + BZB
    T_right(x) = BZC * (10x * exp(-Interval(10) / 3 * x))^(19) + BZB

    return PwMap(
        [T_left_leq_1_8, T_left_geq_1_8, T_right],
        [Interval(0), Interval(1) / 8, Interval(3) / 10, Interval(1)],
        [
            T_left_leq_1_8(Interval(0)) T_left_leq_1_8(Interval(1) / 8)
            T_left_geq_1_8(Interval(1) / 8) T_left_geq_1_8(Interval(3) / 10)
            T_right(Interval(3) / 10) T_right(Interval(1))
        ],
    )
end

# Create the BZ map
D = BZ_fixed()

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
savefig(p, joinpath(@__DIR__, "BZ_noisy_system_plot.png"))

println("Plot saved to examples/BZ_noisy_system_plot.png")
