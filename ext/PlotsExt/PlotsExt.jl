module PlotsExt

using RigorousInvariantMeasures, Plots
using RecipesBase
using LaTeXStrings
using IntervalArithmetic: mid, Interval
#COV_EXCL_START
@recipe f(::Type{ApproxInducedLSV}, D::ApproxInducedLSV) = x -> plottable(D, x)
#COV_EXCL_STOP

#COV_EXCL_START
@userplot PlotC2
@recipe function f(h::PlotC2)
    if length(h.args) != 2 ||
       (typeof(h.args[1]) != C2) ||
       !(typeof(h.args[2]) <: AbstractVector)
        error("Plot C2 needs as an input a C2 Basis and a vector")
    end

    B = h.args[1]
    w = h.args[2]

    seriestype := :path
    collect(B), mid.(w)
end
#COV_EXCL_STOP

#COV_EXCL_START
"""
Plots a function in the Hat basis
"""
@recipe function f(B::Hat, w::AbstractVector)

    legend --> :bottomright

    if eltype(w) <: Interval
        w = mid.(w)
    end

    @series begin
        seriestype --> :path
        label --> L"f_{\delta}"
        ylims --> (0, NaN)
        B.p, vcat(w, w[end])
    end
end

"""
Displays error on a function in the Hat basis
"""
@recipe function f(B::Hat, error::Number, w)

    if eltype(w) <: Interval
        w = mid.(w)
    end

    if isfinite(error)
        @series begin
            seriestype --> :path
            seriesalpha --> 0.5
            fillrange --> vcat(w, w[end]) .- error
            label --> "Error area"
            B.p, vcat(w, w[end]) .+ error
        end
    end
end
#COV_EXCL_STOP

#COV_EXCL_START
"""
Plots a function in the Ulam basis
"""
@recipe function f(B::Ulam, w::AbstractVector)

    legend --> :bottomright

    if eltype(w) <: Interval
        w = mid.(w)
    end

    @series begin
        seriestype --> :steppost
        label --> L"f_{\delta}"
        ylims --> (0, NaN)
        B.p, vcat(w, w[end])
    end
end

"""
Displays error on a function in the Ulam basis

The w argument is unused, but kept for compatibility with other functions
for different bases
"""
@recipe function f(B::Ulam, error::Number, w = nothing)

    if isfinite(error)
        @series begin
            seriestype --> :path
            seriesalpha --> 0.5
            fillrange --> 0
            label --> "L1 Error"
            [0; sqrt(error)], [sqrt(error); sqrt(error)]
        end
    end
end
#COV_EXCL_START
@recipe f(::Type{PM}, D::PM) where {PM<:PwMap} = x -> plottable(D, x)
#COV_EXCL_STOP

#COV_EXCL_START
"""
Plots a function in the HatNP basis
"""
@recipe function f(B::HatNP, w::AbstractVector)

    legend --> :bottomright

    if eltype(w) <: Interval
        w = mid.(w)
    end

    @series begin
        seriestype --> :path
        label --> L"f_{\delta}"
        ylims --> (0, NaN)
        B.p, vcat(w, w[end])
    end
end

"""
Displays error on a function in the Hat basis
"""
@recipe function f(B::HatNP, error::Number, w)

    if eltype(w) <: Interval
        w = mid.(w)
    end

    if isfinite(error)
        @series begin
            seriestype --> :path
            seriesalpha --> 0.5
            fillrange --> vcat(w, w[end]) .- error
            label --> "Error area"
            B.p, vcat(w, w[end]) .+ error
        end
    end
end
#COV_EXCL_STOP

using LinearAlgebra: eigen

import RigorousInvariantMeasures: plot_noisy_system

#COV_EXCL_START
"""
    plot_noisy_system(D::PwMap, K::UniformKernelUlam, w, L; n_samples=1000, n_points=500)

Plot a 3-panel visualization of a dynamical system with uniform noise:

1. **Dynamic with noise**: Point cloud showing x → T(x) + ξₙ where ξₙ is sampled
   uniformly from [-ξ/2, ξ/2] with amplitude ξ = (2l+1)/k extracted from the kernel.

2. **Stationary measure**: The approximated invariant density using the Ulam basis.

3. **Spectrum**: Non-rigorous eigenvalues of the discretized operator with unit circle.

# Arguments
- `D::PwMap`: The piecewise monotonic map
- `K::UniformKernelUlam`: Uniform noise kernel (provides noise amplitude)
- `w`: Coefficients of the stationary measure in the Ulam basis
- `L`: Discretized operator matrix (used for spectrum computation)

# Keyword Arguments
- `n_samples::Int=1000`: Number of noise samples per x-point in the point cloud
- `n_points::Int=500`: Number of x-points to sample for the point cloud
"""
function plot_noisy_system(
    D::PwMap,
    K::UniformKernelUlam{BC},
    w,
    L;
    n_samples::Int = 1000,
    n_points::Int = 500
) where {BC}
    B = K.B
    k = length(B)
    l = K.l
    ξ = (2l + 1) / k  # noise amplitude

    # Extract domain bounds from PwMap
    dom = RigorousInvariantMeasures.domain(D)
    x_min = Float64(mid(Interval(dom[1])))
    x_max = Float64(mid(Interval(dom[2])))

    # Panel 1: Dynamic with noise point cloud
    xs = range(x_min, x_max, length = n_points)
    x_cloud = Float64[]
    y_cloud = Float64[]

    for x in xs
        Tx = plottable(D, x)
        if Tx !== nothing
            for _ in 1:n_samples
                ξₙ = ξ * (rand() - 0.5)  # Uniform in [-ξ/2, ξ/2]
                push!(x_cloud, x)
                push!(y_cloud, Tx + ξₙ)
            end
        end
    end

    p1 = scatter(
        x_cloud, y_cloud,
        markersize = 0.5,
        markeralpha = 0.1,
        markerstrokewidth = 0,
        color = :blue,
        label = nothing,
        xlabel = "x",
        ylabel = "T(x) + ξ",
        title = "Dynamic with noise (ξ = $(round(ξ, digits=4)))"
    )
    # Add the deterministic map
    x_det = range(x_min, x_max, length = 1000)
    y_det = [plottable(D, x) for x in x_det]
    valid_idx = findall(y -> y !== nothing, y_det)
    plot!(p1, x_det[valid_idx], Float64.(y_det[valid_idx]),
          color = :red, linewidth = 1.5, label = "T(x)")

    # Panel 2: Stationary measure
    w_plot = eltype(w) <: Interval ? mid.(w) : w
    p2 = plot(
        B.p, vcat(w_plot, w_plot[end]),
        seriestype = :steppost,
        fillrange = 0,
        fillalpha = 0.3,
        color = :green,
        linewidth = 1.5,
        label = nothing,
        xlabel = "x",
        ylabel = "density",
        title = "Stationary measure",
        ylims = (0, maximum(w_plot) * 1.1)
    )

    # Panel 3: Spectrum with unit circle
    eigenvalues = eigen(Matrix(L)).values

    # Unit circle
    θ = range(0, 2π, length = 200)
    circle_x = cos.(θ)
    circle_y = sin.(θ)

    p3 = plot(
        circle_x, circle_y,
        color = :gray,
        linestyle = :dash,
        linewidth = 1,
        label = "Unit circle",
        aspect_ratio = :equal,
        xlabel = "Re(λ)",
        ylabel = "Im(λ)",
        title = "Spectrum"
    )
    scatter!(p3, real.(eigenvalues), imag.(eigenvalues),
             color = :purple,
             markersize = 3,
             label = "Eigenvalues")

    # Combine into 2x2 layout with empty bottom-right
    l = @layout [a b; c _]
    return plot(p1, p2, p3, layout = l, size = (900, 800))
end
#COV_EXCL_STOP

end
