using InvariantMeasures
using ValidatedNumerics
using Glob
using Plots

using Serialization

ENV["GKSwstype"]="nul" # for headless displays

LorenzMap(θ, α) = PwMap([x->θ*(0.5-x)^α, x->1-θ*(x-0.5)^α],
                    [@interval(0), @interval(0.5), @interval(1)]; infinite_derivative=true)

function f(n)
    D0 = LorenzMap(109/64, 51/64)
    D = D0 ∘ D0 ∘ D0
    B = Ulam(n)
    Q = DiscretizedOperator(B, D)
    return B, D, Q
end

prefix = "Lorenz3"

function save_coarse_data(K)
    compute_coarse_grid_quantities(f, 4; m=24)  # warm-up for precompilation
    for n = 2 .^ K
        print("Coarse+fine $n...")
        F, C = compute_coarse_grid_quantities(f, n; m=24)
        serialize("$prefix-$n-coarse.juliaserialize", C)
        serialize("$prefix-$n-fine.juliaserialize", F)
        println("done.")
    end
end

function save_fine_data(K)
    compute_fine_grid_quantities(f, 4)  # warm-up for precompilation
    for n = 2 .^ K
        print("Fine only $n...")
        F = compute_fine_grid_quantities(f, n)
        serialize("$prefix-$n-fine.juliaserialize", F)
        println("done.")
    end
end

function plot_data()
    Cdict = Dict{Int64, InvariantMeasures.CoarseGridQuantities}()
    Fdict = Dict{Int64, InvariantMeasures.FineGridQuantities}()
    for filename in glob("$prefix-*-coarse.juliaserialize")
        C = deserialize(filename)
        Cdict[length(C.B)] = C
    end
    for filename in glob("$prefix-*-fine.juliaserialize")
        F = deserialize(filename)
        Fdict[length(F.B)] = F
    end

    onegrid_errors = Float64[]
    onegrid_times = Float64[]
    onegrid_n = Plots.PlotText[]
    for (n, C) in sort(Cdict)
        error, time_breakdown = one_grid_estimate(C, Fdict[n])
        push!(onegrid_errors, error)
        push!(onegrid_times, sum(time_breakdown))
        push!(onegrid_n, Plots.PlotText(" $n", font(8)))
    end

    p = plot(
        onegrid_times, onegrid_errors,
        title = "Errors vs. times, experiment $prefix",
        yscale = :log10,
        xscale = :log10,
        color = :red,
        label = "One-grid strategy",
        legend = :outerbottom,
        markershape=:circle,
        xlabel = "Time/s",
        ylabel = "Error",
        size = (700, 700)
    )
    annotate!(onegrid_times/1.08, onegrid_errors, onegrid_n, :right)

    for (n, C) in sort(Cdict)
        twogrid_errors = Float64[]
        twogrid_times = Float64[]
        twogrid_n = Plots.PlotText[]
        for (n_fine, F) in sort(Fdict)
            if n_fine <= n
                continue
            end
            error, time_breakdown = two_grid_estimate(C, F)
            push!(twogrid_errors, error)
            push!(twogrid_times, sum(time_breakdown))
            push!(twogrid_n, Plots.PlotText(" $n_fine", font(8)))
        end
        if all(isinf.(twogrid_errors))
            continue
        end
        p = plot(p,
            twogrid_times, twogrid_errors,
            label = "Two-grid strategy (n_c=$n)",
            markershape=:circle
        )
        annotate!(twogrid_times/1.08, twogrid_errors, twogrid_n, :right)
    end
    savefig("$prefix-time-experiment.pdf")
end
