using InvariantMeasures
using ValidatedNumerics
using Glob
using Plots
using LaTeXStrings

using Serialization

ENV["GKSwstype"]="nul" # for headless displays

LorenzMap(θ, α) = PwMap([x->θ*(0.5-x)^α, x->1-θ*(x-0.5)^α],
                    [@interval(0), @interval(0.5), @interval(1)]; infinite_derivative=true)

function get_experiment(prefix)
    if prefix=="Lorenz3"
        f = n-> begin
            D0 = LorenzMap(109/64, 51/64)
            D = D0 ∘ D0 ∘ D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q    
        end
    elseif prefix=="Lorenz2"
        f = n-> begin
            D0 = LorenzMap(109/64, 51/64)
            D = D0 ∘ D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q    
        end
    elseif prefix=="4x_perturbed_Ulam"
        f = n-> begin
            D = mod1_dynamic(x -> 4*x + 0.01*InvariantMeasures.sinpi(8*x))
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="4x_perturbed_Hat"
        f = n-> begin
            D = mod1_dynamic(x -> 4*x + 0.01*InvariantMeasures.sinpi(8*x))
            B = Hat(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    else
        @error "Unknown experiment"
    end
end

function save_coarse_data(prefix, K)
    f = get_experiment(prefix)
    compute_coarse_grid_quantities(f, 4; m=24)  # warm-up for precompilation
    for n = 2 .^ K
        print("Coarse+fine $n...")
        F, C = compute_coarse_grid_quantities(f, n; m=24)
        serialize("$prefix-$n-coarse.juliaserialize", C)
        serialize("$prefix-$n-fine.juliaserialize", F)
        println("done.")
    end
end

function save_fine_data(prefix, K)
    f = get_experiment(prefix)
    compute_fine_grid_quantities(f, 4)  # warm-up for precompilation
    for n = 2 .^ K
        print("Fine only $n...")
        F = compute_fine_grid_quantities(f, n)
        serialize("$prefix-$n-fine.juliaserialize", F)
        println("done.")
    end
end

function plot_error_time(prefix)
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
        xlabel = "Total CPU Time/s",
        ylabel = "Error bound proved",
        size = (700, 700)
    )
    annotate!(onegrid_times[[1,end]], onegrid_errors[[1,end]]*1.16, onegrid_n[[1,end]], :bottom)

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
        annotate!(twogrid_times[[1,end]], twogrid_errors[[1,end]]*1.16, twogrid_n[[1,end]], :bottom)
    end
    savefig("$prefix-time-experiment.pdf")
end

function plot_norm_bounds_kinds(prefix, n; num_norms=30, num_norms_computed=10)
    C = deserialize("$prefix-$n-coarse.juliaserialize")
    F = deserialize("$prefix-$n-fine.juliaserialize")
    computed_norms = fill(NaN, num_norms)
    # TODO: need to refactor stuff so that norms_of_powers are saved
    actual_computed_norms = norms_of_powers(C.B, weak_norm(C.B), num_norms_computed, Q, integral_covector(C.B))
    computed_norms[1:min(num_norms,num_norms_computed)] = actual_computed_norms[1:min(num_norms,num_norms_computed)]

    trivial_norms = norms_of_powers_trivial(F.normQ, num_norms)
    (dfly_strongs, dfly_norms) = norms_of_powers_dfly(C.B, C.D, num_norms; dfly_coefficients=C.dfly_coefficients)
    norms = min.(trivial_norms, computed_norms, dfly_norms)
    better_norms = refine_norms_of_powers(norms, num_norms)
    p = plot(trivial_norms,
        label = L"$\|Q\|^k$",
        yscale = :log10,
        legend = :bottomleft,
        title = "Available norm bounds, $prefix",
        xlabel = L"$k$",
        ylabel = L"bound to $\|Q\^k|_U\|$",
        markershape=:circle
        )
    plot!(p, computed_norms,
        label = "computational bounds",
        markershape=:x
        )

    plot!(p, dfly_norms,
        label = L"DFLY $2\times 2$ matrix bounds",
        markershape=:diamond,
        )

    plot!(p, better_norms,
        label = "min(previous) + refinement",
        markershape=:+
        )

    savefig(p, "norm-bounds-$prefix-$n.pdf")
end