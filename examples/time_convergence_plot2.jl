using RigorousInvariantMeasures
using IntervalArithmetic

using Glob
using Plots
using StatsPlots
using FastRounding
using LaTeXStrings

using Serialization

ENV["GKSwstype"] = "nul" # for headless displays

function get_experiment(prefix)
    if prefix == "Lorenz3"
        f = n -> begin
            D0 = Lorenz()
            D = D0 ∘ D0 ∘ D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix == "Lorenz2"
        f = n -> begin
            D0 = Lorenz()
            D = D0 ∘ D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix == "4x_perturbed_Ulam"
        f = n -> begin
            D = mod1_dynamic(x -> 4x + RigorousInvariantMeasures.sinpi(8x) / 100)
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix == "4x_perturbed_Hat"
        f = n -> begin
            D = mod1_dynamic(x -> 4x + RigorousInvariantMeasures.sinpi(8x) / 100)
            B = Hat(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix == "Lanford_Ulam"
        f = n -> begin
            D = mod1_dynamic(x -> 2x + 0.5 * x * (1 - x))
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix == "Lanford3_Hat"
        f = n -> begin
            D0 = mod1_dynamic(x -> 2x + 0.5 * x * (1 - x))
            # Taking an iterate is necessary here to get a DFLY inequality with A < 1
            D = D0 ∘ D0 ∘ D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix == "PiecewiseLinear2"
        f =
            n -> begin
                D = PwMap(
                    [
                        x -> 17x / 5,
                        x -> (34 * ((17x - 5) / 17) / 25 + 3) * ((17x - 5) / 17),
                        x -> (34 * ((17x - 10) / 17) / 25 + 3) * ((17x - 10) / 17),
                        x -> 17 * ((17x - 15) / 17) / 5,
                    ],
                    [
                        Interval(0),
                        Interval(5) / 17,
                        Interval(10) / 17,
                        Interval(15) / 17,
                        Interval(1),
                    ],
                    [
                        Interval(0) Interval(1)
                        Interval(0) Interval(1)
                        Interval(0) Interval(1)
                        Interval(0) @interval(0.4)
                    ],
                )
                B = Ulam(n)
                Q = DiscretizedOperator(B, D)
                return B, D, Q
            end
    elseif prefix == "PiecewiseLinear"
        f =
            n -> begin
                D = PwMap(
                    [x -> 2.5x, x -> 4x - 1, x -> 4x - 2, x -> 4x - 3],
                    [
                        @interval(0),
                        @interval(0.25),
                        @interval(0.5),
                        @interval(0.75),
                        @interval(1)
                    ],
                )
                B = Ulam(n)
                Q = DiscretizedOperator(B, D)
                return B, D, Q
            end
    elseif prefix == "175"
        f = n -> begin
            D = mod1_dynamic(x -> 17 // 5 * x)
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix == "175_nonlinear"
        f =
            n -> begin
                D = PwMap(
                    [
                        x -> 17x / 5,
                        x -> 34(x - 5 // 17)^2 / 25 + 3(x - 5 // 17),
                        x -> 34(x - 10 // 17)^2 / 25 + 3(x - 10 // 17),
                        x -> 17(x - 15 // 17) / 5,
                    ],
                    [0, @interval(5 / 17), @interval(10 / 17), @interval(15 / 17), 1],
                    [0 1; 0 1; 0 1; 0 @interval(0.4)],
                )
                B = Ulam(n)
                Q = DiscretizedOperator(B, D)
                return B, D, Q
            end
    else
        @error "Unknown experiment"
    end
end

function save_coarse_data(prefix, K)
    f = get_experiment(prefix)
    compute_coarse_grid_quantities(f, 4; m = 24)  # warm-up for precompilation
    for n in 2 .^ K
        print("Coarse+fine $n...")
        F, C = compute_coarse_grid_quantities(f, n; m = 24)
        serialize("$prefix-$n-coarse.juliaserialize", C)
        serialize("$prefix-$n-fine.juliaserialize", F)
        println("done.")
    end
end

function save_fine_data(prefix, K)
    f = get_experiment(prefix)
    compute_fine_grid_quantities(f, 4)  # warm-up for precompilation
    for n in 2 .^ K
        print("Fine only $n...")
        F = compute_fine_grid_quantities(f, n)
        serialize("$prefix-$n-fine.juliaserialize", F)
        println("done.")
    end
end

"""
Plot error vs. time in a double-log scale
"""
function plot_error_time(prefix)
    pyplot()

    Cdict = Dict{Int64,RigorousInvariantMeasures.CoarseGridQuantities}()
    Fdict = Dict{Int64,RigorousInvariantMeasures.FineGridQuantities}()
    for filename in glob("$prefix-*-coarse.juliaserialize")
        C = deserialize(filename)
        Cdict[length(C.B)] = C
    end
    for filename in glob("$prefix-*-fine.juliaserialize")
        F = deserialize(filename)
        Fdict[length(F.B)] = F
    end
    ks = collect(keys(Cdict) ∪ keys(Fdict))
    logks = Int.(log2.(ks))

    onegrid_errors = Float64[]
    onegrid_times = Float64[]
    onegrid_n = Int[]
    for (n, C) in sort(Cdict)
        error, time_breakdown = one_grid_estimate(C, Fdict[n])
        push!(onegrid_errors, error)
        push!(onegrid_times, sum(time_breakdown))
        push!(onegrid_n, n)
    end
    onegrid_k = Int.(log2.(onegrid_n))

    p = plot(
        onegrid_times,
        onegrid_errors,
        title = "Errors vs. times, experiment $prefix",
        yscale = :log10,
        xscale = :log10,
        color = :red,
        zcolor = onegrid_k,
        colorbar_ticks = (logks, map(x -> "2^$x", logks)),
        label = "One-grid strategy",
        legend = :outerbottom,
        markershape = :square,
        xlabel = "Total CPU Time/s",
        ylabel = "Error bound proved",
        size = (700, 700),
    )

    for (n, C) in sort(Cdict)
        twogrid_errors = Float64[]
        twogrid_times = Float64[]
        twogrid_n = Int[]
        for (n_fine, F) in sort(Fdict)
            if n_fine <= n
                continue
            end
            error, time_breakdown = two_grid_estimate(C, F)
            push!(twogrid_errors, error)
            push!(twogrid_times, sum(time_breakdown))
            push!(twogrid_n, n_fine)
        end
        if all(isinf.(twogrid_errors))
            continue
        end
        twogrid_k = Int.(log2.(twogrid_n))
        k = Int(log2(n))
        p = plot(
            p,
            twogrid_times,
            twogrid_errors,
            label = "Two-grid strategy (n_c=2^$k)",
            markershape = :circle,
            zcolor = twogrid_k,
        )
    end
    # t0 = 1e-2
    # t1 = 4e-2
    # e0 = 1e-2
    # exp = 1
    # C = e0/t0^(-exp)
    # e1 = C*t1^(-exp)
    # p = plot(p, [t0,t1], [e0,e1], legend=nothing, color=:black)
    # annotate!(t1, e1, ("err ∼ t^(-1)", :left, 10))
    # exp = 1/2
    # e0 = e0 * 1.2
    # C = e0/t0^(-exp)
    # e1 = C*t1^(-exp)
    # p = plot(p, [t0,t1], [e0,e1], legend=nothing, color=:black)
    # annotate!(t1, e1, ("err ∼ t^(-1/2)", :left, 10))

    savefig("$prefix-time-experiment.pdf")
    p
end

"""
Variant of refine_norms_of_powers that returns the refined norms even if they are worse
than the best choice
"""
function get_refined_norms(norms::Vector, m)
    refined_norms = fill(NaN, m)
    better_norms = fill(NaN, m)
    refined_norms[1] = norms[1]
    better_norms[1] = norms[1]
    for k = 2:m
        refined_norms[k] = minimum(better_norms[i] ⊗₊ better_norms[k-i] for i = 1:k-1)
        if k <= length(norms)
            better_norms[k] = min(refined_norms[k], norms[k])
        else
            better_norms[k] = refined_norms[k]
        end
    end
    return refined_norms
end

function plot_norm_bounds_kinds(prefix, n; num_norms = 30, num_norms_computed = 10)
    # TODO: these are recomputed, not loaded atm
    f = get_experiment(prefix)
    B, D, Q = f(n)
    normQ = opnormbound(B, weak_norm(B), Q)
    computed_norms = fill(Inf, num_norms)
    actual_computed_norms =
        norms_of_powers(B, weak_norm(B), num_norms_computed, Q, integral_covector(B))
    computed_norms[1:min(num_norms, num_norms_computed)] =
        actual_computed_norms[1:min(num_norms, num_norms_computed)]

    trivial_norms = norms_of_powers_trivial(normQ, num_norms)
    (dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, num_norms)
    norms = min.(trivial_norms, computed_norms, dfly_norms)
    refined_norms = get_refined_norms(norms, num_norms)
    p = plot(
        trivial_norms,
        label = "||Q||^k",
        yscale = :log10,
        legend = :bottomleft,
        title = "Available norm bounds, $prefix",
        xlabel = "k",
        ylabel = "bound to ||Q^k_U||",
        markershape = :circle,
    )
    plot!(p, dfly_norms, label = "DFLY-based", markershape = :diamond)

    plot!(p, refined_norms, label = "submultiplicativity", markershape = :+)
    plot!(p, computed_norms, label = "computational bounds", markershape = :star4)

    savefig(p, "norm-bounds-$prefix-$n.pdf")
end

function plot_norm_bounds_kinds_twogrid(
    prefix,
    n,
    n_fine;
    num_norms = 30,
    num_norms_computed = 10,
)
    # TODO: these are recomputed, not loaded atm
    f = get_experiment(prefix)
    B, D, Q = f(n)
    normQ = opnormbound(B, weak_norm(B), Q)
    B_fine, D_fine, Q_fine = f(n_fine)
    normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
    coarse_norms = powernormbounds(B, D, Q = Q, m = num_norms_computed)
    fine_norms = norms_of_powers_from_coarser_grid(
        B_fine,
        B,
        D,
        refine_norms_of_powers(coarse_norms, num_norms),
        normQ_fine,
    )
    trivial_norms = norms_of_powers_trivial(normQ_fine, num_norms)
    (dfly_strongs, dfly_norms) = norms_of_powers_dfly(B_fine, D_fine, num_norms)
    norms = min.(trivial_norms, fine_norms, dfly_norms)
    refined_norms = get_refined_norms(norms, num_norms)
    p = plot(
        trivial_norms,
        label = "||Q||^k",
        yscale = :log10,
        legend = :bottomleft,
        title = "Available norm bounds\n$prefix, n_f=$n_fine, n=$n",
        xlabel = "k",
        ylabel = "bound to ||Q^k_U||",
        markershape = :circle,
    )
    plot!(p, dfly_norms, label = "DFLY-based", markershape = :diamond)

    plot!(p, refined_norms, label = "submultiplicativity", markershape = :+)
    plot!(p, fine_norms, label = "two-grid bounds", markershape = :star4, color = "blue")

    savefig(p, "norm-bounds-$prefix-$n-$n_fine.pdf")
end

function time_breakdown_plot(
    prefix,
    twogrid_nC;
    filter_coarse = nC -> true,
    filter_fine = nF -> true,
)
    Cdict = Dict{Int64,RigorousInvariantMeasures.CoarseGridQuantities}()
    Fdict = Dict{Int64,RigorousInvariantMeasures.FineGridQuantities}()
    for filename in glob("$prefix-*-coarse.juliaserialize")
        C = deserialize(filename)
        Cdict[length(C.B)] = C
    end
    for filename in glob("$prefix-*-fine.juliaserialize")
        F = deserialize(filename)
        Fdict[length(F.B)] = F
    end

    onegrid_errors = Float64[]
    onegrid_times = zeros(5, 0)
    onegrid_n = Int[]
    onegrid_labels = LaTeXString[]
    for (n, C) in sort(Cdict)
        if !filter_coarse(n)
            continue
        end
        error, time_breakdown = one_grid_estimate(C, Fdict[n])
        push!(onegrid_errors, error)
        onegrid_times = [onegrid_times time_breakdown]
        push!(onegrid_n, n)
        push!(onegrid_labels, LaTeXString("\$2^{$(Int(log2(n)))}\$"))
    end

    twogrid_errors = Float64[]
    twogrid_times = zeros(5, 0)
    twogrid_n = Int[]
    twogrid_labels = LaTeXString[]
    C = Cdict[twogrid_nC]
    for (n_fine, F) in sort(Fdict)
        if n_fine <= twogrid_nC
            continue
        end
        if !filter_fine(n_fine)
            continue
        end
        error, time_breakdown = two_grid_estimate(C, F)
        push!(twogrid_errors, error)
        twogrid_times = [twogrid_times time_breakdown]
        push!(twogrid_n, n_fine)
        push!(twogrid_labels, LaTeXString("\$2^{$(Int(log2(n_fine)))}\$"))
    end



    p1 = groupedbar(
        onegrid_times'[:, end:-1:1],
        bar_position = :stack,
        legend_position = :topleft,
        label = ["dfly coefficients" "matrix assembly" "eigenvalue computation" "norms of powers" "error estimation"][
            :,
            end:-1:1,
        ],
        title = "One-grid",
        ylabel = "CPU Time/s",
        xticks = (1:length(onegrid_n), onegrid_labels),
        link = :y,
        size = (700, 700),
    )

    p2 = groupedbar(
        twogrid_times'[:, end:-1:1],
        bar_position = :stack,
        legend = :topleft,
        label = ["dfly coefficients" "coarse matrix+norms" "matrix assembly" "eigenvalue computation" "error estimation"][
            :,
            end:-1:1,
        ],
        title = "Two-grid, n_c=$(twogrid_nC)",
        xticks = (1:length(twogrid_n), twogrid_labels),
        link = :y,
        size = (700, 700),
    )

    p3 = plot(
        onegrid_n,
        onegrid_errors,
        ylabel = "Error bound",
        mark = :dot,
        yscale = :log10,
        xscale = :log10,
        xticks = (onegrid_n, onegrid_labels),
        #        xticks = nothing,
        label = "One-grid",
        legend = false,
        link = :y,
        size = (700, 700),
    )

    p4 = plot(
        twogrid_n,
        twogrid_errors,
        # ylabel = "Error bound",
        mark = :dot,
        yscale = :log10,
        xscale = :log10,
        color = :red,
        xticks = (twogrid_n, twogrid_labels),
        #        xticks = nothing,
        label = "Two-grid, n_c=$(twogrid_nC)",
        legend = false,
        link = :y,
        size = (700, 700),
    )

    p = plot(p1, p2, p3, p4, plot_title = "Time breakdown, experiment $prefix")
    savefig("$prefix-breakdown.pdf")
    p
end
