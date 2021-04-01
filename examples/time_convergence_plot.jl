using InvariantMeasures
using ValidatedNumerics

using Plots
using LaTeXStrings
using StatsPlots

pgfplotsx()

function onegrid(T, Btype, size)

    time_assembling = @elapsed begin
        D = Mod1Dynamic(T)
        # different backend, a tad slower
        # D = mod1_dynamic(x -> x->2*x+0.5*x*(1-x))
        B = Btype(size)
        Q = DiscretizedOperator(B, D)
    end

    time_norms = @elapsed norms = powernormbounds(B, D; Q=Q)
    time_eigen = @elapsed w = invariant_vector(B, Q)
    time_error = @elapsed error = distance_from_invariant(B, D, Q, w, norms)

    times = [time_error time_eigen time_norms time_assembling 0]

    return times, error, (B, D, norms, time_assembling+time_norms)
end

function twogrid(Btype, size, (B, D, norms, time_coarse))
    time_assembling_fine = @elapsed begin
        B_fine = Btype(size)
        Q_fine = DiscretizedOperator(B_fine, D)
    end

    time_norms_fine = @elapsed norms_fine = finepowernormbounds(B, B_fine, D, norms; Q_fine=Q_fine)
    time_eigen_fine = @elapsed w_fine = invariant_vector(B_fine, Q_fine)
    time_error_fine = @elapsed error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)

    times_fine = [time_error_fine time_eigen_fine time_norms_fine time_assembling_fine time_coarse]
    return times_fine, error_fine
end

function time_convergence_plot(T, Btype, k_onegrid, k_twogrid)
    n_onegrid = 2 .^ k_onegrid

    times_onegrid = fill(NaN, length(n_onegrid), 5)
    errors_onegrid = fill(NaN, size(n_onegrid))
    onegrid(T, Btype, n_onegrid[1]) #warmup
    for i in 1:length(n_onegrid)
        times_onegrid[i,:], errors_onegrid[i], _ = onegrid(T, Btype, n_onegrid[i])
    end

    n_twogrid = 2 .^ k_twogrid
    _, _, coarse_data = onegrid(T, Btype, 1024)
    times_twogrid = fill(NaN, length(n_twogrid), 5)
    errors_twogrid = fill(NaN, size(n_twogrid))

    twogrid(Btype, n_twogrid[1], coarse_data) #warmup
    for i in 1:length(n_twogrid)
        times_twogrid[i,:], errors_twogrid[i] = twogrid(Btype, n_twogrid[i], coarse_data)
    end

    p1 = groupedbar(
        times_onegrid,
        bar_position = :stack,
        legend = false,
        label = ["err" "eigen" "norm" "matrix" "coarse"],
        title = "CPU time breakdown (s)",
        xticks = (1:length(n_onegrid), LaTeXString.(raw"2^{" .* string.(k_onegrid) .* raw"}")),
        link = :y,
    )

    p2 = groupedbar(
        times_twogrid,
        bar_position = :stack,
        legend = false,
        label = ["err" "eigen" "norm" "matrix" "coarse"],
        title = "CPU time breakdown (s)",
        xticks = (1:length(n_twogrid), LaTeXString.(raw"2^{" .* string.(k_twogrid) .* raw"}")),
        link = :y,
    )

    p3 = plot(
        n_onegrid,
        errors_onegrid,
        title = "Error",
        mark = :dot,
        yscale = :log10,
        xscale = :log10,
        xticks = (n_onegrid, LaTeXString.(raw"2^{" .* string.(k_onegrid) .* raw"}" )),
        label = "One-grid strategy",
        legend = :bottomleft,
        link = :y,
    )


    p4 = plot(
        n_twogrid,
        errors_twogrid,
        title = "Error",
        mark = :dot,
        yscale = :log10,
        xscale = :log10,
        color = :red,
        xticks = (n_twogrid, LaTeXString.(raw"2^{" .* string.(k_twogrid) .* raw"}")),
        label = "Two-grid strategy",
        legend = :bottomleft,
        link = :y,
    )

    p = plot(p1, p2, p3, p4)

end

p = time_convergence_plot(x->2x+0.5x*(1-x), Ulam, 8:11, 11:14)

# large-scale version:
# p = time_convergence_plot(x->2x+0.5x*(1-x), Ulam, 8:14, 13:26); savefig(p, "time_convergence_plot.tikz"); savefig(p, "time_convergence_plot.pdf");


# savefig(p, "time_convergence_plot.tikz")
