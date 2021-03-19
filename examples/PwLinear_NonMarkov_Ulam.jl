using InvariantMeasures
using ValidatedNumerics

using Plots
using LaTeXStrings
using StatsPlots

function runExperiment()

    # number of norms to compute and extend
    m = 15
    m_extend = 100

    time_assembling = @elapsed begin

		D = PwMap([x->2.5*x, x->4*x-1, x->4*x-2, x-> 4*x-3],
			[Interval(0), Interval(0.25), Interval(0.5), Interval(0.75), Interval(1)])
        B = Ulam(1024)
        Q = DiscretizedOperator(B, D)
    end

    time_norms = @elapsed norms = powernormbounds(B, D, m, m_extend; Q=Q)
    time_eigen = @elapsed w = invariant_vector(B, Q)
    time_error = @elapsed error = distance_from_invariant(B, D, Q, w, norms)

    time_assembling_fine = @elapsed begin
        B_fine = Ulam(2^16)
        Q_fine = DiscretizedOperator(B_fine, D)
    end

    time_norms_fine = @elapsed norms_fine = finepowernormbounds(B, B_fine, D, norms; Q_fine=Q_fine)
    time_eigen_fine = @elapsed w_fine = invariant_vector(B_fine, Q_fine)
    time_error_fine = @elapsed error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)

    p1 = plot(D, title="Dynamic", label=L"T(x)", legend=:bottomright)
    p2 = plot(B, w, title="Invariant measure (n=$(length(B)))")
    p2 = plot!(p2, B, error, w, label="L1 error $(round(error, sigdigits=2))")

    p3 = plot(B_fine, w_fine, title="Invariant measure (n=$(length(B_fine)))")
    p3 = plot!(p3, B_fine, error_fine, w_fine, label="L1 error $(round(error_fine, sigdigits=2))")

    p4 = groupedbar(
        vcat(
            [time_error time_eigen time_norms time_assembling 0],
            [time_error_fine time_eigen_fine time_norms_fine time_assembling_fine time_assembling+time_norms]
        ),
        bar_position = :stack,
        legend = :topleft,
        label = ["err" "eigen" "norm" "matrix" "coarse"],
        title = "CPU time breakdown",
        xticks = (1:2, ["1-grid estimate", "2-grid estimate"])
    )

    plot(p1, p2, p3, p4, layout=4)

end

runExperiment()
