using RigorousInvariantMeasures
using IntervalArithmetic

using Plots
using LaTeXStrings
using StatsPlots

ENV["GKSwstype"]="nul" # for headless displays

LorenzMap(θ, α) = PwMap([x->θ*(0.5-x)^α, x->1-θ*(x-0.5)^α],
                    [x->-1*θ*α*(0.5-x)^(α-1), x->-θ*α*(x-0.5)^(α-1)],
                    [@interval(0), @interval(0.5), @interval(1)],
                    [θ*(Interval(0.5))^α Interval(0.0);
                    Interval(1.0)  1-θ*(Interval(0.5))^α]; infinite_derivative=true)

include("warmup.jl")

function runExperiment()

    time_assembling = @elapsed begin
        D0 = LorenzMap(109/64, 51/64)
        D = D0 ∘ D0 ∘ D0
        B = Ulam(2^15)
        Q = DiscretizedOperator(B, D)
    end
    time_dfly = @elapsed dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D)
    time_norms = @elapsed norms = powernormbounds(B, D, Q=Q, m=20)
    time_eigen = @elapsed w = invariant_vector(B, Q)
    time_error = @elapsed error = distance_from_invariant(B, D, Q, w, norms; dfly_coefficients=dfly_coefficients)
    time_assembling_fine = @elapsed begin
        B_fine = Ulam(2^18)
        Q_fine = DiscretizedOperator(B_fine, D)
    end

    time_norms_fine = @elapsed begin
        normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
        norms2 = refine_norms_of_powers(norms,400)
        norms_fine = finepowernormbounds(B, B_fine, D, norms2; normQ_fine=normQ_fine, dfly_coefficients=dfly_coefficients)
    end
    time_eigen_fine = @elapsed w_fine = invariant_vector(B_fine, Q_fine)
    time_error_fine = @elapsed error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine; dfly_coefficients=dfly_coefficients)

    A, BB = dfly_coefficients
    p1 = plot(D.E, title="Dynamic (dfly coeffs $(round(A, sigdigits=2)),$(round(BB, sigdigits=2)))", label=L"T(x)", legend=:bottomright)
    p2 = plot(B, w, title="Invariant measure (n=$(length(B)))")
    p2 = plot!(p2, B, error, w, label="L1 error $(round(error, sigdigits=2))")

    p3 = plot(B_fine, w_fine, title="Invariant measure (n=$(length(B_fine)))")
    p3 = plot!(p3, B_fine, error_fine, w_fine, label="L1 error $(round(error_fine, sigdigits=2))")

    p4 = groupedbar(
        vcat(
            [time_dfly time_error time_eigen time_norms time_assembling 0],
            [time_dfly time_error_fine time_eigen_fine time_norms_fine time_assembling_fine time_assembling+time_norms]
        ),
        bar_position = :stack,
        legend = :topleft,
        label = ["dfly" "err" "eigen" "norms" "matrix" "coarse"],
        title = "CPU time breakdown (s)",
        xticks = (1:2, ["1-grid estimate", "2-grid estimate"])
    )

    plot(p1, p2, p3, p4, layout=4)
    savefig("Lorenz.png")
end

runExperiment()
