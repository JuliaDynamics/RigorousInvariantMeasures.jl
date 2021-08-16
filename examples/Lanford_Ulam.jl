using InvariantMeasures
using ValidatedNumerics

<<<<<<< HEAD
using Plots
using LaTeXStrings
using StatsPlots

include("warmup.jl")

function runExperiment()

    time_assembling = @elapsed begin
        # Note that (unlike the experiment in [Galatolo, Nisoli] paper) we do not need
        # to take Iterate(D, 2) here

        D = Mod1Dynamic(x->2*x+0.5*x*(1-x))
        # different backend, a tad slower
        # D = mod1_dynamic(x -> 2*x+0.5*x*(1-x))
        B = Ulam(1024)
        Q = DiscretizedOperator(B, D)
    end

    time_norms = @elapsed norms = powernormbounds(B, D; Q=Q)
    time_eigen = @elapsed w = invariant_vector(B, Q)
    time_error = @elapsed error = distance_from_invariant(B, D, Q, w, norms)

    time_assembling_fine = @elapsed begin
        B_fine = Ulam(2^16)
        Q_fine = DiscretizedOperator(B_fine, D)
    end

    time_norms_fine = @elapsed norms_fine = finepowernormbounds(B, B_fine, D, norms; Q_fine=Q_fine)
    time_eigen_fine = @elapsed w_fine = invariant_vector(B_fine, Q_fine)
    time_error_fine = @elapsed error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)

    A, BB = dfly(strong_norm(B), aux_norm(B), D)
    p1 = plot(D, title="Dynamic (dfly coeffs $(round(A, sigdigits=2)),$(round(BB, sigdigits=2)))", label=L"T(x)", legend=:bottomright)
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
        title = "CPU time breakdown (s)",
        xticks = (1:2, ["1-grid estimate", "2-grid estimate"])
    )

    plot(p1, p2, p3, p4, layout=4)

end

runExperiment()
=======
m = 30
m_extend = 100

D = Mod1Dynamic(x->2*x+0.5*x*(1-x))
B = Ulam(1024)
Q = DiscretizedOperator(B, D)

normQ = opnormbound(weak_norm(B), Q)

trivial_norms = norms_of_powers_trivial(weak_norm(B), Q, m)
computed_norms = norms_of_powers(weak_norm(B), m, Q, integral_covector(B))

(dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, m)

norms = min.(trivial_norms, computed_norms, dfly_norms) # in the current version, dfly_norms are always larger and can be omitted

better_norms = refine_norms_of_powers(norms, m_extend)

w = invariant_vector(B, Q)
@show distance_from_invariant(B, D, Q, w, better_norms)

B_fine = Ulam(2^20)
Q_fine = DiscretizedOperator(B_fine, D)
norm_Q_fine = opnormbound(weak_norm(B_fine), Q_fine)

trivial_norms_fine = norms_of_powers_trivial(weak_norm(B_fine), Q_fine, m_extend)
twogrid_norms_fine = norms_of_powers_from_coarser_grid(B_fine, B, D, better_norms, norm_Q_fine)

(dfly_strongs_fine, dfly_norms_fine) = norms_of_powers_dfly(B_fine, D, m_extend)

norms_fine = min.(trivial_norms_fine, twogrid_norms_fine, dfly_norms_fine)

better_norms_fine = refine_norms_of_powers(norms_fine, m_extend)

w_fine = invariant_vector(B_fine, Q_fine)
@show distance_from_invariant(B_fine, D, Q_fine, w_fine, better_norms_fine)
>>>>>>> origin/noise_contract
