import Pkg
Pkg.activate("../")

using RigorousInvariantMeasures, IntervalArithmetic

Lanford = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
Perturbed4x = mod1_dynamic(x -> 4 * x + 0.01 * RigorousInvariantMeasures.sinpi(8 * x))
NonLinearNonMarkov = PwMap(
    [
        x -> 17 * x / 5,
        x -> (34 * ((17 * x - 5) / 17) / 25 + 3) * ((17 * x - 5) / 17),
        x -> (34 * ((17 * x - 10) / 17) / 25 + 3) * ((17 * x - 10) / 17),
        x -> 17 * ((17 * x - 15) / 17) / 5,
    ],
    [Interval(0), Interval(5) / 17, Interval(10) / 17, Interval(15) / 17, Interval(1)],
    [
        Interval(0) Interval(1)
        Interval(0) Interval(1)
        Interval(0) Interval(1)
        Interval(0) @interval(0.4)
    ],
)

Dynamics = [
    (Lanford, "Lanford"),
    (Perturbed4x, "Perturbed4x"),
    (NonLinearNonMarkov, "NonLinearNonMarkov"),
]

size_coarse = 2^15
size_fine = 2^25

@info "Size coarse: $(size_coarse) Size fine: $(size_fine)"

using JLD

file = jldopen("output_$(size_coarse)_$(size_fine).jld", "w")

for Dyn in Dynamics
    D = Dyn[1]
    @info "Computing $(Dyn[2])"
    total_time = @elapsed begin
        B = Ulam(size_coarse)
        Q = DiscretizedOperator(B, D)
        norms = powernormbounds(B, D; Q = Q)

        B_fine = Ulam(size_fine)
        Q_fine = DiscretizedOperator(B_fine, D)

        normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
        norms_fine = finepowernormbounds(B, B_fine, D, norms; normQ_fine = normQ_fine)
        w_fine = invariant_vector(B_fine, Q_fine)
        error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)
        ϕ = discretizationlogder(B_fine, D, degree = 3)
        lyap = integrateobservable(B_fine, ϕ, w_fine, error_fine)
    end
    @info "Computed $lyap, with precision $(diam(lyap)), in $total_time"
    write(file, Dyn[2], (lyap, total_time))
end

close(file)
