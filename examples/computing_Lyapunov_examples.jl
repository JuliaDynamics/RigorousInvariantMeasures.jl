using RigorousInvariantMeasures, IntervalArithmetic

Lanford = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
NonMarkov = PwMap(
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
Perturbed4x = mod1_dynamic(x -> 4 * x + 0.01 * RigorousInvariantMeasures.sinpi(8 * x))

Dynamics = [(Lanford, "Lanford"); (Perturbed4x, "Perturbed4x")]

function compute_lyapunov(D, size_coarse, size_fine)
    B = Ulam(size_coarse)
    Q = DiscretizedOperator(B, D)
    norms = powernormbounds(B, D, Q = Q)

    B_fine = Ulam(size_fine)
    Q_fine = DiscretizedOperator(B_fine, D)

    normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
    norms_fine = finepowernormbounds(B, B_fine, D, norms; normQ_fine = normQ_fine)
    w_fine = invariant_vector(B_fine, Q_fine)
    error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)

    ϕ = discretizationlogder(B_fine, D; degree = 7)
    lyap = integrateobservable(B_fine, ϕ, w_fine, error_fine)
    return lyap
end

using JLD

coarse_size = 2^15
fine_size = 2^30
file = jldopen("computed_lyapunov_$(coarse_size)_$(fine_size).jld", "w")

for D in Dynamics
    lyap = compute_lyapunov(D[1], coarse_size, fine_size)
    write(file, D[2], lyap)
end
