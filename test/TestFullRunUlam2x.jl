@testset "Full run Ulam 2x mod 1" begin

using RigorousInvariantMeasures
D = mod1_dynamic(x -> 2*x)
B = Ulam(1024)
Q = DiscretizedOperator(B, D)

norms = powernormbounds(B, D, Q=Q)

@test norms[10] < 1e-10

w = invariant_vector(B, Q)

import LinearAlgebra: norm
@test norm(w .- 1.0, Inf) <  1e-10

error = distance_from_invariant(B, D, Q, w, norms)

logder = discretizationlogder(B, D)

lyap = integrateobservable(B, logder, w, error)

@test log(2) ∈ lyap

end