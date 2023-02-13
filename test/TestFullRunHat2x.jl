@testset "Full run Hat 2x mod 1" begin

using RigorousInvariantMeasures
D = mod1_dynamic(x -> 2*x)

B = Hat(1024)
Q = DiscretizedOperator(B, D)
@test 0.5 ∈  Q.L[1,1]
@test 0.25 ∈  Q.L[2,1]

norms = powernormbounds(B, D, Q=Q)
@test norms[10] < 1e-5

import LinearAlgebra: norm
w = invariant_vector(B, Q)
@test norm(w .- 1.0, Inf) <  1e-10

error = distance_from_invariant(B, D, Q, w, norms)
@test error <1e-10

end