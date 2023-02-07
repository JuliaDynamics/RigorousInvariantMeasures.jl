using RigorousInvariantMeasures: preimage_monotonic
using StaticArrays
using IntervalArithmetic

@testset "Contractors" begin

@test preimage_monotonic(2, x -> x^2, 1..2, (0,4); ϵ =  1e-13, max_iter = 100) ≈ sqrt(2)
# to ensure quadratic convergence
@test preimage_monotonic(2, x -> x^2, 1..2, (0,4); ϵ = 1e-13, max_iter = 6) ≈ sqrt(2)

end #testset
