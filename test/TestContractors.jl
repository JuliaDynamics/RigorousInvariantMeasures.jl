using RigorousInvariantMeasures: preimage_monotonic
using StaticArrays
using IntervalArithmetic

@testset "Contractors" begin

    @test in_interval(
        sqrt(2),
        preimage_monotonic(2, x -> x^2, interval(1, 2), (0, 4); ϵ = 1e-13, max_iter = 100),
    )
    # to ensure quadratic convergence
    @test in_interval(
        sqrt(2),
        preimage_monotonic(2, x -> x^2, interval(1, 2), (0, 4); ϵ = 1e-13, max_iter = 6),
    )

end #testset
