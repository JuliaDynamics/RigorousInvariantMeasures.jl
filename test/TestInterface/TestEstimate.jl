@testset "Estimator" begin
    using IntervalArithmetic


    function twoxexperiment(n)
        D = mod1_dynamic(x -> 2 * x)
        B = Ulam(32)
        Q = DiscretizedOperator(B, D)
        return B, D, Q
    end

    CGQ, FGQ = RigorousInvariantMeasures.compute_coarse_grid_quantities(twoxexperiment, 64)

    Q = DiscretizedOperator(CGQ.B, CGQ.D)
    norms = powernormbounds(CGQ.B, CGQ.D; Q = Q)

    @test norms == CGQ.norms

    norms_other = RigorousInvariantMeasures.powernormbounds(CGQ.B, CGQ.D, 8, 16; Q = Q)
    @test norms_other == norms

    error, times = one_grid_estimate(CGQ, FGQ)
    @test error == 0.0

    error, times = two_grid_estimate(CGQ, FGQ)
    @test error == 0.0

    # C = RigorousInvariantMeasures.boundnorm(B, P, 10)

    # @test C == [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

end
