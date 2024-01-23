@testset "Estimator" begin
    using IntervalArithmetic

    D = mod1_dynamic(x -> 2 * x)
    B = Ulam(32)
    P = RigorousInvariantMeasures.assemble(B, D; ϵ = 10^(-15), max_iter = 100, T = Float64)

    # C = RigorousInvariantMeasures.boundnorm(B, P, 10)

    # @test C == [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

end
