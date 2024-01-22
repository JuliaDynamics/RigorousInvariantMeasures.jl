@testset "Lorenz2dUlam" begin
    using RigorousInvariantMeasures

    Ginv = RigorousInvariantMeasures.Lorenz2D.Ginverse(x = (2.0)^-1, r = 5.0, c = 0.5)

    @test Ginv(0.5) == 0.0
    @test Ginv(0.5 + 1 / 1024) == 1
    @test Ginv(0.5 - 1 / 1024) == -1

    Ginv = RigorousInvariantMeasures.Lorenz2D.Ginverse(x = -(2.0)^-1, r = 5.0, c = 0.5)
    @test Ginv(-0.5) == 0.0
    @test Ginv(-0.5 + 1 / 1024) == 1
    @test Ginv(-0.5 - 1 / 1024) == -1

    using IntervalArithmetic
    @test 1 ∈ RigorousInvariantMeasures.Lorenz2D._Lorenz_one_dim_map(Interval(1), 2.0, 3.0)
    @test -1 ∈ RigorousInvariantMeasures.Lorenz2D._Lorenz_right_one_dim_map(
        Interval(0),
        2.0,
        3.0,
    )
    @test 1 ∈
          RigorousInvariantMeasures.Lorenz2D._Lorenz_left_one_dim_map(Interval(0), 2.0, 3.0)

end
