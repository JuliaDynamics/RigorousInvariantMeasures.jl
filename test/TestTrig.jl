@testset "Trig" begin

    using RigorousInvariantMeasures
    using IntervalArithmetic

    @test RigorousInvariantMeasures.sinpi(0) == 0
    @test RigorousInvariantMeasures.sinpi(1) == 0
    @test RigorousInvariantMeasures.sinpi(2) == 0

    @test RigorousInvariantMeasures.cospi(0) == 1
    @test RigorousInvariantMeasures.cospi(1) == -1
    @test RigorousInvariantMeasures.cospi(2) == 1

    @test RigorousInvariantMeasures.find_quadrantspi(0.0) == [0.0, 0.0]
    @test RigorousInvariantMeasures.find_quadrantspi(0.7) == [1.0, 1.0]
    @test RigorousInvariantMeasures.find_quadrantspi(1.1) == [2.0, 2.0]
    @test RigorousInvariantMeasures.find_quadrantspi(1.7) == [3.0, 3.0]
    @test RigorousInvariantMeasures.find_quadrantspi(2.0) == [4.0, 4.0]

    @test RigorousInvariantMeasures.sinpi(interval(0, 2)) == interval(-1, 1)
    @test RigorousInvariantMeasures.sinpi(interval(0, 1)) == interval(0, 1)
    @test RigorousInvariantMeasures.sinpi(interval(1, 2)) == interval(-1, 0)

    @test RigorousInvariantMeasures.cospi(interval(0, 2)) == interval(-1, 1)
    @test RigorousInvariantMeasures.cospi(interval(0, 1)) == interval(-1, 1)
    @test RigorousInvariantMeasures.cospi(interval(0, 0.5)) == interval(0, 1)
    @test RigorousInvariantMeasures.cospi(interval(0.5, 1)) == interval(-1, 0)

    using TaylorSeries

    x = Taylor1([0.0, 1.0], 1)
    y = Taylor1([0.0, π], 1)

    @test RigorousInvariantMeasures.sinpi(x) == y

    y = Taylor1([1.0, 0, 0], 1)
    @test RigorousInvariantMeasures.cospi(x) == y

end
