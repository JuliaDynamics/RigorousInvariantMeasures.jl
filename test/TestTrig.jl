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

    @test RigorousInvariantMeasures.sinpi(Interval(0, 2)) == Interval(-1, 1)
    @test RigorousInvariantMeasures.sinpi(Interval(0, 1)) == Interval(0, 1)
    @test RigorousInvariantMeasures.sinpi(Interval(1, 2)) == Interval(-1, 0)

    @test RigorousInvariantMeasures.cospi(Interval(0, 2)) == Interval(-1, 1)
    @test RigorousInvariantMeasures.cospi(Interval(0, 1)) == Interval(-1, 1)
    @test RigorousInvariantMeasures.cospi(Interval(0, 0.5)) == Interval(0, 1)
    @test RigorousInvariantMeasures.cospi(Interval(0.5, 1)) == Interval(-1, 0)

    x = Taylor1([0.0, 1.0], 1)
    y = Taylor1([0.0, π], 1)

    @test RigorousInvariantMeasures.sinpi(x) == y

    y = Taylor1([1.0, 0,0], 1)
    @test RigorousInvariantMeasures.cospi(x) == y

end