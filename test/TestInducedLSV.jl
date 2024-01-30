@testset "Test InducedLSV" begin
    using RigorousInvariantMeasures


    @test RigorousInvariantMeasures.CoordinateChange(0.0) == -1
    @test RigorousInvariantMeasures.CoordinateChange(1.0) == 1.0
    @test RigorousInvariantMeasures.InvCoordinateChange(-1.0) == 0.0 
    @test RigorousInvariantMeasures.InvCoordinateChange(1.0) == 1.0

    x = rand()
    y = RigorousInvariantMeasures.CoordinateChange(x)
    @test RigorousInvariantMeasures.InvCoordinateChange(y) == x

    D = ApproxInducedLSV(1/8, 10)

    @test D.α == 1/8

    mid = [0.5^(20-i+2) for i in 2:20]
    test = RigorousInvariantMeasures.ShootingLSV(20, 0.5, 0.0)

    @test all(mid .∈ test[2:end])

    D = ApproxInducedLSV(0.0, 10)

    @test RigorousInvariantMeasures.derleft(D, 0.0) == 2.0
    
    domains = RigorousInvariantMeasures.GetDomains(10, 0.0)
    @test domains[end] == Interval(0.5+(1/2^2), 0.5+(1/2)^1)
    @test domains[end-1] == Interval(0.5+(1/2)^3, 0.5+(1/2)^2)
    @test domains[end-2] == Interval(0.5+(1/2)^4, 0.5+(1/2)^3)


end