@testset "Test InducedLSV" begin
    using RigorousInvariantMeasures


    @test RigorousInvariantMeasures.CoordinateChange(0.0) == -1
    @test RigorousInvariantMeasures.CoordinateChange(1.0) == 1.0
    @test RigorousInvariantMeasures.InvCoordinateChange(-1.0) == 0.0
    @test RigorousInvariantMeasures.InvCoordinateChange(1.0) == 1.0

    x = rand()
    y = RigorousInvariantMeasures.CoordinateChange(x)
    @test RigorousInvariantMeasures.InvCoordinateChange(y) == x

    D = ApproxInducedLSV(1 / 8, 10)

    @test D.α == 1 / 8

    mid = [0.5^(20 - i + 2) for i = 2:20]
    test = RigorousInvariantMeasures.ShootingLSV(20, 0.5, 0.0)

    @test in_interval.(all(mid, test[2:end]))

    D = ApproxInducedLSV(0.0, 10)

    @test RigorousInvariantMeasures.derleft(D, 0.0) == 2.0

    domains = RigorousInvariantMeasures.GetDomains(10, 0.0)
    @test isequal_interval(domains[end], interval(0.5 + (1 / 2^2), 0.5 + (1 / 2)^1))
    @test isequal_interval(domains[end-1], interval(0.5 + (1 / 2)^3, 0.5 + (1 / 2)^2))
    @test isequal_interval(domains[end-2], interval(0.5 + (1 / 2)^4, 0.5 + (1 / 2)^3))


end
