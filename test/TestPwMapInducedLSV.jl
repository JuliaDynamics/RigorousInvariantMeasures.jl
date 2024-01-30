@testset "TestPwMapInducedLSV" begin
    using RigorousInvariantMeasures
    T = Float64

    D = RigorousInvariantMeasures.PwDynamicApproxInducedLSV(0.0, 10; T = T)

    for (i, br) in enumerate(reverse(branches(D))[1:end-1])
        @test br.X == (Interval{T}(0.5+0.5^(i+1)), Interval{T}(0.5+0.5^(i)))
    end

    test_br = reverse(branches(D))

    x = (test_br[2].X[1]+test_br[2].X[2])/2
    @test 4*x ⊆ test_br[2].f(x)

    x = (test_br[3].X[1]+test_br[3].X[2])/2
    @test 8*x ⊆ test_br[3].f(x)

    br = test_br[end]
    x = (br.X[1]+br.X[2])/2
    @test 0.75 ∈ br.f(x)

    ϕ = RigorousInvariantMeasures.CoordinateChangePwMap()
    ψ = RigorousInvariantMeasures.InvCoordinateChangePwMap()

    @test ϕ.branches[1].f(0.5) == 0.0
    @test ϕ.branches[1].f(1.0) == 1.0
    @test ψ.branches[1].f(0.0) == 0.5
    @test ψ.branches[1].f(1.0) == 1.0




    D = RigorousInvariantMeasures.RescaledApproxInducedLSV(0.0, 10; T = Float64)

    for (i, br) in enumerate(reverse(branches(D.E))[1:end-1])
        @test br.X == (Interval{T}(0.0+0.5^(i)), Interval{T}(0.0+0.5^(i-1)))
    end


end