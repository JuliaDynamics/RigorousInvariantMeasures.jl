@testset "Ulam basis" begin


    using RigorousInvariantMeasures
    using IntervalArithmetic

    B = Ulam(4)
    @test B.p == LinRange(0.0, 1.0, 5)
    @test length(B) == 4

    @test B[1](0) == 1
    @test B[4](0) == 0

    @test nonzero_on(B, (@interval(0.1), @interval(0.3))) == (1, 2)
    @test nonzero_on(B, (@interval(0), @interval(1))) == (1, 4)
    @test nonzero_on(B, (@interval(0.3), @interval(0.31))) == (2, 2)
    @test nonzero_on(B, (@interval(1), @interval(1))) == (4, 4)

    @test RigorousInvariantMeasures.relative_measure(
        (@interval(0.5), @interval(1)),
        (@interval(0), @interval(1)),
    ) == 0.5

    @test is_refinement(Ulam(8), Ulam(4))
    @test is_refinement(Ulam(8), Ulam(8))
    @test !is_refinement(Ulam(4), Ulam(8))
    @test length(Ulam(4)) == 4

end
