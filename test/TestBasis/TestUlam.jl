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
    @test RigorousInvariantMeasures.aux_normalized_projection_error(Ulam(4)) == 0.0

    B8 = Ulam(8)

    @test one_vector(B8) == ones(8)

    @test RigorousInvariantMeasures.evaluate(B8, 1, 0.0625) == 1.0
    @test RigorousInvariantMeasures.evaluate(B8, 1, 1.0) == 0.0
    @test RigorousInvariantMeasures.evaluate_integral(B8, 1, Float64) == Float64(1) / 8

    @test RigorousInvariantMeasures.strong_weak_bound(B8) == Float64(8)
    @test RigorousInvariantMeasures.aux_weak_bound(B8) == 1.0
    @test RigorousInvariantMeasures.weak_by_strong_and_aux_bound(B8) == (0.0, 1.0)
    @test RigorousInvariantMeasures.bound_weak_norm_from_linalg_norm(B8) == (1.0, 0.0)
    @test RigorousInvariantMeasures.bound_linalg_norm_L1_from_weak(B8) == 1.0
    @test RigorousInvariantMeasures.bound_linalg_norm_L∞_from_weak(B8) == Float64(8)

    D = mod1_dynamic(x -> 2 * x)
    UD = Dual(B8, D; ϵ = 0.00000000001, max_iter = 100)
    @test length(UD) == 16
    @test eltype(UD) == Tuple{Int64,Tuple{Interval,Interval}}

    (i, dual_element) = iterate(UD)

end
