@testset "Ulam basis, 2D, Skew Product" begin

    using RigorousInvariantMeasures
    using RigorousInvariantMeasures.BasisDefinition
    using IntervalArithmetic

    B = Ulam2DSP(4)
    @test B.part_x == LinRange(0.0, 1.0, 5)
    @test B.part_y == LinRange(0.0, 1.0, 5)

    @test length(B) == 16
    @test B[1, 1](0.0, 0.0) == 1
    @test B[1, 1](1.0, 0.0) == 0
    @test B[1, 1](0.0, 1.0) == 0
    @test B[1, 1](1.0, 1.0) == 0

    @test B[4, 4](0.75, 0.75) == 1.0
    @test B[4, 4](1.0, 1.0) == 0.0


    B = Ulam2DSP(4, 5)
    @test length_x(B) == 4
    @test length_y(B) == 5


    @test square_indexes_to_linear(1, 1, 2, 3) == 1
    @test square_indexes_to_linear(2, 1, 2, 3) == 2
    @test square_indexes_to_linear(1, 2, 2, 3) == 3
    @test square_indexes_to_linear(2, 2, 2, 3) == 4
    @test square_indexes_to_linear(1, 3, 2, 3) == 5
    @test square_indexes_to_linear(2, 3, 2, 3) == 6

    @test square_indexes_to_linear(B, 1, 1) == 1
    @test square_indexes_to_linear(B, 1, 2) == 5
    @test square_indexes_to_linear(B, 4, 1) == 4
    @test square_indexes_to_linear(B, 2, 3) == 10
    @test square_indexes_to_linear(B, 4, 5) == 20

    @test linear_indexes_to_square(1, 6, 4) == (1, 1)
    @test linear_indexes_to_square(6, 6, 4) == (6, 1)
    @test linear_indexes_to_square(7, 6, 4) == (1, 2)
    @test linear_indexes_to_square(12, 6, 4) == (6, 2)


    @test linear_indexes_to_square(B, 1) == (1, 1)
    @test linear_indexes_to_square(B, 6) == (2, 2)
    @test linear_indexes_to_square(B, 4) == (4, 1)
    @test linear_indexes_to_square(B, 20) == (4, 5)

    @test getindex_linear(B, 1)(0.0, 0.0) == 1.0
    @test getindex_linear(B, 1)(1.0, 1.0) == 0.0
    @test getindex_linear(B, 6)(0.0, 0.0) == 0.0
    @test getindex_linear(B, 6)(0.25, 0.2) == 1.0

    T = mod1_dynamic(x -> 2 * x)
    D = SkewProductMap(T, [(x, y) -> (y * x) / 8 + 0.25; (x, y) -> (y * x) / 8 + 0.75])

    ϵ = 10^(-14)
    max_iter = 100
    #test_Dual = RigorousInvariantMeasures.Dual(B, D; ϵ, max_iter )
    #@test test_Dual.x == RigorousInvariantMeasures.preimages(B.p_x, D.T; ϵ, max_iter)[1]
    #@test test_Dual.xlabel == RigorousInvariantMeasures.preimages(B.p_x, D.T; ϵ, max_iter)[2]
    #@test length(test_Dual) == 40

    G = D.G[1]
    @test RigorousInvariantMeasures.check_image(B, G, Interval(0.1), Interval(0.2)) ==
          (1, 2)
    B = Ulam2DSP(4, 1024)
    @test RigorousInvariantMeasures.check_image(B, G, Interval(0.1), Interval(0.2)) ==
          (256, 282)


    @test RigorousInvariantMeasures.preimage_fixed_x(D, 1, 0, 0.24, 0.26; ϵ, max_iter) ==
          (Interval(0), Interval(1))
    @test RigorousInvariantMeasures.preimage_fixed_x(D, 1, 0, 0.24, 0.241; ϵ, max_iter) ==
          (Interval(∅), Interval(∅))
    @test RigorousInvariantMeasures.preimage_fixed_x(
        D,
        1,
        0.125,
        0.25,
        0.275;
        ϵ,
        max_iter,
    ) == (Interval(0), Interval(1))
    @test RigorousInvariantMeasures.preimage_fixed_x(D, 1, 0.125, 0.5, 0.6; ϵ, max_iter) ==
          (Interval(∅), Interval(∅))
    @test RigorousInvariantMeasures.preimage_fixed_x(
        D,
        1,
        0.125,
        0.23,
        0.275;
        ϵ,
        max_iter,
    ) == (Interval(0), Interval(1))


    P, err = RigorousInvariantMeasures.rectangle_preimage(
        D,
        1,
        0.1,
        0.2,
        0.25,
        1,
        2;
        ϵ,
        max_iter,
    )

    import Polyhedra as PH
    import GLPK
    lib = PH.DefaultLibrary{Float64}(GLPK.Optimizer)

    A = [
        0.1 9.15934e-16
        0.15 9.4369e-16
        0.2 9.15934e-16
        0.2 1.0
        0.15 1.0
        0.1 1.0
    ]
    @test sum(abs.(P.vrep.V - A)) <= 1e-10

end
