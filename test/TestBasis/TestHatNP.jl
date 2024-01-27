@testset "Non Periodic Hat basis" begin
    using IntervalArithmetic

    using RigorousInvariantMeasures
    using RigorousInvariantMeasures:
        HatFunction, HatFunctionOnTorus, IntervalOnTorus, nonzero_on, is_refinement

    f = HatFunction(1.0, 2, 3)
    @test f(1.5) == 0.5
    @test f(1 .. 1.5) == 0 .. 0.5

    f = HatFunction(0.125, 0.25, 0.375)
    x = Interval(0.1875 .. 0.25)
    @test f(x) == 0.5 .. 1

    f = HatFunction(0, 0.125, 0.25)
    x = Interval(0 .. 0.0625)
    @test f(x) == 0 .. 0.5

    B = HatNP(4)

    @test length(B) == 5
    @test nonzero_on(B, (0.1 .. 0.3, NaN)) == (1, 3)
    @test nonzero_on(B, (0 .. 1, NaN)) == (1, 5)
    @test nonzero_on(B, (0.3 .. 0.31, NaN)) == (2, 3)
    @test nonzero_on(B, (0.9 .. 0.91, NaN)) == (4, 5)

    @test is_refinement(HatNP(8), HatNP(4))
    @test is_refinement(HatNP(8), HatNP(8))
    @test !is_refinement(HatNP(4), HatNP(8))

    ϕ = B[1]
    @test ϕ(0.0) == 1.0
    @test ϕ(0.25) == 0.0
    @test ϕ(-0.25) == 0.0

    ϕ = B[5]
    @test ϕ(1.0) == 1.0
    @test ϕ(0.75) == 0.0
    @test ϕ(1.25) == 0.0

    ϕ = B[3]
    @test ϕ(0.5) == 1.0
    @test ϕ(0.25) == 0.0
    @test ϕ(0.75) == 0.0

    @test integral_covector(B) == 1 / 5 * ones(Interval{Float64}, 5)'
    @test one_vector(B) == ones(5)

    D = mod1_dynamic(x -> 2 * x)

    B = HatNP(4)

    @test RigorousInvariantMeasures.strong_norm(B) == Lipschitz
    @test RigorousInvariantMeasures.weak_norm(B) == Linf
    @test RigorousInvariantMeasures.aux_norm(B) == L1

    BU = Ulam(4)
    BH = HatNP(4)
    v = ones(4)
    @test RigorousInvariantMeasures.change_of_basis(BU, BH, v) == ones(5)


    @test RigorousInvariantMeasures.evaluate_integral(B, 1, Float64) == 0.25

    Q = DiscretizedOperator(B, D)

    @test Q.L[:, 1] == [0.5; 0.25; 0; 0; 0]
    @test Q.L[:, 2] == [0.0; 0.25; 0.5; 0.25; 0]
    @test Q.L[:, 3] == [0.5; 0.25; 0.0; 0.25; 0.5]
    @test Q.L[:, 4] == [0; 0.25; 0.5; 0.25; 0.0]
    @test Q.L[:, 5] == [0; 0; 0; 0.25; 0.5]

    B = HatNP(2^10)
    Q = DiscretizedOperator(B, D)

    norms = powernormbounds(B, D, Q = Q)
    @test norms[10] < 1.0
    @test all(norms[1:9] .>= 1.0)

    @test RigorousInvariantMeasures.invariant_measure_strong_norm_bound(B, D) == 0.0
end
