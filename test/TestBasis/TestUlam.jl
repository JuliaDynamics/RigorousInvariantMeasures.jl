using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Ulam basis" begin

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

    @test RigorousInvariantMeasures.strong_weak_bound(B8) == Float64(16)

    # The constant is sharp, and is attained on a vector of average zero: two cells
    # sharing no boundary carry +1 and -1, so each of the four jumps is counted twice.
    let v = [0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0]
        M1n = RigorousInvariantMeasures.strong_weak_bound(B8)
        L1 = RigorousInvariantMeasures.normbound(B8, RigorousInvariantMeasures.L1, v)
        @test sum(v) == 0
        @test L1 == 2 / 8
        @test RigorousInvariantMeasures.normbound(
            B8, RigorousInvariantMeasures.VarInterval, v) == 4.0
        @test RigorousInvariantMeasures.normbound(
            B8, RigorousInvariantMeasures.VarCircle, v) == 4.0
        @test M1n * L1 == 4.0
    end
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

    # the two variations: VarInterval leaves the jump across 0 out, VarCircle counts it, and
    # both are exact for a piecewise constant function, so the bound equals the sum
    let B = Ulam(8)
        v = [1.0, 1.0, 3.0, 3.0, -1.0, -1.0, 0.0, 0.5]
        int_var = 2.0 + 4.0 + 1.0 + 0.5
        @test normbound(B, VarInterval, v) == int_var
        @test normbound(B, VarCircle, v) == int_var + 0.5
        @test normbound(B, TotalVariation, v) == int_var        # the alias is the interval one
        @test strong_norm(B) == VarInterval
        for trial = 1:50
            w = randn(1 + trial)
            Bw = Ulam(length(w))
            @test normbound(Bw, VarInterval, w) ≈
                  sum(abs(w[i+1] - w[i]) for i = 1:length(w)-1)
            @test normbound(Bw, VarCircle, w) ≈
                  normbound(Bw, VarInterval, w) + abs(w[1] - w[end])
        end
        @test normbound(Ulam(1), VarCircle, [2.0]) == 0.0
        @test normbound(B, VarCircle, ones(8)) == 0.0
    end

end

# A method dispatching on the strong-norm parameter; the two definitions are at top
# level, so the second is a method of the same function as the first.
which_variation(::Ulam) = :interval
which_variation(::Ulam{VarCircle}) = :circle

@testset "Ulam norm kinds as type parameters" begin

    B = Ulam(16)

    # the defaults are the norms the basis carried before the parameters existed
    @test B isa Ulam{VarInterval,L1}
    @test strong_norm(B) == VarInterval
    @test weak_norm(B) == L1
    @test aux_norm(B) == L1

    # and every other accessor returns what it did
    @test B.p == LinRange(0.0, 1.0, 17)
    @test length(B) == 16
    @test one_vector(B) == ones(16)
    @test integral_covector(B) == 1 / 16 * ones(Interval{Float64}, 16)'
    @test is_integral_preserving(B)
    @test is_refinement(Ulam(32), B)
    @test length(AverageZero(B)) == 15
    @test length(collect(AverageZero(B))) == 15
    @test collect(AverageZero(B))[1] ==
          [i == 1 ? 1.0 : (i == 2 ? -1.0 : 0.0) for i = 1:16]
    @test RigorousInvariantMeasures.evaluate(B, 1, 1 / 32) == 1
    @test RigorousInvariantMeasures.evaluate_integral(B, 1, Float64) == 1 / 16
    @test RigorousInvariantMeasures.weak_projection_error(B) == 0.5 / 16
    @test RigorousInvariantMeasures.aux_normalized_projection_error(B) == 0.0
    @test RigorousInvariantMeasures.strong_weak_bound(B) == 32.0
    @test RigorousInvariantMeasures.aux_weak_bound(B) == 1.0
    @test RigorousInvariantMeasures.weak_by_strong_and_aux_bound(B) == (0.0, 1.0)
    @test RigorousInvariantMeasures.bound_weak_norm_from_linalg_norm(B) == (1.0, 0.0)
    @test RigorousInvariantMeasures.bound_linalg_norm_L1_from_weak(B) == 1.0
    @test RigorousInvariantMeasures.bound_linalg_norm_L∞_from_weak(B) == 16.0
    @test RigorousInvariantMeasures.bound_weak_norm_abstract(B) == 1.0

    # the partition may be given instead of its size, and the keywords are optional
    @test Ulam(LinRange(0.0, 1.0, 17)) isa Ulam{VarInterval,L1}
    @test Ulam(LinRange(0.0, 1.0, 17)).p == B.p

    # a basis measured in the variation on the circle; only the strong parameter moves
    BC = Ulam(16; strong = VarCircle)
    @test BC isa Ulam{VarCircle,L1}
    @test strong_norm(BC) == VarCircle
    @test weak_norm(BC) == L1
    @test aux_norm(BC) == L1
    @test BC.p == B.p
    @test length(BC) == 16
    @test Ulam(16; strong = VarCircle, weak = L1) isa Ulam{VarCircle,L1}
    @test Ulam(LinRange(0.0, 1.0, 17); strong = VarCircle) isa Ulam{VarCircle,L1}

    # `normbound` bounds the vector in the kind it is handed, so it does not move with
    # the parameter: the two bases return the same number against each kind
    let v = [1.0, 1.0, 3.0, 3.0, -1.0, -1.0, 0.0, 0.5], A = randn(8, 8)
        B8, B8C = Ulam(8), Ulam(8; strong = VarCircle)
        for N in (VarInterval, VarCircle, TotalVariation, L1)
            @test normbound(B8, N, v) == normbound(B8C, N, v)
        end
        @test normbound(B8, VarInterval, v) == 7.5
        @test normbound(B8, VarCircle, v) == 8.0
        @test normbound(B8, L1, v) == normbound(L1, v)
        @test opnormbound(B8, L1, A) == opnormbound(B8C, L1, A)
    end

    # the parameter is reachable by dispatch
    @test which_variation(B) == :interval
    @test which_variation(BC) == :circle

end
