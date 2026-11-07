using RigorousInvariantMeasures
using IntervalArithmetic
using Random

@testset "DFLY on the circle (VarCircle)" begin
    RIM = RigorousInvariantMeasures
    β = (1 + sqrt(23)) / 2
    D = mod1_dynamic(x -> β * x; full_branch = false)

    @testset "the lift statement" begin
        @test has_monotone_mod1_lift(D)
        c = dfly(VarCircle, L1, D)
        @test strong_kind(c) === VarCircle
        @test weak_kind(c) === L1

        # dist = θ = 0 and one partial branch, so κ = λ_L + λ_R = 2/β and A = B = 2/β,
        # both rounded up
        @test c.A ≈ 2 / β
        @test c.B ≈ 2 / β
        @test c.A >= 2 / β
        @test c.B >= 2 / β
        @test c.B / (1 - c.A) ≈ 2 / (β - 2)

        # the interval pair for the same map pays 2/min|I_i| for a point value
        ci = dfly(VarInterval, L1, D)
        @test ci.B ≈ 2 / (1 - 2 / β)
        @test c.B < ci.B / 9
        # and doubling the interval pair, which is the one conversion that needs no
        # hypothesis, loses A < 1
        @test 2 * ci.A > 1

        # the orientation does not matter: -βx mod 1 has the same pair
        Dm = mod1_dynamic(x -> -β * x)
        @test has_monotone_mod1_lift(Dm)
        cm = dfly(VarCircle, L1, Dm)
        @test cm.A ≈ 2 / β
        @test cm.B ≈ 2 / β

        # a degree-d covering with |T'(0^+)| = |T'(1^-)| has κ = 0, hence the pair of the
        # full-branch interval statement: for x -> 2x that is (1/2, 0), and B/(1-A) = 0 is
        # the variation of Lebesgue
        c2 = dfly(VarCircle, L1, mod1_dynamic(x -> 2 * x))
        @test c2.A == 0.5
        @test c2.B == 0.0

        # same for a nonlinear one, where dist > 0 but the two ends of the lift still carry
        # the same derivative
        D4 = mod1_dynamic(x -> 4 * x + 0.01 * RIM.sinpi(8 * x))
        c4, c4i = dfly(VarCircle, L1, D4), dfly(VarInterval, L1, D4)
        @test c4.A ≈ c4i.A
        @test c4.B ≈ c4i.B

        # a lift whose derivative differs at the two ends: T̃(x) = 2.5x - x²/2 has
        # |T̃'(0)| = 2.5 and |T̃'(1)| = 1.5, no interior jump of 1/|T'| and T̃(1) - T̃(0) = 2,
        # so κ = |λ_L - λ_R|
        D3 = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x), full_branch = true)
        c3 = dfly(VarCircle, L1, D3)
        λ_L, λ_R, κ = 1 / 2.5, 1 / 1.5, abs(1 / 2.5 - 1 / 1.5)
        @test c3.A ≈ λ_R + κ / 2
        @test c3.B ≈ 1 / 1.5^2 + κ
    end

    @testset "the branchwise statement" begin
        # mixed orientations: no monotone lift, and both branches are onto the circle, so
        # A = 2λ + dist·h = 1 and B = 2 dist = 0
        tent = PwMap([x -> 2x, x -> 2 - 2x], [0, 0.5, 1])
        @test !has_monotone_mod1_lift(tent)
        ct = dfly(VarCircle, L1, tent)
        @test ct.A == 1.0
        @test ct.B == 0.0

        # the beta shift read through the statement that assumes nothing: its partial
        # branch costs 3λ, which is what keeps this pair off this map
        cb = RIM._dfly_circle_branchwise(D)
        @test cb.A ≈ 3 / β
        @test cb.B ≈ 2 / (β - 2)
        @test cb.A > 1

        # 4.5x mod 1 has a partial branch of length 1 - 4/4.5 and λ = 1/4.5 < 1/3, so both
        # statements are available there and the lift one is sharper
        D45 = mod1_dynamic(x -> 4.5 * x)
        c45, cb45 = dfly(VarCircle, L1, D45), RIM._dfly_circle_branchwise(D45)
        @test c45.A ≈ 2 / 4.5
        @test c45.B ≈ 2 / 4.5
        @test cb45.A ≈ 3 / 4.5
        @test cb45.B ≈ 2 * (1 / 4.5) / (1 - 4 / 4.5)
        @test c45.B / (1 - c45.A) < cb45.B / (1 - cb45.A)

        # the image endpoints must be exactly 0 and 1: an enclosure around them is not the
        # hypothesis, and the branchwise pair is what comes back
        off = PwMap(
            [x -> 3x + interval(0.1), x -> 3x - interval(0.9), x -> 3x - interval(1.9)],
            [0, 0.3, 1 / 3 + 0.3, 1],
        )
        @test !has_monotone_mod1_lift(off)
        @test dfly(VarCircle, L1, off).A == RIM._dfly_circle_branchwise(off).A

        # and the sharp routine refuses a map that does not satisfy its hypothesis
        @test_throws ArgumentError RIM._dfly_circle_from_lift(tent)
    end

    @testset "which of the two pairs comes back" begin
        # 4x on (0,1/2) continued by 8x on (1/2,1), cut at the integer preimages: a lift,
        # with one interior jump of 1/|T'|, at 1/2, and T(1^-) = 1 = T(0^+) on the circle,
        # so θ = |1/4-1/8| and κ = θ + |λ_L-λ_R| = 1/4
        Dj = PwMap(
            [x -> 4x, x -> 4x - 1, x -> 8x - 4, x -> 8x - 5, x -> 8x - 6, x -> 8x - 7],
            [0, 1 / 4, 1 / 2, 5 / 8, 3 / 4, 7 / 8, 1],
        )
        @test has_monotone_mod1_lift(Dj)
        cl, cb = RIM._dfly_circle_from_lift(Dj), RIM._dfly_circle_branchwise(Dj)
        @test cl.A == 0.25 + 0.125
        @test cl.B == 0.25
        # every branch is onto the circle and the map is piecewise linear, so the
        # branchwise pair has B = 2 dist = 0, which is sharp: Lebesgue is invariant
        @test cb.A == 0.5
        @test cb.B == 0.0
        # and that is the pair that comes back, B/(1-A) ranking them
        c = dfly(VarCircle, L1, Dj)
        @test (c.A, c.B) == (cb.A, cb.B)
        @test RIM.invariant_measure_strong_norm_bound(Ulam(64), Dj; dfly_coefficients = c) ==
              0.0

        # a real gap in B/(1-A) decides, in either argument order, even against a
        # smaller A: 1.0 against 2.0 here
        p1 = DFLYCoefficients(VarCircle, L1, 0.5, 0.5)
        p2 = DFLYCoefficients(VarCircle, L1, 0.25, 1.5)
        @test RIM._better_dfly_pair(p1, p2) === p1
        @test RIM._better_dfly_pair(p2, p1) === p1
        # a pair that does not contract ranks last, and between two such the smaller A wins
        @test RIM._better_dfly_pair(DFLYCoefficients(VarCircle, L1, 1.5, 0.0), p1) === p1
        @test RIM._better_dfly_pair(
            DFLYCoefficients(VarCircle, L1, 1.5, 0.0),
            DFLYCoefficients(VarCircle, L1, 2.5, 0.0),
        ).A == 1.5
        # and B/(1-A) within the tolerance is a tie, settled on A: a residue of 1e-16 in B
        # must not decide A by a factor two
        residue = DFLYCoefficients(VarCircle, L1, 1 / 3, 1.6653345369377348e-16)
        exact = DFLYCoefficients(VarCircle, L1, 2 / 3, 0.0)
        @test RIM._better_dfly_pair(residue, exact) === residue
        @test RIM._better_dfly_pair(exact, residue) === residue
        # equal B/(1-A), here 2.0 both ways, is settled on A as well
        @test RIM._better_dfly_pair(
            DFLYCoefficients(VarCircle, L1, 0.5, 1.0),
            DFLYCoefficients(VarCircle, L1, 0.25, 1.5),
        ).A == 0.25
    end

    @testset "the dyadic and triadic maps" begin
        # Where the two statements disagree only by a rounding residue in B the pair with
        # the smaller A comes back. The lift B is exactly zero for 2x and 4x, where 1/2 and
        # 1/4 are representable, and a few units in the last place for 3x and 5x; see the
        # note in the docstring of _dfly_circle_from_lift.
        for d in (2, 3, 4, 5)
            Dd = mod1_dynamic(x -> d * x)
            cl, cb = RIM._dfly_circle_from_lift(Dd), RIM._dfly_circle_branchwise(Dd)
            c = dfly(VarCircle, L1, Dd)
            @test cl.A ≈ 1 / d
            @test cl.B < 1e-15
            @test cb.A ≈ 2 / d
            @test cb.B == 0.0
            # the lift pair, the one with A = 1/d, is what comes back
            @test (c.A, c.B) == (cl.A, cl.B)
            @test c.A ≈ 1 / d
        end
        # the residue is there and is not zero for 3x: two interior branch points and the
        # two ends of the domain, about one unit in the last place of λ each
        @test 0 < RIM._dfly_circle_from_lift(mod1_dynamic(x -> 3 * x)).B < 1e-15
        @test RIM._dfly_circle_from_lift(mod1_dynamic(x -> 4 * x)).B == 0.0
    end

    @testset "what has no circle statement" begin
        # the variation of VarCircle is the one on [0,1] with 0 and 1 identified
        @test_throws ArgumentError dfly(VarCircle, L1, PwMap([x -> 2x], [0, 0.5]))
        # dfly_inf_der has no circle counterpart
        @test_throws ArgumentError dfly(VarCircle, L1, PwMap([x -> sqrt(x)], [0, 1]))
    end

    @testset "measured against the pair" begin
        # The Ulam projection is the average on the cells, and the cyclic sum of the jumps
        # of the averages is bounded by Var_{S^1}(f), so Var_{S^1}(ΠLφ) ≤ Var_{S^1}(Lφ) and
        # the pair holds for the discretized operator as well. This measures that, over the
        # extremal shapes and their iterates; it checks the implementation and proves
        # nothing.
        function var_circ(w)
            n = length(w)
            s = sum(abs(interval(w[i+1]) - interval(w[i])) for i = 1:n-1)
            return s + abs(interval(w[1]) - interval(w[n]))
        end
        l1(w) = sum(abs(interval(x)) for x in w) / interval(length(w))

        c = dfly(VarCircle, L1, D)
        worst, worst_const = 0.0, 0.0
        for n in (64, 256)
            Bas = Ulam(n)
            Q = DiscretizedOperator(Bas, D)
            xc = [(i - 0.5) / n for i = 1:n]
            j1, j2 = Int(floor(n / β)), Int(floor(2n / β))
            vs = [
                ("constant", ones(n)),
                ("spike_1", [i == 1 ? 1.0 : 0.0 for i = 1:n]),
                ("spike_n", [i == n ? 1.0 : 0.0 for i = 1:n]),
                ("spike_j1", [i == j1 ? 1.0 : 0.0 for i = 1:n]),
                ("spike_j2", [i == j2 ? 1.0 : 0.0 for i = 1:n]),
                ("two_spikes_gap", [i in (j1, j1 + 2) ? 1.0 : 0.0 for i = 1:n]),
                ("pm_spikes_gap", [i == j1 ? 1.0 : (i == j1 + 2 ? -1.0 : 0.0) for i = 1:n]),
                ("spikes_wrap", [i in (1, 3, n - 1) ? 1.0 : 0.0 for i = 1:n]),
                ("alternating", [iseven(i) ? 1.0 : -1.0 for i = 1:n]),
                ("alt_spikes", [iseven(i) ? (isodd(i ÷ 2) ? 1.0 : -1.0) : 0.0 for i = 1:n]),
                ("half_circle", [i <= n ÷ 2 ? 1.0 : 0.0 for i = 1:n]),
                ("smooth_bump", [exp(-50 * (x - 0.5)^2) for x in xc]),
                ("smooth_sin", [1 + 0.5 * sin(2π * x) for x in xc]),
                ("const_plus_spike", [i == j1 ? 11.0 : 1.0 for i = 1:n]),
                ("const_minus_spike", [i == j1 ? -9.0 : 1.0 for i = 1:n]),
                ("step_at_beta_2", [x < β - 2 ? 1.0 : 0.0 for x in xc]),
                ("random", (Random.seed!(20261107); randn(n))),
            ]
            for (tag, v0) in vs
                v = copy(v0)
                for k = 1:5
                    w = Q.L * interval.(v)
                    num = var_circ(w)
                    # the instructed measurement is the sup of that same cyclic sum
                    @test normbound(Bas, VarCircle, w) == sup(num)
                    den = interval(c.A) * var_circ(v) + interval(c.B) * l1(v)
                    # a violation would be the whole enclosure of Var_{S^1}(ΠLφ_v) sitting
                    # above the bound
                    @test inf(num) <= sup(den)
                    if inf(den) > 0
                        worst = max(worst, sup(num) / inf(den))
                        tag == "constant" && k == 1 &&
                            (worst_const = max(worst_const, sup(num) / inf(den)))
                    end
                    v = mid.(w)
                end
            end
        end
        # the only excess over 1 is the width of the assembled operator's enclosure
        @test worst <= 1 + 1e-8
        # and the constant function attains the pair: Var_{S^1}(L1) = 2/β = B‖1‖_1 exactly
        @test worst_const >= 1 - 1e-8
    end
end
