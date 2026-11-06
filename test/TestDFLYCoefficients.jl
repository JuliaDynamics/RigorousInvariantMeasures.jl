@testset "DFLY coefficients carry their norms" begin

    using RigorousInvariantMeasures
    RIM = RigorousInvariantMeasures

    β = (1 + sqrt(23)) / 2
    D = mod1_dynamic(x -> β * x; full_branch = false)
    c = dfly(VarInterval, L1, D)

    # the non-full-branch constants: A = 2/β and B = 2/min|I_i|, the shortest branch
    @test c.A ≈ 2 / β
    @test c.B ≈ 2 / (1 - 2 / β)
    @test strong_kind(c) === VarInterval
    @test weak_kind(c) === L1

    # it stands in for the pair it replaces
    A, B = c
    @test (A, B) == (c.A, c.B)
    @test length(c) == 2
    @test c[1] == c.A && c[2] == c.B
    @test Tuple(c) == (c.A, c.B)
    @test collect(c) == [c.A, c.B]
    @test_throws BoundsError c[3]
    @test occursin("VarInterval", sprint(show, c))

    # a bare tuple is still accepted and records no convention
    @test strong_kind((0.5, 1.0)) === nothing
    @test weak_kind((0.5, 1.0)) === nothing

    # the same operator has one pair on the interval and a sharper one on the circle, and
    # the a priori bound on the density is read in the norm of the pair it came from
    Bas = Ulam(256)
    circle = DFLYCoefficients(VarCircle, L1, 2 / β, 2 / β)
    @test strong_kind(circle) === VarCircle
    us_int = RIM.invariant_measure_strong_norm_bound(Bas, D; dfly_coefficients = c)
    us_cir = RIM.invariant_measure_strong_norm_bound(Bas, D; dfly_coefficients = circle)
    @test us_int ≈ c.B / (1 - c.A)
    @test us_cir ≈ 2 / (β - 2)
    @test us_cir < us_int

    # the other dfly methods report their own norms
    @test strong_kind(dfly(Lipschitz, L1, mod1_dynamic(x -> 2 * x))) === Lipschitz
    @test strong_kind(dfly(W{1,1}, L1, mod1_dynamic(x -> 2 * x))) === W{1,1}
    @test strong_kind(dfly(VarInterval, L1, RIM.UniformKernelUlamPeriodic(Ulam(64), 3))) ===
          VarInterval

    # a candidate measured beside a pair is measured in the pair's own strong norm
    w = [1.0, 2.0, 2.0, 0.5]
    @test normbound(Ulam(4), something(strong_kind(circle), strong_norm(Ulam(4))), w) ==
          normbound(Ulam(4), VarCircle, w)
    @test normbound(Ulam(4), something(strong_kind(c), strong_norm(Ulam(4))), w) ==
          normbound(Ulam(4), VarInterval, w)

end
