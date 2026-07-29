using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Analytic DFLY (Aη and Eρ)" begin
    RIM = RigorousInvariantMeasures

    @testset "strip expansion" begin
        # Ground truth: x -> 2x is z -> z^2, which sends |z| = e^{2πη} to
        # e^{4πη}, so the certified strip must be 2η, approached from below.
        for η in (0.05, 0.1, 0.2)
            η′ = strip_expansion(z -> z^2, η; n = 2048)
            @test η′ <= 2η + 1e-12          # rigorous lower bound
            @test η′ > 0.98 * 2η            # and tight
            @test η′ > η                    # the strip really is enlarged
        end
        # Both boundary circles are certified.
        η_in, η_out = annulus_expansion(z -> z^2, 0.1; n = 2048)
        @test η_in > 0.1 && η_out > 0.1
        # A rotation enlarges nothing.
        @test strip_expansion(z -> z, 0.1; n = 512) <= 0.1 + 1e-12
    end

    @testset "Aη constants" begin
        η, η′ = 0.1, 0.2
        As = [first(analytic_dfly(Aη(η), η′, 1.0, K)) for K in (0, 5, 10, 20)]
        Bs = [last(analytic_dfly(Aη(η), η′, 1.0, K)) for K in (0, 5, 10, 20)]
        @test issorted(As; rev = true)      # A decreases in K
        @test issorted(Bs)                  # B grows in K
        @test all(As .> 0) && all(Bs .>= 1)
        @test As[end] < 1e-4                # geometric decay

        # B(0) sums the single mode k = 0.
        @test last(analytic_dfly(Aη(η), η′, 1.0, 0)) ≈ 1.0 rtol = 1e-12

        K, A, B = analytic_dfly_choose_K(Aη(η), η′, 1.0; target_A = 0.5)
        @test A <= 0.5
        @test first(analytic_dfly(Aη(η), η′, 1.0, K - 1)) > 0.5   # K is minimal

        # The gain must be strictly positive.
        @test_throws ErrorException analytic_dfly(Aη(0.2), 0.2, 1.0, 5)
        @test_throws ErrorException analytic_dfly(Aη(0.2), 0.1, 1.0, 5)
    end

    @testset "Eρ constants" begin
        ρ, ρ′ = 1.5, 2.25                   # ρ′ = ρ^2, the T_2 gain
        As = [first(analytic_dfly(Eρ(ρ), ρ′, 1.0, K)) for K in (0, 5, 10, 20)]
        Bs = [last(analytic_dfly(Eρ(ρ), ρ′, 1.0, K)) for K in (0, 5, 10, 20)]
        @test issorted(As; rev = true)
        @test issorted(Bs)
        @test As[end] < 1e-2

        K, A, B = analytic_dfly_choose_K(Eρ(ρ), ρ′, 1.0; target_A = 0.5)
        @test A <= 0.5
        @test first(analytic_dfly(Eρ(ρ), ρ′, 1.0, K - 1)) > 0.5

        # The L1(mu) operator bound scales B and only B — it plays no part in
        # the contraction, which comes purely from the analytic gain.
        A1, B1 = analytic_dfly(Eρ(ρ), ρ′, 1.0, 8)
        A2, B2 = analytic_dfly(Eρ(ρ), ρ′, 1.0, 8; L1μ_bound = 3.0)
        @test A1 == A2
        @test B2 ≈ 3 * B1 rtol = 1e-12

        @test_throws ErrorException analytic_dfly(Eρ(2.0), 1.5, 1.0, 5)
    end

    @testset "the ellipse gain feeds the DFLY" begin
        # bernstein_expansion certifies ρ′, which is exactly what analytic_dfly
        # consumes: the two halves compose.
        ρ = 1.5
        ρ′ = bernstein_expansion(z -> 2z^2 - 1, ρ; n = 4096)
        @test ρ′ > ρ
        K, A, B = analytic_dfly_choose_K(Eρ(ρ), ρ′, 1.0; target_A = 0.4)
        @test A <= 0.4 && B > 0
    end
end
