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

    @testset "both halves of the DFLY are necessary (doubling map)" begin
        # For x -> 2x the transfer operator acts on Fourier modes as
        # L e_m = e_{m/2} for even m, 0 for odd.  That single fact pins the
        # structure of the analytic DFLY from both sides.
        D = mod1_dynamic(x -> 2 * x)
        B = FourierAnalytic(16, 4096, W{3,1})
        Q = mid.(RIM.assemble(B, D; ϵ = 1e-14, max_iter = 100))
        n, K = length(B), B.k
        idx(m) = m >= 0 ? m + 1 : n + m + 1

        for m in (2, 4, 6, 8)
            col = Q[:, idx(m)]
            @test isapprox(col[idx(m ÷ 2)], 1.0; atol = 1e-12)
            @test maximum(abs.(col[setdiff(1:n, [idx(m ÷ 2)])])) < 1e-14
        end

        # f = e_{2k} lies in L¹ ∩ A_η, with ||f||_{L¹} = 1 and Lf = e_k.
        # ||Lf||_{A_η} = e^{2πηk} is unbounded in k, so there is no constant C
        # with ||Lf||_{A_η} ≤ C||f||_{L¹}: the L¹ bound on the coefficients is
        # uniform in k and cannot carry the weighted sum on its own.
        η = 0.1
        ratios_L1 = [exp(2π * η * k) for k in (2, 4, 8)]
        @test issorted(ratios_L1)
        @test ratios_L1[end] > 100

        # The strong norm, on the other hand, contracts geometrically:
        # ||Lf||_{A_η}/||f||_{A_η} = e^{-2πηk} -> 0. That is the analytic gain,
        # and it is invisible to any real-variable quantity on [0,1].
        ratios_strong = [exp(-2π * η * k) for k in (2, 4, 8)]
        @test issorted(ratios_strong; rev = true)
        @test ratios_strong[end] < 0.01

        # Hence the split: A(K) from the analytic gain, B(K) from the L¹ bound.
        Kc, A, Bc = analytic_dfly_choose_K(Aη(η), 2η, 1.0; target_A = 0.5)
        @test A <= 0.5 && Bc > 1
    end

    @testset "a priori constants from the complex neighbourhood" begin
        # min |4z| over ∂E_ρ is attained on the minor axis, |z| = (ρ-1/ρ)/2,
        # so the certified value must be 2(ρ - 1/ρ).
        for ρ in (1.5, 2.0, 3.0)
            md = min_modulus_on_ellipse(z -> 4 * z, ρ; n = 4096)
            @test md <= 2 * (ρ - 1 / ρ) + 1e-9
            @test md > 0.999 * 2 * (ρ - 1 / ρ)
        end
        @test min_modulus_on_circle(z -> complex(interval(2.0), interval(0.0)), 0.1;
                                    n = 32) ≈ 2.0 rtol = 1e-12

        # C₂ = #branches / min|T'|, from the complex neighbourhood only.
        @test analytic_transfer_bound(2.0, 2) ≈ 1.0 rtol = 1e-12
        @test analytic_transfer_bound(4.0, 2) ≈ 0.5 rtol = 1e-12
        @test_throws ErrorException analytic_transfer_bound(0.0, 2)

        # The degenerate LY is the continuity constant of L on A_η: it holds for
        # every f, carries no auxiliary term, and need not be < 1 — compactness,
        # not contraction, is what drives the certification.
        Tprime = z -> complex(interval(2.0), interval(0.0))
        for η in (0.1, 0.2, 0.3)
            η′ = strip_expansion(z -> z^2, η; n = 2048)
            C₂ = analytic_transfer_bound(min_modulus_on_circle(Tprime, η; n = 64), 2)
            @test C₂ ≈ 1.0 rtol = 1e-12
            A, B = analytic_dfly_degenerate(Aη(η), η′, C₂)
            @test B == 0.0
            @test A >= C₂                       # G(δ) ≥ 1
            @test isfinite(A)
        end

        # Bigger gain ⇒ smaller constant, monotonically.
        As = map((0.1, 0.2, 0.3)) do η
            η′ = strip_expansion(z -> z^2, η; n = 2048)
            first(analytic_dfly_degenerate(Aη(η), η′, 1.0))
        end
        @test issorted(As; rev = true)

        # Chebyshev side, same shape.
        for ρ in (1.5, 2.0, 3.0)
            ρ′ = bernstein_expansion(z -> 2z^2 - 1, ρ; n = 4096)
            C₂ = analytic_transfer_bound(min_modulus_on_ellipse(z -> 4z, ρ; n = 2048), 2)
            A, B = analytic_dfly_degenerate(Eρ(ρ), ρ′, C₂)
            @test B == 0.0
            @test isfinite(A) && A > 0
        end
    end

end
