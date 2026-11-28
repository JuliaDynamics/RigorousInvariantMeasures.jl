using RigorousInvariantMeasures
using IntervalArithmetic
using FFTW  # FFTWExt: the doubling-map testset assembles a FourierAnalytic basis

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
        As = [first(analytic_dfly(Eρ(ρ), ρ′, 1.0, K; L1μ_bound = 1.0)) for K in (0, 5, 10, 20)]
        Bs = [last(analytic_dfly(Eρ(ρ), ρ′, 1.0, K; L1μ_bound = 1.0)) for K in (0, 5, 10, 20)]
        @test issorted(As; rev = true)
        @test issorted(Bs)
        @test As[end] < 1e-2

        K, A, B = analytic_dfly_choose_K(Eρ(ρ), ρ′, 1.0; target_A = 0.5, L1μ_bound = 1.0)
        @test A <= 0.5
        @test first(analytic_dfly(Eρ(ρ), ρ′, 1.0, K - 1; L1μ_bound = 1.0)) > 0.5

        # The L1(mu) operator bound scales B and only B — it plays no part in
        # the contraction, which comes purely from the analytic gain.
        A1, B1 = analytic_dfly(Eρ(ρ), ρ′, 1.0, 8; L1μ_bound = 1.0)
        A2, B2 = analytic_dfly(Eρ(ρ), ρ′, 1.0, 8; L1μ_bound = 3.0)
        @test A1 == A2
        @test B2 ≈ 3 * B1 rtol = 1e-12

        @test_throws ErrorException analytic_dfly(Eρ(2.0), 1.5, 1.0, 5; L1μ_bound = 1.0)
    end

    @testset "the ellipse gain feeds the DFLY" begin
        # bernstein_expansion certifies ρ′, which is exactly what analytic_dfly
        # consumes: the two halves compose.
        ρ = 1.5
        ρ′ = bernstein_expansion(z -> 2z^2 - 1, ρ; n = 4096)
        @test ρ′ > ρ
        K, A, B = analytic_dfly_choose_K(Eρ(ρ), ρ′, 1.0; target_A = 0.4, L1μ_bound = 1.0)
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

    @testset "annulus geometry and its certificate" begin
        # T(x) = 2x + a sin(2πx)/(2π), a = 1/2, degree 2.  On Im x = ±r the exact values are
        # min |Im T| = 2r - a sinh(2πr)/(2π) and min |T'| = 2 - a cosh(2πr) (at Re x = 1/2),
        # and Re T' > 0 on the closed strip iff a cosh(2πr) < 2, i.e. r < acosh(4)/(2π) ≈ 0.3285.
        twoπ = 2 * interval(π)
        a = interval(0.5)
        T(x) = 2x + a * sin(twoπ * x) / twoπ
        Rex(r) = 2r - 0.5 * sinh(2π * r) / (2π)
        Wex(r) = 2 / (2 - 0.5 * cosh(2π * r))
        @test_throws ArgumentError AnnulusGeometry(0.2, 0.1, 0.3)
        G = AnnulusGeometry(0.1, 0.12, 0.14)
        C = certify_compatibility(G, T; degree = 2, n = 512, cover = 16)
        @test C.G === G && C.lift === T && C.degree == 2
        @test Rex(0.1) - 1e-3 < C.R_image <= Rex(0.1)
        @test Wex(0.1) <= C.W < Wex(0.1) + 1e-3
        @test_throws ArgumentError certify_compatibility(G, T; degree = 1, n = 512, cover = 16)
        @test_throws ArgumentError certify_compatibility(AnnulusGeometry(0.1, 0.12, 0.15), T;
                                                         degree = 2, n = 512, cover = 16)
        @test_throws ArgumentError certify_compatibility(AnnulusGeometry(0.35, 0.4, 0.45), T;
                                                         degree = 2, n = 512, cover = 16)
        # the search: r = 0.35 fails the derivative check; R_image / r decreases in r here
        rs = [0.05, 0.1, 0.2, 0.35]
        S = search_annulus_geometry(T, rs; degree = 2, n = 512, cover = 16)
        @test !S.table[end].feasible && S.r == 0.05 && S.cost == S.R_image / S.r
        Sw = search_annulus_geometry(T, rs; degree = 2, cost = (r, R, W) -> R - r,
                                     n = 512, cover = 16)
        @test Sw.cost == maximum(e.cost for e in Sw.table)
        Rb = Sw.r + (Sw.R_image - Sw.r) / 2
        @test certify_compatibility(AnnulusGeometry(Sw.r, (Sw.r + Rb) / 2, Rb), T;
                                    degree = 2, n = 512, cover = 16).R_image == Sw.R_image
        @test_throws ArgumentError search_annulus_geometry(T, [0.35]; degree = 2, n = 512,
                                                           cover = 16)
    end

    @testset "annulus: L¹ Lasota–Yorke, coefficient decay, Galerkin defect" begin
        # lift x -> 2x + sin(2πx)/10: Im(x ± ir) maps to ±(2r + cos(2πx) sinh(2πr)/10), so
        # min |Im| on ∂S_r is 2r - sinh(2πr)/10, 0.2393 at r = 0.2
        twoπ = 2 * interval(π)
        dbl = x -> 2x + sin(twoπ * x) / 10
        G = AnnulusGeometry(0.2, 0.21, 0.23)
        C = certify_compatibility(G, dbl; degree = 2, n = 512, cover = 16)
        # the coefficient bound |(L1)^_j| = |∫ e^{-2πij(2y + sin(2πy)/10)}| <= e^{-2πR|j|}, by a
        # periodic trapezoid rule in Float64 (spectrally accurate here; an independent check)
        q = IntervalArithmetic.sup(annulus_coefficient_decay(C))
        n = 4096
        ys = (0:n-1) ./ n
        for j in 1:6
            c = abs(sum(exp(-2π * im * j * (2y + sin(2π * y) / 10)) for y in ys) / n)
            @test c <= q^j
        end
        # the constants against their formulas, and the full-norm form of analytic_dfly
        for M in (0, 5, 40)
            c = annulus_dfly_l1(C, M)
            qq = exp(-2π * (0.23 - 0.21))
            @test IntervalArithmetic.inf(c.A_M) <= 2qq^(M + 1) / (1 - qq) * (1 + 1e-12)
            @test 2qq^(M + 1) / (1 - qq) <= IntervalArithmetic.sup(c.A_M) * (1 + 1e-12)
            Bp = 2 * sum((exp(2π * 0.21 * j) for j in 1:M); init = 0.0)
            @test IntervalArithmetic.inf(c.Bp_M) <= Bp * (1 + 1e-12) + 1e-300
            @test Bp <= IntervalArithmetic.sup(c.Bp_M) * (1 + 1e-12) + 1e-300
            A, B = analytic_dfly(Aη(0.21), 0.23, 1.0, M)
            @test A == IntervalArithmetic.sup(c.A_M)
            @test B ≈ 1 + IntervalArithmetic.mid(c.Bp_M) rtol = 1e-12
        end
        M, c = annulus_dfly_l1_smallest_M(C)
        @test IntervalArithmetic.sup(c.A_M) < 1 &&
              IntervalArithmetic.sup(annulus_dfly_l1(C, M - 1).A_M) >= 1
        @test_throws ArgumentError annulus_dfly_l1_smallest_M(C; Mmax = 0)
        # the defect: its two terms, decreasing in N, and M_w passed explicitly
        δ = annulus_galerkin_defect(C, 32; M_w = 1.25)
        δex = 2exp(-2π * 0.23 * 33) / (1 - exp(-2π * 0.23)) + 1.25exp(-2π * 0.21 * 33)
        @test IntervalArithmetic.inf(δ) <= δex * (1 + 1e-12) &&
              δex <= IntervalArithmetic.sup(δ) * (1 + 1e-12)
        @test IntervalArithmetic.sup(annulus_galerkin_defect(C, 64)) <
              IntervalArithmetic.inf(annulus_galerkin_defect(C, 32))
        cb = annulus_dfly_l1(C, M; T = BigFloat)
        @test cb.A_M isa Interval{BigFloat}
        @test IntervalArithmetic.sup(cb.A_M) <= IntervalArithmetic.sup(c.A_M)
    end

    @testset "invariant_measure_strong_norm_bound refuses a continuity constant" begin
        D = mod1_dynamic(x -> 2 * x; circle_smoothness = :ω)
        B = FourierAnalytic(16, 64; η = 0.1)
        @test_throws ArgumentError RigorousInvariantMeasures.invariant_measure_strong_norm_bound(B, D)
        # a pair with A < 1 and the L¹ norm: B/(1 - A)
        c = DFLYCoefficients(Aη(0.1), L1, 0.5, 3.0)
        @test RigorousInvariantMeasures.invariant_measure_strong_norm_bound(B, D;
                                                                            dfly_coefficients = c) >= 6.0
    end

end
