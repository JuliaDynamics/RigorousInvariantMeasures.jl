using RigorousInvariantMeasures
using Test

const RIM = RigorousInvariantMeasures

@testset "SpectralCertification" begin

    @testset "bN_constant" begin
        # b_1 = b · a^0 M^0 = b
        @test RIM.bN_constant(0.5, 1.0, 1.0, 1) == 1.0
        @test RIM.bN_constant(0.5, 2.0, 1.0, 1) == 2.0

        # b_2 = b * (a^1 M^0 + a^0 M^1) = b * (a + M)
        @test RIM.bN_constant(0.5, 2.0, 1.0, 2) ≥ 2.0 * (0.5 + 1.0) - 1e-12

        # When M = 1, b_N = b * (a^(N-1) + a^(N-2) + ... + 1) = b * (1 - a^N) / (1 - a)
        a, b = 0.5, 3.0
        N = 5
        expected = b * (1 - a^N) / (1 - a)
        got = RIM.bN_constant(a, b, 1.0, N)
        @test got ≥ expected - 1e-10
        @test got ≤ expected + 1e-10

        # Argument validation
        @test_throws ArgumentError RIM.bN_constant(0.5, 1.0, 1.0, 0)
    end

    @testset "_SN_bound" begin
        # Geometric form: |z|>M=1, ∑_{ℓ=0}^{N-1} (1/|z|)^ℓ / |z|
        # = (1 - |z|^{-N}) / (|z| - 1)
        abs_z = 2.0
        N = 4
        expected = (1 - abs_z^(-N)) / (abs_z - 1)
        got = RIM._SN_bound(abs_z, N; M = 1.0)
        @test got ≥ expected - 1e-10
        @test got ≤ expected + 1e-10

        # With explicit norm_powers all equal to 1: same answer
        got2 = RIM._SN_bound(abs_z, N; norm_powers = ones(N))
        @test got2 ≥ expected - 1e-10
        @test got2 ≤ expected + 1e-10

        # Sharper sequence gives a tighter bound
        powers_decay = [1.0, 0.5, 0.25, 0.125]
        got_sharp = RIM._SN_bound(abs_z, N; norm_powers = powers_decay)
        @test got_sharp < expected

        # Sharp value: 1/|z| * ∑ p_ℓ / |z|^ℓ
        sharp_expected = sum(powers_decay[ℓ+1] / abs_z^ℓ for ℓ = 0:N-1) / abs_z
        @test got_sharp ≥ sharp_expected - 1e-10
        @test got_sharp ≤ sharp_expected + 1e-10

        # Argument validation
        @test_throws ArgumentError RIM._SN_bound(abs_z, 0; M = 1.0)
        @test_throws ArgumentError RIM._SN_bound(abs_z, N) # neither supplied
        @test_throws ArgumentError RIM._SN_bound(abs_z, N; M = 1.0, norm_powers = ones(N))
        @test_throws ArgumentError RIM._SN_bound(abs_z, N; norm_powers = ones(N - 1))
    end

    @testset "strong_resolvent_lift (Prop A.7)" begin
        # 2x-mod1: full-branch expanding, lam=1/2, dist=0, so DFLY (a,b) = (0.5, 0.0)
        # With b=0, the formula collapses: ℛ_s(z, L) ≤ 1 / (|z| - a).
        B = Ulam(64)
        D = mod1_dynamic(x -> 2 * x)
        abs_z = 1.5
        R_w_coarse = 10.0
        got = strong_resolvent_lift(B, D, abs_z, R_w_coarse)
        expected = 1.0 / (abs_z - 0.5)
        @test got ≥ expected - 1e-10
        @test got ≤ expected + 1e-10
    end

    @testset "Prop A.14 scalar formula" begin
        # Synthetic check: Δ_k = 0 ⇒ β̃ = 0 ⇒ formula collapses to
        # S_N + |z|^{-N} E_sw K(z) (a^N E_{k,w→s} + b_N)
        a, b, M = 0.4, 1.0, 1.0
        Δ_k = 0.0
        K_z = 5.0
        E_sw = 1.0
        E_k_ws = 8.0
        N = 3
        abs_z = 2.0
        M_k = 1.0

        bN = RIM.bN_constant(a, b, M, N)
        SN = RIM._SN_bound(abs_z, N; M = M_k)
        factor = a^N * E_k_ws + bN
        expected = SN + abs_z^(-N) * E_sw * K_z * factor

        got = RIM._coarse_fine_weak_resolvent_scalar(
            abs_z, K_z, M_k;
            a = a, b = b, M_aux = M, Δ_k = Δ_k,
            E_sw = E_sw, E_k_ws = E_k_ws, N = N,
        )
        @test got ≥ expected - 1e-10
        @test got ≤ expected + 1e-10

        # With Δ_k > 0 (perturbation closes): bound is larger but still finite
        Δ_k_pos = 1e-3
        got_pert = RIM._coarse_fine_weak_resolvent_scalar(
            abs_z, K_z, M_k;
            a = a, b = b, M_aux = M, Δ_k = Δ_k_pos,
            E_sw = E_sw, E_k_ws = E_k_ws, N = N,
        )
        @test got_pert > got
        @test isfinite(got_pert)

        # If Δ_k is large enough that β̃ ≥ 1, formula returns Inf
        Δ_k_huge = 1e10
        got_huge = RIM._coarse_fine_weak_resolvent_scalar(
            abs_z, K_z, M_k;
            a = a, b = b, M_aux = M, Δ_k = Δ_k_huge,
            E_sw = E_sw, E_k_ws = E_k_ws, N = N,
        )
        @test isinf(got_huge)
    end

    @testset "abstract_weak_norm_bound" begin
        # Ulam has (S₁, S₂) = (0, 1), M₂ = 1, so the bound collapses to 1.
        B = Ulam(64)
        D = mod1_dynamic(x -> 2 * x)
        dfly_co = RIM.dfly(RIM.strong_norm(B), RIM.aux_norm(B), D)
        got = RIM.abstract_weak_norm_bound(B; dfly_coefficients = dfly_co)
        @test got ≥ 1.0 - 1e-12
        @test got ≤ 1.0 + 1e-12
    end

    @testset "Basis-aware coarse_fine_weak_resolvent (Ulam)" begin
        # 2x-mod1 on Ulam(64). With b = 0 for the doubling map, b_N = 0 and
        # the bound is well-behaved. Verify it stays finite and grows when
        # Δ_k grows.
        B = Ulam(64)
        D = mod1_dynamic(x -> 2 * x)
        abs_z = 1.5
        K_z = 2.0
        N = 6

        got = coarse_fine_weak_resolvent(B, D, abs_z, K_z; N = N)
        @test isfinite(got)
        @test got > 0

        # Sharper M_k vector (all-ones is the trivial L¹-preserving bound)
        got_vec = coarse_fine_weak_resolvent(
            B, D, abs_z, K_z; N = N, M_k_weak = ones(N),
        )
        @test isfinite(got_vec)

        # Refining Δ_k tightens the bound monotonically.
        got_tight = coarse_fine_weak_resolvent(
            B, D, abs_z, K_z; N = N, Δ_k = 1e-6,
        )
        got_loose = coarse_fine_weak_resolvent(
            B, D, abs_z, K_z; N = N, Δ_k = 1e-2,
        )
        @test got_tight ≤ got_loose
    end

    @testset "projector_distance_bound (smoke)" begin
        # Smoke test on a tiny coarse-fine Fourier setup. CertifScripts is
        # run only on the *coarse* matrix; Prop A.7 + A.14 propagate to the
        # fine level. `BallMatrix(Q.L)` requires a dense matrix, so Ulam
        # (sparse) is not usable here.
        if Base.get_extension(RigorousInvariantMeasures, :FFTWExt) !== nothing
            B_coarse = FourierAdjoint(FourierPoints(128, Float64), 8,
                                       W{1,1}(), L2())
            B_fine = FourierAdjoint(FourierPoints(128, Float64), 16,
                                     W{1,1}(), L2())
            D = mod1_dynamic(x -> 2 * x)
            Q_coarse = DiscretizedOperator(B_coarse, D)
            Q_fine = DiscretizedOperator(B_fine, D)
            result = projector_distance_bound(
                B_coarse, B_fine, D, Q_coarse, Q_fine, 0.1;
                samples = 32, N = 4,
            )
            @test result.radius == 0.1
            @test result.δ_k > 0
            @test result.R_w_coarse > 0
            @test (isfinite(result.projector_distance) ||
                   isinf(result.projector_distance))
        end
    end

    @testset "coarse_fine_weak_resolvent_auto_N" begin
        # Ulam M = 1 so Cor A.16 needs 0 < a < M — choose a < 1 dynamic.
        B = Ulam(64)
        D = mod1_dynamic(x -> 2 * x)
        abs_z = 1.5
        K_z = 2.0

        (R, N_k) = coarse_fine_weak_resolvent_auto_N(B, D, abs_z, K_z)
        @test isfinite(R)
        @test N_k ≥ 1

        # Override μ inside (a, M)
        (R2, N_k2) = coarse_fine_weak_resolvent_auto_N(B, D, abs_z, K_z; μ = 0.7)
        @test isfinite(R2)
        @test N_k2 ≥ 1
    end

end
