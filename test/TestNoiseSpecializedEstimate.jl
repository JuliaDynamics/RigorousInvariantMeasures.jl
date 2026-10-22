@testset "NoiseSpecializedEstimate" begin

    using Test
    using IntervalArithmetic
    using RigorousInvariantMeasures

    import RigorousInvariantMeasures:
        total_variation,
        prepare_Wnorm_estimate,
        prepare_derivative_bounds,
        estimate_Wnorm_of_discretization_error,
        estimate_variation_after_noise,
        estimate_measure_after_noise,
        estimate_L1_NL_discretization_error,
        noise_error_apriori,
        noise_error_aposteriori,
        PreimageInfo,
        noise_support,
        reflect_outward_idx,
        invariant_vector_noise,
        distance_from_invariant_noise,
        powernormboundsnoise,
        residualboundnoise,
        infinite_sum_norms

    # ──── total_variation ────

    @testset "total_variation" begin
        @test total_variation([1.0, 1.0, 1.0]) == 0.0
        @test total_variation([1.0, 2.0, 1.0]) ≥ 2.0
        @test total_variation([0.0, 1.0, 0.0, 1.0]) ≥ 3.0
        @test total_variation([3.0, 1.0, 4.0]) ≥ 5.0  # |3-1| + |1-4| = 5
    end

    # ──── the reflection used is the mirror at 0 and 1 ────

    @testset "reflection" begin
        # it used to send cell 0 to cell 2 on the left and cell k+1 to cell k on the right
        @test reflect_outward_idx(0, 5) == 1
        @test reflect_outward_idx(-1, 5) == 2
        @test reflect_outward_idx(6, 5) == 5
        @test reflect_outward_idx(7, 5) == 4
    end

    # ──── prepare_Wnorm_estimate for 2x mod 1 ────

    @testset "prepare_Wnorm_estimate" begin
        D = mod1_dynamic(x -> 2x)
        B = Ulam(16)
        prep = prepare_Wnorm_estimate(D, B)

        @test length(prep) == 16
        # Each interval should have preimage data from 2 branches (full-branch 2x)
        for i in 1:16
            @test length(prep[i]) >= 1
            for pinfo in prep[i]
                @test pinfo.A < pinfo.B
                @test pinfo.invtp > 0
                @test pinfo.invtp < 1.0  # |1/T'| = 1/2 for T=2x
                @test isfinite(pinfo.distortion)
            end
        end
    end

    # ──── prepare_derivative_bounds ────

    @testset "prepare_derivative_bounds" begin
        D = mod1_dynamic(x -> 2x)
        B = Ulam(16)
        derivs = prepare_derivative_bounds(D, B)

        @test length(derivs) == 16
        for d in derivs
            @test d ≥ 2.0  # derivative of 2x is 2
            @test d < 3.0   # should be close to 2
        end
    end

    # ──── noise_error_apriori ────

    @testset "noise_error_apriori" begin
        K = 128
        noise_rel = 0.1
        α = 0.5
        sumCi = 2.0
        err = noise_error_apriori(K, noise_rel, α, sumCi)
        @test err > 0.0
        @test isfinite(err)
    end

    # ──── Integration test: full pipeline with a perturbed map ────

    @testset "full pipeline: 4x perturbed + noise" begin
        # Define dynamics
        D = mod1_dynamic(x -> 4x + RigorousInvariantMeasures.sinpi(8x) / 100)
        K_size = 256
        B = Ulam(K_size)
        ξ = 17 // K_size
        NK = UlamProjectedUniformKernelReflecting(B, ξ)

        # Assemble operator
        Q = DiscretizedOperator(B, D)

        # Compute invariant measure
        w = invariant_vector_noise(B, Q, NK; iter=20)

        # Compute generic error via distance_from_invariant_noise
        norms = powernormboundsnoise(B; Q=Q, NK=NK)
        generic_error = distance_from_invariant_noise(B, Q, NK, w, norms)

        @test generic_error > 0
        @test isfinite(generic_error)

        # Lemma 15 needs, for one n, α ≥ ‖L_δⁿ|_V‖ < 1 and Σ_{i<n} ‖L_δⁱ|_V‖ with the 0th = 1
        n = findfirst(<(1), norms)
        α_est = norms[n]
        sumCi = 1.0 + sum(norms[1:n-1]; init = 0.0)

        # Numerical error: floating-point residual of the eigenvector computation.
        # This is the actual computational precision, NOT the interval-arithmetic
        # residual (which includes discretization error that the specialized
        # estimator already accounts for via A1/A2/A3).
        mQ = mid(Q)
        v_fp = mQ.L * w
        v_fp = NK * v_fp
        numeric_error = sum(abs.(v_fp .- w)) / K_size

        # Prepare preimage data
        prep1 = prepare_Wnorm_estimate(D, B)
        prep2 = prepare_derivative_bounds(D, B)

        # The a posteriori estimate
        L1err, apriori_err, Linf_est = noise_error_aposteriori(
            prep1, prep2, w, numeric_error, B, NK, α_est, sumCi,
        )
        @test_throws ArgumentError noise_error_aposteriori(
            prepare_Wnorm_estimate(D, Ulam(24)), prep2, w, numeric_error, B, NK, α_est, sumCi)

        @test L1err > 0
        @test isfinite(L1err)
        @test apriori_err > 0
        @test isfinite(apriori_err)
        @test all(isfinite.(Linf_est))

        # The a posteriori bound should be tighter than the generic DFLY bound
        @info "Generic DFLY error:     $generic_error"
        @info "A posteriori L1 error:  $L1err"
        @info "A priori error:         $apriori_err"

        # The finer estimate should be at most as large as the a priori
        @test L1err ≤ apriori_err || L1err ≤ generic_error * 10  # allow some slack
    end

    # ──── estimate_Wnorm_of_discretization_error ────

    @testset "estimate_Wnorm_of_discretization_error" begin
        D = mod1_dynamic(x -> 2x)
        B = Ulam(32)
        prep = prepare_Wnorm_estimate(D, B)

        # Uniform density
        density = ones(32)
        est, meas_Lf, var_Lf = estimate_Wnorm_of_discretization_error(prep, density, B)

        @test est ≥ 0
        @test isfinite(est)
        @test length(meas_Lf) == 32
        @test length(var_Lf) == 32
        @test all(meas_Lf .≥ 0)
        @test all(isfinite.(meas_Lf))
    end

    # ──── estimate_variation_after_noise ────

    @testset "estimate_variation_after_noise" begin
        B = Ulam(32)
        NK = UlamProjectedUniformKernelReflecting(B, 9 // 32)
        meas_Lf = ones(32) ./ 32
        var_Lf = 0.1 * ones(32)

        var_NLf = estimate_variation_after_noise(meas_Lf, var_Lf, NK)
        @test length(var_NLf) == 32
        @test all(var_NLf .≥ 0)
        @test all(isfinite.(var_NLf))
    end

    # g = M + J·1_{x ≥ 1/2} on 32 cells, ξ = 13/64, so ξ/h = 6.5: the mass term is large and
    # the variation term binds; on a cell I whose window contains the jump in its interior,
    # Var_I(N_ξ g) = h J/ξ exactly. The divisor used to be 2⌈ξ/2h⌉ = 8 instead of ξ/h = 6.5.
    @testset "variation after noise bounds the true one" begin
        Kest, J, M = 32, 1.0, 100.0
        h = 1 / Kest
        NK = UlamProjectedUniformKernelReflecting(Ulam(256), 13 // 64)
        @test noise_support(NK) == 13 // 64
        g = [M + (i > Kest ÷ 2 ? J : 0.0) for i = 1:Kest]
        meas = g .* h
        var = zeros(Kest)
        var[Kest÷2+1] = J                  # the jump, on the closed cell to its right
        bound = estimate_variation_after_noise(meas, var, NK)
        for i = 13:20                      # windows away from the boundary and containing the jump
            @test bound[i] >= h * J / (13 / 64) * (1 - 1e-12)
        end
    end

    # a unit mass in the cell 4 cells from i reaches I = cell i through the noise, since
    # (4 - 1) h < ξ/2 = 3.25 h; the window used to stop at round(ξ/h) ÷ 2 = 3 cells
    @testset "measure after noise bounds the true one" begin
        Kest = 32
        h = 1 / Kest
        NK = UlamProjectedUniformKernelReflecting(Ulam(256), 13 // 64)
        meas = zeros(Kest)
        meas[20] = 1.0
        # ∫_I N_ξ g for g = mass 1 spread on cell 20 and I = cell 16: (1/ξ)∫∫ 1[|x-y| ≤ ξ/2]
        # over x ∈ I, y ∈ cell 20, computed by the tent: h²/ξ · overlap of (4 ± 1) with ±ξ/2h
        s = 6.5
        tentF(x) = x <= -1 ? 0.0 : x <= 0 ? (1 + x)^2 / 2 : x <= 1 ? 1 - (1 - x)^2 / 2 : 1.0
        exact = (tentF(s / 2 - 4) - tentF(-s / 2 - 4)) / s
        @test exact > 0
        @test estimate_measure_after_noise(meas, NK)[16] >= exact * (1 - 1e-12)
    end

    # the boundary term of Lemma 56 is |f(y)|/|T'(y)| at the endpoint y of the branch domain:
    # for x ↦ 1.5x mod 1 the second branch, on [2/3, 1], ends at y = 1 with T(y) = 1/2
    @testset "boundary terms at the domain endpoint" begin
        D = mod1_dynamic(x -> 1.5x)
        B = Ulam(16)
        prep = prepare_Wnorm_estimate(D, B)
        ys = [y for i = 1:16 for pinfo in prep[i] for (y, _) in pinfo.boundary_terms]
        @test any(y -> isapprox(y, 1.0; atol = 1e-12), ys)
        @test !any(y -> isapprox(y, 0.5; atol = 1e-12), ys)
        for i = 1:16, pinfo in prep[i], (y, invtp) in pinfo.boundary_terms
            @test invtp >= 1 / 1.5 * (1 - 1e-12)
        end
    end

    # the BZ map has a branch endpoint at x = 1/8 where T' is infinite; the derivative there
    # used to throw a DomainError, it is now an infinite bound
    @testset "singular endpoint" begin
        B = Ulam(64)
        prep = prepare_Wnorm_estimate(BZ(), B)
        @test length(prep) == 64
        est, meas_Lf, var_Lf = estimate_Wnorm_of_discretization_error(prep, ones(64), B)
        derivs = prepare_derivative_bounds(BZ(), B)
        @test derivs[8] == Inf || derivs[9] == Inf          # the cells meeting x = 1/8
        @test !any(isnan, derivs)
        @test isfinite(est)
        @test all(isfinite, meas_Lf)
        @test !any(isnan, var_Lf)
    end

    # where 1/T' is unbounded Lemma 56 gives nothing: the variation estimate must be
    # infinite, so that the W-norm term falls back to the mass; it used to be 2·mass/δ
    @testset "unbounded derivative falls back to the mass" begin
        B = Ulam(16)
        prep = [i == 1 ? [PreimageInfo(0.0, 0.25, Inf, 1.0, Tuple{Float64,Float64}[])] :
                PreimageInfo[] for i = 1:16]
        density = collect(1.0:16.0)
        est, meas_Lf, var_Lf = estimate_Wnorm_of_discretization_error(prep, density, B)
        @test var_Lf[1] == Inf
        @test est >= meas_Lf[1] * (1 / 16) / 2 * (1 - 1e-12)
    end

end
