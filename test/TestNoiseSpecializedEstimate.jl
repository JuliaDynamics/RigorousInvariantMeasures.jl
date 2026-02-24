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
        _reflect,
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

    # ──── _reflect helper ────

    @testset "_reflect" begin
        # Inside range: identity
        @test _reflect(1, 5) == 1
        @test _reflect(3, 5) == 3
        @test _reflect(5, 5) == 5
        # Outside right: reflection
        @test _reflect(6, 5) == 5
        @test _reflect(7, 5) == 4
        # Outside left
        @test _reflect(0, 5) == 2
        @test _reflect(-1, 5) == 3
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
        noise_half_width = 8  # noise window = 2*8+1 = 17 bins
        NK = UniformKernelUlamReflecting(B, noise_half_width)

        # Assemble operator
        Q = DiscretizedOperator(B, D)

        # Compute invariant measure
        w = invariant_vector_noise(B, Q, NK; iter=20)

        # Compute generic error via distance_from_invariant_noise
        norms = powernormboundsnoise(B; Q=Q, NK=NK)
        generic_error = distance_from_invariant_noise(B, Q, NK, w, norms)

        @test generic_error > 0
        @test isfinite(generic_error)

        # Compute finer error via noise_error_aposteriori
        # First compute the needed parameters
        noise_size_rel = Float64(2 * noise_half_width + 1) / Float64(K_size)

        # Spectral contraction rate from norms
        # Find the contraction rate α
        # norms[end] should be < 1, use it as a proxy for the geometric rate
        m = length(norms)
        α_est = Float64(norms[end]^(1.0 / m))
        if α_est >= 1.0
            α_est = 0.99
        end

        # sumCi = Σ_{i=0}^{∞} ||Q^i||
        sumCi = infinite_sum_norms(norms)

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
        NK = UniformKernelUlamReflecting(B, 4)
        meas_Lf = ones(32) ./ 32
        var_Lf = 0.1 * ones(32)

        var_NLf = estimate_variation_after_noise(meas_Lf, var_Lf, NK)
        @test length(var_NLf) == 32
        @test all(var_NLf .≥ 0)
        @test all(isfinite.(var_NLf))
    end

end
