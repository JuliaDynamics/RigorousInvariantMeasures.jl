@testset "UlamProjectedUniformKernel" begin

    using Test
    using IntervalArithmetic
    using RigorousInvariantMeasures
    import RigorousInvariantMeasures: UlamProjectedUniformKernel, UniformKernelUlamPeriodic,
        opnormbound, opradius, nonzero_per_row, dfly, L1, TotalVariation

    # The entry at offset d of π N_ξ π in cell units, computed from the definition and not
    # from the tent primitive the implementation uses: (1/s) ∫_0^1 |[x+d-s/2, x+d+s/2] ∩ [0,1]| dx,
    # the integrand being piecewise linear in x and integrated exactly between its breakpoints.
    function reference_weight(s::Rational{BigInt}, d::Integer)
        overlap(x) = max(zero(x), min(x + d + s / 2, one(x)) - max(x + d - s / 2, zero(x)))
        pts = sort(unique([zero(s), one(s),
            filter(p -> 0 < p < 1, [-d + s / 2, 1 - d + s / 2, -d - s / 2, 1 - d - s / 2])...]))
        integral = sum((pts[i+1] - pts[i]) * (overlap(pts[i]) + overlap(pts[i+1])) / 2
                       for i = 1:length(pts)-1)
        return integral / s
    end

    # the dense matrix of a kernel, column by column
    dense(K, k) = reduce(hcat, [K * Float64.(1:k .== j) for j = 1:k])

    @testset "row against the definition, s = $s" for s in
            (Rational{BigInt}(5), Rational{BigInt}(5, 2), Rational{BigInt}(13, 4),
             Rational{BigInt}(2), Rational{BigInt}(7), Rational{BigInt}(173, 8))
        k = 64
        B = Ulam(k)
        K = UlamProjectedUniformKernelPeriodic(B, s / k)
        M = dense(K, k)
        i = k ÷ 2                     # an interior row, away from the wrap-around
        for j = 1:k
            w = reference_weight(s, i - j)
            @test abs(M[i, j] - Float64(w)) <= 4eps() * max(Float64(w), 1 / Float64(s))
        end
        @test (2K.m + 1) + 2 * (K.c₁ + K.c₂) == Float64(s)
    end

    @testset "odd support: (1/8, 7/8, 1, …, 1, 7/8, 1/8)/(2l+1)" begin
        for l in (1, 2, 86)
            K = UlamProjectedUniformKernelPeriodic(Ulam(1024), (2l + 1) // 1024)
            @test (K.m, K.c₁, K.c₂) == (l - 1, 7 / 8, 1 / 8)
        end
    end

    @testset "Markov and symmetric, $bc" for bc in (:periodic, :reflecting)
        k = 128
        K = UlamProjectedUniformKernel(Val(bc), Ulam(k), 21 // 4 / k)
        M = dense(K, k)
        @test maximum(abs.(sum(M, dims = 1) .- 1)) <= 1e-14
        @test maximum(abs.(sum(M, dims = 2) .- 1)) <= 1e-14
        @test M ≈ M' atol = 1e-16
        @test all(M .>= 0)
    end

    @testset "distance from the flat window is δ/(2ξ)" begin
        k, l = 256, 4
        P = dense(UlamProjectedUniformKernelPeriodic(Ulam(k), (2l + 1) // k), k)
        W = dense(UniformKernelUlamPeriodic(Ulam(k), l), k)
        @test maximum(sum(abs.(P - W), dims = 1)) ≈ 1 / (2(2l + 1)) rtol = 1e-12
    end

    @testset "interval product encloses the exact one" begin
        k = 256
        s = Rational{BigInt}(37, 4)
        K = UlamProjectedUniformKernelReflecting(Ulam(k), s / k)
        v = [exp(-20 * ((i - 0.5) / k - 0.05)^2) * (1 + Base.sinpi(7i / k)) for i = 1:k]
        # exact product in rationals, the reflected weights folded into the matrix
        exact = zeros(Rational{BigInt}, k)
        D = K.m + 2
        for i = 1:k, d = -D:D
            j = RigorousInvariantMeasures.get_idx(Val(:reflecting), i - d, k)
            exact[i] += reference_weight(s, d) * Rational{BigInt}(v[j])
        end
        w = K * interval.(v)
        @test all(BigFloat(inf(w[i])) <= exact[i] <= BigFloat(sup(w[i])) for i = 1:k)
        @test maximum(radius.(w)) < 1e-14
        wf = K * copy(v)
        @test all(in_interval(wf[i], w[i]) for i = 1:k)
    end

    @testset "constants" begin
        k = 512
        K = UlamProjectedUniformKernelReflecting(Ulam(k), 11 // k)
        @test opnormbound(Ulam(k), L1, K) == 1.0
        @test opradius(L1, K) == 0.0
        @test nonzero_per_row(K) == 2K.m + 5
        @test dfly(TotalVariation, L1, K) == (0.0, 2k / 11)
    end

    @testset "rejected input" begin
        @test_throws ArgumentError UlamProjectedUniformKernelPeriodic(Ulam(64), 1 // 64)
        @test_throws ArgumentError UlamProjectedUniformKernelPeriodic(Ulam(64), 1 // 3)
        @test_throws ArgumentError UlamProjectedUniformKernel(Val(:open), Ulam(64), 5 // 64)
    end
end
