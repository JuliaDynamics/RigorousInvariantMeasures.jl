@testset "NoiseKernel regressions" begin

    using Test
    using IntervalArithmetic
    using RigorousInvariantMeasures
    import RigorousInvariantMeasures:
        UniformKernelUlamPeriodic, UniformKernelUlamReflecting, gamma

    # A density that is not symmetric and piles mass against the left edge, so
    # that a boundary condition applied wrongly shows up in the total mass.
    edge_heavy(k) = [exp(-20 * ((i - 0.5) / k - 0.05)^2) for i = 1:k]

    @testset "reflecting kernel is Markov" begin
        # Reflecting used to throw a BoundsError, and reflecting the input
        # rather than folding the output loses O(1/k) of the mass.
        for k in (64, 256, 1024)
            B = Ulam(k)
            v = edge_heavy(k)
            for K in (UniformNoiseUlam(0.05, B, :reflecting),
                      UniformKernelUlamReflecting(B, max(1, k ÷ 40)))
                w = K * copy(v)
                @test sum(w) ≈ sum(v) rtol = 1e-13
            end
        end
    end

    @testset "reflecting kernel accepts intervals" begin
        # There was no interval method for the reflecting branch at all.
        k = 256
        B = Ulam(k)
        v = edge_heavy(k)
        K = UniformNoiseUlam(0.05, B, :reflecting)
        w = K * copy(v)
        wi = K * interval.(v)
        @test all(in_interval(w[i], wi[i]) for i = 1:k)
    end

    @testset "Kahan kernel enclosure is not vacuous" begin
        # γₖ was hard-coded to 1.0, which made the radius ‖v‖₁/n: about 9.75 on
        # entries of order one at k = 1024. It should sit near the rounding
        # level instead, and it must still contain the float result.
        for k in (1024, 4096)
            B = Ulam(k)
            l = k ÷ 40
            K = UniformKernelUlamPeriodic(B, l)
            v = edge_heavy(k)
            w = K * copy(v)
            wi = K * interval.(v)
            @test maximum(radius.(wi)) < 1e-9
            @test all(in_interval(w[i], wi[i]) for i = 1:k)
            # and it must not be so tight that it stops being an enclosure
            @test maximum(radius.(wi)) >= gamma(Float64, 2l + 2) * sum(abs, v) / (2l + 1)
        end
    end

    @testset "scratch buffers carry the bound type" begin
        # `w`/`z` and the Kahan scratch used to be fixed to Vector{Float64},
        # which boxed every access and would silently narrow a BigFloat input.
        B = Ulam(64)
        @test eltype(UniformNoiseUlam(0.05, B).w) === Float64
        @test eltype(UniformNoiseUlam(BigFloat, 0.05, B).w) === BigFloat
        @test eltype(UniformKernelUlamPeriodic(B, 3).scratch_ext) === Float64
        @test eltype(UniformKernelUlamPeriodic(BigFloat, B, 3).scratch_ext) === BigFloat
        setprecision(256) do
            K = UniformKernelUlamPeriodic(BigFloat, B, 3)
            v = [BigFloat(1) / BigFloat(i) for i = 1:64]
            @test precision((K*copy(v))[1]) == 256
        end
    end

    @testset "old periodic kernel does not allocate per entry" begin
        # `sum(M.v .* h)` allocated a length-n temporary on each of the k
        # iterations: 222 MB per application at k = 16384.
        k = 4096
        B = Ulam(k)
        K = UniformNoiseUlam(0.05, B)
        v = edge_heavy(k)
        K * copy(v)                       # warm up
        allocated = @allocated K * copy(v)
        @test allocated < 20 * k          # a few vectors, not k of them
    end

end
