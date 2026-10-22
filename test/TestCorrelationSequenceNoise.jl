@testset "Noisy Ulam operator: correlations and two-grid norms" begin

    using Test
    using IntervalArithmetic
    using RigorousInvariantMeasures
    using SparseArrays, LinearAlgebra
    import RigorousInvariantMeasures: correlation_sequence_noise, norms_of_powers_noise,
        norms_of_powers_from_coarser_grid_noise, ProjectedFunction, L1, mid

    @testset "correlation_sequence_noise encloses a 256-bit iteration" begin
        k = 256
        B = Ulam(k)
        Q = DiscretizedOperator(B, BZ())
        ξ = 37 // 1024
        NK = UlamProjectedUniformKernelReflecting(B, ξ)
        φ = [(i - 0.5) / k for i = 1:k]
        ϕ = ProjectedFunction(B, φ, 1.0)
        g = [Base.sinpi(3 * (i - 0.5) / k) + 0.1 for i = 1:k]
        nlags = 200
        C = correlation_sequence_noise(ϕ, Q, NK, g, nlags)

        # the same products in 256-bit arithmetic: the midpoint matrix of Q, the kernel
        # built at BigFloat, the pairing in BigFloat
        setprecision(BigFloat, 256) do
            M = BigFloat.(mid.(Q.L))
            NKbig = UlamProjectedUniformKernel(Val(:reflecting), BigFloat, B, ξ)
            v = BigFloat.(g)
            ok = true
            for n = 1:nlags
                c = sum(BigFloat.(φ) .* v) / k
                ok &= BigFloat(inf(C[n])) <= c <= BigFloat(sup(C[n]))
                v = NKbig * (M * v)
            end
            @test ok
        end
        # the radius grows at most linearly: each step adds the same rounding constant
        # times ‖v‖, and ‖v‖ does not grow, the operator being Markov
        @test radius(C[end]) <= 2 * nlags * radius(C[2]) + 1e-15
        @test maximum(radius.(C)) < 1e-8
    end

    @testset "g_error enters through ‖φ‖_∞" begin
        k = 128
        B = Ulam(k)
        Q = DiscretizedOperator(B, BZ())
        NK = UlamProjectedUniformKernelReflecting(B, 9 // 128)
        ϕ = ProjectedFunction(B, ones(k), 2.0)
        g = ones(k)
        C0 = correlation_sequence_noise(ϕ, Q, NK, g, 5)
        C1 = correlation_sequence_noise(ϕ, Q, NK, g, 5; g_error = 1e-3)
        @test all(radius(C1[n]) >= radius(C0[n]) + 2.0 * 1e-3 * (1 - 1e-12) for n = 1:5)
    end

    @testset "two-grid norms bound the true norms of the fine operator" begin
        # coarse 64 and fine 256 cells, the same ξ: s = 21 fine cells, 5.25 coarse ones
        ξ = 21 // 256
        Bc, Bf = Ulam(64), Ulam(256)
        Qc, Qf = DiscretizedOperator(Bc, BZ()), DiscretizedOperator(Bf, BZ())
        NKc = UlamProjectedUniformKernelReflecting(Bc, ξ)
        NKf = UlamProjectedUniformKernelReflecting(Bf, ξ)
        m = 30
        coarse = norms_of_powers_noise(Bc, L1, m, Qc, NKc, integral_covector(Bc))
        fine = norms_of_powers_from_coarser_grid_noise(Bf, Bc, Qf, NKf, coarse)

        # true ‖L̃ᵏ|_V‖_{L¹} of the fine midpoint operator: the unit ball of V is the
        # convex hull of (e_i - e_j)/2, so it is the largest column difference over 2
        n = length(Bf)
        L = reduce(hcat, [NKf * (mid.(Qf.L) * Float64.(1:n .== j)) for j = 1:n])
        P = Matrix(1.0I(n))
        for kk = 1:m+1
            P = L * P
            truenorm = maximum(sum(abs.(P[:, i] .- P[:, j])) / 2 for i = 1:n for j = i+1:n)
            @test truenorm <= fine[kk] + 1e-12
        end
    end
end
