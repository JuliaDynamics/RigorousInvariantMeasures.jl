using RigorousInvariantMeasures

using LinearAlgebra: I, opnorm
using SparseArrays: sparse
using Test
using IntervalArithmetic
@testset "Norms of powers" begin

    import Random
    Random.seed!(1) #ensures tests are deterministic

    n = 9
    m = 10
    M = 0.2 * randn(n, n)
    R = 1e-14 * rand(n, n)

    LL = sparse(interval.(M .- R, M .+ R))

    Q = IntegralPreservingDiscretizedOperator(LL)

    U = [ones(1, n - 1); -Matrix(I, n - 1, n - 1)]

    @test norms_of_powers(Hat(9), Linf, m, Q, zeros(Interval{Float64}, 1)) ≈
          [opnorm(M^k * U, Inf) for k = 1:m]

    @test norms_of_powers(Ulam(9), L1, m, Q, zeros(Interval{Float64}, 1)) ≈
          [opnorm(M^k * U, 1) for k = 1:m]

    e = ones(n)
    f = adjoint(e) / n
    intf = map(interval, f)
    Q = NonIntegralPreservingDiscretizedOperator(LL, e, f)

    @test norms_of_powers(Hat(9), Linf, m, Q, intf) ≈
          [opnorm((M + e * (f - f * M))^k * U, Inf) for k = 1:m]
    @test norms_of_powers(Ulam(9), L1, m, Q, intf) ≈
          [opnorm((M + e * (f - f * M))^k * U, 1) for k = 1:m]

    @test refine_norms_of_powers([0.5, 1, 2, 0.001]) == [0.5, 0.25, 0.125, 0.001]
    @test refine_norms_of_powers([0.5, 1, 2, 1e-3], 8) ==
          [0.5, 0.25, 0.125, 1e-3, 0.5e-3, 0.25e-3, 0.125e-3, 1.0000000000000002e-6] # also tests correct rounding
    @test refine_norms_of_powers([2, 0.2, 0.1], 4) == [2, 0.2, 0.1, 0.04000000000000001]

    D = mod1_dynamic(x -> 2 * x)

    B = Ulam(16)
    Q = DiscretizedOperator(B, D)

    norms = powernormbounds(B, D, Q = Q)

    B_fine = Ulam(1024)
    Q_fine = DiscretizedOperator(B_fine, D)

    normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
    norms_fine = finepowernormbounds(B, B_fine, D, norms; normQ_fine = normQ_fine)

    # These two regressions are proportional to strong_weak_bound(B_fine), which was
    # corrected from n to the sharp 2n on 2026-11-07, so both doubled.
    @test norms_of_powers_from_coarser_grid(B_fine, B, D, norms, normQ_fine)[1:3] ==
          [193.0, 289.0, 337.0]

    @test norms_fine[1:11] == ones(11)

    @test norms_of_powers_dfly(B_fine, D, 3) == ([1024.0, 512.0, 256.0], [1.0, 1.0, 1.0])

    # the recursive two-grid bound: never above the one of Theorem 4.5 with the Lasota-Yorke
    # strong norms, and above the true norms of the fine operator on U⁰, computed densely as
    # the largest column difference over 2 (the unit ball of U⁰ in L¹ is the hull of
    # (e_i - e_j)/2); a nonlinear map, so that the strong norms matter
    @testset "recursive two-grid bound" begin
        D2 = mod1_dynamic(x -> 2x + 0.1 * RigorousInvariantMeasures.sinpi(2x))
        Bc, Bf = Ulam(32), Ulam(128)
        Qc, Qf = DiscretizedOperator(Bc, D2), DiscretizedOperator(Bf, D2)
        coarse = powernormbounds(Bc, D2; Q = Qc)
        nQ = opnormbound(Bf, weak_norm(Bf), Qf)
        m = min(length(coarse), 40)
        old = norms_of_powers_from_coarser_grid(Bf, Bc, D2, coarse[1:m], nQ)
        rec = norms_of_powers_from_coarser_grid_recursive(Bf, Bc, D2, coarse[1:m], nQ)
        @test all(rec .<= old)
        @test rec[end] < old[end]
        n = length(Bf)
        L = mid.(Matrix(Qf.L))
        P = Matrix{Float64}(I, n, n)
        for k = 1:m
            P = L * P
            truenorm = maximum(sum(abs.(P[:, i] .- P[:, j])) / 2 for i = 1:n for j = i+1:n)
            @test truenorm <= rec[k] + 1e-12
        end
        # on this pair of grids no fine norm of a power is below one and finepowernormbounds
        # raises an error rather than returning a vector whose tail is not summable
        @test_throws ErrorException finepowernormbounds(Bc, Bf, D2, coarse[1:m])
        # the default normQ_fine used to read weak_norm(B_fine)DiscretizedOperator(...),
        # a product of a type and an operator, which threw; on the doubling map's grids the
        # bound contracts and the vector is extended until its last entry is below the target
        nf = finepowernormbounds(B, B_fine, D, norms)
        @test length(nf) >= length(norms) && nf[end] < 1e-8
    end

end
