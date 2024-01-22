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

    LL = sparse(interval_from_midpoint_radius.(M, R))

    Q = IntegralPreservingDiscretizedOperator(LL)

    U = [ones(1, n - 1); -Matrix(I, n - 1, n - 1)]

    @test norms_of_powers(Hat(9), Linf, m, Q, zeros(Interval{Float64}, 1)) ≈
          [opnorm(M^k * U, Inf) for k = 1:m]

    @test norms_of_powers(Ulam(9), L1, m, Q, zeros(Interval{Float64}, 1)) ≈
          [opnorm(M^k * U, 1) for k = 1:m]

    e = ones(n)
    f = adjoint(e) / n
    intf = map(Interval, f)
    Q = NonIntegralPreservingDiscretizedOperator(LL, e, f)

    @test norms_of_powers(Hat(9), Linf, m, Q, intf) ≈
          [opnorm((M + e * (f - f * M))^k * U, Inf) for k = 1:m]
    @test norms_of_powers(Ulam(9), L1, m, Q, intf) ≈
          [opnorm((M + e * (f - f * M))^k * U, 1) for k = 1:m]

    @test refine_norms_of_powers([0.5, 1, 2, 0.001]) == [0.5, 0.25, 0.125, 0.001]
    @test refine_norms_of_powers([0.5, 1, 2, 1e-3], 8) ==
          [0.5, 0.25, 0.125, 1e-3, 0.5e-3, 0.25e-3, 0.125e-3, 1.0000000000000002e-6] # also tests correct rounding
    @test refine_norms_of_powers([2, 0.2, 0.1], 4) == [2, 0.2, 0.1, 0.04000000000000001]

end
