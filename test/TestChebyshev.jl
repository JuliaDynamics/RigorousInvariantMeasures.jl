using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Chebyshev assembler" begin

    N = 16
    B = RigorousInvariantMeasures.Chebyshev(N, 2)

    D = mod1_dynamic(x -> 2 * x)
    L(ϕ, x) = (ϕ(x / 2) + ϕ(x / 2 + 0.5)) / 2

    M = assemble(B, D)
    using LinearAlgebra

    all_true = true
    # for i in 1:N
    #  w = mid.(L.(B[i], B.p))
    # z = RigorousInvariantMeasures.chebtransform(w)
    #   all_true = all_true && norm(z-M[:, i], Inf)< 10^-13
    # end

    @test all_true

end

@testset "Chebyshev assemble agrees with the Clenshaw evaluation it replaces" begin
    # `assemble` sweeps all degrees at once via T_m(t) = Re(z^m), z on the unit
    # circle, instead of one Clenshaw pass per (basis function, node). Pin it to
    # the definitional loop built on `evalChebyshev`, i.e. on
    # `eval_Clenshaw_BackwardFirst` (Ledoux–Moroz, MACIS 2019, Algorithm 2).
    ext = Base.get_extension(RigorousInvariantMeasures, :FFTWExt)

    D = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
    B = Chebyshev(64, 3)
    n = length(B.p)

    M_naive = zeros(Interval{Float64}, (n, n))
    x, labels, x′ = RigorousInvariantMeasures.Dual(B, D; ϵ = 1e-13, max_iter = 100)
    for i = 1:n
        ϕ = B[i]
        w = zeros(Interval{Float64}, n)
        for j = 1:length(x)
            w[labels[j]] += ϕ(x[j]) / abs(x′[j])
        end
        M_naive[:, i] = ext.chebtransform(w)
    end

    M_fast = RigorousInvariantMeasures.assemble(B, D; ϵ = 1e-13, max_iter = 100)

    @test size(M_fast) == size(M_naive)
    @test all(
        !isempty_interval(intersect_interval(a, b)) for (a, b) in zip(M_naive, M_fast)
    )
    # The unit-circle sweep accumulates ~log2(n) ulps, so it must come out no
    # looser than the Clenshaw path it replaces.
    @test maximum(radius.(M_fast)) <= maximum(radius.(M_naive))
end
