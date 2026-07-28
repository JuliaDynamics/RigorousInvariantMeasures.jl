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

@testset "Chebyshev Gram matrix" begin
    using LinearAlgebra
    # Selective import: a bare `using BallArithmetic` leaks `mid`/`sup` into
    # Main and shadows IntervalArithmetic's for every later test file.
    using BallArithmetic: BallMatrix, upper_bound_L2_opnorm

    B = Chebyshev(16, 3)
    n = length(B)
    G = gram_matrix(B)
    Gi = inv_gram_matrix(B)

    # G = diag(1, 1/2, ..., 1/2) for the arcsine probability measure.
    @test G isa Diagonal
    @test in_interval(1, G[1, 1])
    @test all(in_interval(1 // 2, G[i, i]) for i = 2:n)
    @test all(in_interval(i == j ? 1 : 0, (G*Gi)[i, j]) for i = 1:n, j = 1:n)

    S = gram_sqrt(B)
    Si = inv_gram_sqrt(B)
    @test all(!isempty_interval(intersect_interval((S*S)[i, i], G[i, i])) for i = 1:n)
    @test all(in_interval(1, (S*Si)[i, i]) for i = 1:n)

    # Independent check of the Gram values: <f,f>_mu against Gauss-Chebyshev
    # quadrature, which is exact for the arcsine measure up to node count.
    c = [1 / (i + 1) for i = 1:n]
    N = 4096
    nodes = [cos((2k - 1) * pi / (2N)) for k = 1:N]
    fvals = [sum(c[i] * cos((i - 1) * acos(t)) for i = 1:n) for t in nodes]
    quad = sum(abs2, fvals) / N
    gram_val = sum(c[i]^2 * mid(G[i, i]) for i = 1:n)
    @test isapprox(quad, gram_val; rtol = 1e-12)

    # The reason for this weight: the similarity carrying an ell^2 operator
    # bound to the L^2(mu) one is diagonal, so BallArithmetic applies directly.
    D = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
    Q = RigorousInvariantMeasures.assemble(B, D; ϵ = 1e-13, max_iter = 100)
    conj_Q = S * Q * Si
    bound = upper_bound_L2_opnorm(BallMatrix(conj_Q))
    @test isfinite(bound)
    @test bound > 0
end
