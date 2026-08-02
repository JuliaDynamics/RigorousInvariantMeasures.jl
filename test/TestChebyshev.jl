using RigorousInvariantMeasures
using IntervalArithmetic

@testset "chebtransform recovers known Chebyshev coefficients" begin
    # This used to be `all_true = true; @test all_true` with the real comparison
    # commented out, which is how a factor-of-2(N-1) normalization error in
    # `chebtransform` survived: it divided by n on top of interval_fft's own
    # 1/(2n).  Ground truth here is a polynomial with known coefficients.
    ext = Base.get_extension(RigorousInvariantMeasures, :FFTWExt)

    for k in (8, 16)
        B = RigorousInvariantMeasures.Chebyshev(k, 2)
        N = length(B)
        t = 2 .* mid.(B.p) .- 1                       # Chebyshev points in [-1,1]

        # f = 3*T_0 - 2*T_1 + 5*T_4
        f(x) = 3 - 2 * x + 5 * cos(4 * acos(clamp(x, -1, 1)))
        w = [interval(f(tt)) for tt in t]

        expected = zeros(N)
        expected[1], expected[2], expected[5] = 3, -2, 5

        got = ext.chebtransform(w)
        @test all(in_interval(expected[i], got[i]) for i = 1:N)
        @test maximum(abs.(mid.(got) .- expected)) < 1e-13
    end
end

@testset "Chebyshev assembler preserves the integral covector" begin
    # The transfer operator satisfies ∫Lf = ∫f, so the integral covector is a
    # left eigenvector of eigenvalue 1.  With the old normalization the
    # eigenvalue came out as 1/(2(N-1)) instead, which is what flagged the bug.
    D = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
    for k in (8, 16)
        B = RigorousInvariantMeasures.Chebyshev(k, 3)
        Q = mid.(assemble(B, D; ϵ = 1e-13, max_iter = 100))
        v = vec(mid.(collect(integral_covector(B))))
        @test maximum(abs.(transpose(Q) * v .- v)) < 1e-2
    end
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

@testset "Lebesgue Gram matrix, inverse and average-zero restriction" begin
    using LinearAlgebra
    using BallArithmetic: BallMatrix

    B = Chebyshev(16, 3)
    n = length(B)

    GL = gram_matrix(B; measure = :lebesgue)
    @test in_interval(1, GL[1, 1])              # ∫ 1 dx
    @test in_interval(1 // 3, GL[2, 2])         # ∫ (2x-1)^2 dx
    # entries with (i-1)+(j-1) odd vanish
    @test all(in_interval(0, GL[i, j]) for i = 1:n, j = 1:n if isodd(i + j))

    # The identity the restriction rests on: the integral covector is the first
    # column of the Lebesgue Gram matrix, because T_0 = 1.
    v = collect(integral_covector(B))
    @test all(!isempty_interval(intersect_interval(v[i], GL[i, 1])) for i = 1:n)

    GLi = inv_gram_matrix(B; measure = :lebesgue)
    P = GL * GLi
    @test all(in_interval(i == j ? 1 : 0, P[i, j]) for i = 1:n, j = 1:n)

    c_leb, C_n = l2_measure_conversion_bounds(B)
    @test c_leb >= sqrt(pi / 2)                 # uniform direction
    @test C_n > 1                               # finite-dimensional direction
    @test C_n < 2 * sqrt(n)

    # Restriction: (n-1)x(n-1) block, and the spectral radius must be stable in n
    # (it is the second eigenvalue of the transfer operator).
    D = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
    ρs = Float64[]
    for k in (12, 24)
        Bk = Chebyshev(k, 3)
        Q = RigorousInvariantMeasures.assemble(Bk, D; ϵ = 1e-13, max_iter = 100)
        blk, chol = gram_restrict_to_average_zero(Bk, BallMatrix(Q))
        @test chol.success
        @test size(blk.c) == (length(Bk) - 1, length(Bk) - 1)
        @test maximum(blk.r) < 1e-10
        push!(ρs, maximum(abs.(eigvals(blk.c))))
    end
    @test abs(ρs[1] - ρs[2]) < 1e-5             # converged second eigenvalue
    @test ρs[1] < 1                             # spectral gap on average-zero
end

@testset "Chebyshev norm parameterization" begin
    RIM = RigorousInvariantMeasures

    # Default is W^{k,1} strong / L2 weak, as for FourierAnalytic.
    B = Chebyshev(16, 3)
    @test typeof(strong_norm(B)) == W{3,1}
    @test weak_norm(B) == L2
    @test aux_norm(B) == L1

    # The Taylor-Crush C1 path is untouched.
    @test weak_norm(Chebyshev(16, 3, C1)) == C1
    @test length(Chebyshev(16, 3)) == length(Chebyshev(16, 3, C1))

    # The L2 projection error comes from the W^{k,1} coefficient decay via
    # Theorems 3.12/3.13; since dx is a probability measure ||.||_L2 <= ||.||_inf
    # with constant 1, so it coincides with the aux (C0) bound.
    for n in (16, 32)
        Bn = Chebyshev(n, 3)
        @test RIM.weak_projection_error(Bn) == RIM.aux_normalized_projection_error(Bn)
    end

    # ...and it decays like n^{-ν}: refining by 2 should gain about 2^3 = 8.
    e16 = RIM.weak_projection_error(Chebyshev(16, 3))
    e32 = RIM.weak_projection_error(Chebyshev(32, 3))
    @test 5 < e16 / e32 < 12

    @test RIM.aux_weak_bound(Chebyshev(16, 3)) == 1.0
    @test RIM.aux_weak_bound(Chebyshev(16, 3, C1)) == 1.0
    @test RIM.bound_weak_norm_from_linalg_norm(Chebyshev(16, 3)) == (1.0, 0.0)
    @test RIM.weak_by_strong_and_aux_bound(Chebyshev(16, 3)) == (1.0, 0.0)

    # All interface constants finite and positive for every weak norm.
    for Bx in (Chebyshev(16, 3), Chebyshev(16, 3, C1))
        for g in (RIM.weak_projection_error, RIM.aux_normalized_projection_error,
                  RIM.strong_weak_bound, RIM.bound_linalg_norm_L1_from_weak,
                  RIM.bound_linalg_norm_L∞_from_weak)
            v = g(Bx)
            @test isfinite(v) && v > 0
        end
    end
end

@testset "Bernstein ellipses" begin
    RIM = RigorousInvariantMeasures

    # bernstein_parameter inverts bernstein_point exactly.
    for ρ in (1.1, 1.3, 2.0, 4.0), θ in (0.0, 0.137, 0.25, 0.5, 0.9)
        z = bernstein_point(interval(ρ), interval(θ))
        @test in_interval(ρ, bernstein_parameter(z))
    end

    # [-1,1] is the degenerate ellipse ρ = 1.
    for x in (-1.0, -0.4, 0.0, 0.7, 1.0)
        @test in_interval(1, bernstein_parameter(complex(interval(x), interval(0.0))))
    end

    # Ground truth: the Chebyshev polynomials map E_ρ onto E_{ρ^m}.
    T2(z) = 2z^2 - 1
    T3(z) = 4z^3 - 3z
    for ρ in (1.2, 1.5, 2.0), (m, f) in ((2, T2), (3, T3))
        got = bernstein_expansion(f, ρ; n = 2048)
        @test got <= ρ^m + 1e-12          # it is a rigorous LOWER bound
        @test got > 0.95 * ρ^m            # ...and a tight one
    end

    # The gap is pure discretization: refining by 4 should shrink it ~4x or more.
    gaps = [2.25 - bernstein_expansion(T2, 1.5; n = n) for n in (256, 1024, 4096)]
    @test all(gaps .> 0)
    @test gaps[2] < gaps[1] / 3
    @test gaps[3] < gaps[2] / 3

    # Expansion test: T_2 expands every ellipse, the identity expands none.
    for ρ in (1.2, 2.0)
        expands, ρ_img = expands_bernstein_ellipse(T2, ρ; n = 2048)
        @test expands
        @test ρ_img > ρ
        @test !first(expands_bernstein_ellipse(identity, ρ; n = 256))
    end

    # The [0,1] -> [-1,1] conjugation: x -> 2x mod 1 becomes t -> 2t+1 on the
    # first branch, and doubling the angle doubles the ellipse parameter.
    f01 = x -> 2 * x
    @test to_symmetric_interval(f01)(interval(0.0)) == interval(1.0)

    # Eρ basis: geometric projection error, decaying like ρ^{-n}.
    for ρ in (1.5, 2.0)
        B16, B32 = Chebyshev(16, Eρ(ρ)), Chebyshev(32, Eρ(ρ))
        @test typeof(strong_norm(B16)) == Eρ
        # B = 0 in the analytic DFLY, so the auxiliary norm is multiplied by zero
        # and the weak norm is free to be the one the approximation error uses.
        @test weak_norm(B16) == L2
        @test aux_norm(B16) == L1
        e16 = RIM.weak_projection_error(B16)
        e32 = RIM.weak_projection_error(B32)
        @test 0 < e32 < e16
        @test isapprox(e16 / e32, ρ^(length(B32) - length(B16)); rtol = 1e-8)
        @test RIM.weak_by_strong_and_aux_bound(B16) == (1.0, 0.0)
        @test isfinite(RIM.strong_weak_bound(B16))
    end
end

@testset "L2(dμ) weak norm with L1(dμ) auxiliary" begin
    RIM = RigorousInvariantMeasures

    Bμ = Chebyshev(16, 3, L2μ)
    @test weak_norm(Bμ) == L2μ
    @test aux_norm(Bμ) == L1μ           # both against μ
    @test aux_norm(Chebyshev(16, 3)) == L1

    # Matching the measures makes aux_weak_bound a plain Cauchy-Schwarz 1;
    # pairing L1(dx) with L2(dμ) instead would cost π/(2√2).
    @test RIM.aux_weak_bound(Bμ) == 1.0

    # Parseval constants: no conversion factor, unlike the Lebesgue weak norm,
    # which pays C_n for the same quantities.
    n = length(Bμ)
    @test RIM.bound_linalg_norm_L1_from_weak(Bμ) ≈ sqrt(2n) rtol = 1e-12
    @test RIM.bound_linalg_norm_L∞_from_weak(Bμ) ≈ sqrt(2) rtol = 1e-12
    @test RIM.bound_linalg_norm_L1_from_weak(Bμ) <
          RIM.bound_linalg_norm_L1_from_weak(Chebyshev(16, 3))
    @test RIM.bound_weak_norm_from_linalg_norm(Bμ) == (1.0, 0.0)
    @test RIM.weak_by_strong_and_aux_bound(Bμ) == (1.0, 0.0)

    # μ is a probability measure, so the sup-norm projection bound carries over
    # unchanged: same number as for L2(dx).
    @test RIM.weak_projection_error(Bμ) == RIM.weak_projection_error(Chebyshev(16, 3))

    # The analytic basis: geometric projection error under either weak norm.
    BE = Chebyshev(16, Eρ(1.5))
    @test weak_norm(BE) == L2
    @test isfinite(RIM.weak_projection_error(BE)) && RIM.weak_projection_error(BE) > 0
end
