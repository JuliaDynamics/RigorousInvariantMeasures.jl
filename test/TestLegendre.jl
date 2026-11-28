using RigorousInvariantMeasures, IntervalArithmetic, BallArithmetic, LinearAlgebra, Test
using FFTW, GenericFFT  # assemble(::Chebyshev, …) in Float64 and BigFloat, for route 2

@testset "Legendre basis" begin
    RIM = RigorousInvariantMeasures
    IA = IntervalArithmetic
    _in(x, b) = abs(x - b.c) <= b.r                       # x ∈ ball b
    inball(M, B::BallMatrix) = all(abs.(M .- B.c) .<= B.r)

    @testset "changes of basis" begin
        N = 40
        L = RIM.legendre_to_chebyshev_rational(N)
        C = RIM.chebyshev_to_legendre_rational(N)
        @test C * L == I && L * C == I
        @test istriu(L) && istriu(C)
        @test all(iszero(L[k+1, n+1]) for n = 0:N, k = 0:N if isodd(n - k))
        # closed forms, Λ(z) = Γ(z+1/2)/Γ(z+1), via Λ(z+1) = Λ(z)(z+1/2)/(z+1)
        Λ(z) = (λ = z == floor(z) ? sqrt(big(π)) : 2 / sqrt(big(π));
                for w in (z == floor(z) ? 0 : 1/2):1:z-1; λ *= (w + 1/2) / (w + 1); end; λ)
        for n = 0:N, k = 0:n
            isodd(n - k) && continue
            l = (2 - (k == 0)) / big(π) * Λ((n - k) / 2) * Λ((n + k) / 2)
            @test Float64(L[k+1, n+1]) ≈ Float64(l) rtol = 1e-12
            if n >= 1
                c = k == n ? sqrt(big(π)) / (2Λ(n)) :
                    -n * (k + 1/2) * Λ((n - k - 2) / 2) * Λ((n + k - 1) / 2) / ((n + k + 1) * (n - k))
                @test Float64(C[k+1, n+1]) ≈ Float64(c) rtol = 1e-12
            end
        end
        Ci = BallMatrix(RIM._rational_interval.(Float64, C))
        Li = BallMatrix(RIM._rational_interval.(Float64, L))
        @test inball(Matrix(1.0I, N + 1, N + 1), Ci * Li)
        setprecision(BigFloat, 512) do
            Cb = BallMatrix(RIM._rational_interval.(BigFloat, C))
            @test maximum(Cb.r) < 1e-140
        end
    end

    @testset "evaluation" begin
        B = Legendre(5)
        for x in (0.0, 0.3, 0.71, 1.0)
            t = 2x - 1
            @test RIM.evaluate(B, 3, x) ≈ sqrt(5) * (3t^2 - 1) / 2
            @test RIM.evaluate(B, 4, x) ≈ sqrt(7) * (5t^3 - 3t) / 2
            @test in_interval(sqrt(7) * (5t^3 - 3t) / 2, RIM.evaluate(B, 4, interval(x)))
        end
        @test integral_covector(B) * one_vector(B) == [1.0]
        @test_throws BoundsError RIM.evaluate(B, 7, 0.5)
    end

    f = x -> 2x + x * (1 - x) / 2
    D = mod1_dynamic(f; full_branch = true)

    @testset "Galerkin matrix of the identity is I" begin
        Did = PwMap([x -> x], [interval(0), interval(1)])
        A = galerkin_matrix(Legendre(10), Did; nsub = 16)
        @test inball(Matrix(1.0I, 11, 11), A)
        @test maximum(A.r) < 1e-12
    end

    @testset "Lanford Galerkin matrix, route 1" begin
        N = 12
        A = galerkin_matrix(Legendre(N), D)
        @test maximum(A.r) < 1e-12
        @test inball(permutedims([1.0; zeros(N)]), BallMatrix(A.c[1:1, :], A.r[1:1, :]))
        # non-rigorous reference: composite Simpson on the two branches
        c = IA.mid(endpoints(D)[2])
        e(k, x) = sqrt(2k + 1) * RIM._legendre_P(k, 2x - 1)
        function simpson(g, a, b; n = 4000)
            h = (b - a) / n
            return h / 3 * (g(a) + g(b) + 4sum(g(a + (2i - 1)h) for i = 1:n÷2) +
                            2sum(g(a + 2i * h) for i = 1:n÷2-1))
        end
        for (j, k) in ((1, 0), (2, 3), (5, 5), (12, 7), (0, 12))
            ref = simpson(x -> e(k, x) * e(j, f(x)), 0.0, c) +
                  simpson(x -> e(k, x) * e(j, f(x) - 1), c, 1.0)
            @test abs(ref - A.c[j+1, k+1]) < A.r[j+1, k+1] + 1e-9
        end
        # the BigFloat output carries the working precision
        # with the branch point enclosed at 512 bits as well
        Ab = setprecision(BigFloat, 512) do
            Db = mod1_dynamic(f; full_branch = true, T = BigFloat, ϵ = big(2.0)^-480)
            galerkin_matrix(Legendre(N; T = BigFloat), Db; precision = 512, order = 120)
        end
        @test maximum(Float64.(Ab.r)) < 1e-60
        @test all(abs.(Float64.(Ab.c) - A.c) .<= A.r .+ 1e-30)
    end

    @testset "geometry, compatibility, route 2, change of basis" begin
        @test_throws ArgumentError EllipseGeometry(3.5, 3.4, 4.0)
        G = EllipseGeometry(3.5, 3.75, 4.0)
        Cert = certify_compatibility(G, D)
        @test Cert.G === G && Cert.dynamic === D && Cert.R_image > 4.03 && Cert.W < 1.9
        # T'(x) = 5/2 - x vanishes at t = 4, of Bernstein parameter ≈ 7.87
        @test_throws ArgumentError certify_compatibility(EllipseGeometry(8.0, 8.5, 9.0), D)
        # the image of ∂E_{3.5} is certified outside E_{4.03} only
        @test_throws ArgumentError certify_compatibility(EllipseGeometry(3.5, 3.75, 4.1), D)
        @test_throws ArgumentError certify_compatibility(G, mod1_dynamic(x -> 1.5x))   # not full branch

        # the search over r: r = 8 fails the derivative check, the best r certifies
        rs = [2.0, 3.0, 3.5, 5.0, 8.0]
        S = search_ellipse_geometry(D, rs; n = 256, cover = 32)
        @test length(S.table) == length(rs) && !S.table[end].feasible
        @test S.cost == maximum(e.cost for e in S.table) && S.R_image > S.r
        @test S.cost == S.R_image / S.r
        Rb = S.r + (S.R_image - S.r) * 2 / 3
        Gb = EllipseGeometry(S.r, (S.r + Rb) / 2, Rb)
        @test certify_compatibility(Gb, D; n = 256, cover = 32).R_image == S.R_image
        # a prescribed cost: the smallest W
        Sw = search_ellipse_geometry(D, rs; cost = (r, R, W) -> -W, n = 256, cover = 32)
        @test Sw.W == minimum(e.W for e in Sw.table if e.feasible)
        @test_throws ArgumentError search_ellipse_geometry(D, [8.0]; n = 256, cover = 32)

        # one call for a family: T_a(x) = 2x + a x(1-x), a ∈ [0.49, 0.51]
        a = interval(0.49, 0.51)
        Dfam = mod1_dynamic(x -> 2x + a * x * (1 - x); full_branch = true)
        Cfam = certify_compatibility(G, Dfam)
        @test Cfam.W >= Cert.W && Cfam.R_image <= Cert.R_image
        for a0 in (0.49, 0.5, 0.51)
            Ca = certify_compatibility(G, mod1_dynamic(x -> 2x + a0 * x * (1 - x); full_branch = true))
            @test Ca.W <= Cfam.W && Ca.R_image >= Cfam.R_image
        end

        N = 12
        B = Legendre(N, G)
        @test B.geometry === G
        A1 = galerkin_matrix(B, D)
        A2 = galerkin_matrix_from_chebyshev(B, D, Cert; Nprime = 64)
        @test all(abs.(A1.c - A2.c) .<= A1.r + A2.r)
        @test maximum(A2.r) < 1e-10
        @test_throws ArgumentError galerkin_matrix_from_chebyshev(B, D, Cert; Nprime = 48)
        # route 2 in BigFloat, on a map built at that precision: it intersects both
        # Float64 routes and is far narrower
        setprecision(BigFloat, 512) do
            Db = mod1_dynamic(f; full_branch = true, T = BigFloat, ϵ = big(2.0)^-460)
            Cb = certify_compatibility(G, Db)
            Ab = galerkin_matrix_from_chebyshev(Legendre(N, G; T = BigFloat), Db, Cb; Nprime = 256)
            @test Ab isa BallMatrix{BigFloat}
            @test maximum(Ab.r) < 1e-130
            @test all(abs.(BigFloat.(A2.c) .- Ab.c) .<= BigFloat.(A2.r) .+ Ab.r)
            @test all(abs.(BigFloat.(A1.c) .- Ab.c) .<= BigFloat.(A1.r) .+ Ab.r)
        end
        D2 = mod1_dynamic(f; full_branch = true)                 # equal map, another object
        @test_throws ArgumentError galerkin_matrix_from_chebyshev(B, D2, Cert; Nprime = 64)
        G2 = EllipseGeometry(3.5, 3.75, 4.0)                     # equal values, another object
        @test_throws ArgumentError galerkin_matrix_from_chebyshev(Legendre(N, G2), D, Cert; Nprime = 64)
        @test_throws ArgumentError galerkin_matrix_from_chebyshev(Legendre(N), D, Cert; Nprime = 64)

        Bc = Chebyshev(16, G)
        @test Bc.geometry === G && Bc.strong.ρ == 3.75
        @test Chebyshev(16, Eρ(3.75)).geometry === nothing    # the old constructors unchanged
        Bl, c2l, l2c = legendre_basis(Bc)
        @test Bl.N == 16 && Bl.geometry === G
        @test inball(Matrix(1.0I, 17, 17), c2l * l2c)
        # T_2(2x-1) = 2(2x-1)^2 - 1 = (4/3)·P_2 - (1/3)·P_0 = (4/(3√5)) e_2 - (1/3) e_0
        v = change_of_basis(Bc, Bl, [0.0, 0.0, 1.0, zeros(14)...])
        @test _in(-1 / 3, v[1]) && _in(4 / (3sqrt(5)), v[3]) && _in(0.0, v[2])
        w = change_of_basis(Bl, Bc, [-1 / 3, 0.0, 4 / (3sqrt(5)), zeros(14)...])
        @test _in(1.0, w[3]) && abs(w[1].c) <= w[1].r + 1e-15
        @test_throws ArgumentError change_of_basis_matrices(Bc, Legendre(16, G2))
        @test_throws ArgumentError change_of_basis_matrices(Bc, Legendre(16))
        @test_throws ArgumentError change_of_basis_matrices(Chebyshev(16, Eρ(3.75)), Bl)
        @test_throws ArgumentError change_of_basis_matrices(Chebyshev(8, G), Bl)
        @test_throws ArgumentError legendre_basis(Chebyshev(16, Eρ(3.75)))

        # constants from the basis and the certificate
        M, cM = legendre_dfly_smallest_M(B, Cert)
        @test IA.sup(cM.A_M) < 1 && IA.sup(legendre_dfly(B, Cert, M - 1).A_M) >= 1
        @test IA.sup(cM.M_w) >= l2_transfer_norm_bound(D)
        @test_throws ArgumentError legendre_dfly(Legendre(N, G2), Cert, 10)
        Ml, cl = legendre_dfly_l1_smallest_M(B, Cert)
        @test Ml == M && IA.sup(cl.A_M) == IA.sup(cM.A_M)
        @test IA.mid(cl.Bp_M) ≈ sum((2j + 1) * 3.75^j for j = 0:M) rtol = 1e-12
        @test_throws ArgumentError legendre_dfly_l1(Legendre(N, G2), Cert, 10)
        @test_throws ArgumentError legendre_dfly_smallest_M(Legendre(N), Cert)
    end

    @testset "L² resolvent on the mean-zero subspace" begin
        A = galerkin_matrix(Legendre(16), D)
        res = l2_resolvent_bound_average_zero(A)
        σ = svdvals(Matrix(1.0I, 16, 16) - A.c[2:end, 2:end])[end]
        @test 1 / σ <= res.bound <= (1 / σ) * (1 + 1e-8)
        @test res.Rbar == max(1, res.bound)
        @test opnorm(A.c[2:end, 2:end]) < l2_transfer_norm_bound(D)
    end

    @testset "constants of the Legendre appendix" begin
        ρ, R, r, W = 3.75, 4.0271, 3.5, 1.8984
        for M in (0, 5, 30)
            CM = sqrt(sum((2j + 1) * ρ^(2j) for j = 0:M))
            @test in_interval(CM, legendre_C_M(ρ, M)) || isapprox(CM, IA.mid(legendre_C_M(ρ, M)); rtol = 1e-14)
            τ2 = sum(((ρ / R)^k + CM * R^(-k))^2 for k = M+1:M+600)   # the series of lem:legendre-tail
            @test sqrt(τ2) ≈ IA.mid(legendre_tau_M(ρ, R, M)) rtol = 1e-10
        end
        @test IA.mid(legendre_l2_tail(ρ, 10)) ≈ ρ^(-11) / sqrt(1 - ρ^(-2))
        @test IA.mid(legendre_l2_tail_geometric(2.0, 0.5, 3)) ≈ 2 * 0.5^4 / sqrt(0.75)
        c = RIM._legendre_dfly(ρ, R, r, W, 60)
        @test IA.mid(c.C₂) ≈ 2W / sqrt(1 - (r / ρ)^2)
        M, cM = RIM._legendre_dfly_smallest_M(ρ, R, r, W)    # the hand values: M = 55
        @test M == 55
        @test IA.sup(cM.A_M) < 1 && IA.sup(RIM._legendre_dfly(ρ, R, r, W, M - 1).A_M) >= 1
        @test_throws ArgumentError legendre_tau_M(4.0, 3.0, 5)
    end
end
