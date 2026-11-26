using Symbolics, SymbolicUtils, RigorousInvariantMeasures, IntervalArithmetic, Test

# Load the extension module
const HDF = Base.get_extension(RigorousInvariantMeasures, :SymbolicsExt)

@testset "Test HigherDFLY" begin

    # ──── SymbL constructor ────

    @testset "SymbL constructor" begin
        @variables z
        P = HDF.SymbL(3, z)
        @test P.n == 3
        @test P.f === z
    end

    # ──── Diff (Eq. 3.2: (L_k h)' = L_{k+1}(h') + k·L_k(h·D_f)) ────

    @testset "Diff single operator" begin
        @variables x
        @variables (f(x))[1:4]
        @variables (DD(x))[1:4]
        @variables Dist(x)

        ∂ = Differential(x)

        # Differentiate L_1(f) → [L_2(f'), L_1(1·f·Dist)]
        P = HDF.SymbL(1, f[1])
        result = HDF.Diff(P, ∂, Dist)
        @test length(result) == 2
        @test result[1].n == 2
        @test isequal(result[1].f, ∂(f[1]))
        @test result[2].n == 1
        @test isequal(result[2].f, f[1] * Dist)

        # Differentiate L_3(f) → [L_4(f'), L_3(3·f·Dist)]
        P = HDF.SymbL(3, f[1])
        result = HDF.Diff(P, ∂, Dist)
        @test result[1].n == 4
        @test result[2].n == 3
        @test isequal(result[2].f, 3 * f[1] * Dist)
    end

    @testset "Diff vector of operators" begin
        @variables x
        @variables (f(x))[1:4]
        @variables Dist(x)
        ∂ = Differential(x)

        v = [
            HDF.SymbL(1, f[1]),
            HDF.SymbL(2, f[2]),
        ]
        result = HDF.Diff(v, ∂, Dist)
        @test length(result) == 4  # 2 terms per input
    end

    # ──── compute_dfly_k_fi_DDi ────

    @testset "compute_dfly_k_fi_DDi k=0" begin
        @variables x
        @variables (f(x))[1:3]

        # k=0: just L_1(f[1])
        v = HDF.compute_dfly_k_fi_DDi(0)
        @test v isa HDF.SymbL
        @test v.n == 1
        @test isequal(v.f, f[1])
    end

    @testset "compute_dfly_k_fi_DDi k=1" begin
        @variables x
        @variables (f(x))[1:3]
        @variables (DD(x))[1:3]

        # k=1: (L_1 f)' = L_2(f') + L_1(f·DD[1])
        v = HDF.compute_dfly_k_fi_DDi(1)
        @test length(v) == 2
        @test v[1].n == 2
        @test isequal(v[1].f, f[2])
        @test v[2].n == 1
        @test isequal(v[2].f, DD[1] * f[1])
    end

    @testset "compute_dfly_k_fi_DDi k=2" begin
        @variables x
        @variables (f(x))[1:4]
        @variables (DD(x))[1:4]

        # k=2: second derivative of L_1(f)
        # (L_1 f)'' should produce terms with n=3, n=2, n=1
        v = HDF.compute_dfly_k_fi_DDi(2)
        @test length(v) >= 3

        ns = [t.n for t in v]
        @test 3 in ns  # L_3(f'') term
        @test 1 in ns  # L_1(...) terms from distortion chain

        # All terms should be nonzero
        for t in v
            @test !isequal(t.f, 0)
        end
    end

    @testset "compute_dfly_k_fi_DDi k=3" begin
        @variables x
        @variables (f(x))[1:5]
        @variables (DD(x))[1:5]

        v = HDF.compute_dfly_k_fi_DDi(3)
        @test length(v) >= 4

        ns = [t.n for t in v]
        @test 4 in ns  # L_4(f''') from the highest chain
        @test 1 in ns  # L_1(...) terms
    end

    # ──── the expansion as integer monomials ────

    # (Lf)'' = L3 f'' + 3 L2(f' D1) + L2(f D2) + L1(f D1^2), D_m = (1/T')^(m)
    @testset "k = 2 expansion" begin
        merged = Dict{Tuple{Int,Int,Vector{Int}},Int}()
        for t in HDF.dfly_monomials(2)
            key = (t.n, t.i, t.p)
            merged[key] = get(merged, key, 0) + t.coeff
        end
        @test merged == Dict(
            (3, 2, [0, 0]) => 1,
            (2, 1, [1, 0]) => 3,
            (2, 0, [0, 1]) => 1,
            (1, 0, [2, 0]) => 1,
        )
        @test HDF.dfly_monomials(0) == [(n = 1, i = 0, coeff = 1, p = Int[])]
    end

    # ──── certified constants on T(x) = 2x + x(1-x)/2 mod 1 ────
    #
    # T'(y) = 5/2 - y on [0,1], so (1/T')^(i)(y) = i!/(5/2 - y)^(i+1), largest at y = 1
    # where 1/T' = 2/3. Hence sup|(1/T')^(i)| |1/T'|^(n-1) = i! (2/3)^(i+n), and
    # c10 = 4/9, c11 = 2/3, c20 = sup 3/(5/2-y)^4 = 48/81, c21 = 3 (2/3)^3 = 8/9, c22 = 4/9.

    D = mod1_dynamic(x -> 2x + x * (1 - x) / 2)
    exact_C = [
        1.0 0 0
        4 / 9 2 / 3 0
        48 / 81 8 / 9 4 / 9
    ]
    encloses(bound, exact) = exact <= bound <= exact * (1 + 1e-3)

    @testset "distortion_vals" begin
        vals = HDF.distortion_vals(D, 2)
        for i = 0:2, n = 1:3
            exact = i == 0 ? (2 / 3)^n : factorial(i) * (2 / 3)^(i + n)
            @test encloses(vals[i+1, n], exact)
        end
    end

    @testset "dfly_coefficient_matrix" begin
        C = HDF.dfly_coefficient_matrix(D, 2)
        for j = 0:2, i = 0:2
            if i > j
                @test C[j+1, i+1] == 0
            else
                @test encloses(C[j+1, i+1], exact_C[j+1, i+1])
            end
        end
    end

    @testset "dfly(W{2,1}, L1, D)" begin
        C = HDF.dfly_coefficient_matrix(D, 2)
        c8 = dfly(W{2,1}, L1, D; A_target = 0.8)
        @test c8.strong isa WeightedSobolev
        @test strong_kind(c8) === c8.strong && weak_kind(c8) == L1
        @test c8.strong.a[2] == 1 && c8.strong.A_target == 0.8
        @test c8.strong.a[1] ≈ C[3, 2] / (0.8 - C[2, 2]) rtol = 1e-14
        @test 0.8 <= c8.A <= 0.8 + 1e-12
        # k = 2: B(A*) = c20 + c10 c21 / (A* - c11), not below its value at the exact c
        @test c8.B ≈ C[3, 1] + C[2, 1] * C[3, 2] / (0.8 - C[2, 2]) rtol = 1e-12
        @test c8.B >= 48 / 81 + (4 / 9) * (8 / 9) / (0.8 - 2 / 3)

        c = dfly(W{2,1}, L1, D)
        Astar = c.strong.A_target
        @test C[2, 2] < Astar < 1
        @test c.A <= Astar + 1e-12 && c.A < 1
        @test c.B / (1 - c.A) <= c8.B / (1 - c8.A)

        @test WeightedSobolev(D, 2; A_target = 0.8).a == c8.strong.a
        @test dfly(W{2,1}, L1, D; weights = c8.strong.a).A == c8.A

        c3 = dfly(W{3,1}, L1, D)
        @test c3.A < 1 && isfinite(c3.B) && c3.strong.k == 3
    end

    @testset "fixed weights on a perturbed map" begin
        norm = WeightedSobolev(D, 2; A_target = 0.8)
        Dp = mod1_dynamic(x -> 2x + x * (1 - x) / 2 + x * (1 - x) / 100)
        cp = dfly(W{2,1}, L1, Dp; weights = norm)
        @test cp.strong === norm
        @test cp.A < 1 && isfinite(cp.B)
    end

    @testset "Chebyshev basis on the weighted seminorm" begin
        norm = WeightedSobolev(D, 2; A_target = 0.8)
        p = RigorousInvariantMeasures.ChebPoints(16, Float64)
        Bc1 = RigorousInvariantMeasures.Chebyshev(p, 2, norm, C1())
        a1 = norm.a[1]
        @test RigorousInvariantMeasures.weak_by_strong_and_aux_bound(Bc1) ==
              (max(2 / a1, 1.0), 1.0)
        Bl2 = RigorousInvariantMeasures.Chebyshev(p, 2, norm, L2())
        S1, S2 = RigorousInvariantMeasures.weak_by_strong_and_aux_bound(Bl2)
        @test S1 >= 1 / a1 && S2 == 1.0
        Bμ = RigorousInvariantMeasures.Chebyshev(p, 2, norm, L2μ())
        @test RigorousInvariantMeasures.weak_by_strong_and_aux_bound(Bμ)[2] >= π / 2

        c = dfly(norm, L1, D)
        @test strong_kind(c) == norm
        us = RigorousInvariantMeasures.invariant_measure_strong_norm_bound(Bc1, D)
        @test us >= c.B / (1 - c.A)
        other = WeightedSobolev(D, 2; A_target = 0.9)
        @test_throws ArgumentError RigorousInvariantMeasures.invariant_measure_strong_norm_bound(
            Bc1,
            D;
            dfly_coefficients = dfly(other, L1, D),
        )

        @test_throws ErrorException RigorousInvariantMeasures.strong_weak_bound(Bc1)
        # projection errors are per unit of ||f^(k)||_1 <= ||f||_a, as for W^{k,1}
        Bw = RigorousInvariantMeasures.Chebyshev(p, 2, W{2,1}(), L2())
        @test RigorousInvariantMeasures.aux_normalized_projection_error(Bl2) ==
              RigorousInvariantMeasures.aux_normalized_projection_error(Bw)
        @test_throws ArgumentError RigorousInvariantMeasures.weak_projection_error(Bc1)
        # a plain W^{2,1} basis refuses constants stated in the weighted seminorm
        @test_throws ArgumentError RigorousInvariantMeasures.invariant_measure_strong_norm_bound(
            Bw,
            D,
        )
        @test_throws ArgumentError RigorousInvariantMeasures.weak_by_strong_and_aux_bound(
            RigorousInvariantMeasures.Chebyshev(p, 3, norm, C1()),
        )
    end

    @testset "Chebyshev(n, k, D)" begin
        RIM = RigorousInvariantMeasures
        for k = 2:4
            B = RIM.Chebyshev(32, k, D)
            @test RIM.strong_norm(B) isa WeightedSobolev && RIM.strong_norm(B).k == k
            c = dfly(RIM.strong_norm(B), RIM.aux_norm(B), D)
            @test strong_kind(c) == RIM.strong_norm(B)
            us = RIM.invariant_measure_strong_norm_bound(B, D)
            @test isfinite(us) && us >= c.B / (1 - c.A)
        end
        B8 = RIM.Chebyshev(32, 2, D; A_target = 0.8, weak = C1)
        @test RIM.strong_norm(B8).A_target == 0.8 && RIM.weak_norm(B8) == C1
    end

    @testset "Fourier bases on the weighted seminorm" begin
        RIM = RigorousInvariantMeasures
        Dc = mod1_dynamic(x -> 2x + RIM.sinpi(2x) / 20; circle_smoothness = 3)
        B = FourierAnalytic(16, 64, Dc, 2)
        @test RIM.strong_norm(B) isa WeightedSobolev && RIM.strong_norm(B).k == 2
        a = RIM.strong_norm(B).a
        S1, S2 = RIM.weak_by_strong_and_aux_bound(B)
        @test 1 / a[1] <= S1 <= nextfloat(1 / a[1]) && S2 == 1.0
        @test RIM.strong_weak_bound(B) >= a[1] * 2π * 16 + (2π * 16)^2
        @test RIM.weak_projection_error(B) ==
              RIM.weak_projection_error(FourierAnalytic(RIM.FourierPoints(64, Float64), 16, W{2,1}(), L2()))
        c = dfly(RIM.strong_norm(B), L1, Dc)
        @test RIM.invariant_measure_strong_norm_bound(B, Dc) >= c.B / (1 - c.A)
        @test FourierAdjoint(16, 64, Dc, 2).strong == RIM.strong_norm(B)
        # the basis needs the map asserted C^s on the circle
        @test_throws ArgumentError FourierAnalytic(16, 64, Dc, 4)
        @test_throws ArgumentError RIM.invariant_measure_strong_norm_bound(B, D)
        @test_throws MethodError FourierAnalytic(16, 64, D, 2)
    end

    @testset "check_dfly_coefficients on every route" begin
        RIM = RigorousInvariantMeasures
        Bw = RIM.Chebyshev(RIM.ChebPoints(16, Float64), 2, W{2,1}(), L2())
        c = dfly(W{2,1}, L1, D)              # stated in WeightedSobolev
        @test_throws ArgumentError RIM.basis_dfly(Bw, D)
        @test_throws ArgumentError RIM.check_dfly_coefficients(Bw, D, c)
        @test_throws ArgumentError RIM.norms_of_powers_dfly(Bw, D, 3)
        @test_throws ArgumentError RIM.strong_resolvent_lift(Bw, D, 0.9, 1.0)
        # a bare tuple has no record of its norm and is not checked
        @test RIM.check_dfly_coefficients(Bw, D, (c.A, c.B)) == (c.A, c.B)
        # the basis built on the norm of the pair passes
        Ba = RIM.Chebyshev(32, 2, D)
        @test RIM.basis_dfly(Ba, D).strong == RIM.strong_norm(Ba)
        # Fourier of Sobolev order 2 needs the map asserted C^2 on the circle
        Dc = mod1_dynamic(x -> 2x + RIM.sinpi(2x) / 20; circle_smoothness = 3)
        Bf = FourierAnalytic(16, 64, Dc, 2)
        @test RIM.basis_dfly(Bf, Dc).strong == RIM.strong_norm(Bf)
        @test_throws ArgumentError RIM.check_dfly_coefficients(Bf, D, dfly(RIM.strong_norm(Bf), L1, Dc))
    end

    @testset "guards" begin
        @test_throws ArgumentError dfly(W{2,2}, L1, D)
        @test_throws ArgumentError dfly(W{2,1}, L1, D; A_target = 0.5)
        @test_throws ArgumentError dfly(W{2,1}, L1, D; A_target = 1.0)
        @test_throws ArgumentError dfly(W{2,1}, L1, D; A_target = 0.8, weights = [1.0, 1.0])
        @test_throws ArgumentError dfly(W{2,1}, L1, D; weights = [1.0, 1.0, 1.0])

        Dβ = mod1_dynamic(x -> 5x / 2)
        @test !is_full_branch(Dβ)
        @test_throws ArgumentError dfly(W{2,1}, L1, Dβ)

        # full branch, T' = 2 + (π/2) cos(2πx) >= 0.42, so sup|1/T'| > 1
        Dslow = mod1_dynamic(x -> 2x + RigorousInvariantMeasures.sinpi(2x) / 4)
        @test is_full_branch(Dslow)
        @test_throws ErrorException dfly(W{2,1}, L1, Dslow)

        @test_throws ArgumentError WeightedSobolev(2, [1.0])
        @test_throws ArgumentError WeightedSobolev(2, [1.0, 2.0])
        @test_throws ArgumentError WeightedSobolev(2, [-1.0, 1.0])
    end

end
