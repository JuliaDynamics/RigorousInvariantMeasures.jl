using Symbolics, SymbolicUtils, RigorousInvariantMeasures

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

    # ──── _optimize_mult ────

    @testset "_optimize_mult basic" begin
        @variables x
        @variables (f(x))[1:3]
        @variables (DD(x))[1:3]

        # Create a simple Mul term: DD[1] * f[2]
        term = (DD[1] * f[2]).val
        if SymbolicUtils.ismul(term)
            # vals[i, j]: vals[2, j] = bound for DD[1] with weight (T')^{j-1}
            vals = ones(3, 3)
            vals[2, 1] = 2.0   # ||DD[1]||_∞ = 2
            vals[2, 2] = 0.5   # ||DD[1]/(T')||_∞ = 0.5

            result = HDF._optimize_mult(2, 2, term, vals)
            # DD[1] has pow=1, n=2
            # bound = vals[2, 2] * vals[2, 1]^(1-1) = 0.5 * 1 = 0.5
            # result = 1 * 0.5 * f[2]
            @test !isequal(result, 0)
        end
    end

    @testset "_optimize_mult exponent correctness" begin
        @variables x
        @variables (f(x))[1:4]
        @variables (DD(x))[1:4]

        # Term with DD[1]^2 * f[1]
        term_expr = 2 * DD[1]^2 * f[1]
        if SymbolicUtils.ismul(term_expr.val)
            vals = ones(4, 4)
            vals[2, 1] = 3.0   # ||DD[1]||_∞ = 3
            vals[2, 2] = 1.5   # ||DD[1]/(T')||_∞ = 1.5

            result = HDF._optimize_mult(3, 2, term_expr.val, vals)
            # DD[1] has pow=2, it's the last (and only) DD term, n=2
            # bound = vals[2, 2] * vals[2, 1]^(2-1) = 1.5 * 3.0 = 4.5
            # result = 2 * 4.5 * f[1] = 9.0 * f[1]
            if SymbolicUtils.ismul(result.val)
                @test result.val.coeff ≈ 9.0
            end
        end
    end

    @testset "_optimize_mult with multiple DD terms" begin
        @variables x
        @variables (f(x))[1:4]
        @variables (DD(x))[1:4]

        # Term: DD[1]^2 * DD[2] * f[1]
        term_expr = DD[1]^2 * DD[2] * f[1]
        if SymbolicUtils.ismul(term_expr.val)
            vals = ones(4, 4)
            vals[2, 1] = 3.0   # ||DD[1]||_∞ = 3
            vals[2, 2] = 1.5   # ||DD[1]/(T')||_∞ = 1.5
            vals[3, 1] = 2.0   # ||DD[2]||_∞ = 2
            vals[3, 2] = 0.8   # ||DD[2]/(T')||_∞ = 0.8

            result = HDF._optimize_mult(3, 2, term_expr.val, vals)
            # Highest DD index is 2 (DD[2]), pow=1, n=2
            # bound = vals[3, 2] * vals[3, 1]^(1-1) = 0.8 * 1 = 0.8
            # Then j=1: DD[1] has pow=2
            # bound *= vals[2, 1]^2 = 3.0^2 = 9.0
            # total bound = 0.8 * 9.0 = 7.2
            # result = 1 * 7.2 * f[1]
            if SymbolicUtils.ismul(result.val)
                @test result.val.coeff ≈ 7.2
            end
        end
    end

    # ──── optimize_coefficients ────

    @testset "optimize_coefficients k=1" begin
        v = HDF.compute_dfly_k_fi_DDi(1)

        vals = ones(3, 3)
        vals[2, 1] = 0.5   # ||DD[1]||_∞
        vals[2, 2] = 0.3   # ||DD[1]/(T')||_∞

        result = HDF.optimize_coefficients(1, v, vals)
        @test !isequal(result, 0)
    end

    @testset "optimize_coefficients k=2" begin
        v = HDF.compute_dfly_k_fi_DDi(2)

        vals = ones(4, 4)
        vals[2, 1] = 0.5
        vals[2, 2] = 0.3
        vals[2, 3] = 0.2
        vals[3, 1] = 0.4
        vals[3, 2] = 0.25

        result = HDF.optimize_coefficients(2, v, vals)
        @test !isequal(result, 0)
    end

    # ──── substitute_values ────

    @testset "substitute_values" begin
        @variables x
        @variables (f(x))[1:3]
        @variables (DD(x))[1:3]

        v = HDF.compute_dfly_k_fi_DDi(1)
        # v = [L_2(f[2]), L_1(DD[1]*f[1])]

        # vals = [λ, DD_1_val, DD_2_val]
        vals = [0.5, 0.3, 0.1]

        w = HDF.substitute_values(2, v, vals)
        @test length(w) == length(v)

        # First term: L_2(f[2]) → λ^(2-1) * f[2] = 0.5 * f[2]
        @test !isequal(w[1], 0)

        # Second term: L_1(DD[1]*f[1]) → λ^(1-1) * substitute(DD[1]*f[1], DD[1]=>0.3)
        #            = 1.0 * 0.3 * f[1]
        @test !isequal(w[2], 0)
    end

end
