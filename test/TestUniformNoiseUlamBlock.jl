@testset "UniformNoiseUlam, block summation" begin

    using Test
    using Random
    using IntervalArithmetic
    using RigorousInvariantMeasures
    import RigorousInvariantMeasures: gamma, wrap_idx, reflect_outward_idx

    # exact window sums in BigFloat, with the same boundary conditions
    function exact_apply(v, l, bc)
        k = length(v)
        n = 2l + 1
        idx(j) = bc === :periodic ? wrap_idx(j, k) : reflect_outward_idx(j, k)
        vb = big.(v)
        [sum(vb[idx(j)] for j = (i-l):(i+l)) / n for i = 1:k]
    end

    setprecision(BigFloat, 256) do
        Random.seed!(20260911)
        u = eps(Float64) / 2
        for bc in (:periodic, :reflecting), k in (1, 7, 64, 1000), l in (0, 1, 3, 10, 40)
            B = Ulam(k)
            Kb = UniformKernelUlam(Val(bc), B, l; summation = :block)
            Ks = UniformKernelUlam(Val(bc), B, l)
            n = 2l + 1
            γ = gamma(Float64, n)
            for gen in (k -> randn(k), k -> exp.(4 .* randn(k)) .* sign.(randn(k)))
                v = gen(k)
                e = exact_apply(v, l, bc)
                # the local bound |fl((Kv)_i) - (Kv)_i| ≤ γₙ (K|v|)_i
                w = Kb * v
                loc = exact_apply(abs.(v), l, bc)
                @test all(abs.(big.(w) .- e) .<= γ .* loc)
                # the two schemes agree to rounding
                @test isapprox(w, Ks * v; rtol = 0, atol = 64 * n * u * sum(abs, v))
                # the interval version encloses the exact image of every point of the input
                r = 1e-9 .* abs.(v) .+ 1e-12
                vi = [interval(v[i] - r[i], v[i] + r[i]) for i = 1:k]
                wi = Kb * vi
                for x in (v, v .+ r, v .- r, v .+ r .* (2 .* rand(k) .- 1))
                    ex = exact_apply(x, l, bc)
                    @test all(inf.(wi) .<= ex .<= sup.(wi))
                end
            end
        end
    end

    # the scheme is part of the type, and the default is unchanged
    B = Ulam(16)
    @test UniformKernelUlamPeriodic(B, 2) isa UniformKernelUlam{:periodic,Float64,:sliding}
    @test UniformKernelUlamPeriodic(B, 2; summation = :block) isa
          UniformKernelUlam{:periodic,Float64,:block}
    @test_throws ArgumentError UniformKernelUlamPeriodic(B, 2; summation = :pairwise)
end
