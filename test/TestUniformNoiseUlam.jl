@testset "UniformNoiseUlam" begin

    using Test
    using IntervalArithmetic
    using RigorousInvariantMeasures

    import RigorousInvariantMeasures: wrap_idx, reflect_outward_idx


    k = 5

    # ---- Periodic wrap_idx tests ----
    @test wrap_idx(1, k) == 1
    @test wrap_idx(5, k) == 5
    @test wrap_idx(6, k) == 1    # wraps around
    @test wrap_idx(0, k) == 5    # negative shift equivalent
    @test wrap_idx(-1, k) == 4

    # ---- Reflecting reflect_outward_idx tests ----
    # Inside the interval: identity
    @test reflect_outward_idx(1, k) == 1
    @test reflect_outward_idx(5, k) == 5
    @test reflect_outward_idx(3, k) == 3

    # Just outside the left boundary
    @test reflect_outward_idx(0, k) == 1
    @test reflect_outward_idx(-1, k) == 2
    @test reflect_outward_idx(-2, k) == 3
    @test reflect_outward_idx(-3, k) == 4
    @test reflect_outward_idx(-4, k) == 5  # folds into 2k=10, mirror back

    # Just outside the right boundary
    @test reflect_outward_idx(6, k) == 5
    @test reflect_outward_idx(7, k) == 4
    @test reflect_outward_idx(8, k) == 3
    @test reflect_outward_idx(9, k) == 2
    @test reflect_outward_idx(10, k) == 1
    @test reflect_outward_idx(11, k) == 1  # period-2k wrap, then mirror

    # --- Tests ---

    # 1. Identity when l = 0
    function test_identity()
        B = Ulam(10)
        v = rand(length(B))
        Kp = UniformKernelUlamPeriodic(B, 0)
        Kr = UniformKernelUlamReflecting(B, 0)
        @test v ≈ (Kp * v)
        @test v ≈ (Kr * v)
    end

    # 2. Periodic kernel averages correctly
    function test_periodic_average()
        B = Ulam(5)
        v = collect(1.0:5.0)  # [1,2,3,4,5]
        K = UniformKernelUlamPeriodic(B, 1)
        w = K * v
        # manual check for first entry (neighbors are 5,1,2)
        expected1 = (5 + 1 + 2) / 3
        @test isapprox(w[1], expected1)
    end

    # 3. Reflecting kernel averages correctly
    function test_reflecting_average()
        B = Ulam(5)
        v = collect(1.0:5.0)  # [1,2,3,4,5]
        K = UniformKernelUlamReflecting(B, 1)
        w = K * v
        # at i=1, neighbors are 2 and reflected neighbor also 2
        expected1 = (2*1 + 2) / 3
        @test isapprox(w[1], expected1)
    end

    # 4. Normalization (stochasticity)
    function test_stochasticity()
        B = Ulam(100)
        v = rand(length(B))
        Kp = UniformKernelUlamPeriodic(B, 3)
        Kr = UniformKernelUlamReflecting(B, 3)
        @test isapprox(sum(Kp * v), sum(v); rtol=1e-12)
        @test isapprox(sum(Kr * v), sum(v); rtol=1e-12)
    end

    # 5. Interval consistency
    function test_interval_enclosure()
        B = Ulam(10)
        v = rand(length(B))
        iv = map(Interval, v)
        K = UniformKernelUlamPeriodic(B, 2)
        w_float = K * v
        w_interval = K * iv
        @test all(w_float[i] ∈ w_interval[i] for i in 1:length(v))
    end

    #import Statistics   
    # 6. Convergence under repeated application
    function test_convergence()
        B = Ulam(20)
        v = rand(length(B))
        K = UniformKernelUlamPeriodic(B, 2)
        w = copy(v)
        for _ in 1:300
            mul!(K, w)   # make sure mul! writes in-place
        end

        mv = sum(v) / length(v)

        @test isapprox(w, fill(mv, length(v)); rtol=1e-8)
    end

    # run all
    test_identity()
    test_periodic_average()
    test_reflecting_average()
    test_stochasticity()
    test_interval_enclosure()
    test_convergence()

end
