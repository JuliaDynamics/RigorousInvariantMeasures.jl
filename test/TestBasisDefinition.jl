@testset "Test BasisDefinition.jl" begin


    using RigorousInvariantMeasures
    
    struct TestBasis <: RigorousInvariantMeasures.Basis end

    B = TestBasis()

    @test_logs (:error, "Not Implemented") length(B)
    @test_logs (:error, "Not Implemented") is_dual_element_empty(B, 1.0)
    @test_logs (:error, "Not Implemented") nonzero_on(B, 1.0)
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.evaluate(B, 1, 1.0)
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.evaluate_integral(B, 1.0)
    @test_logs (:error, "Must be specialized") strong_norm(B)
    @test_logs (:error, "Must be specialized") weak_norm(B)
    @test_logs (:error, "Must be specialized") aux_norm(B)
    @test_logs (:error, "Not Implemented") is_refinement(B, B)
    @test_logs (:error, "Must be specialized") integral_covector(B)
    @test_logs (:error, "Must be specialized") one_vector(B)
    @test is_integral_preserving(B) == false

    U0 = AverageZero(B)
    @test_logs (:error, "Not Implemented") Base.iterate(U0, 1)

    BU = Ulam(1024)
    BU0 = AverageZero(BU)
    @test length(BU0) == 1023

    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.weak_projection_error(B)
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.aux_normalized_projection_error(
        B,
    )
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.strong_weak_bound(B)
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.aux_weak_bound(B)
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.weak_by_strong_and_aux_bound(B)
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.bound_weak_norm_from_linalg_norm(
        B,
    )
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.bound_linalg_norm_L1_from_weak(B)
    @test_logs (:error, "Not Implemented") RigorousInvariantMeasures.bound_linalg_norm_L∞_from_weak(B)

    D = mod1_dynamic(x -> 2 * x)
    @test_logs (:error, "Must be specialized") RigorousInvariantMeasures.invariant_measure_strong_norm_bound(
        B,
        D,
    )
    @test_logs (:error, "Must be specialized") bound_weak_norm_abstract(
        B,
        D,
    )


    @test_logs (:error, "Must be specialized") opnormbound(
        B,
        L1,
        [
            1.0 0.0
            0.0 1.0
        ],
    )
    @test_logs (:error, "Must be specialized") normbound(B, L1, [1.0; 1.0])



end
