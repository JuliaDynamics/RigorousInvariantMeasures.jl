@testset "Test BasisDefinition.jl" begin


    using RigorousInvariantMeasures

    struct TestBasis <: RigorousInvariantMeasures.Basis end

    B = TestBasis()

    # The abstract `Basis` interface is exposed via `function … end`
    # declarations: any subtype that doesn't provide a method gets a clean
    # `MethodError`, instead of the previous `@error "Not Implemented"`
    # which logged but silently returned `nothing`.
    @test_throws MethodError length(B)
    @test_throws MethodError is_dual_element_empty(B, 1.0)
    @test_throws MethodError nonzero_on(B, 1.0)
    @test_throws MethodError RigorousInvariantMeasures.evaluate(B, 1, 1.0)
    @test_throws MethodError RigorousInvariantMeasures.evaluate_integral(B, 1.0)
    @test_throws MethodError strong_norm(B)
    @test_throws MethodError weak_norm(B)
    @test_throws MethodError aux_norm(B)
    @test_throws MethodError is_refinement(B, B)
    @test_throws MethodError integral_covector(B)
    @test_throws MethodError one_vector(B)
    @test is_integral_preserving(B) == false

    U0 = AverageZero(B)
    @test_throws MethodError Base.iterate(U0, 1)

    BU = Ulam(1024)
    BU0 = AverageZero(BU)
    @test length(BU0) == 1023

    @test_throws MethodError RigorousInvariantMeasures.weak_projection_error(B)
    @test_throws MethodError RigorousInvariantMeasures.aux_normalized_projection_error(B)
    @test_throws MethodError RigorousInvariantMeasures.strong_weak_bound(B)
    @test_throws MethodError RigorousInvariantMeasures.aux_weak_bound(B)
    @test_throws MethodError RigorousInvariantMeasures.weak_by_strong_and_aux_bound(B)
    @test_throws MethodError RigorousInvariantMeasures.bound_weak_norm_from_linalg_norm(B)
    @test_throws MethodError RigorousInvariantMeasures.bound_linalg_norm_L1_from_weak(B)
    @test_throws MethodError RigorousInvariantMeasures.bound_linalg_norm_L∞_from_weak(B)

    D = mod1_dynamic(x -> 2 * x)
    @test_throws MethodError RigorousInvariantMeasures.invariant_measure_strong_norm_bound(
        B,
        D,
    )
    @test_throws MethodError bound_weak_norm_abstract(B, D)


    @test_throws MethodError opnormbound(
        B,
        L1,
        [
            1.0 0.0
            0.0 1.0
        ],
    )
    @test_throws MethodError normbound(B, L1, [1.0; 1.0])



end
