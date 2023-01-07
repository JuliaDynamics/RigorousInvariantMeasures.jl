@testset "Test BasisDefinition.jl" begin


using RigorousInvariantMeasures
import RigorousInvariantMeasures.BasisDefinition

struct TestBasis <: RigorousInvariantMeasures.Basis
end

B = TestBasis()

@test_logs (:error,"Not Implemented") length(B)
@test_logs (:error,"Not Implemented") BasisDefinition.is_dual_element_empty(B, 1.0) 
@test_logs (:error,"Not Implemented") BasisDefinition.nonzero_on(B, 1.0) 
@test_logs (:error,"Not Implemented") BasisDefinition.evaluate(B, 1, 1.0) 
@test_logs (:error,"Not Implemented") BasisDefinition.evaluate_integral(B, 1.0)
@test_logs (:error,"Must be specialized") BasisDefinition.strong_norm(B)
@test_logs (:error,"Must be specialized") BasisDefinition.weak_norm(B)
@test_logs (:error,"Must be specialized") BasisDefinition.aux_norm(B)
@test_logs (:error,"Not Implemented") BasisDefinition.is_refinement(B, B) 
@test_logs (:error,"Must be specialized") BasisDefinition.integral_covector(B)
@test_logs (:error,"Must be specialized") BasisDefinition.one_vector(B)
@test BasisDefinition.is_integral_preserving(B) == false

U0 = BasisDefinition.AverageZero(B)
@test_logs (:error,"Not Implemented") Base.iterate(U0, 1) 

BU = Ulam(1024)
BU0 = BasisDefinition.AverageZero(BU)
@test length(BU0) == 1023

@test_logs (:error,"Not Implemented") BasisDefinition.weak_projection_error(B) 
@test_logs (:error,"Not Implemented") BasisDefinition.aux_normalized_projection_error(B)
@test_logs (:error,"Not Implemented") BasisDefinition.strong_weak_bound(B)
@test_logs (:error,"Not Implemented") BasisDefinition.aux_weak_bound(B)
@test_logs (:error,"Not Implemented") BasisDefinition.weak_by_strong_and_aux_bound(B)
@test_logs (:error,"Not Implemented") BasisDefinition.bound_weak_norm_from_linalg_norm(B)
@test_logs (:error,"Not Implemented") BasisDefinition.bound_linalg_norm_L1_from_weak(B)
@test_logs (:error,"Not Implemented") BasisDefinition.bound_linalg_norm_L∞_from_weak(B)

D = mod1_dynamic(x->2*x)
@test_logs (:error,"Must be specialized") BasisDefinition.invariant_measure_strong_norm_bound(B, D)
@test_logs (:error,"Must be specialized") BasisDefinition.bound_weak_norm_abstract(B, D)


# Failing due to a delicate issue, while L1 is a subtype of NormKind
# Type{L1} is not a subtype of Type{NormKind}
# 
# 
# @test_logs (:error,"Must be specialized") BasisDefinition.opnormbound(B, L1, [1.0 0.0;
#                                                                             0.0  1.0])
# @test_logs (:error,"Must be specialized") BasisDefinition.normbound(B, L1, [1.0; 1.0])



end