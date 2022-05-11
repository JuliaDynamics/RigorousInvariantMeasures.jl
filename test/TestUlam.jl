@testset "Ulam basis" begin

using ValidatedNumerics
using InvariantMeasures
using InvariantMeasures.BasisDefinition


B = Ulam(4)
@test nonzero_on(B, (@interval(0.1), @interval(0.3))) == (1,2)
@test nonzero_on(B, (@interval(0), @interval(1))) == (1,4)
@test nonzero_on(B, (@interval(0.3), @interval(0.31))) == (2,2)
@test nonzero_on(B, (@interval(1), @interval(1))) == (4,4)

@test is_refinement(Ulam(8), Ulam(4))
@test is_refinement(Ulam(8), Ulam(8))
@test !is_refinement(Ulam(4), Ulam(8))
@test length(Ulam(4)) == 4

end
