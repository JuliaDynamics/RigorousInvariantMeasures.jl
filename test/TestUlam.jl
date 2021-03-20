@testset "Hat basis" begin

using ValidatedNumerics
using InvariantMeasures
using InvariantMeasures.BasisDefinition


B = Ulam(4)
@test nonzero_on(B, (@interval(0.1), @interval(0.3))) == (1,2)
@test nonzero_on(B, (@interval(0), @interval(1))) == (1,4)
@test nonzero_on(B, (@interval(0.3), @interval(0.31))) == (2,2)
@test nonzero_on(B, (@interval(1), @interval(1))) == (4,4)


end
