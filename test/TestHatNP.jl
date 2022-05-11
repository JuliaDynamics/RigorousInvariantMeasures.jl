@testset "Hat basis" begin

using ValidatedNumerics
using InvariantMeasures
using InvariantMeasures: HatFunction, HatFunctionOnTorus, IntervalOnTorus, nonzero_on, is_refinement

f = HatFunction(1., 2, 3)
@test f(1.5) == 0.5
@test f(1..1.5) == 0..0.5

f = HatFunction(0.125, 0.25, 0.375)
x = Interval(0.1875..0.25)
@test f(x) == 0.5..1

f = HatFunction(0, 0.125, 0.25)
x = Interval(0..0.0625)
@test f(x) == 0..0.5

B = HatNP(4)
@test nonzero_on(B, (0.1..0.3, NaN)) == (1,3)
@test nonzero_on(B, (0..1, NaN)) == (1,5)
@test nonzero_on(B, (0.3..0.31, NaN)) == (2,3)
@test nonzero_on(B, (0.9..0.91, NaN)) == (4,5)

@test is_refinement(Hat(8), Hat(4))
@test is_refinement(Hat(8), Hat(8))
@test !is_refinement(Hat(4), Hat(8))


end
