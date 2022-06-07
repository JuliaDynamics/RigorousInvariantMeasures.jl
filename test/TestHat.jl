@testset "Periodic Hat basis" begin

using ValidatedNumerics
using RigorousInvariantMeasures
using RigorousInvariantMeasures: HatFunction, HatFunctionOnTorus, IntervalOnTorus, nonzero_on, is_refinement

f = HatFunction(1., 2, 3)
@test f(1.5) == 0.5
@test f(1..1.5) == 0..0.5

f = HatFunctionOnTorus(0.125, 0.25, 0.375)
x = IntervalOnTorus(0.375..1.1875)
@test f(x) == 0..0.5

@test f(IntervalOnTorus(0..1)) == 0..1

f = HatFunctionOnTorus(0.125, 0.25, 0.375)
x = IntervalOnTorus(3.1875..3.25)
@test f(x) == 0.5..1

f = HatFunctionOnTorus(0, 0.125, 0.25)
x = IntervalOnTorus(0..0.0625)
@test f(x) == 0..0.5

f = HatFunctionOnTorus(0.875, 0, 0.125)
x = IntervalOnTorus(0..0.0625)
@test f(x) == 0.5..1

f = HatFunctionOnTorus(0.875, 0, 0.125)
x = IntervalOnTorus(0.9375..1)
@test f(x) == 0.5..1

B = Hat(4)
@test nonzero_on(B, (0.1..0.3, NaN)) == (1,3)
@test nonzero_on(B, (0..1, NaN)) == (1,4)
@test nonzero_on(B, (0.3..0.31, NaN)) == (2,3)

@test is_refinement(Hat(8), Hat(4))
@test is_refinement(Hat(8), Hat(8))
@test !is_refinement(Hat(4), Hat(8))


end
