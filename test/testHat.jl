@testset "Hat basis" begin

using ValidatedNumerics
using InvariantMeasures.HatBasis: HatFunction, IntervalOnTorus, Hat, nonzero_on, EquispacedPartition

f = HatFunction(1., 2, 3)
@test f(1.5) == 0.5
@test f(1..1.5) == 0..0.5

f = HatFunctionOnTorus(0.125, 0.25, 0.375)
I = IntervalOnTorus(0.375..1.1875)
@test f(I) == 0..0.5

@test f(0..1) == 0..1

f = HatFunctionOnTorus(0.125, 0.25, 0.375)
I = IntervalOnTorus(3.1875..3.25)
@test f(I) == 0.5..1

f = HatFunctionOnTorus(0, 0.125, 0.25)
I = IntervalOnTorus(0..0.0625)
@test f(I) == 0..0.5

f = HatFunctionOnTorus(0.875, 0, 0.125)
I = IntervalOnTorus(0..0.0625)
@test f(I) == 0.5..1

f = HatFunctionOnTorus(0.875, 0, 0.125)
I = IntervalOnTorus(0.9375..1)
@test f(I) == 0.5..1

B = Hat(EquispacedPartition{Float64}(4))
@test nonzero_on(B, 0.1..0.3) == (1,3)
@test nonzero_on(B, 0..1) == (1,4)
@test nonzero_on(B, 0.3..0.31) == (2,3)

end
