using ValidatedNumerics
using InvariantMeasures.HatBasis: HatFunction, IntervalOnTorus

f = HatFunction(1., 2, 3)
@test f(1.5) == 0.5
@test f(1..1.5) == 0..0.5

f = HatFunction(0.125, 0.25, 0.375)
I = IntervalOnTorus(0.375..1.1875)
@test f(I) == 0..0.5

@test f(0..1) == 0..1

f = HatFunction(0.125, 0.25, 0.375)
I = IntervalOnTorus(0.1875..0.25)
@test f(I) == 0.5..1
