using InvariantMeasures
using ValidatedNumerics

D = Mod1Dynamic(x->2*x)

@test D.T(0.1) == 0.2

@test dfly(Lipschitz, L1, D) == (0.5, 0.0)
@test dfly(TotalVariation, L1, D) == (0.5, 0.0)


D = PwMap([ x-> x^2+0.25, x -> 4*x-2, x -> 4*x-3], [0, 0.5, 0.75, 1])

@test InvariantMeasures.DynamicDefinition.der(D, 0.1..0.1) ≈ 0.2
@test InvariantMeasures.DynamicDefinition.der(D, 0.2..0.3) ≈ 0.4..0.6
@test InvariantMeasures.DynamicDefinition.der(D, 0.4..0.6) ≈ 0.8..4
@test InvariantMeasures.DynamicDefinition.der(D, 0.7..0.8) ≈ 4

D = PwMap([ x-> 4*x, x -> 2*x-0.5, x -> 4*x-3], [0, 0.25, 0.75, 1])

@test InvariantMeasures.DynamicDefinition.der(D, 0.1..0.3) ≈ 2..4
@test InvariantMeasures.DynamicDefinition.der(D, 0.1..0.2) ≈ 4..4
