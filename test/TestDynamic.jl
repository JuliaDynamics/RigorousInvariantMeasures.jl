using InvariantMeasures
using ValidatedNumerics

@testset "Dynamics" begin

D = Mod1Dynamic(x->2*x)

@test D.T(0.1) == 0.2

@test dfly(Lipschitz, L1, D) == (0.5, 0.0)
@test dfly(TotalVariation, L1, D) == (0.5, 0.0)


D = PwMap([ x-> x^2+0.25, x -> 4*x-2, x -> 4*x-3], [0, 0.5, 0.75, 1])

@test InvariantMeasures.DynamicDefinition.derivative(D, 0.1..0.1) ≈ 0.2
@test InvariantMeasures.DynamicDefinition.derivative(D, 0.2..0.3) ≈ 0.4..0.6
@test InvariantMeasures.DynamicDefinition.derivative(D, 0.4..0.6) ≈ 0.8..4
@test InvariantMeasures.DynamicDefinition.derivative(D, 0.7..0.8) ≈ 4

D = PwMap([ x-> 4*x, x -> 2*x-0.5, x -> 4*x-3], [0, 0.25, 0.75, 1])

@test InvariantMeasures.DynamicDefinition.derivative(D, 0.1..0.3) ≈ 2..4
@test InvariantMeasures.DynamicDefinition.derivative(D, 0.1..0.2) ≈ 4..4

D = mod1_dynamic(x -> 3.5x, 0..1)

@test D.endpoints ≈ [0, 2/7, 4/7, 6/7, 1]
@test D.Ts[2](0.5) == 0.75
@test D.is_full == [1, 1, 1, 0]
@test D.orientations == [1,1,1,1]

D = mod1_dynamic(x -> 3.5x + 0.5)
@test D.endpoints ≈ [0, 1/7, 3/7, 5/7, 1]
@test D.Ts[end](1) == 1
@test D.is_full == [0,1,1,1]
@test D.orientations == [1,1,1,1]

D = mod1_dynamic(x -> -3.5x + 0.5)
@test D.endpoints ≈ [0, 1/7, 3/7, 5/7, 1]
@test D.Ts[end](5/7) == 1
@test D.is_full == [0,1,1,1]
@test D.orientations == [-1,-1,-1,-1]


end
