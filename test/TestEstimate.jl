@testset "Estimator" begin


D = Mod1Dynamic(x->2*x)
B = Ulam(32)
P = assemble(B, D)

# C = InvariantMeasures.boundnorm(B, P, 10)

# @test C == [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

end
