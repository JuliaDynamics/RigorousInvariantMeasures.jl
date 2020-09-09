

D = Mod1Dynamic(x->2*x)

@test D.T(0.1) == 0.2

@test dfly(Lipschitz, L1, D) == (0.5, 0.0)
@test dfly(TotalVariation, L1, D) == (0.5, 0.0)
