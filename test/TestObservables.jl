
@testset "Observables" begin
B = Ulam(4)

Obs = Observable(B, x->x)

@test Obs.infbound.hi >= 1

@test 0.125 ∈ Obs.v[1]

end