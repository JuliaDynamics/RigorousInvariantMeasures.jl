
@testset "Observables" begin
using IntervalArithmetic

B = Ulam(4)

Obs = Observable(B, x->x)

@test Obs.infbound.hi >= 1

@test 0.125 ∈ Obs.v[1]

D = mod1_dynamic(x->2*x)

B = Ulam(1024)
logder = RigorousInvariantMeasures.discretizationlogder(B, D)

@test log(2) ∈ logder.v[5]

end