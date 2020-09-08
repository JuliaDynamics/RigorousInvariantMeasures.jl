using InvariantMeasures
using ValidatedNumerics

D = Mod1Dynamic(x->2*x)
B = Hat(EquispacedPartition{Float64}(8))
P = assemble(B, D)

Ptrue = [
        0.5 0.25 0    0    0    0    0   0.25;
        0   0.25 0.5  0.25 0    0    0   0   ;
        0   0    0    0.25 0.5  0.25 0   0   ;
        0   0    0    0    0    0.25 0.5 0.25;
        0.5 0.25 0    0    0    0    0   0.25;
        0   0.25 0.5  0.25 0    0    0   0   ;
        0   0    0    0.25 0.5  0.25 0   0   ;
        0   0    0    0    0    0.25 0.5 0.25;
		]
Ptrue = Ptrue'

@test all(contains_zero.(P-Ptrue))

@test opnormbound(DiscretizedOperator(B, D), L1()) == 1
@test opnormbound(DiscretizedOperator(B, D), Linf()) == 1
