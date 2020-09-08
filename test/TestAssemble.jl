using InvariantMeasures
using ValidatedNumerics

D = Mod1Dynamic(x->2*x)
B = Ulam(8)
P = assemble(B, D)

Ptrue = [
		0.5 0.5 0 0 0 0 0 0;
		0 0 0.5 0.5 0 0 0 0;
		0 0 0 0 0.5 0.5 0 0;
		0 0 0 0 0 0 0.5 0.5;
		0.5 0.5 0 0 0 0 0 0;
		0 0 0.5 0.5 0 0 0 0;
		0 0 0 0 0.5 0.5 0 0;
		0 0 0 0 0 0 0.5 0.5;
		]
Ptrue = Ptrue'


@test all(contains_zero.(P-Ptrue))

@test opnormbound(L1(), DiscretizedOperator(B, D)) == 1
@test opnormbound(Linf(), DiscretizedOperator(B, D)) == 1
