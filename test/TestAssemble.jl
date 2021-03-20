using InvariantMeasures
using ValidatedNumerics

@testset "Ulam assembler" begin

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

# mod1_dynamic with a non-Markov dynamic
D = mod1_dynamic(x->x+0.5)
B = Ulam(8)
P = assemble(B, D)

Ptrue = [
0  0  0  0  1  0  0  0;
0  0  0  0  0  1  0  0;
0  0  0  0  0  0  1  0;
0  0  0  0  0  0  0  1;
1  0  0  0  0  0  0  0;
0  1  0  0  0  0  0  0;
0  0  1  0  0  0  0  0;
0  0  0  1  0  0  0  0;
		]

@test all(contains_zero.(P-Ptrue))

@test opnormbound(L1, DiscretizedOperator(B, D)) == 1
@test opnormbound(Linf, DiscretizedOperator(B, D)) == 1

end
