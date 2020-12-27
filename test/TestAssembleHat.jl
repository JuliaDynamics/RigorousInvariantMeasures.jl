using InvariantMeasures
using ValidatedNumerics
using LinearAlgebra

@testset "Hat assembler" begin

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

@test opnormbound(L1, DiscretizedOperator(B, D)) == 1
@test opnormbound(Linf, DiscretizedOperator(B, D)) == 1

Q = DiscretizedOperator(B, D)
@test size(Q) == (8,8)

n = size(Q)[1]
e = randn(n)
x = similar(e)
x[:] .= 2.

mQ = mid(Q)
@test mul!(x, mQ, e, 1, 0) == mQ*e

fQ = mQ.L + mQ.e * mQ.w
@test fQ * e ≈ mQ * e

end
