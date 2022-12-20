using RigorousInvariantMeasures

using LinearAlgebra
using IntervalArithmetic

@testset "Hat assembler" begin

using RigorousInvariantMeasures: L1, Linf

D = mod1_dynamic(x->2*x)
B = Hat(8)
P = RigorousInvariantMeasures.assemble(B, D; ϵ = 0.0, max_iter = 100, T = Float64)

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

# @test opnormbound(B,L1, DiscretizedOperator(B, D)) == 1 # not defined anymore now that we include B in the signature
@test opnormbound(B, Linf, DiscretizedOperator(B, D)) >= 1

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
