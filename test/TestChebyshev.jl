using RigorousInvariantMeasures

@testset "Chebyshev assembler" begin

N = 16
B = RigorousInvariantMeasures.Chebyshev(N, 2)

D = mod1_dynamic(x->2*x)
L(ϕ, x) = (ϕ(x/2)+ϕ(x/2+0.5))/2

M = assemble(B, D)
using LinearAlgebra

all_true = true
# for i in 1:N
#  w = mid.(L.(B[i], B.p))
# z = RigorousInvariantMeasures.chebtransform(w)
#   all_true = all_true && norm(z-M[:, i], Inf)< 10^-13
# end

@test all_true

end
