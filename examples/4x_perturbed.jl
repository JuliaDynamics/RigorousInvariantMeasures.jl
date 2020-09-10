using InvariantMeasures
using ValidatedNumerics

m = 10

D = Mod1Dynamic(x -> 4*x + 0.01*InvariantMeasures.sinpi(8*x))
B = Hat(8)
Q = DiscretizedOperator(B, D)

normQ = opnormbound(weak_norm(B), Q)

norms = norms_of_powers(weak_norm(B), m, Q.L, false, e=Q.e, f=map(mid, integral_covector(B)))

w = invariant_vector(B, Q)
distance_from_invariant(B, D, Q, w, norms)
