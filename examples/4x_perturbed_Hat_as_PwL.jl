using InvariantMeasures
using ValidatedNumerics
using TaylorSeries

m = 15
m_extend = 100
eps = 1e-12

T = x-> 4*x + 0.01*InvariantMeasures.sinpi(8*x)
Tprime = x-> T(TaylorSeries.Taylor1([x, 1], 1))[1]
a1 = InvariantMeasures.root(x->T(x)-1, Tprime, 0..1, eps)
a2 = InvariantMeasures.root(x->T(x)-2, Tprime, 0..1, eps)
a3 = InvariantMeasures.root(x->T(x)-3, Tprime, 0..1, eps)

D = PwMap([x->T(x), x->T(x)-1, x->T(x)-2, x->T(x)-3],
	[Interval(0), a1, a2, a3, Interval(1)], fill(true, 4))

B = Hat(1024)
Q = DiscretizedOperator(B, D)

normQ = opnormbound(weak_norm(B), Q)

trivial_norms = norms_of_powers_trivial(weak_norm(B), Q, m)
computed_norms = norms_of_powers(weak_norm(B), m, Q, integral_covector(B))

(dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, m)

norms = min.(trivial_norms, computed_norms, dfly_norms) # in the current version, dfly_norms are always larger and can be omitted

better_norms = refine_norms_of_powers(norms, m_extend)

w = invariant_vector(B, Q)
@show distance_from_invariant(B, D, Q, w, better_norms)

B_fine = Hat(2^16)
Q_fine = DiscretizedOperator(B_fine, D)
norm_Q_fine = opnormbound(weak_norm(B_fine), Q_fine)

trivial_norms_fine = norms_of_powers_trivial(weak_norm(B_fine), Q_fine, m_extend)
twogrid_norms_fine = norms_of_powers_from_coarser_grid(B_fine, B, D, better_norms, norm_Q_fine)

(dfly_strongs_fine, dfly_norms_fine) = norms_of_powers_dfly(B_fine, D, m_extend)

norms_fine = min.(trivial_norms_fine, twogrid_norms_fine, dfly_norms_fine)

better_norms_fine = refine_norms_of_powers(norms_fine, m_extend)

w_fine = invariant_vector(B_fine, Q_fine)
@show distance_from_invariant(B_fine, D, Q_fine, w_fine, better_norms_fine)
