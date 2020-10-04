using InvariantMeasures
using ValidatedNumerics

m = 30
m_extend = 100

D = Mod1Dynamic(x->2*x+0.5*x*(1-x))
B = Ulam(1024)
Q = DiscretizedOperator(B, D)

normQ = opnormbound(weak_norm(B), Q)

trivial_norms = norms_of_powers_trivial(weak_norm(B), Q, m)
computed_norms = norms_of_powers(weak_norm(B), m, Q, integral_covector(B))

(dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, m)

norms = min.(trivial_norms, computed_norms, dfly_norms) # in the current version, dfly_norms are always larger and can be omitted

better_norms = refine_norms_of_powers(norms, m_extend)

w = invariant_vector(B, Q)
@show distance_from_invariant(B, D, Q, w, better_norms)

B_fine = Ulam(2^20)
Q_fine = DiscretizedOperator(B_fine, D)
norm_Q_fine = opnormbound(weak_norm(B_fine), Q_fine)

trivial_norms_fine = norms_of_powers_trivial(weak_norm(B_fine), Q_fine, m_extend)
twogrid_norms_fine = norms_of_powers_from_coarser_grid(B_fine, B, D, better_norms, norm_Q_fine)

(dfly_strongs_fine, dfly_norms_fine) = norms_of_powers_dfly(B_fine, D, m_extend)

norms_fine = min.(trivial_norms_fine, twogrid_norms_fine, dfly_norms_fine)

better_norms_fine = refine_norms_of_powers(norms_fine, m_extend)

w_fine = invariant_vector(B_fine, Q_fine)
@show distance_from_invariant(B_fine, D, Q_fine, w_fine, better_norms_fine)
