using InvariantMeasures
using ValidatedNumerics

using Plots
using LaTeXStrings
using StatsPlots

"""
Simple and quick experiment to precompile all needed functions before timing measurements
"""

D0 = mod1_dynamic(x->2*x+0.5*x*(1-x))
D = D0 ∘ D0
B = Ulam(8)
Q = DiscretizedOperator(B, D)

norms = powernormbounds(B, D; Q=Q)

B_fine = Ulam(128)
Q_fine = DiscretizedOperator(B_fine, D)

norms_fine = finepowernormbounds(B, B_fine, D, norms; Q_fine=Q_fine)
w_fine = invariant_vector(B_fine, Q_fine)
error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)

D = mod1_dynamic(x -> 4*x + 0.01*InvariantMeasures.sinpi(8*x))
B = Hat(64)
Q = DiscretizedOperator(B, D)

norms = powernormbounds(B, D; Q=Q)

B_fine = Hat(128)
Q_fine = DiscretizedOperator(B_fine, D)

norms_fine = finepowernormbounds(B, B_fine, D, norms; Q_fine=Q_fine)
w_fine = invariant_vector(B_fine, Q_fine)
error_fine = distance_from_invariant(B_fine, D, Q_fine, w_fine, norms_fine)
