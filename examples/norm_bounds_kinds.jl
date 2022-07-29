using RigorousInvariantMeasures
using IntervalArithmetic

using Plots
using LaTeXStrings

D = mod1_dynamic(x -> 4*x + 0.01*RigorousInvariantMeasures.sinpi(8*x))
#D = Mod1Dynamic(x -> 16*x + 0.01*RigorousInvariantMeasures.sinpi(32*x))
#D = Mod1Dynamic(x->2*x+0.5*x*(1-x))
# D = PwMap(
# [x -> 17x/5,
#  x -> 34(x-5//17)^2/25 + 3(x-5//17),
#  x -> 34(x-10//17)^2/25 + 3(x-10//17),
#  x -> 17(x-15//17)/5
# ],
# [0, @interval(5/17), @interval(10/17), @interval(15/17), 1]
# )

num_norms = 30

B = Hat(1024)
Q = DiscretizedOperator(B, D)

normQ = opnormbound(B, weak_norm(B), Q)
trivial_norms = norms_of_powers_trivial(normQ, num_norms)
computed_norms = norms_of_powers(B, weak_norm(B), num_norms, Q, integral_covector(B))

(dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, num_norms)

norms = min.(trivial_norms, computed_norms, dfly_norms)
better_norms = refine_norms_of_powers(norms, num_norms)

pgfplotsx()
p = plot(trivial_norms,
    label = L"$\|Q\|^k$",
    yscale = :log10,
    legend = :bottomleft,
    title = "Various bounds for norms of powers",
    xlabel = L"$k$",
    ylabel = L"bound to $\|Q^k|_U\|$"
    )
plot!(p, computed_norms,
    label = "computational bounds")

plot!(p, dfly_norms,
    label = L"DFLY $2\times 2$ matrix bounds")

plot!(p, better_norms,
    label = "min(previous) + refinement")

# savefig(p, "norm_bounds_kinds.tikz")
