using InvariantMeasures
using ValidatedNumerics

using Plots
using LaTeXStrings

D = Mod1Dynamic(x -> 4*x + 0.01*InvariantMeasures.sinpi(8*x))
#D = Mod1Dynamic(x -> 16*x + 0.01*InvariantMeasures.sinpi(32*x))
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

trivial_norms = norms_of_powers_trivial(weak_norm(B), Q, num_norms)
computed_norms = norms_of_powers(weak_norm(B), num_norms, Q, integral_covector(B))

(dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, num_norms)

norms = min.(trivial_norms, computed_norms, dfly_norms)
better_norms = refine_norms_of_powers(norms, num_norms)

p = plot(trivial_norms,
    label = "power of ||Q||",
    yscale = :log10,
    legend = :bottomleft,
    title = "Various bounds for norms of powers",
    xlabel = "k",
    ylabel = "bound to ||Q^k|_U||"
    )
plot!(p, computed_norms,
    label = "computational bound")

plot!(p, dfly_norms,
    label = "DFLY 2x2 matrix bound")

plot!(p, better_norms,
    label = "minimum + refinement")
