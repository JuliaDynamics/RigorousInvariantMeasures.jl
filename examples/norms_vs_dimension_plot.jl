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

dims = 2 .^ (1:12)
num_norms = 30

norms = fill(NaN, (num_norms, length(dims)))

for (i, d) in enumerate(dims)
    B = Ulam(d)
    Q = DiscretizedOperator(B, D)
    norms[:,i] = norms_of_powers(weak_norm(B), num_norms, Q, integral_covector(B))
end

plot(norms,
    label = string.(dims'),
    yscale= :log10,
    legend = :topright,
#    ylims = (1e-6, 1)
    xlabel = L"k",
    ylabel = L"$\|Q^k|_U\|$",
    title = "Norm bounds for discretizations of different size"
)
