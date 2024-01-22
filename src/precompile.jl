# following the tutorial on https://julialang.org/blog/2021/01/precompile_tutorial/
# we implement some precompile solutions

precompile(DiscretizedOperator, (Ulam, PwMap))
precompile(DiscretizedOperator, (Hat, PwMap))
precompile(
    powernormbounds,
    (Ulam, PwMap, Int64, Int64, IntegralPreservingDiscretizedOperator),
)
precompile(
    powernormbounds,
    (Hat, PwMap, Int64, Int64, NonIntegralPreservingDiscretizedOperator),
)
precompile(invariant_vector, (Ulam, IntegralPreservingDiscretizedOperator))
precompile(invariant_vector, (Hat, NonIntegralPreservingDiscretizedOperator))
