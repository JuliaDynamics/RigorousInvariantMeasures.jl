using RigorousInvariantMeasures

# Defining the base dynamic

const α = 6.0::Float64
const s = 3.0::Float64

"""
    ContractingLorenz1D(; α , s)

Construct a dynamic representing the one-dimensional map for the contracting Lorenz Flow,
with parameters `α` and `s`, i. e.

``T(x) = -α*(0.5-x)^s+1`` if ``0<x<0.5`` and ``T(x) = α*(x-0.5)^s`` for ``0.5<x<1``
"""
ContractingLorenz1D(; α, s) =
    PwMap([x -> -α * (0.5 - x)^s + 1, x -> α * (x - 0.5)^s], [0, 0.5, 1])
