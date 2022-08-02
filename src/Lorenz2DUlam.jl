module Lorenz2D

const α = 1.75::Float64
const s = 3.0::Float64

"""
    ContractingLorenz1D(; α , s)

Construct a dynamic representing the one-dimensional map for the contracting Lorenz Flow,
with parameters `α` and `s`, i. e.

``T(x) = -α*(-x)^s+1`` if ``-1<x<0`` and ``T(x) = α*(x)^s-1`` for ``0<x<1``
"""

import RigorousInvariantMeasures: PwMap

ContractingLorenz1D(; α , s) = PwMap([x-> -α*(-x)^s+1, x-> α*(x)^s-1], [-1, 0, 1])

T = ContractingLorenz1D(α = α , s = s)

"""
    FiberMapLorenz1D(; x, r, c)

For a fixed ``x`` return the fiber map of the skew product, with parameter ``r`` and ``c``,
i.e., if ``x>0`` we have ``G(x, y) = 2^(-r)y x^r+c``
and for ``x<0`` we have ``G(x, y) = 2^(-r)*y*(-x)^r-c``

This map is represented as a Branch object, where we pass the function, the value of the ``y``-derivative
(which is a constant depending on ``x``)
"""
# to simplify implementation, we use the Julia Polyhedra library
import Polyhedra as PH
import GLPK
lib = PH.DefaultLibrary{Float64}(GLPK.Optimizer)

import IntervalArithmetic
import IntervalArithmetic: Interval, mid

import RigorousInvariantMeasures: preimages

Ginverse(; x, r, c) = x > 0 ?
                v -> ((v-c)*2^r)/(x^r) :
                v -> ((v+c)*2^r)/((-x)^r)




"""
    PreimageRectangleLorenz(; branch, preim_x_left, preim_x_right, y_lower, y_upper, k)


"""
function PreimageRectangleLorenz2(;  preim_x_left::Interval, 
                                    preim_x_right::Interval,
                                    y_lower, 
                                    y_upper, 
                                    k)

    preim_x = [preim_x_left+(preim_x_right-preim_x_left)/k*i for i in 0:k]

    err_x = maximum(IntervalArithmetic.radius.(preim_x))
    err_y = 0.0

    y_s = NTuple{2, Interval}[]

    for x in preim_x
        G_inv = Ginverse(x = x, r = 5.0, c = 0.5)
        preim_y_low = min(max(G_inv(y_lower), -1), 1)
        preim_y_up = max(min(G_inv(y_upper), 1), -1)
        err_y = maximum([err_y, IntervalArithmetic.radius(preim_y_low), IntervalArithmetic.radius(preim_y_up)])
        push!(y_s, (preim_y_low, preim_y_up))
    end
    n = length(preim_x)
    A = Matrix{Float64}(undef, 2*n, 2)
    for (i, x) in enumerate(preim_x)
        A[i, :] = [mid(x) mid(y_s[i][1])]
    end
    
    for (i, x) in enumerate(reverse(preim_x))
        A[n+i, :] = [mid(x) mid(y_s[end-i+1][2])]
    end
    
    P = PH.polyhedron(PH.vrep(A), lib)
    return P, max(err_x, err_y)
end



"""
    PreimageRectangleLorenz(; T, x_left, x_right, y_lower, y_upper, k)

Return a Polyhedra that approximates the preimage of the rectangle ``[x_{left}, x_{right}] × [y_{lower}, y_{upper}]``
through the skew product map ``F(x, y) = (T(x), G(x, y)) and the maximum of the error in the computation 
of the vertices.
"""
function PreimageRectangleLorenz(; T, x_left, x_right, y_lower, y_upper, k)
    x = [Interval(-1); [Interval(x_left)+i*Interval(x_right-x_left)/k for i in 0:k]; Interval(1)]
    
    
    preim_x = preimages(x, T.branches[2])
    @info preim_x
    
    err_x = maximum(IntervalArithmetic.radius.(preim_x[1]))
    y_s = NTuple{2, Interval}[]
    err_y = 0.0
    for x in preim_x[1][2:end]
        G_inv = Ginverse(x = x, r = 5.0, c = 0.5)
        preim_y_low = min(max(G_inv(y_lower), -1), 1)
        preim_y_up = max(min(G_inv(y_upper), 1), -1)
        err_y = maximum([err_y, IntervalArithmetic.radius(preim_y_low), IntervalArithmetic.radius(preim_y_up)])
        push!(y_s, (preim_y_low, preim_y_up))
    end
    n = length(preim_x[1][2:end])
    A = Matrix{Float64}(undef, 2*n, 2)
    for (i, x) in enumerate(preim_x[1][2:end])
        A[i, :] = [mid(x) mid(y_s[i][1])]
    end

    for (i, x) in enumerate(reverse(preim_x[1][2:end]))
        A[n+i, :] = [mid(x) mid(y_s[end-i+1][2])]
    end

    P = PH.polyhedron(PH.vrep(A), lib)
    return P, max(err_x, err_y)
end

end