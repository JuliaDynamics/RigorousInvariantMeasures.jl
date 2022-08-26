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

# these are the coordinate change maps
ϕ = PwMap([x->2*x-1, ], [0, 1])
ψ = PwMap([x-> x/2+1/2, ], [-1, 1])

# this is the map on [0, 1], incredibly it works!!!
D = ψ∘T∘ϕ


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
    PreimageRectangleLorenz(; preim_x_left, preim_x_right, y_lower, y_upper, k)

Return a Polyhedra that approximates the preimage of the rectangle 
``[T(preimx_{left}), T(preimx_{right})] × [y_{lower}, y_{upper}]``
through the skew product map ``F(x, y) = (T(x), G(x, y)) and the maximum of the 
error in the computation of the vertices.
The argument `k` is the number of equispaced linear interpolation points in the
`x` direction.
"""
function PreimageRectangleLorenz(;  preim_x_left::Interval, 
                                    preim_x_right::Interval,
                                    y_lower, 
                                    y_upper, 
                                    k,
                                    r,
                                    c)
    preim_x = [preim_x_left+(preim_x_right-preim_x_left)/k*i for i in 0:k]

    err_x = maximum(IntervalArithmetic.radius.(preim_x))
    err_y = 0.0

    y_s = NTuple{2, Interval}[]

    for x in preim_x
        G_inv = Ginverse(x = x, r = r, c = c)
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

using IntervalArithmetic
function _Lorenz_one_dim_map(x::Interval, α, s)
    x_left = x ∩ @interval -1 0
    x_right = x ∩ @interval 0 1
    return -α*(-x_left)^s+1 ∪ α*(x_right)^s-1
end

_Lorenz_left_one_dim_map(x::Interval, α, s) = -α*(-(x ∩ @interval -1 0))^s+Interval(1)
_Lorenz_right_one_dim_map(x::Interval, α, s) = α*(x ∩ @interval 0 1)^s-Interval(1)
_Lorenz_left_fiber_map(x::Interval, y::Interval, r, c) = 2^(-r)*y*(-x)^r-c
_Lorenz_right_fiber_map(x::Interval, y::Interval, r, c) = 2^(-r)*y*x^r+c


function BoundingRandomAttractor(α, r, s, c, ξ)
    Dict1to2 = Dict{Int64, NTuple{2, Int64}}()
    Dict2to1 = Dict{NTuple{2, Int64}, Int64}()

    left_image_rectangle_x = _Lorenz_left_one_dim_map(Interval(-1,0), α, s)
    left_image_rectangle_y = _Lorenz_left_fiber_map(Interval(-1), Interval(-1, 1), r, c)
    @info left_image_rectangle_x
    @info left_image_rectangle_y
end


end