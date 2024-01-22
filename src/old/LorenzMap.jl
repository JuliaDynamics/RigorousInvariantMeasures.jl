# COV_EXCL_START
module LorenzMapDefinition

using ..PwDynamicDefinition: PwMap




export LorenzMap

struct LorenzMap
    D::PwMap
end

LorenzMap(θ, α) = PwMap(
    [x -> θ * (0.5 - x)^α, x -> 1 - θ * (x - 0.5)^α],
    [@interval(0), @interval(0.5), @interval(1)],
)

import TaylorSeries
function dfly(::Type{TotalVariation}, ::Type{L1}, D::LorenzMap)
    #this function uses TaylorSeries to compute 1/T' and its derivatives
    inv_der_and_derivatives(x) = 1 / TaylorSeries.derivative(D.Ts[1](Taylor1([x, 1], 2)))

    inv_der(x) = inv_der_and_derivatives(x)[0]
    distortion(x) = inv_der_and_derivatives(x)[1]

    zero_inv_der(x) = inv_der(x) + 2 * inv_der(Interval(0))


end


end
# COV_EXCL_STOP