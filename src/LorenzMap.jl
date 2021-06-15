module LorenzMapDefinition

using ..PwDynamicDefinition: PwMap
using ValidatedNumerics



export LorenzMap

struct LorenzMap
    D::PwMap
end

LorenzMap(θ, α) = PwMap([x->θ*abs(x-0.5)^α, x->1-θ*abs(x-0.5)^α],
                    [@interval(0), @interval(0.5), @interval(1)])

import ..DynamicDefinition: derivative



end
