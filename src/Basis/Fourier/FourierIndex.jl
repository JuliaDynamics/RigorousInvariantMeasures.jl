include("FourierCommon.jl")
include("FourierAdjoint.jl")
using .AdjointFourierBasis
export FourierAdjoint

include("FourierAnalytic.jl")
using .AnalyticFourierBasis
export FourierAnalytic
