module InvariantMeasures

abstract type NormKind end
struct L1 <: NormKind end
struct Linf <: NormKind end
struct Lipschitz <: NormKind end
struct TotalVariation <: NormKind end

# the module Contractors does not depend on any submodule
include("Contractors.jl")
include("Partition.jl")

include("DynamicDefinition.jl")
include("BasisDefinition.jl")
include("GenericAssembler.jl")
include("GenericEstimate.jl")

include("Mod1Dynamic.jl")

include("UlamBasis.jl")
include("HatBasis.jl")

using .DynamicDefinition, .BasisDefinition, .Mod1DynamicDefinition, .UlamBasis, .GenericEstimate, .HatBasis, .Contractors

include("Norms.jl")
include("pitrig.jl")
include("NormsOfPowers.jl")

#include("Hat.jl")
#include("GenericEstimate.jl")

export Ulam, Mod1Dynamic, Basis, Dynamic, assemble, preim, Hat, EquispacedPartition, L1, Linf, Lipschitz, TotalVariation, norm_of_powers, sinpi, cospi, dfly, DiscretizedOperator, opnormbound

end
