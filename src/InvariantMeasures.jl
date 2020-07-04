module InvariantMeasures

# the module Contractors does not depend on any submodule
include("Contractors.jl")
include("Partition.jl")

include("DynamicDefinition.jl")
# BasisDefinition uses the abstract Dynamic class
include("BasisDefinition.jl")
# Generic assembler
include("GenericAssembler.jl")
include("GenericEstimate.jl")

include("Mod1Dynamic.jl")

include("UlamBasis.jl")




using .DynamicDefinition, .BasisDefinition, .Mod1DynamicDefinition, .UlamBasis, .GenericAssembler, .GenericEstimate

#include("Hat.jl")
#include("GenericEstimate.jl")

export Ulam, Mod1Dynamic, Basis, Dynamic, assemble, preim

end
