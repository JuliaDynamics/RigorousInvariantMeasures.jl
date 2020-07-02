module InvariantMeasures

# the module Contractors does not depend on any submodule
include("Contractors.jl")


include("DynamicDefinition.jl")
# BasisDefinition uses the abstract Dynamic class
include("BasisDefinition.jl")
include("GenericAssembler.jl")


include("Mod1Dynamic.jl")

include("UlamBasis.jl")

# Generic assembler


using .DynamicDefinition, .BasisDefinition, .Mod1DynamicDefinition, .UlamBasis, .GenericAssembler

#include("Hat.jl")
#include("GenericEstimate.jl")

export Ulam, Mod1Dynamic, Basis, Dynamic, assemble, preim

end
