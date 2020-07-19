using InvariantMeasures
using Test

@testset "InvariantMeasures.jl" begin
    # Write your tests here.

include("TestDynamic.jl")
include("TestAssemble.jl")
include("TestAssembleHat.jl")
include("TestEstimate.jl")


end
