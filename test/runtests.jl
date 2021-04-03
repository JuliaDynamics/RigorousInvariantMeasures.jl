using InvariantMeasures
using Test

@testset "InvariantMeasures.jl" begin

    include("TestContractors.jl")
    include("TestDynamic.jl")
    include("TestHat.jl")
    include("TestUlam.jl")
    include("TestAssemble.jl")
    include("TestAssembleHat.jl")
    include("TestEstimate.jl")
    include("TestNormOfPowers.jl")

end
