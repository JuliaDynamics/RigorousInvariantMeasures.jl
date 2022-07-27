using RigorousInvariantMeasures
using Test

@testset "RigorousInvariantMeasures.jl" begin

    include("TestContractors.jl")
    include("TestDynamic.jl")
    include("TestHat.jl")
    include("TestHatNP.jl")
    include("TestUlam.jl")
    include("TestAssemble.jl")
    include("TestAssembleHat.jl")
    include("TestEstimate.jl")
    include("TestNormOfPowers.jl")
    include("TestAbstractConvergence.jl")
    include("TestPreimages.jl")
    include("TestChebyshev.jl")
end
