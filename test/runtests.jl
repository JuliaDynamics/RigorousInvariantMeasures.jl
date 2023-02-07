using RigorousInvariantMeasures
using Test, Documenter

@testset "RigorousInvariantMeasures.jl" begin

    include("TestDifferentiation.jl")
    include("TestBasisDefinition.jl")
    include("TestContractors.jl")
    include("TestDynamic.jl")
    include("TestHat.jl")
    include("TestHatNP.jl")
    include("TestUlam.jl")
    include("TestAssemble.jl")
    include("TestAssembleHat.jl")
    include("TestEstimate.jl")
    include("TestFullRunHat2x.jl")
    include("TestFullRunUlam2x.jl")
    include("TestNormOfPowers.jl")
    include("TestAbstractConvergence.jl")
    include("TestPreimages.jl")
    include("TestChebyshev.jl")
    include("TestObservables.jl")
    include("TestLorenz2DUlam.jl")
    include("TestHigherDFLY.jl")

    DocMeta.setdocmeta!(RigorousInvariantMeasures, :DocTestSetup, :(using RigorousInvariantMeasures); recursive=true)

    @testset "Doctests" begin 
        if Base.VERSION >= v"1.8"
            # It seems that some output formats have changed from 1.7 to 1.8
            # therefore we use this hack to avoid failing doctests
            @info "The doctest implemented for version $(Base.VERSION)"
            doctest(RigorousInvariantMeasures)
        end
    end

end


