using RigorousInvariantMeasures, IntervalArithmetic
using Test, Documenter

DocMeta.setdocmeta!(
    RigorousInvariantMeasures,
    :DocTestSetup,
    :(using RigorousInvariantMeasures);
    recursive = true,
    warn = false,
)

@testset "RigorousInvariantMeasures.jl" begin


    include("TestDifferentiation.jl")
    include("TestBasisDefinition.jl")
    include("TestContractors.jl")
    include("TestDynamic.jl")
    include("TestRun/TestRunIndex.jl")
    include("TestChebyshev.jl")

    include("TestBasis/TestBasisIndex.jl")
    include("TestInterface/TestInterfaceIndex.jl")

    include("TestInducedLSV.jl")
    include("TestPwMapInducedLSV.jl")

    include("TestTrig.jl")

    #include("TestSkewProductMap.jl")
    #include("TestUlam2DSP.jl")
    #include("TestHigherDFLY.jl")



    # @testset "Doctests" begin 
    #     if Base.VERSION >= v"1.8"
    #         # It seems that some output formats have changed from 1.7 to 1.8
    #         # therefore we use this hack to avoid failing doctests
    #         @info "The doctest implemented for version $(Base.VERSION)"
    #         doctest(RigorousInvariantMeasures)
    #     end
    # end

end
