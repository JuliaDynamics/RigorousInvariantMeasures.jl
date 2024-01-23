using RigorousInvariantMeasures

using RigorousInvariantMeasures: assemble, L1, Linf
using IntervalArithmetic

@testset "Ulam assembler" begin

    D = mod1_dynamic(x -> 2 * x)
    B = Ulam(8)
    P = assemble(B, D; ϵ = 10^(-15), max_iter = 100, T = Float64)

    Ptrue = [
        0.5 0.5 0 0 0 0 0 0
        0 0 0.5 0.5 0 0 0 0
        0 0 0 0 0.5 0.5 0 0
        0 0 0 0 0 0 0.5 0.5
        0.5 0.5 0 0 0 0 0 0
        0 0 0.5 0.5 0 0 0 0
        0 0 0 0 0.5 0.5 0 0
        0 0 0 0 0 0 0.5 0.5
    ]
    Ptrue = Ptrue'

    @test all(contains_zero.(P - Ptrue))

    # mod1_dynamic with a non-Markov dynamic
    D = mod1_dynamic(x -> x + 0.5)
    B = Ulam(8)
    P = assemble(B, D; ϵ = 1e-13, max_iter = 100, T = Float64)

    Ptrue = [
        0 0 0 0 1 0 0 0
        0 0 0 0 0 1 0 0
        0 0 0 0 0 0 1 0
        0 0 0 0 0 0 0 1
        1 0 0 0 0 0 0 0
        0 1 0 0 0 0 0 0
        0 0 1 0 0 0 0 0
        0 0 0 1 0 0 0 0
    ]

    @test all(contains_zero.(P - Ptrue))

    @test opnormbound(B, L1, DiscretizedOperator(B, D)) >= 1
    # not defined anymore now that we include the basis
    # @test opnormbound(B, Linf, DiscretizedOperator(B, D)) == 1 

end
