using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Fourier assembler: Common" begin
    coeff = [0.5; 0.0; 0.5]

    @test evalFourier(coeff, 0.0) == 1.0
    @test real(evalFourier(coeff, 0.5)) == 0.0

    Bc = RigorousInvariantMeasures.FourierAdjoint(128, 1024)
    Bf = RigorousInvariantMeasures.FourierAdjoint(256, 1024)

    @test RigorousInvariantMeasures.is_refinement(Bc, Bf) == true
    @test RigorousInvariantMeasures.integral_covector(Bc) == [1.0; zeros(256)]'
    @test RigorousInvariantMeasures.one_vector(Bc) == [1.0; zeros(256)]

    S = AverageZero(Bc)

    @test length(S) == 256

    (v, state) = iterate(S)

    w = zeros(257)
    w[2] = 1.0
    @test v == w

    (v, state) = iterate(S, state)

    w = zeros(257)
    w[3] = 1.0
    @test v == w


end
