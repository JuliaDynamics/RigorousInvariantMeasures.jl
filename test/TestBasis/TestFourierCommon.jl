using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Fourier assembler: Common" begin
    coeff = [0.5; 0.0; 0.5]

    @test evalFourier(coeff, 0.0) == 1.0
    @test real(evalFourier(coeff, 0.5)) == 0.0

    Bc = RigorousInvariantMeasures.FourierAdjoint(128, 1024)
    Bf = RigorousInvariantMeasures.FourierAdjoint(256, 1024)

    @test RigorousInvariantMeasures.is_refinement(Bf, Bc) == true
    @test RigorousInvariantMeasures.is_refinement(Bc, Bf) == false
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

@testset "Fourier projection errors" begin
    RIM = RigorousInvariantMeasures
    # per unit of ||f^(s)||_1: (2π)^{-s} (2K^{1-2s}/(2s-1) + (2K+1)(2s/(s-1))^2 (n-K)^{-2s})^{1/2}
    for (K, n, s) in ((64, 256, 2), (16, 64, 4))
        B = FourierAnalytic(FourierPoints(n, Float64), K, W{s,1}(), L2())
        e = (2π)^(-s) * sqrt(2 * float(K)^(1 - 2s) / (2s - 1) + (2K + 1) * (2s / (s - 1))^2 * float(n - K)^(-2s))
        @test e <= RIM.weak_projection_error(B) <= e * (1 + 1e-12)
        @test RIM.aux_normalized_projection_error(B) == RIM.weak_projection_error(B)
        Bj = FourierAdjoint(FourierPoints(n, Float64), K, W{s,1}(), L2())
        @test RIM.weak_projection_error(Bj) == RIM.weak_projection_error(B)
    end
    @test_throws ArgumentError RIM.weak_projection_error(
        FourierAnalytic(FourierPoints(256, Float64), 64, W{1,1}(), L2()),
    )
    @test_throws ArgumentError RIM.weak_projection_error(
        FourierAnalytic(FourierPoints(64, Float64), 64, W{2,1}(), L2()),
    )
    # Cω: the tail 2q^{K+1}/(1-q), q = e^{-2π/10}, decreases in K
    e8 = RIM.weak_projection_error(FourierAdjoint(FourierPoints(64, Float64), 8))
    e16 = RIM.weak_projection_error(FourierAdjoint(FourierPoints(64, Float64), 16))
    q = exp(-2π / 10)
    @test 2q^9 / (1 - q) <= e8 <= 2q^9 / (1 - q) * (1 + 1e-12) && e16 < e8
end
