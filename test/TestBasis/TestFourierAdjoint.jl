using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Fourier assembler: Adjoint" begin
    B = RigorousInvariantMeasures.FourierAdjoint(128, 1024)
    @test length(B) == 257
    @test lastindex(B) == 257

    v = zeros(128)

    G(x) = 2 * x

    P = RigorousInvariantMeasures.assemble(B, G)

    real_P = real.(P)

    M = zeros(257, 257)
    M[1, 1] = 1.0

    for i = 2:129
        if (i - 1) % 2 == 0
            M[(i-1)÷2+1, i] = 1.0
        end
    end
    for i = 1:128
        if i % 2 == 0
            M[257-(i÷2)+1, 257-i+1] = 1.0
        end
    end

    @test all(M .∈ real_P)

    D = mod1_dynamic(x -> 2 * x)

    P = RigorousInvariantMeasures.assemble(B, D)

    @test all(M .∈ real_P)

end
