using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Fourier assembler: Analytic" begin
    B = RigorousInvariantMeasures.AnalyticFourierBasis.FourierAnalytic(128, 1024)
    length(B) == 128

    v = zeros(128)

    D = mod1_dynamic(x->2*x)

    P = RigorousInvariantMeasures.AnalyticFourierBasis.assemble_standard(B, D)

    real_P = real.(P)

    M = zeros(257, 257)
    M[1, 1] = 1.0

    for i in 2:129
        if (i-1) % 2 ==0 
            M[(i-1)÷2+1, i] = 1.0
        end
    end
    for i in 1:128
        if i % 2 ==0 
            M[257-(i÷2)+1, 257-i+1] = 1.0
        end
    end

    
    @test all(M .∈ real_P)




    
end
