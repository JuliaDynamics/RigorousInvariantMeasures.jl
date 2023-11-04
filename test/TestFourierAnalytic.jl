using RigorousInvariantMeasures
using IntervalArithmetic

@testset "Fourier assembler" begin
    B = RigorousInvariantMeasures.AnalyticFourierBasis.FourierAnalytic(128, 1024)
    length(B) == 128

    v = zeros(128)
    
end
