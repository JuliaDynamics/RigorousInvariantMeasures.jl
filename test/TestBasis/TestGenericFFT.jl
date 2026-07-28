using RigorousInvariantMeasures
using IntervalArithmetic
using GenericFFT

@testset "Arbitrary-precision interval_fft and assembly (GenericFFTExt)" begin
    @test Base.get_extension(RigorousInvariantMeasures, :GenericFFTExt) !== nothing

    setprecision(BigFloat, 256) do
        # FFT of the all-ones vector: N at the DC bin, 0 elsewhere.
        N = 8
        v = Complex{Interval{BigFloat}}[
            complex(interval(BigFloat, 1), interval(BigFloat, 0)) for _ = 1:N
        ]
        y = RigorousInvariantMeasures.interval_fft(v)
        @test in_interval(1, real(y[1]))
        @test all(in_interval(0, real(y[i])) for i = 2:N)
        # The point of the arbitrary-precision path: the enclosure is far
        # tighter than anything Float64 could give.
        @test radius(real(y[1])) < 1e-60

        # The Fourier assembler is generic, so it carries BigFloat through.
        D = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
        B = FourierAnalytic(8, 64, W{3,1}; T = BigFloat)
        M = RigorousInvariantMeasures.assemble(B, D; ϵ = 1e-20, max_iter = 1000)
        @test eltype(M) == Complex{Interval{BigFloat}}
        @test size(M) == (17, 17)
    end
end
