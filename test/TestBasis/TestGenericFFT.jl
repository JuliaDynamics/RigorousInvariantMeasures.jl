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

        # End-to-end arbitrary precision: with a BigFloat dynamic, the
        # preimage bisection, the frequency sweep and the FFT all stay in
        # BigFloat, so the assembled enclosure goes far below Float64 eps.
        # (This is what `unbounded_like` in Preimages.jl unblocked; the
        # placeholder used to be Interval{Float64} and pinned the whole
        # pipeline to double precision.)
        f(x) = 2 * x + interval(BigFloat, 1) / 2 * x * (1 - x)
        xb = (interval(BigFloat, 5) - sqrt(interval(BigFloat, 17))) / 2
        z, o = interval(BigFloat, 0), interval(BigFloat, 1)
        D = RigorousInvariantMeasures.PwMap(
            [x -> f(x), x -> f(x) - 1],
            [z, xb, o],
            [z o; z o];
            full_branch = true,
        )
        B = FourierAnalytic(8, 64, W{3,1}; T = BigFloat)

        cd = RigorousInvariantMeasures.Dual(B, D; ϵ = 1e-30, max_iter = 2000)
        @test cd.x[3] isa Interval{BigFloat}
        @test maximum(radius.(cd.x)) < 1e-28

        M = RigorousInvariantMeasures.assemble(B, D; ϵ = 1e-30, max_iter = 2000)
        @test eltype(M) == Complex{Interval{BigFloat}}
        @test size(M) == (17, 17)
        @test maximum(max(radius(real(w)), radius(imag(w))) for w in M) < 1e-28
    end
end
