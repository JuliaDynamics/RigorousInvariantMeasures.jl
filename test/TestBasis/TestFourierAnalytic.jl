using RigorousInvariantMeasures
using IntervalArithmetic
@testset "Fourier assembler: Analytic" begin
    B = RigorousInvariantMeasures.FourierAnalytic(128, 1024)
    @test length(B) == 257

    v = zeros(128)

    D = mod1_dynamic(x -> 2 * x)

    P = RigorousInvariantMeasures.assemble(
        B,
        D;
        ϵ = 0.00000001,
        max_iter = 100,
        T = Float64,
    )

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


    @test all(in_interval.(M, real_P))

end

@testset "FourierAnalytic assemble agrees with the naive per-basis-function loop" begin
    # `assemble_common` sweeps all frequencies at once, using ϕ_m(x) = z^m and
    # the m ↦ -m conjugate symmetry, instead of evaluating a rigorous complex
    # exponential per (basis function, dual node) pair. This pins the fast path
    # to the definition it replaces, on a nonlinear map (so the dual nodes are
    # genuinely irrational).
    ext = Base.get_extension(RigorousInvariantMeasures, :FFTWExt)
    interval_fft = ext.interval_fft

    D = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))
    B = FourierAnalytic(24, 128, W{3,1})
    n = length(B)
    k = (n - 1) ÷ 2

    # Definitional assembler: column i is the FFT of eval_on_dual(B, ·, B[i]).
    M_naive = zeros(Complex{Interval{Float64}}, (n, n))
    cd = RigorousInvariantMeasures.Dual(B, D; ϵ = 1e-13, max_iter = 100)
    for i = 1:n
        F = interval_fft(RigorousInvariantMeasures.eval_on_dual(B, cd, B[i]))
        M_naive[:, i] = [F[1:k+1]; F[end-k+1:end]]
    end

    M_fast = RigorousInvariantMeasures.assemble(B, D; ϵ = 1e-13, max_iter = 100)

    @test size(M_fast) == size(M_naive)
    # Both are rigorous enclosures of the same matrix, so every entry must
    # intersect; neither need contain the other.
    @test all(
        !isempty_interval(intersect_interval(real(a), real(b))) &&
        !isempty_interval(intersect_interval(imag(a), imag(b))) for
        (a, b) in zip(M_naive, M_fast)
    )
    # ... and the enclosures must actually be tight, not vacuously overlapping.
    @test maximum(max(radius(real(z)), radius(imag(z))) for z in M_fast) < 1e-11
end
