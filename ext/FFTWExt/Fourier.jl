# Assemble methods for the Fourier basis family. Depend on `interval_fft`
# from IntervalFFT.jl.

function RigorousInvariantMeasures.assemble_common(
    B::Fourier,
    D;
    ϵ = 0.0,
    max_iter = 100,
    T = Float64,
)
    n = length(B)
    @info n
    k = (n - 1) ÷ 2
    @info k

    M = zeros(Complex{Interval{Float64}}, (n, n))
    computed_dual = RigorousInvariantMeasures.Dual(B, D; ϵ, max_iter)
    for i = 1:n
        ϕ = B[i]
        w = RigorousInvariantMeasures.eval_on_dual(B, computed_dual, ϕ)
        FFTw = interval_fft(w)
        M[:, i] = [FFTw[1:k+1]; FFTw[end-k+1:end]]
    end
    return M
end

function RigorousInvariantMeasures.assemble(
    B::FourierAdjoint,
    D;
    ϵ = 0.0,
    max_iter = 100,
    T = Float64,
)
    return RigorousInvariantMeasures.assemble_common(
        B,
        D;
        ϵ = 0.0,
        max_iter = 100,
        T = Float64,
    )'
end

function RigorousInvariantMeasures.assemble(
    B::FourierAnalytic,
    D::Dynamic;
    ϵ = 0.0,
    max_iter = 100,
    T = Float64,
)
    return RigorousInvariantMeasures.assemble_common(B, D; ϵ, max_iter, T)
end
