# Assemble for the Chebyshev basis. Uses `interval_fft` from IntervalFFT.jl
# via the `chebtransform` helper (DCT-via-FFT round-trip).

function chebtransform(w)
    n = length(w) - 1
    z = interval_fft([reverse(w); w[2:end-1]]) / n
    t = real.(z[1:length(w)])
    t[1] /= 2
    t[end] /= 2
    return Interval.(t)
end

function RigorousInvariantMeasures.assemble(
    B::Chebyshev,
    D::Dynamic;
    ϵ = 1e-13,
    max_iter = 100,
    T = Float64,
)
    n = length(B.p)
    M = zeros(Interval{T}, (n, n))
    x, labels, x′ = RigorousInvariantMeasures.Dual(B, D; ϵ, max_iter)
    for i = 1:n
        ϕ = B[i]
        w = zeros(Interval{Float64}, n)
        for j = 1:length(x)
            w[labels[j]] += ϕ(x[j]) / abs(x′[j])
        end
        M[:, i] = chebtransform(w)
    end
    return M
end
