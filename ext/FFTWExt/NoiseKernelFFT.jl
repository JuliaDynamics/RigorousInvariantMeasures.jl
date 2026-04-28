# DiscretizedNoiseKernelFFT — the FFT-based uniform-noise kernel. Currently
# unused (no callers in src/, test/, or examples/), but kept here so the main
# module no longer pulls FFTW just to declare the type.

struct DiscretizedNoiseKernelFFT{S<:AbstractVector,T<:AbstractVector} <:
       RigorousInvariantMeasures.NoiseKernel
    v::S
    Mfft::T
    rad::Any
    P::Any
end

DiscretizedNoiseKernelFFT(v::Vector{Real}) =
    DiscretizedNoiseKernelFFT(v, FFTW.fft(v), 0, FFTW.plan_fft(mid.(v)))
DiscretizedNoiseKernelFFT(v::Vector{Interval{T}}) where {T} = DiscretizedNoiseKernelFFT(
    v,
    FFTW.fft(mid.(v)),
    opnormbound(L2, radius.(v)),
    FFTW.plan_fft(mid.(v)),
)
Mfft(Q::DiscretizedNoiseKernelFFT) = Q.Mfft

# Ulam discretization of uniform noise of size ξ on a partition of size k.
function UniformNoiseFFT(ξ, k, boundarycondition = :periodic)
    n = Int64(floor(ξ * k))
    if boundarycondition == :periodic
        v = Interval.([ones(n); zeros(k - n)]) * (1 / (2 * ξ))
        v += Interval.([zeros(k - n); ones(n)]) * (1 / (2 * ξ))
        v[n+1] += @interval (ξ - n * 1 / k) * k / (2 * ξ)
        v[k-n-1] += @interval (ξ - n * 1 / k) * k / (2 * ξ)
    end
    return DiscretizedNoiseKernelFFT(v)
end

function Base.:*(M::DiscretizedNoiseKernelFFT, v)
    P = M.P
    w = P * v
    @info w
    w = M.Mfft .* w
    @info w
    return real.(P \ w) / length(v)
end
