# `Interval{Float64}` transform. The midpoint/radius split and the conversion
# back to rectangular enclosures live in `src/IntervalFFTCommon.jl`; the error
# bound comes from `BallArithmetic`'s `fft(::BallVector)`. `GenericFFTExt`
# supplies the arbitrary-precision counterpart.

@doc raw"""
    interval_fft(v::Vector{Complex{Interval{Float64}}})
    interval_fft(v::Vector{Interval{Float64}})

Rigorously enclosing FFT of `v`, normalized by `length(v)`, using FFTW for the
floating-point transform of the midpoints.

FFTW computes its twiddle factors to full accuracy, so the Brisebarre–Muller–
Picot bound applies directly.
"""
interval_fft(v::Vector{Complex{Interval{Float64}}}) =
    interval_fft_via_ball(FFTW.fft, v)

interval_fft(v::Vector{Interval{Float64}}) = interval_fft(v .+ 0im)
