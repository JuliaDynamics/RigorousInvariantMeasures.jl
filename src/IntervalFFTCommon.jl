using IntervalArithmetic
# Selective import: a bare `using BallArithmetic` would shadow
# IntervalArithmetic's `sup` inside this module.
using BallArithmetic: BallVector, add_up, mul_up, sqrt_up

export interval_fft

@doc raw"""
    interval_fft(v)

Rigorously enclosing FFT of `v`, normalized by `length(v)`.

The generic function lives in the package so that more than one extension can
supply the underlying floating-point transform:

- `FFTWExt` (`using FFTW`) handles `Interval{Float64}`.
- `GenericFFTExt` (`using GenericFFT`) handles `Interval{BigFloat}`.

Without one of those loaded, callers hit a `MethodError`.

The error bound itself is not computed here: both methods route through
[`interval_fft_via_ball`](@ref), which defers to `BallArithmetic`'s
`fft(::BallVector)` and its Brisebarre–Muller–Picot (ARITH 2023) a-priori
bound.
"""
function interval_fft end

@doc raw"""
    midradius_complex_interval(v) -> (midpoints, radii)

Per-element midpoint/radius split of a `Vector{Complex{Interval{T}}}`. The
radius is widened up so the disk it describes contains the rectangular
`Complex{Interval}` enclosure — which is the form `BallArithmetic` wants.
"""
function midradius_complex_interval(v::AbstractVector{Complex{Interval{T}}}) where {T}
    n = length(v)
    mid_vec = zeros(Complex{T}, n)
    rad_vec = zeros(T, n)
    @inbounds for i = 1:n
        re_m, re_r = mid(real(v[i])), radius(real(v[i]))
        im_m, im_r = mid(imag(v[i])), radius(imag(v[i]))
        mid_vec[i] = re_m + im * im_m
        rad_vec[i] = sqrt_up(add_up(mul_up(re_r, re_r), mul_up(im_r, im_r)))
    end
    return mid_vec, rad_vec
end

@doc raw"""
    interval_fft_via_ball(fftfun, v)

Shared body of [`interval_fft`](@ref): split `v` into midpoints and radii, hand
them to `BallArithmetic` as a `BallVector`, transform, and read the resulting
balls back as rectangular `Complex{Interval}` enclosures, normalized by
`length(v)`.

`fftfun` is `FFTW.fft` — a `BallVector` method for it is supplied by
`BallArithmetic`'s own `FFTExt`. It is passed in rather than called directly
because this file must not depend on FFTW. Which floating-point kernel that
method ends up using depends on what is loaded: FFTW for `Float64`, GenericFFT
for `BigFloat`.

The enclosure is the Brisebarre–Muller–Picot (ARITH 2023) bound, which assumes
a radix-2 transform with **correctly rounded** twiddle factors. FFTW satisfies
that; see `GenericFFTExt` for why the arbitrary-precision path computes at
raised precision to stay inside the same assumption.
"""
function interval_fft_via_ball(
    fftfun,
    v::AbstractVector{Complex{Interval{T}}},
) where {T}
    N = length(v)
    mid_vec, rad_vec = midradius_complex_interval(v)
    y = fftfun(BallVector(mid_vec, rad_vec))
    return [
        (
            interval(real(y.c[i]) - y.r[i], real(y.c[i]) + y.r[i]) +
            im * interval(imag(y.c[i]) - y.r[i], imag(y.c[i]) + y.r[i])
        ) / N for i = 1:N
    ]
end
