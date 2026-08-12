"""
    GenericFFTExt

Triggered by `using GenericFFT`. Supplies the arbitrary-precision
(`Interval{BigFloat}`) method of `interval_fft`, so the Fourier and Chebyshev
assemblers can run at precisions FFTW cannot reach.

`GenericFFT` reexports `FFTW`, so loading it also brings up `FFTWExt` and the
`Float64` path.
"""
module GenericFFTExt

using RigorousInvariantMeasures
using IntervalArithmetic
using BallArithmetic
using GenericFFT

import RigorousInvariantMeasures: interval_fft, interval_fft_via_ball

@doc raw"""
    interval_fft(v::Vector{Complex{Interval{BigFloat}}})
    interval_fft(v::Vector{Interval{BigFloat}})

Rigorously enclosing FFT at arbitrary precision. The floating-point transform
of the midpoints is GenericFFT's; the enclosure is `BallArithmetic`'s
Brisebarre–Muller–Picot bound, which is generic in the float type.

**Why the precision is raised.** The BMP bound assumes a radix-2 transform
whose twiddle factors are correctly rounded (relative error `≤ u`).
GenericFFT's `generic_fft_pow2!` instead generates them with the classical
Numerical-Recipes trigonometric recurrence, whose drift over the `m` steps of a
butterfly row is *not* covered by that assumption — applying the bound at the
working precision would therefore not be rigorous.

The transform is consequently run at twice the working precision, while the
bound is evaluated at the working precision. Any twiddle drift is then bounded
by the recurrence length times the raised-precision unit roundoff, i.e. by
roughly `N² · u²` in terms of the working `u` — many orders of magnitude below
the reported radius, so the reported bound is valid with room to spare. (At
`N = 2^14` and 256 bits: drift `≲ 1e-150` against a radius `~1e-72`.)

Non-power-of-two lengths route through GenericFFT's Bluestein path, which is
not radix-2 at all; `BallArithmetic` warns, and the result should not be
treated as certified.
"""
function interval_fft(v::Vector{Complex{Interval{BigFloat}}})
    return setprecision(BigFloat, 2 * precision(BigFloat)) do
        interval_fft_via_ball(GenericFFT.AbstractFFTs.fft, v)
    end
end

interval_fft(v::Vector{Interval{BigFloat}}) = interval_fft(v .+ 0im)

end  # module GenericFFTExt
