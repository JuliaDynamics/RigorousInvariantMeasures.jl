"""
    FFTWExt

Triggered by `using FFTW`. Provides the FFT-based assembly methods for the
Fourier and Chebyshev bases, and the FFT-based uniform-noise kernel
(`DiscretizedNoiseKernelFFT`). Without this extension loaded, the basis
types are still constructible and most introspection methods (`length`,
`getindex`, `weak_norm`, …) work, but `assemble(B, D)` calls — and any
pipeline that goes through them — raise a `MethodError`.

The rigorously-enclosed FFT (`interval_fft`, internal) uses FFTW for the
floating-point transform of the per-element midpoints and adds a
Higham 1996 a-priori bound on the per-entry error.
"""
module FFTWExt

using RigorousInvariantMeasures
using IntervalArithmetic
using FastRounding
using FFTW

import RigorousInvariantMeasures: assemble, assemble_common, opnormbound, L2, W,
    Fourier, FourierAdjoint, FourierAnalytic, Chebyshev, Dynamic,
    NoiseKernel, NormKind, Observable, ProjectedFunction

include("IntervalFFT.jl")
include("Fourier.jl")
include("FourierObservables.jl")
include("Chebyshev.jl")
include("NoiseKernelFFT.jl")

end  # module FFTWExt
