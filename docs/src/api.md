# Library

---

```@meta
CurrentModule = RigorousInvariantMeasures
```

## Contractors
These are interval contractors implemented here; in the future they may be substituted by established libraries, as [IntervalRootFinding.jl]()
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["Contractors.jl"]
```

## Preimages
Optimized methods to compute preimages of ``1`` dimensional maps.

```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["Preimages.jl"]
```

# Differentiation interface
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["differentiation_interface.jl"]
```

# Rigorously enclosed FFT

The rigorous error bound is **not** computed here: it comes from
`BallArithmetic`'s `fft(::BallVector)`, which implements the a-priori
Brisebarre–Muller–Picot (ARITH 2023) bound. What lives in this package is the
adapter — the `Complex{Interval}` ↔ `BallVector` conversion, the `1/N`
normalization, and a generic function so that each extension can supply the
floating-point kernel: `FFTWExt` (`using FFTW`) for `Interval{Float64}`,
`GenericFFTExt` (`using GenericFFT`) for `Interval{BigFloat}`.

```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["IntervalFFTCommon.jl"]
```



