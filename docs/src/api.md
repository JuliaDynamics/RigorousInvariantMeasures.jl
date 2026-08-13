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
Backend-independent parts of the interval FFT. The transforms themselves live
in the `FFTWExt` and `GenericFFTExt` extensions, loaded by `using FFTW` and
`using GenericFFT` respectively.

```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["IntervalFFTCommon.jl"]
```



