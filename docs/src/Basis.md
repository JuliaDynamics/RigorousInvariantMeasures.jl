We present some of the basis already implemented in the package

# Abstract basis
```@autodocs
Modules = [Base, 
        RigorousInvariantMeasures]
Pages = ["BasisDefinition.jl"]
```

# The Ulam basis
The Ulam basis associated to a finite partition
`\{0 = x_0, \ldots, x_N = 1\}`, is given by 
the collection of characteristic functions
`\chi_{[x_i,x_{i+1})}` for `i` in `0, \ldots, N`.

The regularity seminorm associated to this basis
is the Variation seminorm, 
`
\textrm{Var}(\phi) = \sup_{\mathcal{P}} \sum |\phi(z_{i+1})-\phi(z_i)|
`

```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["UlamBasis.jl"]
```

# The hat basis on ``S^1``
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["CircleHatBasis.jl"]
```

# The hat basis on ``[0,1]``
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["IntervalHatBasis.jl"]
```

# The C2 basis
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["C2Basis.jl"]
```

# Fourier and Chebyshev bases (require `using FFTW`)

The Fourier (`Fourier`, `FourierAdjoint`, `FourierAnalytic`) and Chebyshev
(`Chebyshev`) basis types live in the main package, so they can be constructed
and inspected without `FFTW`. Their `assemble(B, D)` methods, the
`interval_fft` helper, and the optional FFT-based uniform-noise kernel
(`DiscretizedNoiseKernelFFT`) are provided by the `FFTWExt` extension and
become available once you load `using FFTW`. Calling `assemble` on a Fourier
or Chebyshev basis without `FFTW` loaded raises a `MethodError`.

## The Chebyshev basis
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["NewChebyshev.jl"]
```

## Common Fourier interface
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["FourierCommon.jl"]
```

## The Fourier Adjoint basis 
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["FourierAdjoint.jl"]
```

## The Fourier Analytic basis 
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["FourierAnalytic.jl"]
```
