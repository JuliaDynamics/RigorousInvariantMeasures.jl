# Noise kernel and noise-specialized estimates

The `UniformKernelUlam` family discretizes additive uniform noise on the
Ulam basis. Two boundary policies are exposed: `UniformKernelUlamPeriodic`
(circle) and `UniformKernelUlamReflecting` (interval with reflecting
boundary).

```@autodocs
Modules = [RigorousInvariantMeasures]
Pages = ["NoiseKernel.jl", "UniformNoiseUlam.jl"]
```

## A posteriori noise error estimates

The routines below produce data-dependent error bounds for systems with
additive uniform noise, exploiting the actually-computed density vector
to obtain tighter bounds than the generic DFLY-based
`distance_from_invariant_noise`.

```@autodocs
Modules = [RigorousInvariantMeasures]
Pages = ["NoiseSpecializedEstimate.jl"]
```

## Visualisation (extension)

```@docs
RigorousInvariantMeasures.plot_noisy_system
```
