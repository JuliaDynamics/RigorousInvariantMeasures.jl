# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`RigorousInvariantMeasures.jl` is a Julia package (module name `RigorousInvariantMeasures`, not `InvariantMeasures` as some older README snippets suggest) that computes **rigorously certified** approximations of absolutely continuous invariant measures for 1D dynamical systems. "Rigorous" here means every bound is a validated interval-arithmetic enclosure: the output is a numerical answer paired with a proven upper bound on its error. Any change that loses that guarantee (e.g. replacing `Interval`/`BallMatrix` math with `Float64` on a load-bearing path) is a correctness bug, not a performance choice.

Julia compat is `julia = "1.9"`; CI runs on Julia 1.10 (Ubuntu, x64).

## Commands

All commands run from the repo root unless noted.

```bash
# Run the full test suite (what CI runs)
julia --project=. -e 'using Pkg; Pkg.test()'

# REPL with the package active (preferred dev loop — avoids recompilation)
julia --project=.

# Inside the REPL:
#   ] test                            # full test suite
#   include("test/TestUlam.jl")       # run a single test file (after `using RigorousInvariantMeasures, Test`)
#   include("examples/175_ulam.jl")   # run an example (most require Plots)

# Format (CI enforces SciML style via .JuliaFormatter.Toml)
julia -e 'using JuliaFormatter; format(".")'

# Build docs locally
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

To run a subset of tests, edit `test/runtests.jl` to comment out the `include(...)` lines you don't need — several tests are slow (minutes) because they run full rigorous pipelines.

`ENV["PROGRESS_BARS"] = "false"` silences ProgressMeter output in long computations (checked once at module load in `src/RigorousInvariantMeasures.jl` as `SHOW_PROGRESS_BARS`).

## Platform caveat

`src/pitrig.jl` re-exports `sinpi`/`cospi` backed by `CRlibm.jl`, which is **Linux-only**. Any dynamic defined with trig (including most `mod1_dynamic` examples) will fail to load on macOS/Windows. Don't rewrite these to use `Base.sinpi`: the CRlibm version is load-bearing for correctness — it guarantees `f(1) == 4` exactly in the canonical `4x + 0.01 sinpi(8x)` example, which is required for the Markov property of the discretized operator.

## Architecture

The whole package is one computational pipeline. A typical run composes these pieces in order:

```
Dynamic  ──►  Basis  ──►  DiscretizedOperator (Q)  ──►  powernormbounds  ──►  invariant_vector + distance_from_invariant
```

### Core abstractions (all defined as `abstract type` + multiple dispatch)

- **`Dynamic`** (`src/AbstractDynamicDefinition.jl`) — a 1D map `T: [0,1] → [0,1]`. Concrete subtypes live in `src/PwDynamic.jl` (`PwMap`, `MonotonicBranch`, `mod1_dynamic`) and `src/InducedLSV.jl`. Each subtype must implement `branches`, `nbranches`, `preim`, `endpoints`, `max_inverse_derivative`, `max_distortion`, etc.
- **`Basis`** (`src/Basis/BasisDefinition.jl`) — a finite-dimensional approximation space for densities. Concrete bases live in `src/Basis/`: `Ulam` (piecewise constant, `UlamBasis.jl`), `Hat`/`HatNP` (piecewise linear, `CircleHatBasis.jl`/`IntervalHatBasis.jl`), `C2Basis`, `Fourier`/`FourierAnalytic`/`FourierAdjoint` (`Basis/Fourier/`), `Chebyshev` (`Basis/NewChebyshev.jl`). Each basis declares its `weak_norm`, `strong_norm`, `aux_norm` (subtypes of `NormKind` from `src/Norms.jl` — `L1`, `L2`, `Linf`, `TotalVariation`, `Lipschitz`, `C1`, `W{k,l}`, `Aη`, `Cω`).
- **`DiscretizedOperator`** (`src/GenericAssembler.jl`) — the finite-rank projection of the transfer operator. Comes in `IntegralPreservingDiscretizedOperator` (a matrix) and `NonIntegralPreservingDiscretizedOperator` (sparse + rank-1 correction `L + e*w`). Assembled via `DiscretizedOperator(B, D)`.
- **`DFLY`** (`src/DFLY.jl`) — Doeblin-Fortet / Lasota-Yorke inequality. `dfly(strong_norm, aux_norm, D)` returns `(A, B)` with `‖Lf‖_s ≤ A‖f‖_s + B‖f‖_aux`. The method table is the main per-basis × per-dynamic dispatch surface; missing methods return `@error "Not implemented"`.

### Estimation pipeline (`src/GenericEstimate.jl`, `src/NormsOfPowers.jl`)

- `powernormbounds(B, D; Q)` computes `‖Q^k|_U‖` for `k = 1..k_max` (a posteriori mixing-time estimate; **this is the novel contribution of the underlying papers — don't bypass it with a priori bounds**).
- `finepowernormbounds(B_coarse, B_fine, D, coarse_norms; normQ_fine)` implements the coarse-fine scheme from ref [1] in the README.
- `invariant_vector(B, Q)` returns a non-rigorous eigenvector via Arpack's `eigs`.
- `distance_from_invariant(B, D, Q, w, norms)` is the rigorous L¹ error bound that closes the loop.

### Noise variant (separate pipeline)

`src/NoiseKernel.jl`, `src/UniformNoiseUlam.jl`, `src/NormsOfPowersNoise.jl`, `src/NoiseSpecializedEstimate.jl` implement the same flow for systems with additive uniform noise (refs [2], Nonlinearity 2020). Entry points: `UniformKernelUlam{Periodic,Reflecting}`, `powernormboundsnoise`, `invariant_vector_noise`, `distance_from_invariant_noise`.

### Weak-dependency extensions (`ext/`)

Four extensions loaded via `Project.toml` `[weakdeps]` + `[extensions]`:
- `PlotsExt` — plot recipes for bases, dynamics, noisy systems (`plot_noisy_system` stub declared in main module).
- `CUDAExt` — GPU-accelerated noise kernel.
- `SymbolicsExt` — higher-order DFLY via `Symbolics`/`SymbolicUtils` (`dfly(W{k,1}, L1, D)` for `k ≥ 2`).
- `TaylorModelsExt` — Taylor-model-based norm bounds.

Do **not** add a weakdep to `[deps]` just to make it load unconditionally — Julia's extension system requires trigger packages to be in `[weakdeps]` only. If an extension isn't loading, check that its triggers are all in `[weakdeps]` (see `BugReport.md` for an incident where `Symbolics` in `[deps]` silently broke `SymbolicsExt`).

## Tests

Tests are organized by feature area, not mirroring `src/` one-to-one:
- `test/TestBasis/` — per-basis unit tests (assembly, projection, norm bounds).
- `test/TestInterface/` — interface-level tests (assemble, estimate, preimages, observables).
- `test/TestRun/` — end-to-end rigorous pipeline runs (slow).
- Top-level `test/Test*.jl` — noise, induced-LSV, Chebyshev, differentiation, trig, higher DFLY.

A few `include(...)` lines in `test/runtests.jl` are commented out (Skew-product / Ulam2D); assume they are known-broken, not something to re-enable casually.

## Conventions that matter

- **Interval arithmetic is mandatory on rigorous paths.** `Float64` arithmetic is fine for the non-rigorous eigenvector step (`invariant_vector` calls `eigs` on `mid(Q)`) but the residual bound and norm bounds must stay in interval/ball arithmetic.
- **Directed rounding helpers**: `⊕₊`, `⊗₊`, `square_round(..., RoundUp)`, `Float64(x, RoundUp)` from `FastRounding` appear throughout — these preserve upper-bound validity when converting intervals back to floats. Don't replace them with `+`/`*`.
- **SciML formatter style** (`.JuliaFormatter.Toml`). A `Format.yml` workflow checks this on PRs.
- Long computations use `ProgressMeter`; gate new progress bars on `SHOW_PROGRESS_BARS`.
