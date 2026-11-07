# Attic

Code salvaged from branches that were deleted, kept so it can be read without
digging through the reflog. Nothing here is part of the package: no file in this
directory is `include`d from `src/RigorousInvariantMeasures.jl`, and none of it
has been checked against today's interfaces. Each file records the branch and
the commit it came from in its first lines.

## `FourierToby.Fourier2.jl`, `FourierToby.Fourier.jl`

Toby's Fourier basis, in the two versions that existed. `Fourier2` is the
earlier one, from `ibis/FourierToby` at commit `5308cdd`, which was the only ref
that carried it; `Fourier` is the later one, from `ibis/Lorenz2DRebased` at
commit `cff183c`, which is an ancestor of `main`. The later version renames the
basis, adds a two-dimensional basis `Fourier2D` with its evaluation array, adds
`ConvolutionMatrix`, `FourierCoefficientsError`, `ProjectorError` and
`OperatorError`, and repairs a set of methods that had been left dispatching on
`Chebyshev` after a copy. Both are kept because three things live only in the
earlier version: the half-weights `0.5*exp(±2πikx)` in the evaluation matrix,
the `Fourier2Dual` struct, and an `assemble` that builds the operator on a basis
of twice the size and folds the four corner blocks back, which the later version
comments out in favour of calling `assemble_standard` directly.

The package's own Fourier basis is `src/Basis/Fourier/`, and it is unrelated to
these files.

## `expAlexToby.jl`, `TobyAlex.ipynb`

Spectral-flatness experiments on the noise operator, from the GitHub branch
`spectralflatness` at commit `056f8cb`, which was the only ref that carried
them, in any form.

## `dfly_inf_der_function_barrier.jl`

The function-barrier split of `dfly_inf_der`, from the GitHub branch `master` at
commit `41c661d`. `master` is kept, so this is a copy for reading rather than a
rescue. It was not carried into `src/PwDynamic.jl`, which keeps the loop inline:
the split also changes the interface, scanning the thirteen candidate widths at
a coarse tolerance and recomputing only the best one at a fine tolerance, and
`main` has since moved the dispatch from `TotalVariation` to `VarInterval` and
the return value to a `DFLYCoefficients`. Reconciling the two means comparing
the coefficients the current routine returns against the ones it returned
before, on every map the tests cover, which is a run and not a transcription.
