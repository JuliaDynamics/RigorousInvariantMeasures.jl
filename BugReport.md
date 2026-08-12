# Bug Report: SymbolicsExt extension does not load

## Summary

The `SymbolicsExt` extension (which provides `dfly(W{k,1}, L1, D)` for k >= 2) never loads because `Symbolics` is listed in `[deps]` but not in `[weakdeps]` in `Project.toml`.

## Details

In `Project.toml`, the extension is declared as:
```toml
[extensions]
SymbolicsExt = ["Symbolics", "SymbolicUtils"]
```

However, `Symbolics` was only in `[deps]` (line 30) and `SymbolicUtils` was in both `[deps]` (line 29) and `[weakdeps]` (line 41). Julia's extension mechanism requires **all** trigger packages to be in `[weakdeps]`. Since `Symbolics` was missing from `[weakdeps]`, the extension never triggered.

Calling `dfly(W{3,1}, L1, D)` would always fall through to the generic fallback at `src/DFLY.jl:7` which returns `@error "Not implemented"`.

## Steps to Reproduce

```julia
using RigorousInvariantMeasures
D = mod1_dynamic(x -> 2x + 0.01 * sinpi(2x); full_branch=true)
A, B = dfly(W{3,1}, L1, D)
# Error: Not implemented
```

## Fix Applied

Moved `Symbolics` from `[deps]` only to `[weakdeps]` (it was already in `[extras]` for tests). The main module `src/RigorousInvariantMeasures.jl` does not use Symbolics directly — it's only used in the extension.

Additionally, several other packages that are only used in extensions (CUDA, Adapt, Plots, RecipesBase, LaTeXStrings, StatsPlots, TaylorModels) are also listed in both `[deps]` and `[weakdeps]`. These should probably be reviewed to determine if they can be removed from `[deps]` as well, since being in `[deps]` means they're always installed even when the extension isn't needed. This bloats the dependency tree unnecessarily.

---

# Bug Report: `Float64(Symbolics.Num)` fails in HigherDFLY.jl

## Summary

At `ext/SymbolicsExt/HigherDFLY.jl:200`, the code calls `Float64(Symbolics.substitute(...))` which returns a `Symbolics.Num`. There is no `Float64(::Num)` method, causing a `MethodError`.

## Steps to Reproduce

```julia
using RigorousInvariantMeasures, Symbolics
D = mod1_dynamic(x -> 2x + 0.01 * sinpi(2x); full_branch=true)
A, B = dfly(W{3,1}, L1, D)
# MethodError: no method matching Float64(::Symbolics.Num)
```

## Fix Applied

Changed line 200 from:
```julia
B_val = Float64(Symbolics.substitute(opt, Dict(...)))
```
to:
```julia
B_symbolic = Symbolics.substitute(opt, Dict(...))
B_val = Float64(Symbolics.unwrap(B_symbolic))
```

`Symbolics.unwrap` extracts the underlying numeric value from the `Num` wrapper, which can then be converted to `Float64`.

## Affected Versions

- RigorousInvariantMeasures v0.2.3
- Symbolics v6.58.0
- Julia 1.9+
