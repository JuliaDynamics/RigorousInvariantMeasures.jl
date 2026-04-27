# Working with `IntervalArithmetic` 1.0

`RigorousInvariantMeasures` is built on `IntervalArithmetic.jl` 1.0+. Several
idioms that worked under `IntervalArithmetic` 0.20 were intentionally removed
upstream and are no longer available. This page is a reference for users
porting code (or for anyone reading our source who's used to the older API).

If a `RigorousInvariantMeasures` function returns an interval — for instance
`distance_from_invariant`, `dfly`, or any of the `*norm*` helpers — you almost
always want to **inspect** rather than compare it. The replacements below
cover the common patterns.

## Constructors

| Old (IA 0.20) | New (IA 1.0) | Notes |
|---|---|---|
| `Interval(a, b)` | `interval(a, b)` | The 2-argument capital constructor was removed. |
| `Interval{T}(a, b)` | `interval(T, a, b)` | Single-argument `Interval{T}(x::Real)` (typed conversion of a number) still works. |
| `a..b` | `interval(a, b)` | The `..` operator is gone. (Module-path syntax `using ..ParentModule` is unrelated.) |
| `IntervalArithmetic.atomic(Interval{T}, x)` | `interval(T, x)` | |
| `interval_from_midpoint_radius(m, r)` | `interval(m - r, m + r)` | The named helper was dropped. |

`Interval{T}` as a type annotation is unchanged — it still names the type used
throughout the package.

## Accessors

| Old | New | Notes |
|---|---|---|
| `x.lo`, `x.hi` | `inf(x)`, `sup(x)` | Field access is gone — `Interval` now stores `bareinterval`/`decoration`/`isguaranteed`. |
| `midpoint_radius(x)` | `midradius(x)` | Returns the same `(mid, rad)` tuple. |

A subtle gotcha: `inf(interval(0, 1))` returns `-0.0` (negative zero) in IA 1.0.
This matters whenever the result is used as a search key in a sorted array of
positive `Float64` partition points — `searchsortedlast([0.0, …], -0.0)`
returns `0`, not `1`. The package handles this internally with an
`nz(x) = iszero(x) ? zero(x) : x` guard wherever it calls `searchsortedlast`
on partition points; if you write similar code, do the same.

## Comparisons (the big change)

In IA 1.0, comparisons that involve interval-as-set semantics throw
`InconclusiveBooleanOperation` rather than returning `false`. The replacement
is always a named `*_interval` predicate that makes the intent explicit:

| Old | New |
|---|---|
| `a == b`, `a != b` | `isequal_interval(a, b)`, `!isequal_interval(a, b)` |
| `x ∈ I`, `x ∉ I`, `in(x, I)` | `in_interval(x, I)`, `!in_interval(x, I)` |
| `a ∩ b`, `intersect(a, b)` | `intersect_interval(a, b)` |
| `a ∪ b`, `union(a, b)` | `hull(a, b)` (convex hull) |
| `isempty(I)` | `isempty_interval(I)` |
| `a ⊂ b`, `issubset(a, b)` | `issubset_interval(a, b)` |
| `isfinite(I)` | `isbounded(I)` |

Mixed comparisons between an interval and a real (`x <= 0.5` with `x` an
`Interval{Float64}`) also throw — they fall back to `==` internally. If you
need to test "is the interval entirely below 0.5", extract the bound: write
`sup(x) <= 0.5` (or `inf(x) >= 0.5` for the converse). For functions whose
argument may be either a real or an interval (piecewise definitions, Newton
iterations), the package uses a small dispatch helper:

```julia
_ub(z::Real)     = z
_ub(z::Interval) = sup(z)
```

so the chained comparison `_ub(x) <= 0.5 ? branch_a(x) : branch_b(x)` resolves
unambiguously for both kinds of argument. See `src/InducedLSV.jl` for an
example.

### Why no shim?

You might be tempted to write

```julia
Base.intersect(a::Interval, b::Interval) = IntervalArithmetic.intersect_interval(a, b)
```

and recover the old syntax in one line. **Don't do this.** Upstream removed
the operator deliberately, and a global shim hides the choice from readers,
collides with future upstream definitions, and tends to outlive the migration.
Hand-rewrite each site to the named function — it's a one-time cost and keeps
the call sites self-documenting.

## Empty intervals and the `∅` symbol

The `∅` symbol moved to `IntervalArithmetic.Symbols` in IA 1.0 and is no
longer auto-exported. Use `emptyinterval()` (or `emptyinterval(typeof(x))`
when the element type matters) to construct an empty interval, and
`isempty_interval(x)` to test for one. Comparisons such as `x == ∅`
correspondingly become `isempty_interval(x)`.

## Decorations and `BareInterval`

In IA 1.0 every `Interval{T}` carries a *decoration* — `_com` (common, the
default), `_dac` (defined-and-continuous), `_def` (defined), `_trv` (trivial),
`_ill` (ill-formed). Decorations propagate through arithmetic and through
operations like `intersect_interval`, which conservatively returns a `_trv`
result. For most uses this is invisible, but two things to watch:

- The macro `@round(T, ex_lo, ex_hi)` returns a **`BareInterval{T}`**, not a
  decorated `Interval{T}`. If you need a decorated interval (e.g. so the
  result can flow into the rest of the package), wrap with `interval(...)`:

  ```julia
  return interval(@round(T, expr_lo, expr_hi))
  ```

  See `src/pitrig.jl` for the full pattern across the `sinpi`/`cospi`
  rewrites.

- Decorated equality is what `isequal_interval` checks; comparing a `_com`
  result to a `_trv` result for "the same" set returns `true` only when the
  decorations also match. Use `isequal_interval(bareinterval(a), bareinterval(b))`
  if you want set-only equality without decoration semantics.

## Square root with directed rounding

`sqrt(x, RoundUp)` was removed in IA 1.0. `RigorousInvariantMeasures` uses
`FastRounding.sqrt_round(x, RoundUp)` internally; do the same in any new
rigorous numerics you add to the package.

## Test idioms

When asserting something *contains* a known value, use `in_interval`:

```julia
@test in_interval(0.5, distance_from_invariant(B, D, Q, w, norms))
```

When asserting two computed intervals are bitwise equal (e.g. unit tests for
basis-construction code), use `isequal_interval`. When asserting a vector of
intervals contains a vector of reals, broadcast:

```julia
@test all(in_interval.(expected_floats, computed_intervals))
```

The `≈` operator on intervals also reduces to `==` and will throw — replace
with explicit `in_interval` enclosure checks or with bound-by-bound `≈`:

```julia
inf(a) ≈ inf(b) && sup(a) ≈ sup(b)
```
