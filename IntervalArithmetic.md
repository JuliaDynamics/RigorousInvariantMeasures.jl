# IntervalArithmetic 1.0 — known issues and workarounds

This file tracks open issues with the IntervalArithmetic 1.0 API as they affect this package.
Update it when upstream releases a fix or when the workaround changes.

---

## Scalar multiplication marks results as `_NG` (not guaranteed)

**Status:** open — to be raised with IA maintainers.

### What happens

In IA 1.0, multiplying an interval by a plain Julia scalar (`Int` or `Float64`) marks the
result with the `NG` ("not guaranteed") flag, even when the multiplication is exact:

```julia
julia> using IntervalArithmetic
julia> 2 * interval(0.5)
[1.0, 1.0]_com_NG   # NG despite 2 being an exact integer
julia> interval(2) * interval(0.5)
[1.0, 1.0]_com      # _com when both operands are intervals
julia> interval(0.5) + interval(0.5)
[1.0, 1.0]_com      # _com via addition
```

The `NG` flag propagates: `_com_NG` through `hull()` degrades to `_trv_NG`, and so on.

### Why it matters here

User-supplied branch functions like `x -> 2x` or `x -> 3x - 1` are evaluated on intervals
via plain Julia arithmetic.  The branch endpoint images `f(interval(a))`, `f(interval(b))` end
up `_com_NG`, and this propagates into every preimage interval computed by `Contractors.jl`.

The numerical enclosures produced by the Newton/Krawczyk contractor in `Contractors.jl` **are**
rigorous — the NG flag is a decoration conservatism, not a correctness failure.

### Workaround

We call

```julia
setdisplay(:infsup; decorations = false, ng_flag = false)
```

in the module initialisation (`src/RigorousInvariantMeasures.jl`) to suppress decoration
suffixes from all interval display.  This keeps doctests and REPL output readable until the
upstream issue is resolved.

The call is placed right after the `IntervalArithmetic` import.  Remove it (and update the
doctests) once IA fixes scalar-multiplication decoration.

### Reference

- Upstream issue: (link to be added after filing)
- The same `setdisplay` knob is documented in IA's own `display.jl` under `DisplayOptions`.
