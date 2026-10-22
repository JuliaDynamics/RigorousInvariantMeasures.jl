"""
    struct UniformKernelUlam{BC,T,S} <: NoiseKernel

Ulam discretization of the uniform noise kernel with half-width `l` on a partition
`B::Ulam`. The type parameter `BC` specifies the boundary condition:

- `:periodic`   → periodic wrap-around
- `:reflecting` → mirror reflection (projection π)

`T` is the type of the bounds, `Float64` by default. The scratch buffers carry
it, so applying a kernel to a `Vector{Interval{T′}}` with `T′` wider than `T`
would narrow the midpoints on the way in, which the error term added afterwards
does not account for; build the kernel at the type you mean to use.

`S` is the summation scheme of the window sums, `:sliding` by default:

- `:sliding` → one running sum, updated by adding the entering entry and
  subtracting the leaving one, with Kahan compensation. The error of an entry is
  that of a compensated sum of every term processed so far, so it is bounded in
  terms of the whole vector and not of the window.
- `:block`   → the block scheme of van Herk (Pattern Recognit. Lett. 13 (1992)
  517–521) and Gil and Werman (IEEE TPAMI 15 (1993) 504–507), written there for
  running maxima and valid for any associative operation. The extended vector is
  cut into blocks of length `n = 2l+1`, sums are accumulated inside each block
  from the right and from the left, and every window, which meets at most two
  blocks, is the sum of one right-accumulated and one left-accumulated partial
  sum. No entry is subtracted, and each window sum is a floating-point sum of its
  own `n` terms, so the error of an entry is at most `γₙ` times the kernel applied
  to `|v|` at that entry. The cost is three additions per entry, as for `:sliding`.

# Fields
- `B::Ulam`                 : the Ulam partition
- `l::Int`                  : half-width of the averaging window
- `scratch_ext::Vector{T}`  : workspace of length `k+2l` for building extended vector
- `scratch_sum::Vector{T}`  : workspace of length `k` for storing window sums
- `scratch_suf::Vector{T}`  : workspace of length `k+2l` for the right-accumulated
  block sums of the `:block` scheme (empty for `:sliding`)
"""
struct UniformKernelUlam{BC,T<:Real,S} <: NoiseKernel
    B::Ulam
    l::Int
    scratch_ext::Vector{T}
    scratch_sum::Vector{T}
    scratch_suf::Vector{T}
end

"""
    UniformKernelUlam(::Val{BC}, [T,] B::Ulam, l::Int; summation = :sliding)

Constructor for a `UniformKernelUlam{BC,T,S}` kernel on the partition `B` with
half-width `l`, bounds of type `T` (`Float64` by default) and summation scheme
`S = summation`, either `:sliding` or `:block`.
Allocates the necessary scratch buffers of sizes `k+2l` and `k` (`k = length(B)`).
"""
UniformKernelUlam(bc::Val, B::Ulam, l::Int; summation::Symbol = :sliding) =
    UniformKernelUlam(bc, Float64, B, l; summation)

function UniformKernelUlam(
    ::Val{BC},
    ::Type{T},
    B::Ulam,
    l::Int;
    summation::Symbol = :sliding,
) where {BC,T<:Real}
    summation in (:sliding, :block) ||
        throw(ArgumentError("summation must be :sliding or :block, got :$summation"))
    k = length(B)
    suf = summation === :block ? zeros(T, k + 2l) : T[]
    UniformKernelUlam{BC,T,summation}(B, l, zeros(T, k + 2l), zeros(T, k), suf)
end

"""
    opnormbound(B::Ulam, ::Type{L1}, M::UniformKernelUlam)

Return an upper bound on the operator norm in L¹.  
For the uniform noise kernel, the L¹ operator norm is exactly 1.
"""
opnormbound(B::Ulam, ::Type{L1}, M::UniformKernelUlam{BC}) where {BC} = 1.0

"""
    opradius(::Type{L1}, M::UniformKernelUlam)

Return the "radius" term in the L¹ Lasota–Yorke inequality.  
For the uniform kernel this is zero, but can be extended for interval-aware variants.
"""
opradius(::Type{L1}, M::UniformKernelUlam{BC}) where {BC} = 0.0

"""
    nonzero_per_row(M::UniformKernelUlam)

Return the number of nonzeros per row of the transition matrix associated
with the discretized uniform noise kernel. This equals the window size `2l+1`.
"""
nonzero_per_row(M::UniformKernelUlam{BC}) where {BC} = 2*M.l + 1

"""
    dfly(::Type{TotalVariation}, ::Type{L1}, N::UniformKernelUlam)

Return the Doeblin–Fortet–Lasota–Yorke inequality coefficients for the operator
acting from L¹ to bounded variation (Total Variation).  
For the uniform kernel on `k` bins with half-width `l`, the effective noise size is
ξ = (2l+1)/k, and the inequality is bounded by (0, 1/(2ξ)).
"""
dfly(::Type{TotalVariation}, ::Type{L1}, N::UniformKernelUlam{BC}) where {BC} = begin
    k = length(N.B)
    # `l` is the half-width in cells, so the support is 2l+1 cells wide and the
    # half-width is ξ = (2l+1)/(2k); this read (2l+1)/k, the full width, which
    # halved the constant a second time. With Var(ρ_ξ) = 1/ξ (Lemma 47 of
    # Galatolo-Monge-Nisoli) the bound is 2k/(2l+1); it used to be a quarter of
    # that, and so was not an upper bound.
    ξ = (2N.l + 1) / (2k)
    (0.0, 1 / ξ)
end

"""
    UniformKernelUlamPeriodic(B::Ulam, l::Int)

Construct a **periodic uniform Ulam kernel** on the partition `B` with
half-width `l` (window size = 2l+1).  

This operator acts as a stochastic convolution with uniform noise, where
indices outside `[1,k]` wrap around periodically.
"""
UniformKernelUlamPeriodic(B::Ulam, l::Int; summation::Symbol = :sliding) =
    UniformKernelUlam(Val(:periodic), B, l; summation)
UniformKernelUlamPeriodic(::Type{T}, B::Ulam, l::Int; summation::Symbol = :sliding) where {T<:Real} =
    UniformKernelUlam(Val(:periodic), T, B, l; summation)

"""
    UniformKernelUlamReflecting(B::Ulam, l::Int)

Construct a **reflecting uniform Ulam kernel** on the partition `B` with
half-width `l` (window size = 2l+1).  

This operator acts as a stochastic convolution with uniform noise, where
indices outside `[1,k]` are mapped back into `[1,k]` by the reflecting
projection π (period-2 mirror).
"""
UniformKernelUlamReflecting(B::Ulam, l::Int; summation::Symbol = :sliding) =
    UniformKernelUlam(Val(:reflecting), B, l; summation)
UniformKernelUlamReflecting(::Type{T}, B::Ulam, l::Int; summation::Symbol = :sliding) where {T<:Real} =
    UniformKernelUlam(Val(:reflecting), T, B, l; summation)

"""
    *(K::UniformKernelUlam, v::AbstractVector)

Non-mutating application of the kernel to vector `v`.  
This makes a copy of `v` internally and dispatches to the appropriate `mul!` method.
"""
function Base.:*(K::UniformKernelUlam{BC}, v::AbstractVector) where {BC}
    mul!(K, copy(v))
end

"""
    wrap_idx(i, k)

Periodic index mapping: wraps `i` into the range `1:k`.
"""
@inline wrap_idx(i::Int, k::Int) = (mod(i - 1, k) + 1)

"""
    reflect_outward_idx(i, k)

Reflecting index mapping (projection π):  
maps an index `i` on ℤ into `1:k` by period-2k reflection symmetry.
"""
@inline function reflect_outward_idx(i::Int, k::Int)
    r = mod(i - 1, 2k) + 1
    return r <= k ? r : (2k - r + 1)
end

"""
    get_idx(::Val{:periodic}, i, k)

Return the periodic index corresponding to `i` in `1:k`.
"""
@inline function get_idx(::Val{:periodic}, i::Int, k::Int)
    return wrap_idx(i, k)
end

"""
    get_idx(::Val{:reflecting}, i, k)

Return the reflecting index corresponding to `i` in `1:k`,
according to the projection π.
"""
@inline function get_idx(::Val{:reflecting}, i::Int, k::Int)
    return reflect_outward_idx(i, k)
end

"""
    mul!(K::UniformKernelUlam{BC,T}, v::Vector{T})

In-place application of the uniform kernel to a real vector `v`.  
Uses a preallocated scratch extension vector and sliding-window sum
with Kahan summation for numerical stability.
"""
function mul!(K::UniformKernelUlam{BC,T,:sliding}, v::Vector{T}) where {BC,T}
    k = length(v)
    l = K.l
    n = 2l + 1
    sums = K.scratch_sum
    v_ext = K.scratch_ext

    # build extension with chosen boundary condition
    @inbounds for j = 1:(k+2l)
        idx = get_idx(Val(BC), j - l, k)
        v_ext[j] = v[idx]
    end

    # initial sum
    s = sum(@view v_ext[1:n])
    c = zero(T)
    sums[1] = s

    # sliding window with Kahan
    @inbounds for i = 2:k
        δ = v_ext[i+n-1] - v_ext[i-1]
        y = δ - c
        t = s + y
        c = (t - s) - y
        s = t
        sums[i] = s
    end

    # normalize into v
    @inbounds for i = 1:k
        v[i] = sums[i] / n
    end

    return v
end

"""
    mul!(K::UniformKernelUlam, v::Vector{Interval})

In-place application of the uniform kernel to a vector of intervals.  
The operation is performed on midpoints with sliding window sums,
and then a uniform interval error bound is added to account for radii.
"""
function mul!(K::UniformKernelUlam{BC,T,:sliding}, v::Vector{Interval{T}}) where {BC,T}
    k = length(v)
    l = K.l
    n = 2l + 1
    sums = K.scratch_sum
    v_ext = K.scratch_ext

    midv = mid.(v)
    radv = radius.(v)

    # norms
    nrmv = sum(abs, midv)
    nrmrad = sum(abs, radv)

    # build extended midpoints
    @inbounds for j = 1:(k+2l)
        idx = get_idx(Val(BC), j - l, k)
        v_ext[j] = midv[idx]
    end

    # initial sum
    s = sum(@view v_ext[1:n])
    c = zero(T)
    sums[1] = s

    # sliding window with Kahan
    @inbounds for i = 2:k
        δ = v_ext[i+n-1] - v_ext[i-1]
        y = δ - c
        t = s + y
        c = (t - s) - y
        s = t
        sums[i] = s
    end

    # Error bound. The window sum is accumulated in floating point and then
    # divided by n, so the error on one entry is bounded by
    #
    #     γ_{n+1} · Σ_j |mid(v_j)| / n  +  Σ_j radius(v_j) / n ,
    #
    # the n additions of the sum together with the final division on the left,
    # the propagated input radii on the right. Kahan summation makes the true
    # error much smaller than γ_{n+1} allows, so this is an over-estimate.
    #
    # γₖ was hard-coded to 1.0 here, which made ϵ equal to ‖v‖₁/n and left the
    # enclosure vacuous; at k = 1024, ξ = 0.05 that gave a radius of 9.75 on
    # entries of order one, against 2.4e-11 for the kernel in NoiseKernel.jl.
    # One unit beyond the n additions covers the final division by n.
    δₖ = zero(T)  # the weights are the exact rationals 1/n; no matrix radius
    γₖ = gamma(T, n + 1)
    normMK = one(T)  # ‖N‖_{L¹→L¹} = 1 exactly, the kernel being Markov
    nT = T(n, RoundDown)
    ϵ = ((γₖ ⊗₊ normMK) ⊗₊ nrmv) ⊘₊ nT ⊕₊ ((normMK ⊗₊ nrmrad) ⊘₊ nT)

    # normalize into intervals
    @inbounds for i = 1:k
        v[i] = interval(T, sums[i] / n) + interval(T, -ϵ, ϵ)
    end

    return v
end

# ---------------------------------------------------------------------------
# The block scheme (van Herk; Gil and Werman)
# ---------------------------------------------------------------------------

"""
    _block_window_sums!(sums, ext, suf, k, n)

Write into `sums[i]`, for `i = 1:k`, the floating-point sum of `ext[i:i+n-1]`,
computed without subtraction. The extended vector `ext` (length `k+n-1`) is cut
into blocks `1:n`, `n+1:2n`, …; `suf[t]` is the sum from `t` to the end of its
block, accumulated from the right. A window starting at a block start is that
block, `suf[i]`; any other window starting in block `b` ends in block `b+1`, and
is `suf[i]` plus the sum from the start of block `b+1` to the end of the window,
accumulated from the left in `p`. Each window sum therefore uses `n-1` additions
of its own `n` terms and nothing else.
"""
function _block_window_sums!(sums::AbstractVector{T}, ext::AbstractVector{T},
                             suf::AbstractVector{T}, k::Int, n::Int) where {T}
    N = k + n - 1
    @inbounds for bstart = 1:n:N
        bend = min(bstart + n - 1, N)
        s = ext[bend]
        suf[bend] = s
        for t = (bend-1):-1:bstart
            s = ext[t] + s
            suf[t] = s
        end
    end
    p = zero(T)
    @inbounds for i = 1:k
        e = i + n - 1
        # `p` is reset at every block start, and it is read only for windows that
        # do not start at a block start, whose end `e` lies in the next block and
        # was reached from that block's start; the value it holds when `i` is a
        # block start is never read.
        p = (e - 1) % n == 0 ? ext[e] : p + ext[e]
        sums[i] = (i - 1) % n == 0 ? suf[i] : suf[i] + p
    end
    return sums
end

function _fill_extension!(ext::AbstractVector{T}, v::AbstractVector, ::Val{BC},
                          k::Int, l::Int) where {T,BC}
    @inbounds for j = 1:(k+2l)
        ext[j] = v[get_idx(Val(BC), j - l, k)]
    end
    return ext
end

"""
    mul!(K::UniformKernelUlam{BC,T,:block}, v::Vector{T})

In-place application of the kernel with the block scheme. Every entry of the
result is a window sum of `n = 2l+1` terms, formed with `n-1` additions and no
subtraction, divided by `n`; by the bound for summation in any order
[Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., (4.3)],
together with the final division,

    |fl((Kv)_i) - (Kv)_i| ≤ γₙ (K|v|)_i ,

so that the error is local and, in any monotone norm, at most `γₙ ‖K‖ ‖v‖`.
"""
function mul!(K::UniformKernelUlam{BC,T,:block}, v::Vector{T}) where {BC,T}
    k = length(v)
    l = K.l
    n = 2l + 1
    _fill_extension!(K.scratch_ext, v, Val(BC), k, l)
    _block_window_sums!(K.scratch_sum, K.scratch_ext, K.scratch_suf, k, n)
    @inbounds for i = 1:k
        v[i] = K.scratch_sum[i] / n
    end
    return v
end

"""
    mul!(K::UniformKernelUlam{BC,T,:block}, v::Vector{Interval{T}})

In-place application of the kernel to a vector of intervals with the block
scheme. With `m` and `r` the midpoints and radii of `v`, the result at entry `i`
is centred at the computed `(Km)_i` with radius

    γₙ (K|m|)_i + (K r)_i ,

the first term the rounding of the centre, the second the propagated radii
(`K` is entrywise nonnegative). Both `(K|m|)_i` and `(Kr)_i` are sums of
nonnegative terms, which the block scheme computes with relative error at most
`γₙ`, so each is bounded above by its computed value divided by `1 - γₙ`; the
radius is formed with upward rounding from these. The radius is local to the
window, unlike the uniform one of the `:sliding` scheme.
"""
function mul!(K::UniformKernelUlam{BC,T,:block}, v::Vector{Interval{T}}) where {BC,T}
    k = length(v)
    l = K.l
    n = 2l + 1
    ext, sums, suf = K.scratch_ext, K.scratch_sum, K.scratch_suf

    midv = mid.(v)
    absm = abs.(midv)
    radv = radius.(v)

    _fill_extension!(ext, midv, Val(BC), k, l)
    _block_window_sums!(sums, ext, suf, k, n)
    centre = sums ./ n

    _fill_extension!(ext, absm, Val(BC), k, l)
    _block_window_sums!(sums, ext, suf, k, n)
    kabsm = sums ./ n

    _fill_extension!(ext, radv, Val(BC), k, l)
    _block_window_sums!(sums, ext, suf, k, n)
    krad = sums ./ n

    γₙ = gamma(T, n)
    den = one(T) ⊖₋ γₙ
    @inbounds for i = 1:k
        ϵ = ((γₙ ⊗₊ kabsm[i]) ⊕₊ krad[i]) ⊘₊ den
        v[i] = interval(T, centre[i]) + interval(T, -ϵ, ϵ)
    end
    return v
end

# ---------------------------------------------------------------------------
# The Ulam projection of the uniform kernel
# ---------------------------------------------------------------------------

"""
    UlamProjectedUniformKernel{BC,T} <: NoiseKernel

The operator `π N_ξ π` on the Ulam partition `B` of `k` intervals of width `δ = 1/k`:
`N_ξ` is the convolution with the uniform density `ρ_ξ = ξ⁻¹ 1_{[-ξ/2, ξ/2]}` followed by
the boundary condition `BC` (`:periodic` or `:reflecting`), and `π` the Ulam projection.
It is the operator `π_δ N_ξ π_δ` of Galatolo, Monge, Nisoli, Nonlinearity 33 (2020) 4237,
§2, for which the two-grid estimate of their Lemma 29 is written; in particular the same
`ξ` can be used on a partition and on any refinement of it.

The flat window [`UniformKernelUlam`](@ref) is a different operator: for `ξ = (2l+1)δ` its
distance from this one in `L¹` is `δ/(2ξ)`, and on a dyadic refinement no flat window
represents the same `ξ`, since `2l′+1 = 2ʲ(2l+1)` has no solution.

With `s = ξ/δ`, the entry of the matrix at offset `d = i - j`, before the boundary
condition, is

    w_d = (F(s/2 - d) - F(-s/2 - d)) / s,

`F` the primitive of the unit tent `Λ(t) = max(0, 1 - |t|)`, which is the density of the
difference of two independent uniform points of a cell. So `w_d = 1/s` for
`|d| ≤ m := ⌊s/2⌋ - 1`, the offsets `±(m+1)` and `±(m+2)` carry `c₁/s` and `c₂/s` with
`c₁, c₂ ∈ [0, 1]`, and all the others vanish; `(2m+1) + 2(c₁ + c₂) = s`. For `s = 2l+1`
the row is `(1/8, 7/8, 1, …, 1, 7/8, 1/8)/s`.

The product is the flat window of `2m+1` cells, by the block scheme of
[`UniformKernelUlam`](@ref), plus the two boundary corrections on each side, divided by
`s`. Every entry is then a floating-point sum of `2m+5` nonnegative multiples of the
input, each passing through at most `P = max(2m+3, 5)` roundings, so

    |fl((Kv)_i) - (Kv)_i| ≤ γ_P (K|v|)_i ,

[Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §3.1 and (4.3)]; the
error is local, and at most `γ_P ‖v‖₁` in `L¹`.

`ξ k` must be exactly representable in `T`, and so must `c₁` and `c₂`, which are then exact
dyadic numbers when `s` is; the constructor checks both, so pass `ξ` as a rational, e.g.
`(2l+1)//k`. The support must cover at least two cells, `s ≥ 2`.
"""
struct UlamProjectedUniformKernel{BC,T<:Real} <: NoiseKernel
    B::Ulam
    s::T
    m::Int
    c₁::T
    c₂::T
    scratch_ext::Vector{T}
    scratch_sum::Vector{T}
    scratch_suf::Vector{T}
end

"Primitive of the unit tent Λ(t) = max(0, 1 - |t|), vanishing at -∞."
_tent_primitive(x) =
    x <= -1 ? zero(x) : x <= 0 ? (1 + x)^2 / 2 : x <= 1 ? 1 - (1 - x)^2 / 2 : one(x)

function UlamProjectedUniformKernel(::Val{BC}, ::Type{T}, B::Ulam, ξ::Real) where {BC,T<:Real}
    BC in (:periodic, :reflecting) ||
        throw(ArgumentError("boundary condition must be :periodic or :reflecting, got :$BC"))
    k = length(B)
    sq = Rational{BigInt}(ξ) * k
    s = T(sq)
    Rational{BigInt}(s) == sq ||
        throw(ArgumentError("ξk = $sq is not representable in $T; pass ξ as a rational with ξk representable"))
    sq >= 2 || throw(ArgumentError("the support must cover at least two cells, got ξk = $sq"))
    m = floor(Int, sq / 2) - 1
    cq(d) = _tent_primitive(sq / 2 - d) - _tent_primitive(-sq / 2 - d)
    c₁q, c₂q = cq(m + 1), cq(m + 2)
    @assert (2m + 1) + 2 * (c₁q + c₂q) == sq
    c₁, c₂ = T(c₁q), T(c₂q)
    Rational{BigInt}(c₁) == c₁q && Rational{BigInt}(c₂) == c₂q ||
        throw(ArgumentError("the boundary weights $c₁q, $c₂q are not representable in $T"))
    D = m + 2
    return UlamProjectedUniformKernel{BC,T}(B, s, m, c₁, c₂,
        zeros(T, k + 2D), zeros(T, k), zeros(T, k + 2D))
end
UlamProjectedUniformKernel(bc::Val, B::Ulam, ξ::Real) = UlamProjectedUniformKernel(bc, Float64, B, ξ)

UlamProjectedUniformKernelPeriodic(B::Ulam, ξ::Real) = UlamProjectedUniformKernel(Val(:periodic), B, ξ)
UlamProjectedUniformKernelReflecting(B::Ulam, ξ::Real) = UlamProjectedUniformKernel(Val(:reflecting), B, ξ)

"`L¹` operator norm: the matrix is nonnegative with unit row and column sums."
opnormbound(B::Ulam, ::Type{L1}, M::UlamProjectedUniformKernel) = 1.0

"The weights are stored exactly, so the matrix has no radius."
opradius(::Type{L1}, M::UlamProjectedUniformKernel) = 0.0

"""
Number of terms of one entry, `2m+5`. Callers take `gamma(T, nonzero_per_row(K))` as the
relative rounding of an entry, and `γ_{2m+5}` is at least the `γ_P` of the type docstring.
"""
nonzero_per_row(M::UlamProjectedUniformKernel) = 2M.m + 5

"""
    dfly(::Type{TotalVariation}, ::Type{L1}, N::UlamProjectedUniformKernel)

`Var(π N_ξ π f) ≤ Var(N_ξ π f) ≤ ξ⁻¹ Var(ρ) ‖π f‖₁ ≤ (2/ξ) ‖f‖₁`, the projection not
increasing the variation and `Var(ρ) = 2` for the uniform mother kernel on `[-1/2, 1/2]`
(Galatolo, Monge, Nisoli, Nonlinearity 33 (2020) 4237, Lemma 47); `2/ξ = 2k/s`.
"""
dfly(::Type{TotalVariation}, ::Type{L1}, N::UlamProjectedUniformKernel) =
    (0.0, 2 * length(N.B) / N.s)

_rounding_depth(K::UlamProjectedUniformKernel) = max(2K.m + 3, 5)

"Write `fl(K x)` into `out`; `out` may be `x`."
function _projected_window!(out::AbstractVector{T}, K::UlamProjectedUniformKernel{BC,T},
                            x::AbstractVector{T}) where {BC,T}
    k = length(x)
    m = K.m
    D = m + 2
    ext, sums = K.scratch_ext, K.scratch_sum
    _fill_extension!(ext, x, Val(BC), k, D)
    # the flat part: the window of 2m+1 cells centred at ext[i+D] starts at ext[i+2]
    _block_window_sums!(sums, view(ext, 3:length(ext)), K.scratch_suf, k, 2m + 1)
    @inbounds for i = 1:k
        r = (sums[i] + K.c₁ * (ext[i+1] + ext[i+2D-1])) + K.c₂ * (ext[i] + ext[i+2D])
        out[i] = r / K.s
    end
    return out
end

mul!(K::UlamProjectedUniformKernel{BC,T}, v::Vector{T}) where {BC,T} = _projected_window!(v, K, v)

"""
    mul!(K::UlamProjectedUniformKernel, v::Vector{Interval})

With `m` and `r` the midpoints and radii of `v`, the result at entry `i` is centred at the
computed `(Km)_i` with radius `γ_P (K|m|)_i + (Kr)_i`; the two sums are of nonnegative
terms, computed with relative error at most `γ_P`, and are divided by `1 - γ_P` and
combined with upward rounding, as for the `:block` flat window.
"""
function mul!(K::UlamProjectedUniformKernel{BC,T}, v::Vector{Interval{T}}) where {BC,T}
    midv = mid.(v)
    centre = _projected_window!(similar(midv), K, midv)
    kabsm = _projected_window!(similar(midv), K, abs.(midv))
    krad = _projected_window!(similar(midv), K, radius.(v))
    γ = gamma(T, _rounding_depth(K))
    den = one(T) ⊖₋ γ
    @inbounds for i in eachindex(v)
        ϵ = ((γ ⊗₊ kabsm[i]) ⊕₊ krad[i]) ⊘₊ den
        v[i] = interval(T, centre[i]) + interval(T, -ϵ, ϵ)
    end
    return v
end

Base.:*(K::UlamProjectedUniformKernel, v::AbstractVector) = mul!(K, copy(v))
