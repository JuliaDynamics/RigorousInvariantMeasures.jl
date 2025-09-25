"""
    struct UniformKernelUlam{BC} <: NoiseKernel

Ulam discretization of the uniform noise kernel with half-width `l` on a partition
`B::Ulam`. The type parameter `BC` specifies the boundary condition:

- `:periodic`   → periodic wrap-around
- `:reflecting` → mirror reflection (projection π)

# Fields
- `B::Ulam`              : the Ulam partition
- `l::Int`               : half-width of the averaging window
- `scratch_ext::Vector`  : workspace of length `k+2l` for building extended vector
- `scratch_sum::Vector`  : workspace of length `k` for storing window sums
"""
struct UniformKernelUlam{BC} <: NoiseKernel
    B::Ulam
    l::Int
    scratch_ext::Vector{Float64}
    scratch_sum::Vector{Float64}
end

"""
    UniformKernelUlam(::Val{BC}, B::Ulam, l::Int)

Constructor for a `UniformKernelUlam{BC}` kernel on the partition `B` with half-width `l`.
Allocates the necessary scratch buffers of sizes `k+2l` and `k` (`k = length(B)`).
"""
function UniformKernelUlam(::Val{BC}, B::Ulam, l::Int) where {BC}
    k = length(B)
    UniformKernelUlam{BC}(B, l, zeros(k + 2l), zeros(k))
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
    ξ = (2N.l + 1) / k
    (0.0, 1 / (2ξ))
end

"""
    UniformKernelUlamPeriodic(B::Ulam, l::Int)

Construct a **periodic uniform Ulam kernel** on the partition `B` with
half-width `l` (window size = 2l+1).  

This operator acts as a stochastic convolution with uniform noise, where
indices outside `[1,k]` wrap around periodically.
"""
UniformKernelUlamPeriodic(B::Ulam, l::Int) = UniformKernelUlam(Val(:periodic), B, l)

"""
    UniformKernelUlamReflecting(B::Ulam, l::Int)

Construct a **reflecting uniform Ulam kernel** on the partition `B` with
half-width `l` (window size = 2l+1).  

This operator acts as a stochastic convolution with uniform noise, where
indices outside `[1,k]` are mapped back into `[1,k]` by the reflecting
projection π (period-2 mirror).
"""
UniformKernelUlamReflecting(B::Ulam, l::Int) = UniformKernelUlam(Val(:reflecting), B, l)

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
    mul!(K::UniformKernelUlam, v::Vector{Float64})

In-place application of the uniform kernel to a real vector `v`.  
Uses a preallocated scratch extension vector and sliding-window sum
with Kahan summation for numerical stability.
"""
function mul!(K::UniformKernelUlam{BC}, v::Vector{Float64}) where {BC}
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
    c = 0.0
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
function mul!(K::UniformKernelUlam{BC}, v::Vector{Interval{T}}) where {BC,T}
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
    c = 0.0
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

    # crude but safe error bound
    δₖ = 0.0      # uniform kernel has exact opnorm = 1
    γₖ = 1.0
    normMK = 1.0
    ϵ = (γₖ ⊗₊ normMK) ⊗₊ nrmv / n ⊕₊ (normMK ⊗₊ nrmrad) / n

    # normalize into intervals
    @inbounds for i = 1:k
        v[i] = Interval(sums[i] / n) + Interval(-ϵ, ϵ)
    end

    return v
end
