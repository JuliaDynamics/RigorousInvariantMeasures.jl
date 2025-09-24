struct UniformKernelUlam{BC} <: NoiseKernel
    B::Ulam
    l::Int
    scratch_ext::Vector{Float64}
    scratch_sum::Vector{Float64}
end

function UniformKernelUlam(::Val{BC}, B::Ulam, l::Int) where {BC}
    k = length(B)
    UniformKernelUlam{BC}(B, l, zeros(k + 2l), zeros(k))
end

UniformKernelUlamPeriodic(B::Ulam, l::Int) = UniformKernelUlam(Val(:periodic), B, l)
UniformKernelUlamReflecting(B::Ulam, l::Int) = UniformKernelUlam(Val(:reflecting), B, l)

function Base.:*(K::UniformKernelUlam{BC}, v::AbstractVector) where {BC}
    mul!(K, copy(v))  # non-mutating, returns new vector
end

@inline wrap_idx(i::Int, k::Int) = (mod(i - 1, k) + 1)

@inline function reflect_outward_idx(i::Int, k::Int)
    r = mod(i - 1, 2k) + 1      # fold into 1..2k
    return r <= k ? r : (2k - r + 1)
end

@inline function get_idx(::Val{:periodic}, i::Int, k::Int)
    return wrap_idx(i, k)
end

@inline function get_idx(::Val{:reflecting}, i::Int, k::Int)
    return reflect_outward_idx(i, k)
end

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
