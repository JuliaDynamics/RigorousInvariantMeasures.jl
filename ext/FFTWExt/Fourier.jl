# Assemble methods for the Fourier basis family. Depend on `interval_fft`
# from IntervalFFT.jl.

# Lift a bound or an existing interval into `Interval{T}` without relying on
# implicit Real → Interval conversion, which IntervalArithmetic rejects.
_to_interval(::Type{T}, v::Interval) where {T} = interval(T, inf(v), sup(v))
_to_interval(::Type{T}, v::Real) where {T} = interval(T, v)

# Working precision of a basis: the number type behind its grid points.
_basis_numtype(B::Fourier) = _numtype_of(eltype(B.p))
_numtype_of(::Type{Interval{T}}) where {T} = T
_numtype_of(::Type) = Float64   # abstract/heterogeneous grids fall back

@doc raw"""
    assemble_common(B::Fourier, D; ϵ, max_iter, T)

Assemble the discretized transfer operator on a Fourier basis.

Column `i` is the FFT of the dual evaluation of the basis function
``ϕ_{m_i}(x) = e^{2πi m_i x}``. Written naively that is one rigorous complex
exponential per (basis function, dual node) pair — `n × nd` of them, and a
`Complex{Interval}` `exp` costs ~3 µs, which put the `k = 1024, npts = 16384`
assembly at the multi-day scale.

Two structural facts remove almost all of it:

1. ``ϕ_m(x) = z^m`` with ``z = e^{2πix}``, so only `nd` exponentials are
   needed; the powers follow by binary exponentiation. The `z^(2^p)` table is
   built with direct `exp` calls rather than by repeated squaring — `log2(k)`
   extra exponentials per node is nothing next to the sweep, and it keeps each
   factor at ~1 ulp instead of compounding the squaring error through the
   table.
2. The dual weights are real (see [`dual_nodes`](@ref)), so the dual vector at
   frequency `-m` is the elementwise conjugate of the one at `+m`. Only
   `m = 0..k` is swept.

Everything stays in interval arithmetic, so the result is an enclosure exactly
as before — measurably wider per entry (~2× at small `k`, reaching ~1e-11 at
`k = 1024`) because the powering chain replaces a single `exp`.

`T` is the working precision, taken from the basis grid by default. Note that
[`interval_fft`](@ref) currently only implements `Float64`.
"""
function RigorousInvariantMeasures.assemble_common(
    B::Fourier,
    D;
    ϵ = 0.0,
    max_iter = 100,
    T = _basis_numtype(B),
)
    computed_dual = RigorousInvariantMeasures.Dual(B, D; ϵ, max_iter)
    x, labels, weights = RigorousInvariantMeasures.dual_nodes(B, computed_dual)
    # The duals store `Vector{Interval}` (abstract eltype); narrow here so the
    # kernel below is type-stable.
    return _assemble_fourier(
        T,
        Interval{T}[_to_interval(T, xj) for xj in x],
        collect(Int, labels),
        Interval{T}[_to_interval(T, w) for w in weights],
        length(B),
        max(length(B.p), Int(maximum(labels))),
    )
end

# Function barrier: `T` arrives as a runtime value from the keyword argument,
# so the whole sweep would be dynamically dispatched if it were inlined into
# `assemble_common`. Splitting it out costs one dispatch and buys a ~2×
# speedup on the inner loops.
function _assemble_fourier(
    ::Type{T},
    x::Vector{Interval{T}},
    labels::Vector{Int},
    weights::Vector{Interval{T}},
    n::Int,
    nw::Int,
) where {T}
    k = (n - 1) ÷ 2
    CT = Complex{Interval{T}}
    nd = length(x)

    twopi_i = 2 * interval(T, π) * im

    # z_j^(2^p) for p = 0 … ⌊log2 k⌋
    nbits = k == 0 ? 1 : (floor(Int, log2(k)) + 1)
    pow2 = Matrix{CT}(undef, nd, nbits)
    @inbounds for p = 1:nbits
        s = interval(T, 2)^(p - 1)
        for j = 1:nd
            pow2[j, p] = exp(twopi_i * s * x[j])
        end
    end

    c0 = CT[w + 0im for w in weights]

    M = zeros(CT, (n, n))
    w = Vector{CT}(undef, nw)
    acc = Vector{CT}(undef, nd)
    zero_CT = zero(CT)

    for m = 0:k
        # acc_j = weight_j · z_j^m
        copyto!(acc, c0)
        mm, p = m, 1
        while mm > 0
            if isodd(mm)
                @inbounds for j = 1:nd
                    acc[j] = acc[j] * pow2[j, p]
                end
            end
            mm >>= 1
            p += 1
        end

        fill!(w, zero_CT)
        @inbounds for j = 1:nd
            v = acc[j]
            if !isempty_interval(real(v)) && !isempty_interval(imag(v))
                w[labels[j]] += v
            end
        end

        _store_column!(M, interval_fft(w), m + 1, k, n)
        if m > 0                      # frequency -m: weights are real, so conjugate
            _store_column!(M, interval_fft(conj.(w)), n - m + 1, k, n)
        end
    end
    return M
end

# Keep the [0:k; -k:-1] layout of the basis: head and tail of the transform.
function _store_column!(M, F, col, k, n)
    @inbounds @views begin
        M[1:k+1, col] .= F[1:k+1]
        M[k+2:n, col] .= F[end-k+1:end]
    end
    return M
end

function RigorousInvariantMeasures.assemble(
    B::FourierAdjoint,
    D;
    ϵ = 0.0,
    max_iter = 100,
    T = _basis_numtype(B),
)
    # `'` produces a lazy `LinearAlgebra.Adjoint`; materialize to a plain
    # matrix so downstream consumers (e.g. `BallMatrix(Q.L)` in
    # `norms_of_powers`) can wrap it.
    return Matrix(
        RigorousInvariantMeasures.assemble_common(
            B,
            D;
            ϵ,
            max_iter,
            T,
        )',
    )
end

function RigorousInvariantMeasures.assemble(
    B::FourierAnalytic,
    D::Dynamic;
    ϵ = 0.0,
    max_iter = 100,
    T = _basis_numtype(B),
)
    return RigorousInvariantMeasures.assemble_common(B, D; ϵ, max_iter, T)
end
