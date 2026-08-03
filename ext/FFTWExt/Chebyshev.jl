# Assemble for the Chebyshev basis. Uses `interval_fft` via the `chebtransform`
# helper (DCT-via-FFT round-trip).

@doc raw"""
    chebtransform(w)

Chebyshev coefficients of the function taking the values `w` at the Chebyshev
points, computed as a DCT via the FFT of the mirrored sequence.

For values at the `N = n+1` Chebyshev points the coefficients are
``a_k = \frac1n \Re (F_k)`` with `F` the **unnormalized** DFT of
`[reverse(w); w[2:end-1]]` (length `2n`), and the two endpoint coefficients
halved.

`interval_fft` already divides by its own length `2n`, so the `/n` that used to
sit here made the total normalization `2n²` instead of `n` — every Chebyshev
assembly came out scaled by `1/(2n) = 1/(2(N-1))`. Multiplying by 2 undoes
`interval_fft`'s normalization down to the `1/n` the transform actually wants.
"""
function chebtransform(w)
    # The mirrored sequence has length 2(N-1) = 2n. The certified FFT error
    # estimate is only valid for power-of-two lengths — otherwise `interval_fft`
    # merely WARNS ("The rigorous error estimate works for power of two sizes")
    # once per column, which is easy to lose in the noise of a long run and
    # leaves the operator enclosure unjustified. A run at n = 200 (length 400)
    # produced a "certified" bound of 1e-134 on that basis. Fail instead.
    m = 2 * (length(w) - 1)
    ispow2(m) || throw(ArgumentError(
        "Chebyshev transform length $m is not a power of two (basis has " *
        "$(length(w)) points, degree $(length(w)-1)); the rigorous FFT error " *
        "estimate does not apply. Use a basis whose degree is a power of two."))
    z = 2 * interval_fft([reverse(w); w[2:end-1]])
    t = real.(z[1:length(w)])
    t[1] /= 2
    t[end] /= 2
    return Interval.(t)
end

@doc raw"""
    assemble(B::Chebyshev, D::Dynamic; ϵ, max_iter, T)

Assemble the discretized transfer operator on a Chebyshev basis.

Column `i` is the Chebyshev transform of the dual evaluation of basis function
`i`, i.e. of ``T_{i-1}`` at the dual nodes. Evaluating each basis function
separately means one Clenshaw pass — `O(n)` work — per (basis function, node)
pair, so `O(n²·nd)` overall; that is cubic in the basis size, and it dominated
the assembly (91 % of runtime at `n = 512`).

Instead the degrees are swept in one pass, via ``T_m(\cos θ) = \cos(mθ)``:
with ``z = e^{iθ}`` on the unit circle, ``T_m(t) = \Re(z^m)``, so all degrees
follow from binary exponentiation of one complex number per node. That is
`O(log n)` per (degree, node) pair and brings the assembly down to
`O(n·nd·log n)`.

Note that the obvious alternative — running the three-term recurrence
``T_{m+1} = 2t\,T_m - T_{m-1}`` directly in interval arithmetic — is
**unusable** here. Interval arithmetic cannot see the cancellation between the
two terms, so the radii obey ``r_{m+1} = 2|t|\,r_m + r_{m-1}``: Fibonacci-like
growth, exponential in the degree even for ``|t| \le 1``. This is exactly
Example 1 of Ledoux–Moroz; measured here, it reaches a radius of `4e77` at
`n = 256`. Powers of a unit-modulus number have no such wrapping, and the
`z^(2^p)` table is built by direct `exp` rather than by repeated squaring, so
each factor stays at ~1 ulp and the accumulated radius grows like
``\log_2 n`` ulps.

Rigour here comes from the interval operations themselves — `acos`, `cos`,
`sin` and complex multiplication are all outward rounded — not from a separate
error analysis. In particular this path does **not** use
`eval_Clenshaw_BackwardFirst`, which implements Algorithm 2 (backward error
analysis, Lemma 3) of Ledoux–Moroz: that routine evaluates one polynomial in
`O(n)`, so reusing it would keep the assembly cubic. Its analysis targets
evaluation on intervals of nonzero radius `r`; the nodes here reproduce
`evalChebyshev`, which evaluates at the *midpoint* of `2x - 1`, so `r = 0` and
the two approaches are directly comparable — measured, this one is tighter
(`1.9e-16` against `4.2e-16` at `n = 512`).

# References

* V. Ledoux, G. Moroz, *Evaluation of Chebyshev Polynomials on Intervals and
  Application to Root Finding*, MACIS 2019, LNCS 11989, Springer, 2020.
  [doi:10.1007/978-3-030-43120-4_4](https://doi.org/10.1007/978-3-030-43120-4_4),
  [arXiv:1912.05843](https://arxiv.org/abs/1912.05843).
"""
function RigorousInvariantMeasures.assemble(
    B::Chebyshev,
    D::Dynamic;
    ϵ = 1e-13,
    max_iter = 100,
    T = Float64,
)
    n = length(B.p)
    x, labels, x′ = RigorousInvariantMeasures.Dual(B, D; ϵ, max_iter)
    return _assemble_chebyshev(
        T,
        # evalChebyshev evaluates at mid(2x - 1); reproduce that exactly.
        Interval{T}[interval(T, mid(2 * _to_interval(T, xj) - 1)) for xj in x],
        collect(Int, labels),
        Interval{T}[1 / abs(_to_interval(T, d)) for d in x′],
        n,
    )
end

# Function barrier, as for the Fourier sweep: `T` is a runtime keyword value,
# so the recurrence has to sit behind a type-parameterized call to stay
# type-stable.
function _assemble_chebyshev(
    ::Type{T},
    t::Vector{Interval{T}},
    labels::Vector{Int},
    weights::Vector{Interval{T}},
    n::Int,
) where {T}
    nd = length(t)
    CT = Complex{Interval{T}}
    unit = interval(T, -1, 1)

    # θ_j = acos(t_j). The nodes are preimages of grid points in [0,1] pushed
    # through 2x-1, so they belong to [-1,1]. Insist on it rather than clamping:
    # T_m is defined outside [-1,1] too, and silently projecting a node back
    # onto the interval would return a wrong enclosure, not merely a wider one.
    θ = Vector{Interval{T}}(undef, nd)
    @inbounds for j = 1:nd
        issubset_interval(t[j], unit) ||
            error("Chebyshev node $(t[j]) is not contained in [-1, 1]")
        θ[j] = acos(t[j])
    end

    # z_j^(2^p) for p = 0 … ⌊log2(n-1)⌋, built by direct exp (see docstring).
    nbits = n <= 2 ? 1 : (floor(Int, log2(n - 1)) + 1)
    pow2 = Matrix{CT}(undef, nd, nbits)
    @inbounds for p = 1:nbits
        s = interval(T, 2)^(p - 1)
        for j = 1:nd
            sθ = s * θ[j]
            pow2[j, p] = complex(cos(sθ), sin(sθ))
        end
    end

    M = zeros(Interval{T}, (n, n))
    w = Vector{Interval{T}}(undef, n)
    acc = Vector{CT}(undef, nd)
    one_CT = complex(interval(T, 1), interval(T, 0))
    zero_I = interval(T, 0)

    for m = 0:n-1
        # acc_j = z_j^m, so T_m(t_j) = real(acc_j)
        fill!(acc, one_CT)
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

        fill!(w, zero_I)
        @inbounds for j = 1:nd
            w[labels[j]] += weights[j] * real(acc[j])
        end
        M[:, m+1] = chebtransform(w)
    end
    return M
end
