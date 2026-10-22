"""
Functions to estimate Q|_{U^0}. See our paper for details.
"""

using LinearAlgebra
using SparseArrays

export norms_of_powers_noise, correlation_sequence_noise, noise_step, NoiseStepBound

"""
    NoiseStepBound(map_rounding, map_norm, noise_rounding, noise_norm)

The four constants of the a priori rounding bound of one step `v ← MK(Mv)`, with `M` the
midpoint of an integral-preserving `Q` (Galatolo, Nisoli, Numerical Remark 1 in §8.3 of
"An elementary approach to rigorous approximation of invariant measures", SIAM J. Appl.
Dyn. Syst. 13 (2014)): if `v` is within `ϵ` of the exact vector, the computed `M v` is within
`map_rounding ‖v‖ + map_norm ϵ` of the exact `Q v`, where `map_rounding = γ_z ‖M‖ + ‖rad Q‖`
and `map_norm ≥ ‖Q‖`; the noise step likewise, with `noise_rounding = γ_k ‖MK‖ + δ_k` and
`noise_norm ≥ ‖MK‖`. The norms are those `opnormbound` returns for the basis.
"""
struct NoiseStepBound{T}
    map_rounding::T
    map_norm::T
    noise_rounding::T
    noise_norm::T
end

function NoiseStepBound(B::Basis, N::Type{<:NormKind}, Q::DiscretizedOperator, MK::NoiseKernel)
    is_integral_preserving(Q) ||
        throw(ArgumentError("the step bound is written for an integral-preserving operator"))
    T = typeof(sup(zero(eltype(Q.L))))
    M = mid.(Q.L)
    δ = opnormbound(B, N, radius.(Q.L))
    γz = gamma(T, max_nonzeros_per_row(Q.L))
    nrmM = opnormbound(B, N, M)
    δₖ = opradius(N, MK)
    γₖ = gamma(T, nonzero_per_row(MK))
    nrm_MK = opnormbound(B, N, MK)
    return NoiseStepBound(γz ⊗₊ nrmM ⊕₊ δ, nrmM ⊕₊ δ, γₖ ⊗₊ nrm_MK ⊕₊ δₖ, nrm_MK ⊕₊ δₖ)
end

"""
    noise_step(B, N, M, MK, v, ϵ, nrmv, c::NoiseStepBound) -> (v, ϵ, nrmv)

One step `v ← MK(Mv)` in floating point, `nrmv ≥ ‖v‖` on entry, and the bound `ϵ` on the
distance from the exact iterate carried through it as in [`NoiseStepBound`](@ref); returns
the new vector, its error bound and an upper bound on its norm.
"""
function noise_step(B, N, M, MK, v, ϵ, nrmv, c::NoiseStepBound)
    v = M * v
    ϵ = c.map_rounding ⊗₊ nrmv ⊕₊ c.map_norm ⊗₊ ϵ
    nrmv = opnormbound(B, N, v)
    v = MK * v
    ϵ = c.noise_rounding ⊗₊ nrmv ⊕₊ c.noise_norm ⊗₊ ϵ
    return v, ϵ, opnormbound(B, N, v)
end

@doc raw"""
    correlation_sequence_noise(ϕ::ProjectedFunction{<:Ulam}, Q, MK, g, nlags; g_error = 0.0)

Enclosures of ``\int ϕ\, \tilde L^n g\, dm`` for ``n = 0, …, `nlags`-1``, where
``\tilde L = MK ∘ Q`` is the discretised annealed operator on the Ulam basis of `ϕ` and `Q`
is integral preserving. The iterates are computed in floating point on the midpoint of
`Q`, each within the bound of [`noise_step`](@ref) of the exact ``\tilde L^n g``, starting
from `g_error ≥ \|g_{\mathrm{exact}} - g\|_{L^1}`, and are paired with `ϕ` by
[`integral_pairing`](@ref), which charges ``\|ϕ\|_{L^∞}`` times that bound.

What is enclosed is the correlation of the *discretised* operator: the distance between
``\tilde L^n`` and the operator it discretises is not part of the result.
"""
function correlation_sequence_noise(ϕ::ProjectedFunction{<:Ulam}, Q::DiscretizedOperator,
                                    MK::NoiseKernel, g::AbstractVector, nlags::Integer;
                                    g_error::Real = 0.0)
    B = ϕ.B
    k = length(B)
    c = NoiseStepBound(B, L1, Q, MK)
    M = mid.(Q.L)
    v = collect(float.(g))
    # opnormbound counts Σ|vᵢ|, integral_pairing takes ∫|·| dm = Σ|vᵢ|/k
    ϵ = Float64(g_error, RoundUp) ⊗₊ Float64(k, RoundUp)
    nrmv = opnormbound(B, L1, v)
    C = Vector{typeof(interval(0.0))}(undef, nlags)
    for n = 1:nlags
        C[n] = integral_pairing(ϕ, v, ϵ ⊘₊ Float64(k, RoundDown))
        n < nlags && ((v, ϵ, nrmv) = noise_step(B, L1, M, MK, v, ϵ, nrmv, c))
    end
    return C
end

"""
    norms_of_powers_noise(B::Ulam, N, m, Q, MK, f; normv0, normQ)

Upper bounds on `‖(MK Q)^k|_{U⁰}‖`, `k = 1, …, m`, in the norm `N` as `opnormbound` counts
it, for an integral-preserving `Q` (every Ulam operator is). Each basis vector
`e₁ - e_{j+1}` of the zero-average space is iterated in floating point by
[`noise_step`](@ref), its rounding carried by [`NoiseStepBound`](@ref); the result is
the largest column norm plus its error, which bounds the norm on `U⁰` up to the factor
`2` lost in writing an element of `U⁰` in this basis. Callers should take the minimum
with the trivial bounds, as [`powernormboundsnoise`](@ref) does.

`f`, the integral covector, must be an interval vector with equal entries.
`normv0` and `normQ` may be passed to skip computing `‖e₁ - e_{j+1}‖` and `‖Q‖`.
"""
function norms_of_powers_noise(
    B::Ulam,
    N::Type{L1},
    m::Integer,
    Q::DiscretizedOperator,
    MK::NoiseKernel,
    f::AbstractArray;
    normv0::Real = -1.0, #used as "missing" value
    normQ::Real = -1.0,
)
    @assert eltype(f) <: Interval
    is_integral_preserving(Q) ||
        throw(ArgumentError("norms_of_powers_noise is written for an integral-preserving Q"))
    T = typeof(sup(zero(eltype(Q.L)))) # gets "Float64" from Q.L
    n = size(Q.L, 1)
    M = mid.(Q.L)

    c = NoiseStepBound(B, N, Q, MK)
    if normQ != -1.0
        c = NoiseStepBound(c.map_rounding, normQ, c.noise_rounding, c.noise_norm)
    end

    normcachers = [NormCacher{N}(B, n) for j = 1:m]
    v = zeros(T, n)
    for j = 1:n-1
        v .= zero(T)
        v[1] = one(T)
        v[j+1] = -one(T)
        nrmv = normv0 == -1.0 ? opnormbound(B, N, v) : normv0
        ϵ = zero(T)
        for k = 1:m
            v, ϵ, nrmv = noise_step(B, N, M, MK, v, ϵ, nrmv, c)
            add_column!(normcachers[k], v, ϵ)
        end
    end
    return map(get_norm, normcachers)
end

"""
    norms_of_powers_sequence_noise(B, N, Q, MK, f, shifts)

Like [`norms_of_powers_noise`](@ref), but for a *nonautonomous* (sequential)
cocycle in which a different forcing acts at each step.  The k-th step is

    v  ←  Shift(shifts[k]) ∘ MK ∘ Q  v ,

where `Shift(s)` is the (fractional) circular shift by `s` Ulam bins on the torus.
A fractional shift `Shift(s) = (1−f)·circshift(·,n₀) + f·circshift(·,n₀+1)`
(`n₀ = ⌊s⌋`, `f = s − n₀`) is a convex combination of two basis permutations:
mass-preserving, with L¹ operator norm 1.  An **integer** shift (`f = 0`) is an
exact permutation and introduces **no** floating-point error; a fractional shift
costs one extra rounded mul-add per entry, tracked in the accumulated radius `ϵ`
exactly like the noise step (operator norm 1).  Either way the rigorous error
propagation of `norms_of_powers_noise` carries over.

Using *fractional* shifts lets the forcing track the real kick
`s_n = A·k·sin(2π θ_n)` smoothly in the driving phase `θ_n = θ₀ + nω` (the
irrational rotation), with no rounding of the kick.

`shifts` has length `m` (one shift per step); the returned vector `norms[k]`
bounds `‖ L_{shifts[k]} ⋯ L_{shifts[1]} |_{U⁰} ‖`, the norm of the length-`k`
forced cocycle on the zero-average subspace.  With `shifts .= 0` it reproduces
`norms_of_powers_noise` exactly.
"""
function norms_of_powers_sequence_noise(
    B::Ulam,
    N::Type{L1},
    Q::DiscretizedOperator,
    MK::NoiseKernel,
    f::AbstractArray,
    shifts::AbstractVector{<:Real};
    normv0::Real = -1.0,
    normQ::Real = -1.0,
    normE::Real = -1.0,
    normEF::Real = -1.0,
    normIEF::Real = -1.0,
    normN::Real = -1.0,
    normρ::Real = -1.0,
)
    @assert eltype(f) <: Interval
    m = length(shifts)
    T = typeof(sup(zero(eltype(Q.L))))
    n = size(Q.L, 1)
    M = mid.(Q.L)

    R = radius.(Q.L)
    δ = opnormbound(B, N, R)

    γz = gamma(T, max_nonzeros_per_row(Q.L))
    γn = gamma(T, n + 3)

    ϵ = zero(T)
    nrmM = opnormbound(B, N, M)
    δₖ = opradius(N, MK)
    γₖ = gamma(T, nonzero_per_row(MK))
    nrm_MK = opnormbound(B, N, MK)
    normMK = nrm_MK ⊕₊ δₖ

    # a fractional shift touches 2 entries per output row (one rounded mul-add)
    γₛ = gamma(T, 2)

    if !is_integral_preserving(Q)
        if normE == -1.0
            normE = opnormbound(B, N, Q.e)
        end
        if normEF == -1.0
            normEF = opnormbound(B, N, Q.e * f)
        end
        if normIEF == -1.0
            normIEF = opnormbound(B, N, [Matrix(UniformScaling{Float64}(1), n, n) Q.e * f])
        end
        if normN == -1.0
            normN = opnormbound(B, N, Matrix(UniformScaling{Float64}(1), n, n) - Q.e * f)
        end
    end

    if normQ == -1.0
        if is_integral_preserving(Q)
            normQ = nrmM ⊕₊ δ
        else
            defect = opnormbound(B, N, Q.w)
            normQ = nrmM ⊕₊ δ ⊕₊ normE ⊗₊ defect
        end
    end

    normcachers = [NormCacher{N}(B, n) for j = 1:m]
    midf = map(mid, f)

    v = zeros(T, n)

    for j = 1:n-1
        v .= zero(T)
        v[1] = one(T)
        v[j+1] = -one(T)
        if normv0 == -1.0
            nrmv = opnormbound(B, N, v)
        else
            nrmv = normv0
        end
        ϵ = 0.0
        nrmw = nrmv
        for k = 1:m
            w = M * v
            if is_integral_preserving(Q)
                v = w
                ϵ = (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
            else
                v = w - Q.e * (midf * w)
                new_nrmw = opnormbound(B, N, w)
                ϵ =
                    γn ⊗₊ normIEF ⊗₊ (new_nrmw ⊕₊ normEF ⊗₊ nrmw) ⊕₊
                    normN ⊗₊ (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
                nrmw = new_nrmw
            end

            # the noise step
            nrmv = opnormbound(B, N, v)
            v = w
            w = MK * v
            v = w
            ϵ = (γₖ ⊗₊ nrm_MK ⊕₊ δₖ) ⊗₊ nrmv ⊕₊ normMK ⊗₊ ϵ
            nrmv = opnormbound(B, N, v)   # post-noise norm: input to the shift

            # the forcing step: (fractional) circular shift Shift(shifts[k]).
            # Convex combination of two permutations ⇒ mass-preserving, opnorm 1.
            s = float(shifts[k])
            n0 = floor(Int, s)
            fr = s - n0
            if fr == 0
                v = circshift(v, n0)                                   # exact permutation
            else
                v = (1 - fr) .* circshift(v, n0) .+ fr .* circshift(v, n0 + 1)
                ϵ = γₛ ⊗₊ nrmv ⊕₊ ϵ        # opnorm(shift)=1, radius 0; one extra mul-add/row
            end

            nrmv = opnormbound(B, N, v)
            add_column!(normcachers[k], v, ϵ)
        end
    end
    return map(get_norm, normcachers)
end

"""
Array of "trivial" bounds for the powers of a DiscretizedOperator (on the whole space)
coming from from ||Q^k|| ≤ ||Q||^k
"""
function norms_of_powers_trivial_noise(
    B::Basis,
    N::Type{<:NormKind},
    Q::DiscretizedOperator,
    MK::NoiseKernel,
    m::Integer,
)
    norms = fill(NaN, m)
    δₖ = opradius(N, MK)
    nrm_MK = opnormbound(B, N, MK)
    norms[1] = opnormbound(B, N, Q) ⊗₊ nrm_MK ⊕₊ δₖ
    for i = 2:m
        norms[i] = norms[i-1] ⊗₀ norms[1]
    end
    return norms
end

"""
Arrays of bounds to ||Q^k||_{w → s} = sup_{||f||_w=1} ||Q^k f||_s
and to ||Q^k||_{w}
coming theoretically from iterated DFLY inequalities (the "small matrix method").

Returns two arrays (strongs, norms) of length m:
strongs[k] bounds ||Q^k f||_s, norms[k] bounds ||Q^k f||)
"""
function norms_of_powers_abstract_noise(Bas::Basis, N::NoiseKernel, m)
    A, B = dfly(strong_norm(Bas), aux_norm(Bas), N)
    Eh = aux_normalized_projection_error(Bas)
    M₁n = strong_weak_bound(Bas)
    M₂ = aux_weak_bound(Bas)
    S₁, S₂ = weak_by_strong_and_aux_bound(Bas)

    norms = fill(NaN, m)
    strongs = fill(NaN, m)

    v = Array{Float64}([M₁n; M₂])
    # We evaluate [S₁ S₂] * ([1 0; Eh 1]*[A B; 0 1])^k * [M₁n; M₂] (with correct rounding)
    for k = 1:m
        # invariant: v[1] bounds ||Q^kf||_s for ||f||_w=1
        # v[2] bounds |||Q^kf||| for ||f||_w=1
        v[1] = A ⊗₊ v[1] ⊕₊ B ⊗₊ v[2]
        v[2] = Eh ⊗₊ v[1] ⊕₊ v[2]
        strongs[k] = v[1]
        norms[k] = S₁ ⊗₊ v[1] ⊕₊ S₂ ⊗₊ v[2]
    end
    return strongs, norms
end

"""
Estimate norms of powers from those on a coarser grid (see paper for details)
TODO: Check if it works for other basis types
"""
function norms_of_powers_from_coarser_grid_noise(
    fine_basis::Ulam,
    coarse_basis::Ulam,
    Q::DiscretizedOperator,
    NK::NoiseKernel,
    coarse_norms::Vector,
)
    if !is_refinement(fine_basis, coarse_basis)
        @error "The fine basis is not a refinement of the coarse basis"
    end
    m = length(coarse_norms)
    fine_norms = fill(NaN, m + 1)
    trivial_norms =
        norms_of_powers_trivial_noise(fine_basis, weak_norm(fine_basis), Q, NK, m + 1)

    A, B = dfly(strong_norm(fine_basis), aux_norm(fine_basis), NK)

    # adds a 0th element to strongs
    trivial_norms0(k::Integer) = k == 0 ? 1.0 : trivial_norms[k]
    coarse_norms0(k::Integer) = k == 0 ? 1.0 : coarse_norms[k]

    Kh = weak_projection_error(coarse_basis)

    fine_norms[1] = trivial_norms0(1)

    for k = 1:m
        temp = 0.0
        for i = 0:k-1
            temp =
                temp ⊕₊
                coarse_norms0(i) ⊗₊ (trivial_norms0(k - i) ⊕₊ trivial_norms0(k - i - 1))
        end
        fine_norms[k+1] = coarse_norms0(k) ⊕₊ B ⊗₊ Kh ⊗₊ (temp ⊕₊ 1.0)
    end
    return fine_norms
end

function norms_of_powers_from_coarser_grid_noise_abstract(
    coarse_basis::Ulam,
    NK::NoiseKernel,
    coarse_norms::Vector,
)

    m = length(coarse_norms)
    abstract_norms = fill(NaN, m + 1)
    trivial_norms = ones(m + 1)

    A, B = dfly(strong_norm(coarse_basis), aux_norm(coarse_basis), NK)

    # adds a 0th element to strongs
    trivial_norms0(k::Integer) = k == 0 ? 1.0 : trivial_norms[k]
    coarse_norms0(k::Integer) = k == 0 ? 1.0 : coarse_norms[k]

    Kh = weak_projection_error(coarse_basis)

    fine_norms = fill(NaN, m + 1)
    fine_norms[1] = trivial_norms0(1)

    for k = 1:m
        temp = 0.0
        for i = 0:k-1
            temp =
                temp ⊕₊
                coarse_norms0(i) ⊗₊ (trivial_norms0(k - i) ⊕₊ trivial_norms0(k - i - 1))
        end
        fine_norms[k+1] = coarse_norms0(k) ⊕₊ B ⊗₊ Kh ⊗₊ (temp ⊕₊ 1.0)
    end
    return fine_norms
end


"""
    powernormboundsnoise(B; Q, NK, threshold = 0.1)

Bounds on `‖(NK Q)^k|_{U⁰}‖` for `k = 1, 2, …`: [`norms_of_powers_noise`](@ref) with `m`
doubled from 8 until one bound is below `threshold`, the minimum with the trivial bounds
`(‖Q‖‖NK‖ + δ)^k`, and the extension by `C_{i+j} ≤ C_i C_j`
([`refine_norms_of_powers`](@ref)) until the last bound is below `1e-8`.
"""
function powernormboundsnoise(
    B;
    Q::DiscretizedOperator,
    NK::NoiseKernel,
    threshold = 0.1,
)
    m = 8
    computed_norms = []
    while true
        computed_norms =
            norms_of_powers_noise(B, weak_norm(B), m, Q, NK, integral_covector(B))
        if any(computed_norms .< threshold)
            break
        end
        m = 2 * m
    end
    trivial_norms = norms_of_powers_trivial_noise(B, weak_norm(B), Q, NK, m)
    norms = min.(trivial_norms, computed_norms)

    m_extend = 2 * m
    better_norms = []
    while true
        better_norms = refine_norms_of_powers(norms, m_extend)
        if better_norms[end] < 1e-8
            break
        end
        m_extend = 2 * m_extend
    end

    return better_norms

end

"""
Uses power norm bounds already computed for a coarse operator to estimate
the same norms for a finer operator
"""
function finepowernormboundsnoise(
    B,
    Bfine,
    coarse_norms;
    Qfine::DiscretizedOperator,
    NKfine::NoiseKernel,
)
    m = length(coarse_norms)

    trivial_norms_fine =
        norms_of_powers_trivial_noise(Bfine, weak_norm(Bfine), Qfine, NKfine, m + 1)
    twogrid_norms_fine =
        norms_of_powers_from_coarser_grid_noise(Bfine, B, Qfine, NKfine, coarse_norms)

    norms_fine = min.(trivial_norms_fine, twogrid_norms_fine)

    better_norms_fine = refine_norms_of_powers(norms_fine, m + 1)
    return better_norms_fine
end

function abstractpowernormboundsnoise(B, NK, coarse_norms; m = length(coarse_norms))
    trivial_norms_fine = ones(m + 1)
    abstract_norms_from_coarse =
        norms_of_powers_from_coarser_grid_noise_abstract(B, NK, coarse_norms)

    norms_fine = min.(trivial_norms_fine, abstract_norms_from_coarse)
    better_norms_abstract = refine_norms_of_powers(norms_fine, m + 1)
    return better_norms_abstract
end

"""
Return a numerical approximation to the (hopefully unique) invariant vector
of the dynamic with discretized operator Q.

The vector is normalized so that integral_covector(B)*w ≈ 1
"""
function invariant_vector_noise(
    B::Basis,
    Q::DiscretizedOperator,
    NK::NoiseKernel;
    tol = 0.0,
    iter = 10,
)
    mQ = mid(Q)
    v = one_vector(B)
    for _ = 1:iter
        v = mQ * v
        v = NK * v
    end
    v = v ./ (mid.(integral_covector(B)) * v)
    return v
end

"""
Return an upper bound to Q_h*w - w in the given norm
"""
function residualboundnoise(
    B::Basis,
    N::Type{<:NormKind},
    Q::DiscretizedOperator,
    NK::NoiseKernel,
    w::AbstractVector,
)
    v = Q * w
    v = NK * v # actually, a priori estimate... TODO modify mult NK*v
    return normbound(B, N, v - w)
end

"""
Bounds rigorously the distance of w from the fixed point of Q (normalized with integral = 1),
using a vector of bounds norms[k] ≥ ||Q_h^k|_{U_h^0}||.
"""
function distance_from_invariant_noise(
    B::Basis,
    Q::DiscretizedOperator,
    NK::NoiseKernel,
    w::AbstractVector,
    norms::Vector;
    ε₁::Float64 = residualboundnoise(B, weak_norm(B), Q, NK, w),
    ε₂::Float64 = mag(integral_covector(B) * w - 1),
    normQ::Float64 = opnormbound(B, weak_norm(B), Q),
)
    if ε₂ > 1e-8
        @error "w does not seem normalized correctly"
    end
    _, us = dfly(strong_norm(B), aux_norm(B), NK)
    Cs = infinite_sum_norms(norms)
    Kh = weak_projection_error(B)
    normw = normbound(B, weak_norm(B), w)

    return Cs ⊗₊ (2.0 ⊗₊ Kh ⊗₊ (1.0 ⊕₊ normQ) ⊗₊ us ⊕₊ ε₁ ⊘₊ (1.0 ⊖₋ ε₂)) ⊕₊
           ε₂ ⊘₊ (1.0 ⊖₋ ε₂) ⊗₊ normw
end
