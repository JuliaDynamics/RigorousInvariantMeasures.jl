using LinearAlgebra, Arpack, FastRounding, IntervalArithmetic

export invariant_vector,
    finepowernormbounds,
    powernormbounds,
    distance_from_invariant,
    compute_coarse_grid_quantities,
    compute_fine_grid_quantities,
    one_grid_estimate,
    two_grid_estimate

"""
Return a numerical approximation to the (hopefully unique) invariant vector
of the dynamic with discretized operator Q.

The vector is normalized so that integral_covector(B)*w ≈ 1
"""
function invariant_vector(B::Basis, Q::DiscretizedOperator; tol = 0.0)
    mQ = mid(Q)
    n = size(Q)[1]
    # setting a larger nev seems to slow things down
    Te = eltype(mQ)
    F = eigs(mQ; tol = tol, nev = 1, ritzvec = true, v0 = ones(Te, n))
    w = F[2][:, 1]
    if Te <: Real
        @assert imag(w) ≈ zeros(n)
        w = real(w) # safe for real-matrix bases (Ulam, Hat)
    end
    w = w ./ (mid.(integral_covector(B)) * w) #normalization
    return w
end

"""
Return an upper bound to Q_h*w - w in the given norm
"""
function residualbound(
    B::Basis,
    N::Type{<:NormKind},
    Q::DiscretizedOperator,
    w::AbstractVector,
)
    return normbound(B, N, Q * w - w)
end

"""
Bounds rigorously the distance of w from the fixed point of Q (normalized with integral = 1),
using a vector of bounds norms[k] ≥ ||Q_h^k|_{U_h^0}||.
If ε₁ and normQ are given, then Q can be omitted
"""
function distance_from_invariant(
    B::Basis,
    D::Dynamic,
    Q::Union{DiscretizedOperator,Nothing},
    w::AbstractVector,
    norms::Vector;
    ε₁::Float64 = residualbound(B, weak_norm(B), Q, w),
    ε₂::Float64 = mag(integral_covector(B) * w - 1),
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
)
    if ε₂ > 1e-8
        @error "w does not seem normalized correctly"
    end
    us = invariant_measure_strong_norm_bound(B, D; dfly_coefficients = dfly_coefficients)
    Cs = infinite_sum_norms(norms)
    Kh = weak_projection_error(B)
    normw = normbound(B, weak_norm(B), w)
    normL = bound_weak_norm_abstract(B, D; dfly_coefficients = dfly_coefficients)

    return Cs ⊗₊ (2.0 ⊗₊ Kh ⊗₊ (1.0 ⊕₊ normL) ⊗₊ us ⊕₊ ε₁ ⊘₊ (1.0 ⊖₋ ε₂)) ⊕₊
           ε₂ ⊘₊ (1.0 ⊖₋ ε₂) ⊗₊ normw
end

@doc raw"""
    projection_defect_coefficients(B::Basis, D::Dynamic;
        dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
        normL = bound_weak_norm_abstract(B, D; dfly_coefficients = dfly_coefficients))

Certified constants ``(c_s, c_w)`` such that the discretization defect of the
(integral-preserving) discretized operator ``Q_h`` satisfies, for every ``v``
in the strong space,

```math
\|(L - Q_h)\,v\|_w \;\le\; K_h \,\big( c_s \|v\|_s + c_w \|v\|_{L^1} \big),
\qquad K_h = \texttt{weak\_projection\_error}(B).
```

Generic method, valid for every *compatible discretization* in the sense of
[Galatolo–Monge–Nisoli–Poloni, Chaos Solitons & Fractals 170 (2023) 113329,
Definition 2.7]: by Lemma 3.5 of that paper,
``\|(Q_h - L)f\|_w \le 2K_h(\|L\|_w \|f\|_s + \|Lf\|_s)``, and the one-step
Lasota–Yorke inequality (5) gives ``\|Lf\|_s \le A\|f\|_s + B\|f\|_{L^1}``,
whence

```math
(c_s, c_w) = \big(2(\|L\|_w + A),\; 2B\big).
```

Sharper basis-specific methods may be provided (see the `Ulam` method, where
``(c_s, c_w) = (1+A, B)``).
"""
function projection_defect_coefficients(
    B::Basis,
    D::Dynamic;
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
    normL = bound_weak_norm_abstract(B, D; dfly_coefficients = dfly_coefficients),
)
    A = interval(dfly_coefficients[1])
    Bd = interval(dfly_coefficients[2])
    return (2 * (interval(normL) + A), 2 * Bd)
end

@doc raw"""
    distance_from_invariant_residual(B::Basis, D::Dynamic, Q, w, norms;
        ε₁ = residualbound(B, weak_norm(B), Q, w),
        ε₂ = mag(integral_covector(B) * w - 1),
        dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
        defect_coefficients = projection_defect_coefficients(B, D; dfly_coefficients))

A posteriori (residual-based) rigorous upper bound for the weak-norm distance
``\|h - w\|_w`` between the invariant density ``h`` of the *abstract* transfer
operator ``L`` of the dynamic `D` and a computed candidate `w` in the
approximating space of `B`.

This is the "exchanged" version of [Galatolo–Monge–Nisoli–Poloni, Chaos
Solitons & Fractals 170 (2023) 113329, Theorem 3.4] anticipated in Remark 3.9
of that paper: there, the error ``v = h - w`` is expanded in powers of the
*discretized* operator ``Q_h``, and one pays the a priori strong-norm bound
``\|h\|_s \le B/(1-A)`` (Corollary 2.6) for the unknown density — that is
[`distance_from_invariant`](@ref).  Here the error is instead expanded in
powers of the *abstract* operator ``L``, driven by the computed residual of
the candidate; Remark 3.9 notes this requires summing ``\|L^k r\|_w``, "a
difficult task" a priori — the closed-form two-norm bound below performs that
summation using only the computed norms ``C_k`` of the discretized operator,
the one-step Lasota–Yorke inequality, and the discretization-defect constants.
It is sharper than the a priori bound whenever ``\|w\|_s \ll B/(1-A)``, i.e.
precisely when the Lasota–Yorke constants are poor (``A`` near ``1``, small
branches) — the candidate's strong norm is *computed*, not estimated.

# Framework and assumptions (notation of the cited paper)

* ``L`` preserves the integral, ``\|L\|_{L^1}\le 1``, and satisfies the
  one-step Lasota–Yorke inequality ``\|Lf\|_s \le A\|f\|_s + B\|f\|_{L^1}``
  with ``A < 1`` (`dfly_coefficients`, cf. [`dfly`](@ref)); the norms satisfy
  ``\|\cdot\|_{L^1} \le \|\cdot\|_w`` (Assumption 2.1(4) of the paper).
* `Q` is the rigorously assembled, integral-preserving discretized operator
  ``Q_h`` of a compatible discretization; `norms[k]` ``\ge
  \|Q_h^k|_{V_h^0}\|_w`` are certified (cf. [`powernormbounds`](@ref),
  [`finepowernormbounds`](@ref)); ``K_h`` = [`weak_projection_error`](@ref).
* `defect_coefficients` ``= (c_s, c_w)`` satisfy
  ``\|(L-Q_h)v\|_w \le K_h(c_s\|v\|_s + c_w\|v\|_{L^1})``
  (see [`projection_defect_coefficients`](@ref)).
* ``L`` has a **unique** invariant probability density ``h``; uniqueness
  identifies the Neumann series below with ``\bar w - h`` (any integral-zero
  fixed point of ``L`` vanishes by ergodic decomposition).  In practice this
  is certified from the same `norms` by the small-matrix method
  ([`convergencerateabstract`](@ref), Galatolo–Nisoli–Saussol).
* The strong and ``L^1`` norms of the candidate are bounded rigorously via
  `normbound(B, strong_norm(B), w)` and `normbound(B, aux_norm(B), w)`; the
  basis must provide these methods.

# Derivation

Let ``\bar w = w/\int w`` (normalization defect ``\varepsilon_2`` accounted at
the end) and ``r = L\bar w - \bar w`` the abstract residual, ``i(r) = 0``.
Then

```math
\bar w - h = -(I-L)^{-1}\big|_{V^0}\, r = -\sum_{k\ge0} L^k r,
\qquad \|h - \bar w\|_w \le \sum_{k\ge0} W_k,
```

with ``W_k := \|L^k r\|_w``, ``S_k := \|L^k r\|_s``.

**(1) Both norms of the residual are computable.**  Writing
``r = (L - Q_h)\bar w + (Q_h\bar w - \bar w)``, the second term is the
computed eigen-residual (``\le \varepsilon_1/(1-\varepsilon_2)``) and the
first pays the defect of the *explicit* candidate:

```math
W_0 \le \frac{\varepsilon_1 + K_h (c_s \|w\|_s + c_w \|w\|_{L^1})}{1-\varepsilon_2},
\qquad
S_0 \le \frac{(1+A)\|w\|_s + B\|w\|_{L^1}}{1-\varepsilon_2},
```

the latter from one Lasota–Yorke step applied to `w` and the triangle
inequality.

**(2) Two-channel recursion.**  Telescoping ``L^k = Q_h^k + \sum_{j<k} Q_h^j
(L - Q_h) L^{k-1-j}`` (all vectors have zero integral, so the restricted
norms apply) and iterating the Lasota–Yorke inequality:

```math
S_k \le A^k S_0 + B\sum_{i<k} A^i W_{k-1-i},\qquad
W_k \le C_k W_0 + K_h \sum_{j<k} C_j\,\big(c_s S_{k-1-j} + c_w W_{k-1-j}\big).
```

**(3) Sound closed-form summation.**  Summing over ``k \ge 0`` (Fubini for
nonnegative series), with ``\Sigma_C := \sum_{k\ge0} C_k`` bounded by
[`infinite_sum_norms`](@ref):

```math
\Sigma_S \le \frac{S_0}{1-A} + \frac{B}{1-A}\Sigma_W,\qquad
\Sigma_W \le \Sigma_C W_0 + K_h\Sigma_C\big(c_s\Sigma_S + c_w\Sigma_W\big),
```

whence, **provided the denominator below is positive** (checked; an error is
thrown otherwise),

```math
\|h - \bar w\|_w \;\le\; \Sigma_W \;\le\;
\frac{\Sigma_C W_0 + K_h \Sigma_C\, c_s\, S_0/(1-A)}
     {1 - K_h \Sigma_C \big(c_s\, B/(1-A) + c_w\big)} .
```

The same bound applies to every partial sum, so the Neumann series converges
absolutely and the identity above is justified.  Finally
``\|h - w\|_w \le \Sigma_W + \tfrac{\varepsilon_2}{1-\varepsilon_2}\|w\|_w``.

All arithmetic is carried out in interval arithmetic; the returned value is a
rigorous `Float64` upper bound.

!!! warning
    Do not replace the closed-form total with a truncated recursion and a
    heuristic geometric tail (e.g. capping an observed ratio): observed ratios
    can exceed any a priori cap, and such tails are *not* rigorous.
"""
function distance_from_invariant_residual(
    B::Basis,
    D::Dynamic,
    Q::DiscretizedOperator,
    w::AbstractVector,
    norms::Vector;
    ε₁::Float64 = residualbound(B, weak_norm(B), Q, w),
    ε₂::Float64 = mag(integral_covector(B) * w - 1),
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
    defect_coefficients = projection_defect_coefficients(
        B,
        D;
        dfly_coefficients = dfly_coefficients,
    ),
)
    if ε₂ > 1e-8
        @error "w does not seem normalized correctly"
    end
    A = interval(dfly_coefficients[1])
    Bd = interval(dfly_coefficients[2])
    sup(A) < 1 || error("distance_from_invariant_residual: DFLY A ≥ 1; use an iterate")
    cs = interval(defect_coefficients[1])
    cw = interval(defect_coefficients[2])

    Kh = interval(weak_projection_error(B))
    ΣC = interval(infinite_sum_norms(norms))
    ns = interval(normbound(B, strong_norm(B), w))    # ‖w‖_s (computed candidate)
    na = interval(normbound(B, aux_norm(B), w))       # ‖w‖_{L¹}
    nw = interval(normbound(B, weak_norm(B), w))      # ‖w‖_w

    scale = 1 / (1 - interval(ε₂))
    W₀ = (interval(ε₁) + Kh * (cs * ns + cw * na)) * scale     # ‖r‖_w
    S₀ = ((1 + A) * ns + Bd * na) * scale                      # ‖r‖_s

    den = 1 - Kh * ΣC * (cs * Bd / (1 - A) + cw)
    inf(den) > 0 || error(
        "distance_from_invariant_residual: closed-form denominator not positive; refine the grid or the norms",
    )
    ΣW = (ΣC * W₀ + Kh * ΣC * cs * S₀ / (1 - A)) / den

    return sup(ΣW + interval(ε₂) * scale * nw)
end

# """
# This function returns a sequence of Cᵢ, \\tilde{C}ᵢ for a matrix P
# on a subspace V such that ||P^i|_V||_1\\leq C_i and
# ||P^i|_V||_{\\infty}\\leq \\tilde{C}_i, with respect to the
# """
# function contractmatrix(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
# 	# vector of the Cᵢ for P
# 	C = zeros(m)
# 	S = zeros((length(B), m))
# 	tilde_C = zeros(m)
#
# 	PP = mid.(P)
#
# 	for v in BasisDefinition.AverageZero(B)
# 		λ₁, λ₂ = norm(v, 1), norm(v, Inf)
# 		for i in 1:m
# 			v = PP*v
# 			C[i] = max(C[i], norm(v, 1)/λ₁)
# 			S[:, i]+=abs.(v)/λ₂
# 		end
# 	end
#
#
# 	for k in 1:m
# 		tilde_C[k] = maximum(S[:,k])
# 	end
#
# 	# we keep track of the error due to the basis we have chosen
# 	η₁ = BasisDefinition.spaceconstant(B, Val(:L1))
# 	η₂ = BasisDefinition.spaceconstant(B, Val(:L∞))
#
#  	return η₁*C, η₂*tilde_C
# end
#
# @deprecate contractmatrix norms_of_powers

# """
# This function returns the bound on the weak norm of the discretized operator
# """
# function boundnorm(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
#     W₁, W₂ = bound_weak_norm_from_linalg_norm(B)
#     α₁ = bound_linalg_norm_L1_from_weak(B)
#     α₂ = bound_linalg_norm_L∞_from_weak(B)
#     C, tilde_C = contractmatrix(B, P, m)
#     return (W₁ / α₁) * C + (W₂ / α₂) * tilde_C
# end
# 
# @deprecate boundnorm

"""
Uses different strategies to compute power norm bounds.

If specified, `m` norms of powers are estimated computationally, and then
`m_extend` norms are obtained with a cheaper refinement process. Otherwise
these numbers are selected automatically.

A vector of length m_extend is returned, such that norms[k] ≥ ||Q_h^k|_{U_h^0}||
"""
function powernormbounds(B, D, m, m_extend; Q = DiscretizedOperator(B, D))
    normQ = opnormbound(B, weak_norm(B), Q)
    trivial_norms = norms_of_powers_trivial(normQ, m)
    computed_norms = norms_of_powers(B, weak_norm(B), m, Q, integral_covector(B))

    # not interesting at the moment
    #(dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, m)
    # in the current version, dfly_norms seem to be always larger and could be omitted
    # however they do not cost much to compute
    #norms = min.(trivial_norms, computed_norms, dfly_norms)
    norms = min.(trivial_norms, computed_norms)

    better_norms = refine_norms_of_powers(norms, m_extend)

    return better_norms
end

"""
Computes bounds for norms of powers, taking (optionally) minimum values for the number of norms to compute
"""
function powernormbounds(B, D; Q = DiscretizedOperator(B, D), m = 8, threshold = 0.1)
    computed_norms = []
    while true
        computed_norms = norms_of_powers(B, weak_norm(B), m, Q, integral_covector(B))
        if any(computed_norms .< threshold)
            break
        end
        m = 2 * m
    end
    normQ = opnormbound(B, weak_norm(B), Q)
    trivial_norms = norms_of_powers_trivial(normQ, m)
    # (dfly_strongs, dfly_norms) = norms_of_powers_dfly(B, D, m)
    # in the current version, dfly_norms seem to be always larger and could be omitted
    # however they do not cost much to compute
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
function finepowernormbounds(
    B,
    B_fine,
    D,
    coarse_norms;
    normQ_fine = opnormbound(B_fine, weak_norm(B_fine)DiscretizedOperator(B_fine, D)),
    dfly_coefficients = dfly(strong_norm(B_fine), aux_norm(B_fine), D),
)
    m = length(coarse_norms)

    trivial_norms_fine = norms_of_powers_trivial(normQ_fine, m)
    twogrid_norms_fine = norms_of_powers_from_coarser_grid(
        B_fine,
        B,
        D,
        coarse_norms,
        normQ_fine;
        dfly_coefficients = dfly_coefficients,
    )

    (dfly_strongs_fine, dfly_norms_fine) =
        norms_of_powers_dfly(B_fine, D, m; dfly_coefficients = dfly_coefficients)

    norms_fine = min.(trivial_norms_fine, twogrid_norms_fine, dfly_norms_fine)

    better_norms_fine = refine_norms_of_powers(norms_fine, m)
    return better_norms_fine
end

"""
Struct that encapsulates all the quantities computed from the fine basis that are needed in the two-grid estimate.
It is meant as an intermediate quantity that can be saved on the disk to avoid recomputing Q all the times
"""
struct FineGridQuantities
    B::Basis
    D::Dynamic
    normQ::Any
    w::Any
    ε₁::Any
    ε₂::Any
    time_assembling::Any # time to compute B, D, Q, normQ
    time_eigen::Any      # time to compute w, ε₁, ε₂
end
"""
Struct that encapsulates the additional quantities needed on the coarse basis for a two-grid estimate,
or on the (only) basis for a one-grid estimate. 
It is meant as an intermediate quantity that can be saved on the disk to avoid recomputing Q all the times.
"""
struct CoarseGridQuantities
    B::Basis
    D::Dynamic
    norms::Any
    dfly_coefficients::Any
    time_assembling::Any # time to compute B, D, Q (without normQ, needed in the two-grid estimate)
    time_norms::Any      # time to compute norms
    time_dfly::Any       # time to compute dfly_coefficients
end

"""
Compute FineGridQuantities, given a function f(n) that computes B, D, Q = f(n)
"""
function compute_fine_grid_quantities(f, n)
    time_assembling1 = @elapsed B, D, Q = f(n)
    time_assembling2 = @elapsed normQ = opnormbound(B, weak_norm(B), Q)
    time_eigen1 = @elapsed w = invariant_vector(B, Q)
    time_eigen2 = @elapsed ε₁, ε₂ =
        residualbound(B, weak_norm(B), Q, w), mag(integral_covector(B) * w - 1)
    return FineGridQuantities(
        B,
        D,
        normQ,
        w,
        ε₁,
        ε₂,
        time_assembling1 + time_assembling2,
        time_eigen1 + time_eigen2,
    )
end

"""
Compute FineGridQuantities _and_ CoarseGridQuantities, given a function f(n) that computes B, D, Q = f(n)
"""
function compute_coarse_grid_quantities(f, n; m = 8)
    time_assembling1 = @elapsed B, D, Q = f(n)
    time_assembling2 = @elapsed normQ = opnormbound(B, weak_norm(B), Q)
    time_eigen1 = @elapsed w = invariant_vector(B, Q)
    time_eigen2 = @elapsed ε₁, ε₂ =
        residualbound(B, weak_norm(B), Q, w), mag(integral_covector(B) * w - 1)
    time_norms = @elapsed norms = powernormbounds(B, D, Q = Q, m = m)
    time_dfly = @elapsed dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D)
    return CoarseGridQuantities(
        B,
        D,
        norms,
        dfly_coefficients,
        time_assembling1,
        time_norms,
        time_dfly,
    ),
    FineGridQuantities(
        B,
        D,
        normQ,
        w,
        ε₁,
        ε₂,
        time_assembling1 + time_assembling2,
        time_eigen1 + time_eigen2,
    )
end

"""
Compute a one-grid error estimate.

The first return argument is the error, the second is the time breakdown according to ["dfly", "assembling", "eigen", "norms", "error"]. 
(The sum of that vector is the total time taken)
"""
function one_grid_estimate(C::CoarseGridQuantities, F::FineGridQuantities)
    @assert(C.B == F.B)
    # Dynamics don't compare unfortunately
    time_error = @elapsed error = distance_from_invariant(
        F.B,
        F.D,
        nothing,
        F.w,
        C.norms;
        dfly_coefficients = C.dfly_coefficients,
        ε₁ = F.ε₁,
        ε₂ = F.ε₂,
    )
    return error, [C.time_dfly, F.time_assembling, F.time_eigen, C.time_norms, time_error]
end

"""
Compute a two-grid error estimate.

The first return argument is the error, the second is the time breakdown according to ["dfly", "coarse", "assembling", "eigen", "norms", "error"]. 
(The sum of that vector is the total time taken)
"""
function two_grid_estimate(C::CoarseGridQuantities, F::FineGridQuantities; m_extend = 400)
    @assert(is_refinement(F.B, C.B))
    # Dynamics don't compare unfortunately
    time_error_fine = @elapsed error_fine = distance_from_invariant(
        F.B,
        F.D,
        nothing,
        F.w,
        finepowernormbounds(
            C.B,
            F.B,
            F.D,
            refine_norms_of_powers(C.norms, m_extend);
            normQ_fine = F.normQ,
            dfly_coefficients = C.dfly_coefficients,
        );
        dfly_coefficients = C.dfly_coefficients,
        ε₁ = F.ε₁,
        ε₂ = F.ε₂,
    )
    return error_fine,
    [
        C.time_dfly,
        C.time_assembling + C.time_norms,
        F.time_assembling,
        F.time_eigen,
        time_error_fine,
    ]
end
