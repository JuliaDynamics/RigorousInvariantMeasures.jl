
"""
(A, B) = dfly(strongnorm, auxnorm, dynamic)

Constants (A, B) such that ||Lf||_s ≦ A||f||_s + B||f||_aux
"""
dfly(::Type{<:NormKind}, ::Type{<:NormKind}, ::Dynamic) = @error "Not implemented"
# Instance dispatch fallback: forward to type dispatch for simple norms
dfly(n1::NormKind, n2::Type{<:NormKind}, D::Dynamic) = dfly(typeof(n1), n2, D)

# I don't think this is used in production anymore
# function dfly(N1::Type{VarInterval}, N2::Type{L1}, D::Dynamic)
#     dist = max_distortion(D)
#     lam = max_inverse_derivative(D)

#     if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
#         @error "The function is not expanding"
#     end

#     if is_full_branch(D)
#         return sup(lam), sup(dist)
#     else
#         if !(abs(lam) < 0.5)
#             @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
#         end
#         endpts = endpoints(D)
#         min_width = minimum([endpts[i+1] - endpts[i] for i = 1:length(endpts)-1])
#         return sup(lam), sup(dist) ⊕₊ sup(2 / min_width)
#     end
# end

function dfly(N1::Type{VarInterval}, N2::Type{L1}, D::PwMap)
    if has_infinite_derivative_at_endpoints(D)
        return dfly_inf_der(N1, N2, D, 10^-3)
    end

    dist = max_distortion(D)
    lam = max_inverse_derivative(D)
    vec = endpoints(D)
    disc = maximum(2 / abs(vec[i] - vec[i+1]) for i = 1:nbranches(D))

    if is_full_branch(D)
        if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
            @error "The function is not expanding"
        end
        return DFLYCoefficients(VarInterval, L1, sup(lam), sup(dist))
    else
        if !(abs(2 * lam) < 1)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        return DFLYCoefficients(VarInterval, L1, sup(2 * lam), sup(dist + disc))
    end
end

# ---------------------------------------------------------------------------------------
# VarCircle: the same operator, the variation that counts the jump across 0 ~ 1.
# ---------------------------------------------------------------------------------------

function _assert_unit_circle_domain(D::PwMap)
    a, b = domain(D)
    if !(isequal_interval(a, interval(0)) && isequal_interval(b, interval(1)))
        throw(ArgumentError(
            "VarCircle is the variation on [0,1] with 0 and 1 identified, " *
            "and this PwMap has domain ($a, $b)",
        ))
    end
    return nothing
end

# Whether T(1^-) and T(0^+) are the same point of the circle, 0 and 1 being identified;
# `false` whenever the stored enclosures do not settle it, which is the side that pays more.
function _same_circle_point(y1, y2)
    zeroY, oneY = interval(0), interval(1)
    isequal_interval(y1, y2) && return true
    isequal_interval(y1, zeroY) && isequal_interval(y2, oneY) && return true
    isequal_interval(y1, oneY) && isequal_interval(y2, zeroY) && return true
    return false
end

@doc raw"""
    dfly(VarCircle, L1, D::PwMap) -> DFLYCoefficients

Constants ``(A,B)`` such that the transfer operator of `D`, read as a map of the circle
``[0,1]/\{0\sim 1\}``, satisfies

```math
\mathrm{Var}_{S^1}(Lf) \le A\,\mathrm{Var}_{S^1}(f) + B\,\|f\|_{L^1},
\qquad
\mathrm{Var}_{S^1}(f) = \mathrm{Var}_{(0,1)}(f) + |f(1^-)-f(0^+)| ,
```

the norm of [`VarCircle`](@ref).  Two statements are proved below:
[`_dfly_circle_from_lift`](@ref), under the hypothesis of
[`has_monotone_mod1_lift`](@ref), which every map built by [`mod1_dynamic`](@ref)
satisfies, and [`_dfly_circle_branchwise`](@ref), which assumes nothing beyond
expansivity.  The second is returned when the first does not apply; where both apply the
pair returned is the one with the smaller ``B/(1-A)``, the a priori bound on the strong
norm of the invariant density that [`invariant_measure_strong_norm_bound`](@ref) reads off
the pair, a pair failing ``A<1`` ranking last and a tie going to the lift.

Two terms separate these from `dfly(VarInterval, L1, D)`, which is derived from the same
ingredients -- the minimum of ``|T'|``, the branch lengths and a bound on ``T''/(T')^2``.
The interval statement bounds a point value of ``f`` on the ``i``-th branch by
``\inf_{I_i}|f| + \mathrm{Var}_{I_i}(f)`` and so pays ``\|f\|_{L^1(I_i)}/|I_i|``, which is
where its ``B`` picks up the ``2/\min_i|I_i|`` of the non-full-branch case; on the circle
the two arcs joining two points carry between them the whole variation, so a point value
costs ``\|f\|_{L^1} + \frac12\mathrm{Var}_{S^1}(f)`` instead and the branch length leaves
``B`` altogether.  What appears in its place is the variation ``Lf`` acquires across ``0``,
which the interval norm does not see; for a monotone mod 1 lift that variation is paid for
by the jumps of ``f`` itself except at the two ends of the domain, and the whole of the
circle statement's ``B`` comes from those two ends.

On the beta shift ``T_0(x) = \beta x \bmod 1`` with ``\beta = (1+\sqrt{23})/2`` the pair is
``A = B = 2/\beta = 0.690151\ldots``, against ``A = 0.690151\ldots``,
``B = 2/(1-2/\beta) = 6.45476\ldots`` on the interval, so the a priori bound ``B/(1-A)`` on
the invariant density is ``2/(\beta-2) = 2.22738\ldots`` in a norm larger than the one in
which the interval pair gives ``20.8320\ldots``.

Neither statement dominates the other: the lift pays ``\theta`` and
``|\lambda_L-\lambda_R|`` for the jumps of ``1/|T'|`` at the branch points, which the
branchwise statement charges to the distortion instead, so a piecewise linear map all of
whose branches are onto the circle is sharper read branchwise -- for
``4x`` on ``(0,1/2)`` continued by ``8x`` on ``(1/2,1)`` the branchwise pair is
``(1/2, 0)``, and ``B=0`` is sharp because Lebesgue is invariant, against ``(3/8, 1/4)``
from the lift.

Throws an `ArgumentError` for a `PwMap` whose domain is not ``[0,1]``, and for one with an
infinite derivative at an endpoint: [`dfly_inf_der`](@ref) covers that case on the interval
and has no circle counterpart here.
"""
function dfly(::Type{VarCircle}, ::Type{L1}, D::PwMap)
    _assert_unit_circle_domain(D)
    if has_infinite_derivative_at_endpoints(D)
        throw(ArgumentError(
            "no VarCircle statement is proved for a map with an infinite derivative at " *
            "a branch endpoint; dfly_inf_der covers that case on the interval",
        ))
    end
    c = _dfly_circle_branchwise(D)
    if has_monotone_mod1_lift(D)
        c = _better_dfly_pair(_dfly_circle_from_lift(D), c)
    end
    if !(c.A < 1.0)
        @error "Expansivity is insufficient for a DFLY on the circle; try an iterate." c
    end
    return c
end

# Both pairs are proved, so either may be returned; `B/(1-A)` ranks them because that is
# what the pair is used for downstream, and a pair that does not contract comes last.
_dfly_apriori_bound(c) = c.A < 1.0 ? (c.B ⊘₊ (1.0 ⊖₋ c.A)) : Inf
function _better_dfly_pair(c1, c2)
    h1, h2 = _dfly_apriori_bound(c1), _dfly_apriori_bound(c2)
    if isinf(h1) && isinf(h2)
        return c1.A <= c2.A ? c1 : c2
    end
    return h1 <= h2 ? c1 : c2
end

@doc raw"""
    _dfly_circle_from_lift(D::PwMap) -> DFLYCoefficients

The pair of the following lemma, which needs [`has_monotone_mod1_lift`](@ref).

**Lemma.** Let ``T:[0,1]\to[0,1]`` have branches ``T_i`` on ``I_i=(a_{i-1},a_i)``,
``0=a_0<\dots<a_m=1``, with ``T_i\in C^2(\bar I_i)``,

```math
\lambda := \max_i \sup_{\bar I_i}|1/T_i'| < 1,
\qquad
\mathrm{dist} := \max_i \sup_{\bar I_i}|T_i''/(T_i')^2| < \infty,
```

and let ``T`` be the reduction mod ``1`` of a strictly monotone ``\tilde T:[0,1]\to\mathbb R``,
that is, ``T(a_i^-)=1`` and ``T(a_i^+)=0`` at every interior branch point for increasing
branches, the two exchanged for decreasing ones.  Put

```math
\lambda_L := \frac{1}{|T'(0^+)|},\qquad
\lambda_R := \frac{1}{|T'(1^-)|},\qquad
\theta := \sum_{i=1}^{m-1}\left| \frac{1}{|T'(a_i^-)|} - \frac{1}{|T'(a_i^+)|}\right| ,
```

```math
\kappa := \theta +
\begin{cases}
|\lambda_L-\lambda_R| & \text{if } T(1^-) \text{ and } T(0^+) \text{ are the same point of } S^1,\\
\lambda_L+\lambda_R & \text{otherwise.}
\end{cases}
```

Then, for every ``f\in BV(S^1)``,

```math
\mathrm{Var}_{S^1}(Lf) \le \left(\lambda+\frac{\kappa}{2}\right)\mathrm{Var}_{S^1}(f)
 + (\mathrm{dist}+\kappa)\,\|f\|_{L^1} .
```

*Proof.* Four steps.

(1) *Point values.* For ``x,y`` points of approximate continuity of ``f``, the two arcs of
the circle joining them carry between them ``\mathrm{Var}_{S^1}(f)``, and ``|f(x)-f(y)|``
is at most the variation along either of them, so ``|f(x)-f(y)|\le\frac12\mathrm{Var}_{S^1}(f)``;
integrating in ``y`` over the circle, whose length is ``1``,

```math
\|f\|_\infty \le \|f\|_{L^1} + \tfrac12\mathrm{Var}_{S^1}(f).
```

(2) *``Lf`` is a periodisation.*  Let ``G=\tilde T^{-1}:(\alpha_-,\alpha_+)\to(0,1)``, with
``\{\alpha_-,\alpha_+\}=\{\tilde T(0),\tilde T(1)\}``, let ``u = (f\circ G)\,|G'|`` on
``(\alpha_-,\alpha_+)`` and let ``u`` be extended by ``0`` to ``\mathbb R``.  The preimages
of ``y\in S^1`` are the ``G(y+n)`` with ``y+n\in(\alpha_-,\alpha_+)``, so
``Lf(y) = \sum_{n\in\mathbb Z} u(y+n)``.

(3) *The variation of a periodisation.*  For ``w\in BV(\mathbb R)`` of compact support and
``Pw=\sum_n w(\cdot+n)``, the points of a cyclic partition of the circle, translated by all
integers, form an increasing bi-infinite sequence in ``\mathbb R`` along which the cyclic
sum for ``Pw`` is bounded term by term by the variation sum for ``w``; hence
``\mathrm{Var}_{S^1}(Pw)\le\mathrm{Var}_{\mathbb R}(w)
=\mathrm{Var}_{(\alpha_-,\alpha_+)}(w)+|w(\alpha_-^+)|+|w(\alpha_+^-)|``.  When
``\alpha_+-\alpha_-\in\mathbb N``, which happens exactly when ``T(1^-)`` and ``T(0^+)`` are
the same point of ``S^1``, the two end jumps fall on one point of the circle and only their
difference survives:
``\mathrm{Var}_{S^1}(Pw)\le\mathrm{Var}_{(\alpha_-,\alpha_+)}(w)+|w(\alpha_-^+)-w(\alpha_+^-)|``.

(4) *The two pieces.*  Let ``t_i`` be the integer with ``G(t_i)=a_i``.  On each
``(t_{i-1},t_i)`` the factor ``|G'|`` is absolutely continuous, and
``\mathrm{Var}_J(\phi\psi)\le\sup_J|\psi|\,\mathrm{Var}_J(\phi)+\int_J|\phi|\,|\psi'|``
with ``\phi=f\circ G``, ``\psi=|G'|`` gives
``\lambda\,\mathrm{Var}_{I_i}(f)+\int_{I_i}|f|\,|T''/(T')^2|\,dx`` after the change of
variables ``y=T(x)``, since ``|G''/G'|=|T''/(T')^2|`` and ``\mathrm{Var}(f\circ G)=
\mathrm{Var}(f)`` for monotone ``G``.  At an interior ``t_i`` the jump of ``u`` is at most
``\lambda|f(a_i^+)-f(a_i^-)|`` plus ``\|f\|_\infty`` times the jump of ``1/|T'|`` there.
Summing, and recognising ``\sum_i\mathrm{Var}_{I_i}(f)`` plus the interior jumps of ``f`` as
``\mathrm{Var}_{(0,1)}(f)``,

```math
\mathrm{Var}_{(\alpha_-,\alpha_+)}(u) \le \lambda\,\mathrm{Var}_{(0,1)}(f)
 + \mathrm{dist}\,\|f\|_{L^1} + \theta\,\|f\|_\infty .
```

For the end terms, ``|u(\alpha_-^+)| = \lambda_L|f(0^+)|`` and
``|u(\alpha_+^-)| = \lambda_R|f(1^-)|`` for an increasing lift, the two exchanged for a
decreasing one, and writing ``\lambda_L|a|+\lambda_R|b| =
\lambda_R(|b|-|a|)+(\lambda_L+\lambda_R)|a|`` for ``|a|\le|b|`` and symmetrically otherwise,

```math
|u(\alpha_-^+)|+|u(\alpha_+^-)| \le \lambda\,|f(1^-)-f(0^+)| + (\lambda_L+\lambda_R)\|f\|_\infty ,
```

using ``\max(\lambda_L,\lambda_R)\le\lambda``; in the case of step (3) where only the
difference survives, ``\lambda_L f(0^+)-\lambda_R f(1^-) = \lambda_L(f(0^+)-f(1^-)) +
(\lambda_L-\lambda_R)f(1^-)``, and symmetrically, give the same bound with
``|\lambda_L-\lambda_R|`` in place of ``\lambda_L+\lambda_R``.  In both cases the first
summand is ``\lambda`` times the jump of ``f`` across ``0``, which completes
``\lambda\,\mathrm{Var}_{(0,1)}(f)`` into ``\lambda\,\mathrm{Var}_{S^1}(f)``, and what is
left is ``\kappa\|f\|_\infty``.  Step (1) turns that into
``\kappa\|f\|_{L^1}+\frac{\kappa}{2}\mathrm{Var}_{S^1}(f)``. ``\square``

A map with no partial branch, a ``C^1`` lift, so that ``\theta=0``, and
``|T'(0^+)|=|T'(1^-)|`` has ``\kappa=0``, hence ``A=\lambda`` and ``B=\mathrm{dist}``,
which is the pair of the full-branch case of `dfly(VarInterval, L1, D)`: for
``x\mapsto 2x`` it is ``A=1/2``, ``B=0``, and ``B/(1-A)=0`` is the variation of the
invariant density, Lebesgue.  The beta shift ``\beta x \bmod 1`` has a linear lift, so
``\mathrm{dist}=\theta=0``, and one partial branch, so
``\kappa=\lambda_L+\lambda_R=2/\beta``, giving ``A=B=2/\beta``.
"""
function _dfly_circle_from_lift(D::PwMap)
    _assert_unit_circle_domain(D)
    if !has_monotone_mod1_lift(D)
        throw(ArgumentError(
            "this statement needs D to be the reduction mod 1 of a monotone map, and " *
            "has_monotone_mod1_lift(D) is false",
        ))
    end
    brs = branches(D)
    m = length(brs)

    lam = max_inverse_derivative(D)
    dist = max_distortion(D)
    if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
        @error "The function is not expanding"
    end

    λ_L = abs(inverse_derivative(brs[1].f, brs[1].X[1]))
    λ_R = abs(inverse_derivative(brs[end].f, brs[end].X[2]))
    θ = sum(
        abs(
            abs(inverse_derivative(brs[i].f, brs[i].X[2])) -
            abs(inverse_derivative(brs[i+1].f, brs[i+1].X[1])),
        ) for i = 1:m-1;
        init = interval(0),
    )
    closes = _same_circle_point(brs[1].Y[1], brs[end].Y[2])
    κ = θ + (closes ? abs(λ_L - λ_R) : λ_L + λ_R)

    # `max_inverse_derivative` maximises over the closed branches, so it bounds λ_L and
    # λ_R as well; the `max` only guards against its enclosure being looser than the two
    # endpoint evaluations, which would make A smaller than the λ of the lemma.
    λ = max(sup(lam), sup(λ_L), sup(λ_R))
    A = λ ⊕₊ (sup(κ) ⊘₊ 2.0)
    B = sup(dist) ⊕₊ sup(κ)
    return DFLYCoefficients(VarCircle, L1, A, B)
end

@doc raw"""
    _dfly_circle_branchwise(D::PwMap) -> DFLYCoefficients

The pair of the following lemma, which assumes no relation between the branches; it is what
`dfly(VarCircle, L1, D)` returns when [`has_monotone_mod1_lift`](@ref) fails, and also when
it holds but this pair has the smaller ``B/(1-A)``.

**Lemma.** Let ``T:[0,1]\to[0,1]`` have branches ``T_i\in C^2(\bar I_i)`` on
``I_i=(a_{i-1},a_i)``, ``0=a_0<\dots<a_m=1``, and put ``h_i=|I_i|``,
``\lambda_i=\sup_{\bar I_i}|1/T_i'|``, ``\mathrm{dist}_i=\sup_{\bar I_i}|T_i''/(T_i')^2|``,

```math
A_i = \begin{cases} 2\lambda_i + \mathrm{dist}_i\,h_i, & T(I_i) = (0,1),\\
3\lambda_i, & \text{otherwise,}\end{cases}
\qquad
B_i = \begin{cases} 2\,\mathrm{dist}_i, & T(I_i) = (0,1),\\
\mathrm{dist}_i + 2\lambda_i/h_i, & \text{otherwise.}\end{cases}
```

Then ``\mathrm{Var}_{S^1}(Lf) \le (\max_i A_i)\,\mathrm{Var}_{S^1}(f) +
(\max_i B_i)\,\|f\|_{L^1}`` for every ``f \in BV(S^1)``.

*Proof.*  Write ``Lf=\sum_i\psi_i`` with ``\psi_i=(f\circ g_i)|g_i'|`` on ``J_i=T(I_i)``,
extended by ``0``, and ``V_i=\mathrm{Var}_{I_i}(f)``.  On ``J_i``, the product inequality
``\mathrm{Var}_J(\phi\psi)\le\sup_J|\psi|\mathrm{Var}_J(\phi)+\int_J|\phi||\psi'|`` and the
change of variables ``y=T(x)`` give
``\mathrm{Var}_{J_i}(\psi_i)\le\lambda_iV_i+\mathrm{dist}_i\|f\|_{L^1(I_i)}``.

A point value of ``f`` on ``I_i`` is at most ``\inf_{I_i}|f| + V_i \le
\|f\|_{L^1(I_i)}/h_i + V_i``, and ``|f(a_i^-)-f(a_{i-1}^+)|\le V_i``.

If ``J_i`` is a proper arc, ``\mathrm{Var}_{S^1}(\psi_i)=\mathrm{Var}_{J_i}(\psi_i)
+|\psi_i(\text{left end})|+|\psi_i(\text{right end})|``, the two end terms being at most
``\lambda_i`` times a point value of ``f`` on ``I_i``; this is the ``3\lambda_i`` and the
``\mathrm{dist}_i+2\lambda_i/h_i``.  If instead ``J_i=(0,1)``, so that ``J_i`` is the whole
circle, ``\mathrm{Var}_{S^1}(\psi_i)=\mathrm{Var}_{J_i}(\psi_i)+
|\psi_i(1^-)-\psi_i(0^+)|``, and with ``\delta_i=\mathrm{Var}_{J_i}(|g_i'|)
=\int_{I_i}|T''/(T')^2|\,dx \le \mathrm{dist}_i h_i`` that last jump is at most
``\lambda_i|f(a_i^-)-f(a_{i-1}^+)| + \delta_i\min(|f(a_i^-)|,|f(a_{i-1}^+)|)
\le \lambda_iV_i+\delta_iV_i+\delta_i\|f\|_{L^1(I_i)}/h_i``, where
``\delta_i/h_i\le\mathrm{dist}_i``; this is the ``2\lambda_i+\mathrm{dist}_ih_i`` and the
``2\,\mathrm{dist}_i``.  Summing over the branches, with
``\sum_iV_i\le\mathrm{Var}_{(0,1)}(f)\le\mathrm{Var}_{S^1}(f)`` and
``\sum_i\|f\|_{L^1(I_i)}=\|f\|_{L^1}``, gives the statement. ``\square``

The ``3\lambda_i`` of a partial branch is what keeps this statement off the beta shift,
where it returns ``A=3/\beta=1.0352\ldots``.  The other pair that needs no hypothesis is
the conversion of the interval one by ``\mathrm{Var}_{S^1}\le2\,\mathrm{Var}_{(0,1)}``,
which doubles both constants, and neither dominates the other: on the beta shift the
doubling gives ``A=4/\beta=1.3803\ldots``, worse than the ``3/\beta`` here, while for a
map all of whose branches are onto the circle it gives ``(2\lambda, 2\,\mathrm{dist})``
and so beats this pair by the term ``\mathrm{dist}_ih_i`` in ``A``.  Only this one is
computed; the doubling is a line of arithmetic on `dfly(VarInterval, L1, D)` for a caller
who wants it.
"""
function _dfly_circle_branchwise(D::PwMap)
    _assert_unit_circle_domain(D)
    X = (interval(0), interval(1))
    A, B = 0.0, 0.0
    for br in branches(D)
        # the tolerance the rest of the package maximises with, as in `max_distortion`
        λ_i = branch_inverse_derivative(br, 1e-3)
        dist_i = bound_branch_distortion(br; tol = 1e-3)
        h_i = abs(br.X[2] - br.X[1])
        if is_full_branch(br, X)
            A_i = interval(2) * λ_i + dist_i * h_i
            B_i = interval(2) * dist_i
        else
            A_i = interval(3) * λ_i
            B_i = dist_i + interval(2) * λ_i / h_i
        end
        A = max(A, sup(A_i))
        B = max(B, sup(B_i))
    end
    return DFLYCoefficients(VarCircle, L1, A, B)
end

# W{1,1} → VarInterval bridge: the W^{1,1} seminorm is the total variation on the interval,
# so the constants are the same numbers and only the name of the norm changes.
# Accept both type and instance dispatch
_as_w11(c) = DFLYCoefficients(W{1,1}, L1, c.A, c.B)
dfly(::Type{W{1,1}}, N2::Type{L1}, D::PwMap) = _as_w11(dfly(VarInterval, L1, D))
dfly(::W{1,1}, N2::Type{L1}, D::PwMap) = _as_w11(dfly(VarInterval, L1, D))

# Higher-order W{k,1} (k≥2) is provided by the SymbolicsExt extension.
# Without Symbolics loaded, the generic dfly fallback will @error "Not implemented".

# Aη DFLY: analytic strip norm. Implemented in AnalyticDFLY.jl, which is
# included later — it needs the complex-neighbourhood machinery.

function dfly(::Type{Lipschitz}, ::Type{L1}, D::Dynamic)
    # TODO: should assert that D is globally C2 instead, but we don't have that kind of infrastructure yet.
    @assert is_full_branch(D)

    dist = max_distortion(D)
    #@info dist
    lam = max_inverse_derivative(D)
    #@info lam

    return DFLYCoefficients(Lipschitz, L1, sup(lam * (2 * dist + 1)),
                            sup(dist * (dist + 1)))
end
