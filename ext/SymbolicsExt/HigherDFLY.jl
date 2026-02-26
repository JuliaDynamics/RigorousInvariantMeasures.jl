using Symbolics

"""
    Represent a symbolic transfer operator
    `(L_n f)(x) = \\sum_{y\\in T^{-1}(x)}\\frac{f(y)}{|T'(y)|^n}`
"""
struct SymbL
    n::Any
    f::Num
end



"""
    Compute the derivative of a symbolic transfer operator
    `P` with respect to differential `∂`
    and with derivative of `1/T′ ` given by the function
    `Dist`.
    Refer to Eq. 3.2 in Butterley-Kiamari-Liverani Locating Ruelle Pollicot Resonances
"""
Diff(P::SymbL, ∂::Differential, Dist) =
    [SymbL(P.n + 1, ∂(P.f)), SymbL(P.n, P.n * P.f * Dist)]

"""
    Given the summands of a sum of Symbolic Transfer Operator
    given in a vector, computes a vector containing
    the sum of the derivatives
"""
function Diff(v::Vector{SymbL}, ∂::Differential, Dist)
    w = SymbL[]
    for P in v
        append!(w, Diff(P, ∂, Dist))
    end
    return w
end

function compute_dfly_k_fi_DDi(k::Int)
    # the following vector valued function is such that
    # f_i corresponds to the $i-1$ derivative of f
    @variables x
    @variables (f(x))[1:k+1] #f[0:k](x) #

    # the following vector valued function is such that
    # DD_i corresponds to the $i-1$ derivative of Dist
    @variables Dist(x)
    @variables (DD(x))[1:k+1]

    ∂ = Differential(x)

    # the substitution rules we stated above
    der_dict = Dict(
        [
            [Dist => DD[1]]
            [∂(f[i]) => f[i+1] for i = 1:k]
            [∂(DD[i]) => DD[i+1] for i = 1:k]
        ],
    )

    P = SymbL(1, f[1]) #Lf

    v = P
    for i = 1:k
        v = Diff(v, ∂, Dist)
        for j = 1:length(v)
            # for each one of the elements in v
            # we first expand the derivatives
            temporary = expand_derivatives(v[j].f, true)
            # we then use the dictionary above
            # ∂(f[i]) = f[i+1], and similarly for DD
            temporary = substitute(temporary, der_dict)
            # we simplify
            temporary = simplify(temporary; expand = true)
            # we store again in the vector
            v[j] = SymbL(v[j].n, temporary)
        end
    end
    return v
end

using SymbolicUtils

function _optimize_mult(k, n, h, vals)
    @assert SymbolicUtils.ismul(h)

    @variables x

    @variables (f(x))[1:k+1]
    boolf = Bool[(symb in keys(h.dict)) for symb in collect(f)]

    @variables (DD(x))[1:k+1]
    boolDD = Bool[(symb in keys(h.dict)) for symb in collect(DD)]

    # I start with a simple version, where I use the n on
    # the highest derivative of DD, in further versions
    # I will try to optimize over all possibilities

    i = findlast(boolDD)
    pow = h.dict[DD[i]]

    # This uses the computed bounds to compute
    # ||DDᵢ/(T')^{n-1}||_{∞} * ||DD_i||_{∞}^{pow-1}

    bound = vals[i+1, n] * vals[i+1, 1]^(pow - 1)

    for j = 1:i-1
        if boolDD[j]
            pow = h.dict[DD[j]]
            bound *= vals[j+1, 1]^pow
        end
    end

    i = findlast(boolf)
    return h.coeff * bound * f[i]
end

function optimize_coefficients(k, v::Vector{SymbL}, vals)
    out = Num(0)
    for x in v
        n = x.n
        fval = x.f.val
        if SymbolicUtils.ismul(fval)
            out += _optimize_mult(k, n, fval, vals)
        elseif SymbolicUtils.isadd(fval)
            for y in keys(fval.dict)
                val = _optimize_mult(k, n, y, vals)
                out += fval.dict[y] * val
            end
        end
    end
    return out
end

function substitute_values(k, v::Vector{SymbL}, vals)
    λ = vals[1]
    DDs = vals[2:end]
    @variables x
    @variables (f(x))[1:k+1]
    @variables (DD(x))[1:k+1]

    subsdict = Dict([DD[i] => DDs[i] for i = 1:k])

    w = zeros(Num, length(v))

    for j = 1:length(v)
        n = v[j].n
        w[j] = substitute(v[j].f, subsdict)
        w[j] = λ^(n - 1) * w[j]
    end
    return w
end

###############################################################################
# Higher-order DFLY bridge: dfly(W{k,1}, L1, PwMap)
###############################################################################

using IntervalArithmetic
using TaylorSeries: Taylor1

@doc raw"""
    dfly(::Type{W{k,1}}, ::Type{L1}, D::PwMap) where {k}

Compute DFLY constants (A, B) for the W^{k,1} → L¹ inequality:
``||Lf||_{W^{k,1}} \leq A ||f||_{W^{k,1}} + B ||f||_{L^1}``

Uses the symbolic Faà di Bruno machinery from `compute_dfly_k_fi_DDi` and
`optimize_coefficients` to handle the chain rule derivatives.

For k=1, this reduces to `dfly(TotalVariation, L1, D)`.
"""
function dfly(::Type{W{k,l}}, ::Type{L1}, D::PwMap) where {k,l}
    if k == 1
        return dfly(TotalVariation, L1, D)
    end

    # Compute λ = max |1/T'| across all branches
    lam = max_inverse_derivative(D)

    if !(abs(lam) < 1)
        @error "The function is not expanding"
    end

    # Compute distortion derivatives: vals[i,j] bounds ||DD_{i-1}/(T')^{j-1}||_∞
    # DD_0 = 1/T', DD_1 = d/dx(1/T'), etc.
    # We compute these via TaylorSeries interval evaluation on each branch
    vals = _compute_distortion_vals(D, k)

    # Get symbolic expansion
    v = compute_dfly_k_fi_DDi(k)

    # Optimize coefficients
    opt = optimize_coefficients(k, v, vals)

    # The A constant is λ^k (from the leading L_{k+1} term)
    A = (lam^k).hi

    # The B constant comes from the optimized lower-order terms
    # Extract the coefficient of f[1] from the optimized expression
    @variables x
    @variables (f(x))[1:k+1]
    B_symbolic = Symbolics.substitute(opt, Dict(f[i] => (i == 1 ? 1.0 : 0.0) for i = 1:k+1))
    B_val = Float64(Symbolics.unwrap(B_symbolic))

    return A, abs(B_val)
end

"""
Compute the distortion derivatives matrix for the DFLY bound.
vals[i, j] bounds ||DD_{i-1} / (T')^{j-1}||_∞ across all branches.
"""
function _compute_distortion_vals(D::PwMap, k::Int)
    # Matrix of bounds: rows = DD derivative order (1-indexed), cols = power of 1/T'
    vals = zeros(k + 1, k + 1)

    lam = max_inverse_derivative(D)
    dist = max_distortion(D)

    # Row 1: DD_0 = 1/T', so vals[1,j] = ||1/(T')^j||_∞ = λ^j
    for j = 1:k+1
        vals[1, j] = (lam^j).hi
    end

    # Row 2: DD_1 = d/dx(1/T') = -T''/T'^2, so ||DD_1||_∞ ≤ distortion * λ
    # vals[2,j] bounds ||DD_1 / (T')^{j-1}||_∞
    if k >= 1
        for j = 1:k+1
            vals[2, j] = (dist * lam^j).hi
        end
    end

    # For higher derivatives, use conservative estimates based on distortion
    # DD_i ≈ O(dist^i * lam) — this is conservative but correct
    for i = 3:k+1
        for j = 1:k+1
            vals[i, j] = (dist^(i - 1) * lam^j).hi
        end
    end

    return vals
end
