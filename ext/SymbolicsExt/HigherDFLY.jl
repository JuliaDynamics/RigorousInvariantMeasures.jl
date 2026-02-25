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
    boolf = [(symb in keys(h.dict)) for symb in f]

    @variables (DD(x))[1:k+1]
    boolDD = [(symb in keys(h.dict)) for symb in DD]

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
