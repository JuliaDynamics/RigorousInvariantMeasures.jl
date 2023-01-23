using Symbolics

""" 
    Represent a symbolic transfer operator
    `(L_n f)(x) = \\sum_{y\\in T^{-1}(x)}\\frac{f(y)}{|T'(y)|^n}`
""" 
struct SymbL
    n
    f::Num
end



"""
    Compute the derivative of a symbolic transfer operator  
    `P` with respect to differential `∂`
    and with derivative of `1/T′ ` given by the function 
    `Dist`.
    Refer to Eq. 3.2 in Butterley-Kiamari-Liverani Locating Ruelle Pollicot Resonances
"""
Diff(P::SymbL, ∂::Differential, Dist) = [SymbL(P.n+1, ∂(P.f)), SymbL(P.n, P.n*P.f*Dist)]

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
    der_dict = Dict([[Dist => DD[1] ]; 
            [∂(f[i]) => f[i+1] for i in 1:k]; 
            [∂(DD[i]) => DD[i+1] for i in 1:k]])

    P = SymbL(1, f[1]) #Lf

    v = P   
    for i in 1:k
        v = Diff(v, ∂, Dist)
        for k in 1:length(v)
            # for each one of the elements in v
            # we first expand the derivatives
            temporary = expand_derivatives(v[k].f, true)
            # we know use the dictionary above 
            # ∂(f[i]) = f[i+1], and similarly for DD
            temporary = substitute(temporary, der_dict)
            # we simplify
            temporary = simplify(temporary; expand = true)
            # we store again in the vector
            v[k] = SymbL(v[k].n, temporary)
        end
    end
    return v
end

using SymbolicUtils

function _optimize_mult(k, n, h::SymbolicUtils.Mul, vals)
    
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
    # ||DDᵢ/(T')^n-1||_{∞}*||DD_i||_{∞}^{pow-1}
    
    bound = vals[i+1, n]*vals[i+1, 1]^(pow-1)

    for j in 1:i-1
        if boolDD[j]!= false
            pow = h.dict[DD[j]]
            bound *= vals[j+1, 1]*pow
        end
    end

    i = findlast(boolf)
    return h.coeff*bound*f[i]
end

function optimize_coefficients(k, v::Vector{SymbL}, vals)
    out = Num(0)
    for x in v
        n = x.n
        #@info "x" x
        if typeof(x.f.val) <: SymbolicUtils.Mul
            y = x.f.val
            val = _optimize_mult(k, n, y, vals)
            out+= _optimize_mult(k, n, x.f.val, vals)
        elseif typeof(x.f.val) <: SymbolicUtils.Add
            for y in keys(x.f.val.dict)                
                val = _optimize_mult(k, n, y, vals)
                out+= x.f.val.dict[y]*val
            end
        end
    end
    return out
end

function substitute_values(k, v::Vector{SymbL}, vals)
    λ = vals[1]
    DDs = vals[2:end] 
    @variables x 
    @variables f(x)[1:k] 
    @variables DD(x)[1:k]

    subsdict = Dict([DD[i]=>DDs[i] for i in 1:k])

    w = zeros(Num, length(v))

    for k in 1:length(v)
        n = v[k].n
        f = v[k].f
        #w[k] = substitute(f, subsdict)
        #w[k] = λ^(n-1)*w[k]
    end
    return w
end



