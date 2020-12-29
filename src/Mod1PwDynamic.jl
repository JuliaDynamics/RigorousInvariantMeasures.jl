using ValidatedNumerics
using .DynamicDefinition, .PwDynamicDefinition
using .Contractors

using .DynamicDefinition: derivative

"""
Reimplementation of Mod1Dynamic as a PwMap instance.
We assume that T is monotonic and differentiable, for now (this is not restrictive, for our purposes)
TODO: having the domain as an interval isn't a good idea; it is better to use a pair of intervals for its endpoints
"""
function mod1_dynamic(T::Function, domain = Interval(0,1), ε = 1e-15)
    Tprime = x -> derivative(T, x)

    T0 = T(Interval(domain.lo))
    T1 = T(Interval(domain.hi))
    orientation = sign(T1 - T0)
    @assert(isinteger(orientation))
    orientation = mid(orientation)

    # check monotonicity
    @assert minimise(x -> Tprime(x)*orientation, domain)[1] > 0

    if !isthin(floor(T0))
        @error "T(I.lo) does not have a unique integer part; we did not implement this more complicated case"
    end

    first_branch_integer_part = floor(T0).lo
    if isinteger(T0) && orientation < 0
        first_branch_integer_part -= 1
    end

    last_branch_integer_part = floor(T1).hi
    if isinteger(T1) && orientation > 0
        last_branch_integer_part -= 1
    end

    if orientation > 0
        Ts = [x->T(x)-k for k in first_branch_integer_part:last_branch_integer_part]
        endpoints = [Interval(domain.lo); [root(x->T(x)-k, Tprime, domain, ε) for k in first_branch_integer_part+1:last_branch_integer_part]; Interval(domain.hi)]
        if length(Ts) == 1
            is_full = [isinteger(T0) && isinteger(T1)]
        else
            is_full = [isinteger(T0); fill(true, length(Ts)-2); isinteger(T1)]
        end
    else
        Ts = [x->T(x)-k for k in first_branch_integer_part:-1:last_branch_integer_part]
        endpoints = [Interval(domain.lo); [root(x->T(x)-k, Tprime, domain, ε) for k in first_branch_integer_part:-1:last_branch_integer_part+1]; Interval(domain.hi)]
        if length(Ts) == 1
            is_full = [isinteger(T0) && isinteger(T1)]
        else
            is_full = [isinteger(T0); fill(true, length(Ts)-2); isinteger(T1)]
        end
    end
    return PwMap(Ts, endpoints, is_full, fill(orientation, length(Ts)))
end
