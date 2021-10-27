using ValidatedNumerics
using .DynamicDefinition, .PwDynamicDefinition
using .Contractors

using .DynamicDefinition: derivative

"""
Reimplementation of Mod1Dynamic as a PwMap instance.
We assume that f is monotonic and differentiable, for now (this is not restrictive, for our purposes)
"""
function mod1_dynamic(f::Function, X = (0.,1.), ε = 0.0)
    br = Branch(f, X)

    # check monotonicity
    fprime = x -> derivative(f, x)
    @assert minimise(x -> fprime(x) * (br.increasing ? 1 : -1), hull(Interval.(X)...))[1] > 0

    Yhull = hull(Interval.(br.Y)...)
    possible_integer_parts = floor(Int, Yhull.lo):ceil(Int, Yhull.hi)

    x, integer_parts = preimages(possible_integer_parts, br, possible_integer_parts)

    ep = [x; X[end]]
    Ts = [x->f(x)-k for k in integer_parts]

    n = Base.length(x)
    if br.increasing
        y_endpoints::Matrix{Interval{Float64}} = hcat(fill(0., n), fill(1., n))
    else
        y_endpoints = hcat(fill(1., n), fill(0., n))
    end
    y_endpoints[1, 1] = br.Y[begin] - integer_parts[begin]
    if y_endpoints[1, 1] == 0.
        y_endpoints[1, 1] = 0. # hack to get rid of -0..0 intervals
    end
    y_endpoints[end, end] = br.Y[end] - integer_parts[end]
    if y_endpoints[end, end] == 0.
        y_endpoints[end, end] = 0. # hack to get rid of -0..0 intervals
    end


    # if !isthin(floor(T0))
    #     @error "T(I.lo) does not have a unique integer part; we did not implement this more complicated case"
    # end

    # first_branch_integer_part = floor(T0).lo
    # if isinteger(T0) && orientation < 0
    #     first_branch_integer_part -= 1
    # end

    # last_branch_integer_part = floor(T1).hi
    # if isinteger(T1) && orientation > 0
    #     last_branch_integer_part -= 1
    # end

    # if orientation > 0
    #     Ts = [x->T(x)-k for k in first_branch_integer_part:last_branch_integer_part]
    #     endpoints = [Interval(domain.lo); [preimage(k, T, domain, ε) for k in first_branch_integer_part+1:last_branch_integer_part]; Interval(domain.hi)]
    #     if length(Ts) == 1
    #         y_endpoints = [(isinteger(T0) ? 0 : T0-floor(T0), isinteger(T1) ? 1 : T1-floor(T1))]
    #     else
    #         y_endpoints = [(isinteger(T0) ? 0 : T0-floor(T0)) , 1; fill((0,1), length(Ts)-2); (0, isinteger(T1) ? 1 : T1-floor(T1))]
    #     end
    # else
    #     Ts = [x->T(x)-k for k in first_branch_integer_part:-1:last_branch_integer_part]
    #     endpoints = [Interval(domain.lo); [preimage(k, T, domain, ε) for k in first_branch_integer_part:-1:last_branch_integer_part+1]; Interval(domain.hi)]
    #     if length(Ts) == 1
    #         y_endpoints = [(isinteger(T0) ? 1 : T0-floor(T0), isinteger(T1) ? 0 : T1-floor(T1))]
    #     else
    #         y_endpoints = [(isinteger(T0) ? 1 : T0-floor(T0), 0); fill((1,0), length(Ts)-2); (1, isinteger(T1) ? 0 : T1-floor(T1))]
    #     end
    # end
    return PwMap(Ts, ep, y_endpoints, fill(br.increasing, Base.length(Ts)))
end
