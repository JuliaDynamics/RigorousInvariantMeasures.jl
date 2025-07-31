
## the function branches is after the module

using ..RigorousInvariantMeasures:
    derivative, distortion, inverse_derivative, unique_increasing
using TaylorSeries: Taylor1
using IntervalArithmetic, IntervalOptimisation

export PwMap,
    preim,
    nbranches,
    plottable,
    branches,
    MonotonicBranch,
    mod1_dynamic,
    dfly_inf_der,
    composedPwMap,
    equal_up_to_order,
    has_infinite_derivative_at_endpoints,
    intersect_domain,
    intersect_domain_bool

@doc raw"""
Type used to represent a "branch" of a dynamic. The branch is represented by a map `f` with domain `X=(a,b)`. X[1] and X[2] are interval enclosures of a,b.

The map must be monotonic on [a,b]. Note that this is not the same thing as being monotonic on hull(X[1], X[2]): 
for instance, take the map x → (x-√2)^2 on [√2, 1]: the left endpoint X[1] will be prevfloat(√2)..nextfloat(√2), 
but then the map is not monotonic on the whole hull(X[1], X[2]) because it also contains points that lie left of √2.
This is a tricky case that must be dealt with.

Enclosures `Y[1], Y[2]` for f(a), f(b) and `increasing` may be provided (for instance if we know that `Y=(0,1)`), otherwise they are computed automatically.
"""
struct MonotonicBranch{T<:Function,S<:Interval}
    f::T
    X::Tuple{S,S}
    Y::Tuple{S,S}
    increasing::Bool
end

MonotonicBranch(
    f::Function,
    X,
    Y = (f(X[1]), f(X[2])),
    increasing = unique_increasing(Y[1], Y[2]),
) = MonotonicBranch{typeof(f),typeof(interval(X[1]))}(f, X, Y, increasing)

function equal_up_to_order(X, Y)
    @assert length(X) == 2 && length(Y) == 2
    if eltype(X) <: Interval && !all(isthin.(X))
        return false
    end
    if X[1] == Y[1] && X[2] == Y[2]
        return true
    end
    if X[1] == Y[2] && X[2] == Y[1]
        return true
    end
    return false
end

is_full_branch(b::MonotonicBranch{T,S}, X) where {T,S} = equal_up_to_order(b.Y, X)
is_increasing(b::MonotonicBranch) = b.increasing

import Base.reverse

@doc raw"""
    Base.reverse(br:MonotonicBranch)

"Reverses" the x-axis of a branch: given ``f:[a,b] -> R``, creates a branch with the function ``g:[-b,-a] -> R`` defined as ``g(x) = f(-x)``
"""
Base.reverse(br::MonotonicBranch) = MonotonicBranch(
    x -> br.f(-x),
    (-br.X[2], -br.X[1]),
    (br.Y[2], br.Y[1]),
    !is_increasing(br),
)

@doc raw"""
Dynamic based on a piecewise monotonic map.

The map is defined as T(x) = Ts\[k\](x) if ``x \in [endpoints[k], endpoints[k+1])``.

`y_endpoints` (``k \times 2`` matrix) contains the result of applying Ts to the endpoints of each interval. These can be filled in automatically from `endpoints`,
but sometimes they are known to higher accuracy, for instance for `x -> mod(3x, 1)` we know that it is full-branch exactly.
It is assumed that the map will send its domain hull(endpoints[begin],endpoints[end]) into itself.

the array `branches` is guaranteed to satisfy branches[i].X[end]==branches[i+1].X[begin]
"""
struct PwMap <: Dynamic
    branches::Array{MonotonicBranch,1}
    full_branch::Bool
    function PwMap(branches::Array{<:MonotonicBranch,1}; full_branch = false)
        new(branches, full_branch)
    end
end

function PwMap(
    Ts,
    endpoints,
    y_endpoints_in = hcat(
        [Ts[k](Interval(endpoints[k])) for k = 1:length(Ts)],
        [Ts[k](Interval(endpoints[k+1])) for k = 1:length(Ts)],
    );
    full_branch = false,
)
    branches = MonotonicBranch[]
    for k = 1:length(endpoints)-1
        y_endpoints = (y_endpoints_in[k, 1], y_endpoints_in[k, 2])
        push!(branches, MonotonicBranch(Ts[k], (endpoints[k], endpoints[k+1]), y_endpoints))
    end
    X = (endpoints[begin], endpoints[end])
    full_branch_detected = full_branch || all(is_full_branch(b, X) for b in branches)
    return PwMap(branches; full_branch = full_branch_detected)
end

function intersect_domain(D::PwMap, x)
    return [(Interval(x) ∩ hull(br.X[1], br.X[2])) for br in D.branches]
end

function intersect_domain_bool(D::PwMap, x)
    return [y != ∅ for y in intersect_domain(D, x)]
end

Base.show(io::IO, D::PwMap) =
    print(io, "Piecewise-defined dynamic with $(nbranches(D)) branches")

domain(D::PwMap) = (D.branches[1].X[1], D.branches[end].X[2])
#intersect_domain(T::PwMap, x) = [Interval(x) ∩ hull(br.X[1], br.X[2]) for br in T.branches]
#intersect_domain_bool(T::PwMap, x) = [!(∅ == y) for y in intersect_domain(T, x)]


function is_increasing(D::PwMap)
    inc = is_increasing.(D.branches)
    if all(inc)
        return true
    elseif all(.!inc)
        return false
    else
        error("The given dynamic has branches with different orientations")
    end
end

Base.getindex(D::PwMap, k::Int64) = D.branches[k]

nbranches(D::PwMap) = length(D.branches)
endpoints(D::PwMap) = [[br.X[1] for br in branches(D)]; branches(D)[end].X[2]]
branches(D::PwMap) = D.branches

is_full_branch(D::PwMap) = D.full_branch

"""
Deprecated, but still used in C2Basis
"""
function preim(D::PwMap, k, y, ϵ = 1e-15)
    return preimage(y, D.branches[k]; ϵ = ϵ, max_iter = 100)
end

"""
Intersect an Interval or TaylorSeries with I
"""
restrict(I, x) = I ∩ x
restrict(I, x::Taylor1) = Taylor1([I ∩ x[0]; x[1:end]], x.order)

"""
function that evaluates the k-th branch of a dynamic on a point x
	(assuming it's in its domain, otherwise ∅)
"""
function branch(D::PwMap, k)
    return x -> D[k].f(restrict(hull(D[k].X[1], D[k].X[2]), x))
end

# Unused as of now
# """
# hull of an iterable of intervals
# """
# common_hull(S) = interval(minimum(x.lo for x in S), maximum(x.hi for x in S))


# Rather than defining derivatives of a PwMap, we define Taylor1 expansions directly
@doc raw"""
Function call, and Taylor expansion, of a PwMap.
Note that this ignores discontinuities; users are free to shoot themselves
in the foot and call this on a non-smooth piecewise map. No better solutions for now.
"""
function (D::PwMap)(x::Taylor1)
    fx = fill(∅, x.order + 1)
    x_restricted = deepcopy(x)
    for i = 1:length(D.branches)
        x_restricted[0] = x[0] ∩ hull(D[i].X[1], D[i].X[2])
        if !isempty(x_restricted[0])
            fx_restricted = D[i].f(x_restricted)
            fx = fx .∪ fx_restricted.coeffs
        end
    end
    @debug "Piecewise f($(x)) = $(Taylor1(fx, x.order))"
    return Taylor1(fx, x.order)
end

"""
    branch_inverse_derivative(br::MonotonicBranch; tol = 0.01)

Compute a rigorous bound for the inverse_derivative of a branch
"""
function branch_inverse_derivative(br::MonotonicBranch, tol = 0.01)
    I = hull(br.X[1], br.X[2])
    h = x -> abs(inverse_derivative(br.f, x))
    val, listofboxes = maximise(h, I, tol = tol)
    @debug val, listofboxes
    return val
end

"""
    inverse_derivative(D::PwMap; tol=1e-3)

Compute a rigorous bound for the inverse_derivative of a PwMap
"""
function max_inverse_derivative(D::PwMap, tol = 1e-3)
    #for br in branches(D)
    #    val = maximise(x -> abs(1/derivative(br.f, x)), hull(br.X[1], br.X[2]), tol=tol)[1]
    #        @info val
    #    max_exp = max(val, max_exp)
    #end
    return maximum([branch_inverse_derivative(br, tol) for br in D.branches])
end

"""
    bound_branch_distortion(br::MonotonicBranch; tol = 0.01)

Compute a rigorous bound for the distortion of a branch
on an interval I, defaults to the domain of the branch
"""
function bound_branch_distortion(
    br::MonotonicBranch,
    I = hull(br.X[1], br.X[2]);
    tol = 0.01,
)
    h = x -> abs(distortion(br.f, x))
    val, listofboxes = maximise(h, I, tol = tol)
    @debug val, listofboxes
    return val
end

@doc raw"""
    max_distortion(D::PwMap; tol=1e-3)

Compute a rigorous bound for the distortion of a PwMap

# Example

```jldoctest
julia> using RigorousInvariantMeasures;

julia> D0 = mod1_dynamic(x->2*x+0.5*x*(1-x), full_branch = true)
Piecewise-defined dynamic with 2 branches

julia> max_distortion(D0)
[0.444268, 0.444445]
```
"""
function max_distortion(D::PwMap, tol = 1e-3)
    # max_dist = Interval(0.0)
    # for br in branches(D)
    #     val = maximise(x -> abs(distortion(br.f, x)), hull(br.X[1], br.X[2]), tol=tol)[1]
    #     max_dist = max(val, max_dist)
    # end
    # return max_dist
    return maximum([bound_branch_distortion(br; tol = tol) for br in D.branches])
end

function plottable(D::PwMap, x)
    @assert 0 <= x <= 1
    for k = 1:nbranches(D)
        domain = hull(D[k].X[1], D[k].X[2])
        if x in domain
            return mid(Interval(D[k].f(x)))
        end
    end
end
plottable(D::PwMap) = x -> plottable(D, x)

"""
    has_infinite_derivative_at_endpoints(b::MonotonicBranch)

Returns a pair (left::Bool, right::Bool) that tells if a branch has infinite derivative at any of its endpoints
"""
function has_infinite_derivative_at_endpoints(branch::MonotonicBranch)
    # the Interval(0, 1e-15) summand is there because for some reason TaylorSeries fails on point intervals of singularity but not on larger intervals containing them:
    # derivative(x->x^(6/10), 0..0) # fails
    # derivative(x->x^(6/10), 0..1e-15) # succeeds

    left = !isfinite(derivative(branch.f, branch.X[1] + Interval(0, 1e-15)))
    right = !isfinite(derivative(branch.f, branch.X[2] - Interval(0, 1e-15)))
    return (left, right)
end

"""
has_infinite_derivative_at_endpoints(D::PwMap)

Returns a single bool to tell whether the dynamic has infinite derivative at any of its endpoint
"""
function has_infinite_derivative_at_endpoints(D::PwMap)
    return any(any(has_infinite_derivative_at_endpoints(b)) for b in D.branches)
end

@doc raw"""
    mod1_dynamic(f::Function, ε = 0.0; full_branch = false)

Utility constructor for dynamics Mod 1 on the torus [0,1].
We assume that f is monotonic and differentiable, for now (this is not restrictive, for our purposes)

# Example

```jldoctest
julia> using RigorousInvariantMeasures;

julia> D0 = mod1_dynamic(x->2*x+0.5*x*(1-x), full_branch = true)
Piecewise-defined dynamic with 2 branches
```
"""
function mod1_dynamic(f::Function; ϵ = 0.0, max_iter = 100, full_branch = false)
    X = (0 .. 0, 1 .. 1)
    br = MonotonicBranch(f, X)
    @debug "Auxiliary branch" br

    # check monotonicity

    try
        @assert minimise(
            x -> derivative(f, x) * (is_increasing(br) ? 1 : -1),
            hull(Interval.(X)...),
        )[1] > 0
    catch exc
        if isa(exc, MethodError)
            @warn "It probably is complaining that the result of derivative(f) is not an interval"
            rethrow(exc)
        else
            rethrow(exc)
        end
    end

    Yhull = hull(Interval.(br.Y)...)
    possible_integer_parts = floor(Int, Yhull.lo):ceil(Int, Yhull.hi)
    @debug "Possible integer parts" possible_integer_parts

    x, integer_parts =
        preimages(possible_integer_parts, br, possible_integer_parts; ϵ, max_iter)

    ep = [x; X[end]]
    Ts = [x -> f(x) - k for k in integer_parts]

    n = Base.length(x)
    if is_increasing(br)
        y_endpoints::Matrix{Interval{Float64}} = hcat(fill(0.0, n), fill(1.0, n))
    else
        y_endpoints = hcat(fill(1.0, n), fill(0.0, n))
    end
    y_endpoints[1, 1] = br.Y[begin] - integer_parts[begin]
    if y_endpoints[1, 1] == 0.0
        y_endpoints[1, 1] = 0.0 # hack to get rid of -0..0 intervals
    end
    y_endpoints[end, end] = br.Y[end] - integer_parts[end]
    if y_endpoints[end, end] == 0.0
        y_endpoints[end, end] = 0.0 # hack to get rid of -0..0 intervals
    end
    # not needed, since the check is moved into the PwMap() constructor
    # full_branch_detected = full_branch || all(equal_up_to_order(X, y_endpoints[i,:]) for i in 1:n)

    return PwMap(Ts, ep, y_endpoints; full_branch = full_branch)
end

"""
    Create explicitly D1 ∘ D2 as a PwMap; remark that the endpoints 
    of D1 must be ordered with respect to the canonical order on R.
"""
function composedPwMap(D1::PwMap, D2::PwMap)
    y_endpoints = endpoints(D1)
    x, xlabel = preimages_and_branches(y_endpoints, D2; ϵ = 1e-13, max_iter = 100)
    push!(x, D2.branches[end].X[2])
    branches = [
        MonotonicBranch(
            D1.branches[brid1].f ∘ D2.branches[brid2].f,
            (x[i], x[i+1]),
            D1.branches[brid1].Y,
        ) for (i, (brid1, brid2)) in enumerate(xlabel)
    ]

    # full_branch is defined conservatively; it is possible that non-full-branch dynamics
    # composed could yield a full-branch composed dynamic
    full_branch = is_full_branch(D1) && is_full_branch(D2)

    return PwMap(branches; full_branch = full_branch)
end

using ProgressMeter

@doc raw"""
dfly inequality for maps with infinite derivatives. 
    
The strategy to compute it follows a variant of Lemma 9.1 in the GMNP paper: 
* we find a "problematic set" I by taking a small interval of size radius(branch domain)/2^3 around each endpoint with infinite derivative; 
* we find l such that T >= l for each point in I
* we compute the dfly coefficients as in the lemma.
* we repeat the computation replacing 2^3 with 2^4, 2^5, ... 2^15 and take the best estimate among these.
"""
function dfly_inf_der(::Type{TotalVariation}, ::Type{L1}, D::PwMap, tol = 1e-3)
    leftrightsingularity = Tuple{Bool,Bool}[]
    A = +∞
    B = +∞

    # for each branch, we check if the derivative is infinite at any of the endpoints:
    for br in branches(D)
        (left, right) = has_infinite_derivative_at_endpoints(br)

        push!(leftrightsingularity, (left, right))
    end
    est = +∞
    #@showprogress enabled=SHOW_PROGRESS_BARS  1 "Computing infinite-derivative DFLY..." 
    for i = 3:15
        val = 0.0
        val_summand = Interval(0.0)
        l = 0.0

        for (j, br) in enumerate(branches(D))
            rad = radius(hull(br.X[1], br.X[2]))
            tol = rad / 2^(i + 1)
            left, right = leftrightsingularity[j]
            @debug "branch $j, left_singularity=$left, right_singularity=$right"
            f = inverse_derivative(br.f)
            g = distortion(br.f)
            left_endpoint = br.X[1]
            right_endpoint = br.X[2]
            if left
                left_endpoint += rad / 2^i
            end
            if right
                right_endpoint -= rad / 2^i
            end
            I = hull(left_endpoint, right_endpoint)
            val_br = 2 * maximise(x -> abs(f(x)), I, tol = tol)[1]
            @debug "maximise on $I: $val_br"
            val = max(val, val_br.hi)

            if left
                # we work on the interval [br.X[1], left_endpoint] =: [a, b].
                # we assume here that f(a) = 0, g(a)=∞, and that g is monotonically decreasing on [a,b]
                l_left = abs(g(left_endpoint))
                @debug "left endpoint: g($left_endpoint) = $l_left"
                # we use the fact that the primitive of the distortion g is f=1/T', and compute
                # int_a^b (distorsion) = f(b) - f(a) = f(b)
                val_summand += abs(f(left_endpoint) / 2)
                @debug "abs(f($left_endpoint)/2) = $val_summand"
                l = max(l, abs(l_left).hi)
            end
            if right
                # same reasoning as above but on the right endpoint
                l_right = abs(g(right_endpoint))
                @debug "right endpoint: g($right_endpoint) = $l_right"
                val_summand += abs(f(right_endpoint) / 2)
                @debug "abs(f($right_endpoint)/2) = $val_summand"
                l = max(l, abs(l_right).hi)
            end
        end
        val = val ⊕₊ val_summand.hi
        @debug "i=$i, A=$val, B=$l"

        if val < 1.0 && l ⊘₊ (1.0 ⊖₋ val) < est
            est = l ⊘₊ (1.0 ⊖₋ val)
            A = val
            B = l
            @debug "improving on previous best"
        end
    end
    endpts = endpoints(D)
    min_width = minimum([endpts[i+1] - endpts[i] for i = 1:length(endpts)-1])

    return A, B ⊕₊ (2 / min_width).hi
end
#    aux_der = (1-max_exp)/2
#  @info "aux_der", aux_der

#    for br in branches(D)
#        if is_increasing(br)
#            rr = preimage(aux_der, x -> 1/(2*derivative(br.f, x)), hull(br.X[1], br.X[2]), 10^-13)
#            @info rr
#        else
#            rr = preimage(-aux_der, x -> -1/(2*derivative(br.f, x)), hull(br.X[1], br.X[2]), 10^-13)
#            @info rr
#        end
#    end

#dfly_inf_der(D::ComposedDynamic, tol=1e-3) = dfly_inf_der(D.E, tol)
