"""
Compute preimages of monotonic sequences
"""

using IntervalArithmetic
using .Contractors
using .DynamicDefinition: is_increasing

## Moved the definition of MonotonicBranch in PwDynamicDefinition, with the objective
## of transforming PwMap into an Array of MonotonicBranch


"""
    first_overlapping(y, a)

Smallest possible i such that a is in the semi-open interval [y[i], y[i+1]).

This should work properly even if `a, y` are intervals; in this case it returns the *smallest* possible value of i over all possible "assignments" of a, y inside those intervals.
Assumes y is sorted, i.e., map(y, x->Interval(x).lo) and map(y, x->Interval(x).hi) are sorted.
"""
function first_overlapping(y, a)
    if iszero(a) # avoids -0 crap
        a = zero(a)
    end
    searchsortedlast(y, Interval(a).lo, by=x->Interval(x).hi)
end

"""
    last_overlapping(y, a)

Largest possible j such that a-ε is in the semi-open interval [y[j], y[j+1]).

This should work properly even if `a, y` are intervals; in this case it returns the *largest* possible value of i over all possible "assignments" of a, y inside those intervals.
Assumes y is sorted, i.e., map(y, x->Interval(x).lo) and map(y, x->Interval(x).hi) are sorted.
"""
function last_overlapping(y, a)
    if iszero(a) # avoids -0 crap
        a = zero(a)
    end
    searchsortedfirst(y, Interval(a).hi, by=x->Interval(x).lo) - 1
end

"""
Utility function that estimates the range of a monotone function
"""
range_estimate_monotone(f, X) = hull(f(Interval(X.lo)), f(Interval(X.hi)))

"""
Compute the preimage f⁻¹(y), knowing that it lies inside `search_interval`.
"""
function preimage(y, br::MonotonicBranch, search_interval=hull(br.X...); ϵ, max_iter)
    # Since the branch is monotonic, we can compute the preimages of y.lo and y.hi separately
    # This should give slightly thinner intervals.

    Y = intersect(range_estimate_monotone(br.f, search_interval), Interval(y))

    if isempty(Y)
        return Interval(∅)
    end

    xlo = preimage_monotonic(Y.lo, br.f, search_interval, br.Y; ϵ, max_iter)

    if isthin(Y)
        @debug "preimage of $y on $search_interval: $xlo"
        return xlo
    end

    xhi = preimage_monotonic(Y.hi, br.f, search_interval, br.Y; ϵ, max_iter)

    @debug "preimage of $y on $search_interval: $(hull(xlo, xhi))"
    return hull(xlo, xhi) ∩ search_interval
end


"""
    preimages(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter)

Construct preimages of a partition y under a monotonic branch defined on X = (a, b), propagating additional labels `ylabel`

Construct preimages of a partition `y` under a monotonic branch `br` defined on X = (a, b), propagating additional labels `ylabel`

It is assumed that it cannot happen that ``f(x) < y[1]``.

# Extended help
The sequence y subdivides the y-axis into semi-open intervals [y[l], y[l+1]); each of them is identified by the label `ylabel[l]`. We construct an increasing sequence 
x that splits X (in the x-axis) into semi-open intervals, each of them with f([x[k], x[k+1]) ⊂ [y[l], y[l+1]) for a certain l. 
We set xlabel[k] = ylabel[l], and return the pair (x, xlabel).

It is assumed that it cannot happen that f(x) < y[1].

In the simplest case where D is full-branch, the points in x are preimages of the points in y, but in the general case they can also include D.endpoints:
in general, there may be a certain number of points in y that have no preimage at the beginning and the end of the sequence, because 
they fall out of the range R = [f(a), f(b)]. In the worst case, no point has a preimage, because y[i] < R < y[i+1] for some 
i (or vice versa with orientations), and in this case we just return the 1-element vectors x = [branch.X[1]] and xlabel = [i].

x[begin] always coincides with branch.X[1], while branch.X[2] is "the point after x[end]", and is not stored explicitly in x, for easier composing.
In this way x and xlabel have the same length.

This function fills the array by using a bisection strategy to save computations: if y ∈ [a,b], then f⁻¹(y) ∈ [f⁻¹(a),f⁻¹(b)] (paying attention to orientation).
So we can fill v by filling in first entries `v[k+1]` with higher dyadic valuation of k.

For a dynamic with multiple branches, preimages(y, D) is simply the concatenation of x, xlabel for b in all branches. These values still form an increasing sequence that
splits X into intervals, each of which is mapped into a different semi-open interval [y[k], y[k+1]).

# Example 
```jldoctest
julia> using RigorousInvariantMeasures

julia> D = mod1_dynamic(x->2x)
Piecewise-defined dynamic with 2 branches

julia> RigorousInvariantMeasures.preimages(0:0.1:1, D.branches[1]; ϵ = 10^(-15), max_iter = 100)
(Interval{Float64}[[0, 0], [0.05, 0.0500001], [0.1, 0.100001], [0.149999, 0.15], [0.2, 0.200001], [0.25, 0.25], [0.299999, 0.3], [0.349999, 0.35], [0.4, 0.400001], [0.45, 0.450001]], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
```

"""
function preimages(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter)
    # TODO: there is some efficiency to be gained if we keep track of the value of f on each computed point
    # so that we can pass pairs (x,f(x)) to the underlying preimage() functions. Currently it is recomputed many times,
    # for instance in range_estimate_monotone

    # TODO: consider separate types (maybe via a parametric type) for increasing and decreasing branches.

    if is_increasing(br)
        i = first_overlapping(y, br.Y[1])  # smallest possible i such that a = br.Y[1] is in the semi-open interval [y[i], y[i+1]).
        j = last_overlapping(y, br.Y[2]) # largest possible j such that b-ε, where b = br.Y[2] is in the semi-open interval [y[j], y[j+1]).
        n = j - i + 1
        x = fill((-∞..∞)::typeof(br.X[1]), n)
        xlabel = collect(ylabel[i:j]) # we collect to avoid potential type instability, since this may be an UnitRange while in the other branch we could have a StepRange
        x[1] = br.X[1]
        if n == 1
            return (x, xlabel)
        end
        # the bisection strategy: fill the array in "strides" of length `stride`, then halve the stride and repeat
        # for instance, if the array is 1..13 (with x[1] filled in already), we first take stride=8 and fill in x[9],
        # then stride=4 and fill in x[5], x[13], (note that the distance is `2stride`, since x[9], and in general all the even multiples of `stride`, is already filled in)
        # then stride=2 and fill in x[3], x[7], x[11],
        # then stride=1 and fill in x[2], x[4], x[6], x[8], x[10], x[12]
        # at each step we have bracketed the preimage in a "search range" given by already-computed preimages x[k-stride] and x[k+stride].
        stride = prevpow(2, n-1)
        while stride >= 1
            # fill in v[i] using x[i-stride].lo and x[i+stride].hi as range for the preimage search
            for k = 1+stride:2*stride:n
                search_range = Interval(x[k-stride].lo, (k+stride <= n ? x[k+stride] : br.X[2]).hi)
                x[k] = preimage(y[i-1+k], br, search_range; ϵ, max_iter)
            end
            stride = stride ÷ 2
        end
        return x, xlabel
    else # decreasing branch
        # there is some asymmetry in how the Y axis is treated, since we assume that the function never goes below y[1] but can go above y[end].
        # hence to treat this case it is better to "reverse" the X axis: we switch to [br.X[1], x[1], ... x[d]], with an implied br.X[2] at the end 
        # as the last endpoint, to [-br.X[2], -x[d], -x[d-1], ..., -x[1]], with an implied -br.X[1] at the end.

        revx, revxlabel = preimages(y, reverse(br), ylabel; ϵ, max_iter)
        x = .-reverse(revx[2:end])
        insert!(x, 1, br.X[1])
        xlabel = collect(reverse(revxlabel))
        return x, xlabel
    end
end

"""
    preimages(y, D::Dynamic, ylabel = 1:length(y), ϵ = 0.0; progress = true)

    Construct preimages of an increasing array y under a dynamic, propagating additional labels `ylabel`
"""
function preimages(y, D::Dynamic, ylabel = 1:length(y); ϵ, max_iter)
    results = @showprogress 1 "Computing preimages..." [preimages(y, b, ylabel; ϵ, max_iter) for b in branches(D)]
    x = reduce(vcat, result[1] for result in results)
    xlabel = reduce(vcat, result[2] for result in results)
    return x, xlabel
end

"""
    preimages_and_branches(y, D::Dynamic, ylabel = 1:length(y); ϵ, max_iter)

Works like preimages(y, D) but the second return argument is a vector of pairs (ylabel, branchnumber) that specifies
which branch of D is the source of each interval in the preimages.
"""
# TODO: there is probably a very clever general construction that can reduce `preimages`, `preimages_and_branches`, 
# `preimages_and_derivatives` to three instances of the same method that keeps track of an "abstract" quantity over preimages.
function preimages_and_branches(y, D::PwMap, ylabel = 1:length(y); ϵ, max_iter)
    results = [preimages(y, b, ylabel; ϵ, max_iter) for b in D.branches]
    x = reduce(vcat, result[1] for result in results)
    xlabel = reduce(vcat, [(yid, brid) for yid in result[2]] for (brid, result) in enumerate(results))
    return x, xlabel
end

"""
    preimages_and_derivatives(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, maxiter, left=true)

Compute preimages of D *and* the derivatives f'(x) in each point.

Returns: x, xlabel, x′

We combine the two computations in the same function because in this way it can be implemented more efficiently for composed dynamics.

This assumes that (1) the dynamic is full-branch, and (2) all branches have the same orientation.
This is not restrictive because we'll need it only for the Hat assembler (at the moment).

If `left==true`, `x′` contains the derivatives `f'.(x)`. If right==true, `x′` contains derivatives in the right endpoint of each interval
of the partition, i.e., for each branch, `[f'(x[2]), f'(x[3]), ... f'(x[end]), f'(br.X[2])]`. In any case, `length(x) == length(x′)`.
"""
function preimages_and_derivatives(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter, left=true)
    x, xlabel = preimages(y, br, ylabel; ϵ, max_iter)
    f′ = derivative(br.f)
    x′ = similar(x)
    if left
        x′[1:end] = f′.(x)
    else
        x′ = similar(x)
        x′[1:end-1] = f′.(x[2:end])
        x′[end] = f′(br.X[2])
    end
    return x, xlabel, x′
end
function preimages_and_derivatives(y, D::Dynamic, ylabel = 1:length(y); ϵ, max_iter, left=true)
    @assert is_full_branch(D)
    results = @showprogress 1 "Computing preimages and derivatives..." [preimages_and_derivatives(y, b, ylabel; ϵ, max_iter, left) for b in branches(D)]
    x = reduce(vcat, result[1] for result in results)
    xlabel = reduce(vcat, result[2] for result in results)
    x′ = reduce(vcat, result[3] for result in results)
    return x, xlabel, x′
end

#PwOrComposed = Union{PwMap, ComposedDynamic}
"""
Composed map D1 ∘ D2 ∘ D3. We store with [D1, D2, D3] in this order.

We overwrite ∘ in base, so one can simply write D1 ∘ D2 or ∘(D1, D2, D3) to construct them.
"""
struct ComposedDynamic <: Dynamic
    dyns::Tuple{S, T} where {S, T <: Dynamic}
    E::PwMap
end
Base.:∘(D1::PwMap, D2::PwMap) = ComposedDynamic((D1, D2), composedPwMap(D1, D2))
Base.:∘(D1::PwMap, D2::ComposedDynamic) = ComposedDynamic((D1, D2), composedPwMap(D1, D2.E))
Base.:∘(D1::ComposedDynamic, D2::PwMap) = ComposedDynamic((D1, D2), composedPwMap(D1.E, D2))
Base.:∘(D1::ComposedDynamic, D2::ComposedDynamic) = ComposedDynamic((D1, D2), composedPwMap(D1.E, D2.E))
dfly(N1::Type{<:NormKind}, N2::Type{<:NormKind}, D::ComposedDynamic) = dfly(N1, N2, D.E)
dfly(N1::Type{RigorousInvariantMeasures.TotalVariation}, N2::Type{RigorousInvariantMeasures.L1}, D::ComposedDynamic) = dfly(N1, N2, D.E)
dfly(N1::Type{RigorousInvariantMeasures.Lipschitz}, N2::Type{RigorousInvariantMeasures.L1}, D::ComposedDynamic) = dfly(N1, N2, D.E)


DynamicDefinition.domain(D::ComposedDynamic) = DynamicDefinition.domain(D.dyns[end])
DynamicDefinition.is_increasing(D::ComposedDynamic) = iseven(sum(.!is_increasing.(D.dyns)))

function preimages(z, Ds::ComposedDynamic, zlabel = 1:length(z); ϵ, max_iter)
    @assert length(Ds.dyns)==2
    f = Ds.dyns[1]
    g = Ds.dyns[2]
    y, ylabel = preimages(z, f, zlabel; ϵ, max_iter)
    x, xlabel = preimages(y, g, ylabel; ϵ, max_iter)
    return x, xlabel
end

function preimages_and_derivatives(z, Ds::ComposedDynamic, zlabel = 1:length(z); ϵ, max_iter, left=true)
    @assert length(Ds.dyns)==2
    f = Ds.dyns[1]
    g = Ds.dyns[2]

    @assert is_full_branch(f) && is_full_branch(g)

    if DynamicDefinition.is_increasing(g)
        f_left = left
    else
        f_left = !left
    end

    y, ylabel, dy = preimages_and_derivatives(z, f, zlabel; ϵ, max_iter, left=f_left)
    x, yindex, dx = preimages_and_derivatives(y, g, 1:length(y); ϵ, max_iter, left)
    xlabel = ylabel[yindex]
    dx = dx .* dy[yindex]
    
    return x, xlabel, dx
end

function (D::ComposedDynamic)(x::Taylor1)
    for f in reverse(D.dyns)
        x = f(x)
    end
    return x
end

function DynamicDefinition.endpoints(D::ComposedDynamic)
    v = endpoints(D.dyns[1])
    auxDyn = D.dyns[2]
    x, xlabel = preimages(v, auxDyn; ϵ = 10^-14, max_iter = 100)
    x = [x; domain(auxDyn)[2]]
    return sort!(x)
end

DynamicDefinition.nbranches(D::ComposedDynamic) = length(endpoints(D))-1


# We need a better way to explicit this, at the moment we suppose everything 
# is full branch
DynamicDefinition.is_full_branch(D::ComposedDynamic) = all([is_full_branch(D.dyns[1]);is_full_branch(D.dyns[2])])

## Moved the definition of the abstract type Dual to BasisDefinition.jl

## Moved UlamDual to UlamBasis.jl

## Moved HatDual to HatBasis.jl

## Moved the new assembler to GenericAssembler, renamed the old assembler and DiscretizedOperator
## to _legacy
