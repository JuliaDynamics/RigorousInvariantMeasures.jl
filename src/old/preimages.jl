"""
Compute preimages of monotonic sequences
"""

using IntervalArithmetic
using .Contractors


"""
Strict-precedes on intervals, extended to non-intervals
"""
lt(a::Interval, b::Interval) = strictprecedes(a, b)
lt(a, b::Interval) = a < b.lo
lt(a::Interval, b) = a.hi < b
lt(a, b) = a < b

gt(a, b) = lt(b, a)

"""
Computes the pair (skip, last), where `skip` is the number of elements
at the beginning of `seq.v` that do not intersect with `X`, and `last` is the
last element after `skip` to intersect with `X`.

If `X` falls entirely between `v[i]` and `v[i+1]`, then (i, i) is returned.
"""
function skipandlast(seq, X)
    cmp = ifelse(seq.increasing, lt, gt)
    return (searchsortedfirst(seq.v, X, lt = cmp) - 1, searchsortedlast(seq.v, X, lt = cmp))
end

"""
Type used to store monotonic sequences of preimages. `skip` is a certain number of initial elements that are skipped 
with respect to an original reference array: for instance, preimages of a certain vector `y = [y(1), y(2), y(3), y(4), y(5)]` 
may only contain f^{-1}(y(3)) and f^{-1}(y(4)), so we set skip=2 and construct a v of length 2.
`increasing` tells if the monotonic sequence `v` is increasing (true) or decreasing (false)
"""
struct PointSequence{T<:AbstractVector}
    v::T
    skip::Int
    increasing::Bool
end
PointSequence(v, skip = 0, increasing = unique_increasing(v[begin], v[end])) =
    PointSequence{typeof(v)}(v, skip, increasing)

"""
Type used to represent a "branch" of a dynamic. The branch is represented by a monotonic map `f` with domain `X=(a,b)` with a≤b (where typically a,b are intervals). 
`Y=(f(a),f(b))` and `increasing` may be provided (for instance if we know that `Y=(0,1)`), otherwise they are computed automatically.
"""
struct MonotonicBranch{T,S}
    f::T
    X::Tuple{S,S}
    Y::Tuple{S,S}
    increasing::Bool
end
MonotonicBranch(
    f,
    X,
    Y = (f(Interval(X[1])), f(Interval(X[2]))),
    increasing = unique_increasing(Y[1], Y[2]),
) = MonotonicBranch{typeof(f),typeof(X[1])}(f, X, Y, increasing)

"""
Construct preimages of a monotonic array y under a monotonic function f in a domain X.

In general, there may be a certain number of points in y that have no preimage at the beginning and the end of the sequence, because 
they fall out of the interval R = [f(X.lo), f(X.hi)]. In the worst case, no point has a preimage, because y[i] < R < y[i+1] for some 
i (or vice versa with orientations).

We return a PointSequence where `v` is a vector of preimages, and `skip` is the number of elements of y that have no preimage
at the *beginning* of the interval. So for instance if f and y are increasing but y[i] < R < y[i+1], we return skip=i and v = [].
The number of intervals that have no image at the *end* of the interval is then `length(y) - length(v) - skip`

Rationale: this is more efficient than returning a vector with a lot of empty intervals at the beginning/end, and callers 
need to know how many of the empty values come from the *beginning* of the array rather than its end.

Fills the array by using a bisection strategy to save computations: if y ∈ [a,b], then f⁻¹(y) ∈ [f⁻¹(a),f⁻¹(b)] (paying attention to orientation).
So we can fill v by filling in first entries `v[k+1]` with higher dyadic valuation of k.

Currently this works only for 1-based 1-dimensional arrays y.

TODO: Idea: replace this with a function that returns an *increasing* vector v that is guaranteed to start/end with v = (a,...,b), and a sequence of indices `locs` like k:m or m:-1:k
such that v[i], v[i+1] is inside `locs[i]` for all `i`. This constructs "slices" more explicitly, and composes much better. Unclear: how to keep 'full-branched' info in this framework?
"""
function preimages(seq, branch, ϵ = 0.0)
    v_increasing = !(seq.increasing ⊻ branch.increasing)
    (skip, last) = skipandlast(seq, hull(branch.Y...))

    n = last - skip

    v = fill((-∞ .. ∞)::typeof(Interval(branch.X[1])), n)
    if n == 0
        return PointSequence(v, seq.skip + skip, v_increasing)
    end
    v[1] = preimage(seq.v[skip+1], branch.f, hull(branch.X...), ϵ)
    if n == 1
        return PointSequence(v, seq.skip + skip, v_increasing)
    end
    v[end] = preimage(seq.v[skip+n], branch.f, hull(branch.X...), ϵ)
    stride = prevpow(2, n - 1)
    while stride >= 1
        # fill in v[i] using v[i-stride] and v[i+stride]
        for i = 1+stride:2*stride:n-1
            X = hull(v[i-stride], v[min(i + stride, n)]) #TODO: this hull() could be replaced with the proper [a.lo..b.hi], since we know orientations
            v[i] = preimage(seq.v[skip+i], branch.f, X, ϵ)
        end
        stride = stride ÷ 2
    end
    return PointSequence(v, seq.skip + skip, v_increasing)
end
