"""
Compute preimages of monotonic sequences
"""

using IntervalArithmetic
using .Contractors


# Comparison operator

lt(a::Interval, b::Interval) = strictprecedes(a, b)
lt(a, b::Interval) = a < b.lo
lt(a::Interval, b) = a.hi < b
lt(a, b) = a < b

gt(a, b) = lt(b, a)

"""
Given a monotonic vector `a`, finds how many elements at its beginning are 
strictly smaller (if `a_increasing==true`) or strictly larger (if `a_increasing==false`) than `x`

Must work correctly also if `a` and `x` are intervals.
"""
function skip_beginning(a, x, a_increasing)
    cmp = ifelse(a_increasing, lt, gt)
    return searchsortedfirst(a, x, lt=cmp) - 1
end
"""
Companion to `skip_beginning`: finds the index of the last element of `a` that is *not* strictly larger (if `a_increasing==true`)
or strictly smaller (if `a_increasing==false`) than `x` (returns 0 if all elements are strictly below/above x)
"""
function last_end(a, x, a_increasing)
    cmp = ifelse(a_increasing, lt, gt)
    return searchsortedlast(a, x, lt=cmp)
end


"""
Construct preimages of a monotonic array y under a monotonic function f in a domain X.

In general, there may be a certain number of points in y that have no preimage at the beginning and the end of the sequence, because 
they fall out of the interval R = [f(X.lo), f(X.hi)]. In the worst case, no point has a preimage, because y[i] < R < y[i+1] for some 
i (or vice versa with orientations).

We return the pair `(v, skip)`, where `v` is a vector of preimages, and `skip` is the number of elements of y that have no preimage
at the *beginning* of the interval. So for instance if f and y are increasing but y[i] < R < y[i+1], we return skip=i and v = [].
The number of intervals that have no image at the *end* of the interval is then `length(y) - length(v) - skip`

Rationale: this is more efficient than returning a vector with a lot of empty intervals at the beginning/end, and callers 
need to know how many of the empty values come from the *beginning* of the array rather than its end.

Fills the array by using a bisection strategy to save computations: if y ∈ [a,b], then f⁻¹(y) ∈ [f⁻¹(a),f⁻¹(b)] (paying attention to orientation).
So we can fill v by filling in first entries `v[k+1]` with higher dyadic valuation of k.

Currently this works only for 1-based 1-dimensional arrays y.
"""
function preimages(y, f, X, ϵ = 0.0)
    fa = preimage(X.lo, f, X, ϵ)
    fb = preimage(X.hi, f, X, ϵ)
    f_increasing = unique_increasing(fa, fb) #TODO: we will want to compute these bools outside, I guess, to handle special cases, e.g., length(y)==1
    y_increasing = unique_increasing(y[begin], y[end])
    v_increasing = !(f_increasing ⊻ y_increasing) # unused, but maybe we'll want to return it as well?

    skip = skip_beginning(y, ifelse(f_increasing ⊻ y_increasing, fb, fa), y_increasing)
    last = last_end(y, ifelse(f_increasing ⊻ y_increasing, fa, fb), y_increasing)

    n = last - skip

    v = fill((-∞..∞)::typeof(X), n)
    v[1] = preimage(y[skip+1], f, X, ϵ)
    if n == 1
        return (v, skip)
    end
    v[end] = preimage(y[skip+n], f, X, ϵ)
    stride = prevpow(2, n-1)
    while stride > 1
        # fill in v[i] using v[i-stride] and v[i+stride]
        for i = 1+stride:2*stride:n-1
            # this hull() could be replaced with the proper [a.lo..b.hi], if we know orientations
            X = hull(v[i-stride], v[min(i+stride, n)])
            v[i] = preimage(y[skip+i], f, X, ϵ)
        end
        stride = stride ÷ 2
    end
    return (v, skip)
end
