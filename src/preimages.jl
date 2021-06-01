"""
Compute preimages of monotonic sequences
"""

using IntervalArithmetic
using .Contractors

"""
Type used to represent a "branch" of a dynamic. The branch is represented by a monotonic map `f` with domain `X=(a,b)` with a≤b (where typically a,b are intervals). 
`Y=(f(a),f(b))` and `increasing` may be provided (for instance if we know that `Y=(0,1)`), otherwise they are computed automatically.
"""
struct Branch{T,S}
    f::T
    X::Tuple{S, S}
    Y::Tuple{S, S}
    increasing::Bool
end
Branch(f, X, Y=(f(Interval(X[1])), f(Interval(X[2]))), increasing=unique_increasing(Y[1], Y[2])) = Branch{typeof(f), typeof(X[1])}(f, X, Y, increasing)

"""
Return Branches for a given dynamic, in an iterable
"""
function branches(D::PwMap) # TODO: should probably be rewritten at the interface level, if we go for this approach. D.is_full is not the correct thing for decreasing branches.
    domain = hull(Interval(D.endpoints[begin]), Interval(D.endpoints[end]))
    return [Branch(D.Ts[k], (D.endpoints[k], D.endpoints[k+1]), D.is_full[k] ? (Interval(0),Interval(1)) : (D.Ts[k](Interval(D.endpoints[k])) ∩ domain, D.Ts[k](Interval(D.endpoints[k+1])) ∩ domain), D.orientations[k]==1) for k in 1:length(D.Ts)]
end

"""
Smallest possible i such that a is in the semi-open interval [y[i], y[i+1]).

This should work properly even if `a, y` contain intervals.
Assumes y is sorted, i.e., map(y, x->Interval(x).lo) and map(y, x->Interval(x).hi) are sorted.
"""
function first_overlapping(y, a)
    searchsortedlast(y, Interval(a).lo, by=x->Interval(x).hi)
end

"""
Largest possible j such that a-ε is in the semi-open interval [y[j], y[j+1]).

This should work properly even if `a, y` contain intervals.
Assumes y is sorted, i.e., map(y, x->Interval(x).lo) and map(y, x->Interval(x).hi) are sorted.
"""
function last_overlapping(y, a)
    searchsortedfirst(y, Interval(a).hi, by=x->Interval(x).lo) - 1
end

"""
Construct preimages of an increasing array y under a monotonic branch defined on X = (a, b), propagating additional labels `ylabel`

In general, there may be a certain number of points in y that have no preimage at the beginning and the end of the sequence, because 
they fall out of the range R = [f(a), f(b)]. In the worst case, no point has a preimage, because y[i] < R < y[i+1] for some 
i (or vice versa with orientations).

The sequence y identifies semi-open intervals [y[l], y[l+1]); each of them is identified by label `ylabel[l]`. We construct a sequence x that splits X
into semi-open intervals, each of them with f([x[k], x[k+1]) ⊂ [y[l], y[l+1]) for a different l. We set xlabel[k] = ylabel[l] in this case, and return the pair (x, xlabel).

x[begin] always coincides with branch.X[1], while branch.X[2] is "the point after x[end]", and is not stored explicitly in x, for easier composing.

In this way x and xlabel have the same length.

This function fills the array by using a bisection strategy to save computations: if y ∈ [a,b], then f⁻¹(y) ∈ [f⁻¹(a),f⁻¹(b)] (paying attention to orientation).
So we can fill v by filling in first entries `v[k+1]` with higher dyadic valuation of k.

Currently this works only for 1-based 1-dimensional arrays y.
"""
function preimages(y, br::Branch, ylabel = 1:length(y), ϵ = 0.0)

    if br.increasing
        i = first_overlapping(y, br.Y[1])  # smallest possible i such that a is in the semi-open interval [y[i], y[i+1]).
        j = last_overlapping(y, br.Y[2]) # largest possible j such that b-ε is in the semi-open interval [y[j], y[j+1]).
        n = j - i + 1
        x = fill((-∞..∞)::typeof(Interval(br.X[1])), n)
        xlabel = collect(ylabel[i:j]) # we collect to avoid potential type instability, since this may be an UnitRange while in the other branch we could have a StepRange
        x[1] = br.X[1]
        if n == 1
            return (x, xlabel)
        end
        stride = prevpow(2, n-1)
        while stride >= 1
            # fill in v[i] using x[i-stride].lo and x[i+stride].hi as range for the preimage search
            for k = 1+stride:2*stride:n
                search_range = Interval(x[k-stride].lo, (k+stride <= n ? x[k+stride] : Interval(br.X[2])).hi)
                x[k] = preimage(y[i-1+k], br.f, search_range, ϵ)
            end
            stride = stride ÷ 2
        end
    else # branch decreasing
        i = last_overlapping(y, br.Y[1]) # largest possible j such that b-ε is in the semi-open interval [y[j], y[j+1]).
        j = first_overlapping(y, br.Y[2]) # smallest possible i such that a is in the semi-open interval [y[i], y[i+1]).
        n = i - j + 1
        x = fill((-∞..∞)::typeof(Interval(br.X[1])), n)
        xlabel = collect(ylabel[i:-1:j])
        x[1] = br.X[1]
        if n == 1
            return (x, xlabel)
        end
        stride = prevpow(2, n-1)
        while stride >= 1
            # fill in v[i] using x[i-stride].lo and x[i+stride].hi as range for the preimage search
            for k = 1+stride:2*stride:n
                search_range = Interval(x[k-stride].lo, (k+stride <= n ? x[k+stride] : Interval(br.X[2])).hi)
                x[k] = preimage(y[i+2-k], br.f, search_range, ϵ)
            end
            stride = stride ÷ 2
        end
    end
    return (x, xlabel)
end

function preimages(y, D::Dynamic, ylabel = 1:length(y), ϵ = 0.0)
    results = collect(preimages(y, b, ylabel, ϵ) for b in branches(D))
    x = vcat((result[1] for result in results)...)
    xlabel = vcat((result[2] for result in results)...)
    return x, xlabel
end

function preimages(z, D::ComposedFunction, zlabel = 1:length(z), ϵ = 0.0)
    y, ylabel = preimages(z, D.outer, zlabel, ϵ)
    return preimages(y, D.outer, ylabel, ϵ)
end

struct Dual{Ulam}
    x::Vector{Interval} #TODO: more generic type needed in future
    xlabel::Vector{Int}
    lastpoint::Interval
end

Dual(B, D, ϵ) = Dual{typeof(B)}(preimages(B.p, D, 1:length(B.p)-1, ϵ)..., endpoints(D)[end])

function iterate(dual::Dual, state = 1)
    n = length(dual.x)
    if state < n
        return (dual.xlabel[state], (dual.x[state], dual.x[state+1])), state+1
    elseif state == n
        return (dual.xlabel[n], (dual.x[n], dual.lastpoint)), state+1
    else
        return nothing
    end
end
Base.length(dual::Dual{<:Ulam}) = length(dual.x)
Base.eltype(dual::Dual{<:Ulam}) = Tuple{eltype(dual.xlabel), Tuple{eltype(dual.x), eltype(dual.x)}}

function assemble2(B::Basis, D::Dynamic, ϵ=0.0; T = Float64)
	I = Int64[]
	J = Int64[]
	nzvals = Interval{T}[]
	n = length(B)

	# TODO: reasonable size hint?

	for (i, dual_element) in Dual(B, D, ϵ)
		if !is_dual_element_empty(B, dual_element)
			for (j, x) in ProjectDualElement(B, dual_element)
				push!(I, i)
				push!(J, mod(j,1:n))
				push!(nzvals, x)
			end
		end
	end

	return sparse(I, J, nzvals, n, n)
end

function DiscretizedOperator2(B, D, ϵ=0.0; T = Float64)
	L = assemble2(B, D, ϵ; T)
	if is_integral_preserving(B)
		return IntegralPreservingDiscretizedOperator(L)
	else
		f = integral_covector(B)
		e = one_vector(B)
		w = f - f*L #will use interval arithmetic when L is an interval matrix
		return NonIntegralPreservingDiscretizedOperator(L, e, w)
	end
end
