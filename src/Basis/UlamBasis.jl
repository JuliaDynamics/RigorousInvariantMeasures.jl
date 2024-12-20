using IntervalArithmetic, LinearAlgebra

"""
	Ulam
Ulam basis on [0,1] associated to the partition
``p = \\{x_0 = 0, x_1, \\ldots, x_n=1\\}``
"""
struct Ulam{T<:AbstractVector} <: Basis
    p::T
    # TODO: check in constructor that p is sorted, starts with 0 and ends with 1
end

"""
	Ulam(n::Integer)
	
Equispaced Ulam basis on [0,1] of size n
"""
Ulam(n::Integer) = Ulam(LinRange(0.0, 1.0, n + 1))

@doc raw"""
	Base.length(B::Ulam)
Returns the size of the Ulam basis (the size of the underlying vector -1)
"""
Base.length(B::Ulam) = length(B.p) - 1

@doc raw"""
	Base.getindex(B::Ulam, i::Int)
Returns the i-th element of the Ulam basis as a function.

# Example

```jldoctest
julia> using RigorousInvariantMeasures

julia> B = Ulam(16)
Ulam{LinRange{Float64, Int64}}(LinRange{Float64}(0.0, 1.0, 17))

julia> B[1](1/32)
1

julia> B[2](1/32)
0
```

"""
function Base.getindex(B::Ulam, i::Int)
    return x -> (B.p[i] <= x < B.p[i+1] ? 1 : 0)
end

function is_dual_element_empty(::Ulam, d)
    return isempty(d[1]) || isempty(d[2])
end

#Base.length(S::DualComposedWithDynamic{<:Ulam, <:Dynamic}) = length(S.basis) * nbranches(S.dynamic)

# """
# Returns dual elements which are pairs (i, (a,b))
# i is an interval index, and (a,b) are the endpoints of its preimage
# """
# function Base.iterate(S::DualComposedWithDynamic{<:Ulam, <:Dynamic}, state = (1, 1, nothing))
# 	# a = preim(S.dynamic, k, S.basis.p[i], S.ϵ) is cached in the state (when i \neq 1)
# 	i, k, a = state

# 	if k == nbranches(S.dynamic)+1
# 		return nothing
# 	end

# 	if a == nothing
# 		@assert i==1
# 		a = preim(S.dynamic, k, S.basis.p[i], S.ϵ)
# 	end

# 	# a = preim(S.dynamic, k, S.basis.p[i], S.ϵ) # moved into state
# 	b = preim(S.dynamic, k, S.basis.p[i+1], S.ϵ)

# 	if isempty(a) && !isempty(b)
# 		ep = endpoints(S.dynamic)
# 		if orientation(S.dynamic, k) > 0
# 			a = convert(typeof(b), ep[k])
# 		else
# 			a = convert(typeof(b), ep[k+1])
# 		end
# 	elseif isempty(b) && !isempty(a)
# 		ep = endpoints(S.dynamic)
# 		if orientation(S.dynamic, k) > 0
# 			b = convert(typeof(a), ep[k+1])
# 		else
# 			b = convert(typeof(a), ep[k])
# 		end
# 	end

# 	if i == length(S.basis)
# 		return ((i, (a, b)), (1, k+1, nothing))
# 	else
# 		return ((i, (a, b)), (i+1, k, b))
# 	end
# end

# Base.eltype(f::DualComposedWithDynamic{<:Ulam, <:Dynamic}) = Tuple{Int64,Tuple{Interval{Float64},Interval{Float64}}}

"""
	nonzero_on(B::Ulam, (a, b))

Returns the indices of the elements of the Ulam basis that intersect with the interval y
We do not assume an order of a and b; this should not matter unless
the preimages are computed with very low accuracy.
We assume, though, that y comes from the (possibly inexact) numerical approximation
of an interval in ``[0,1]``, i.e., we restrict to ``y \\cap [0,1]``
"""
function nonzero_on(B::Ulam, (a, b))
    y = hull(a, b)

    # finds in which semi-open interval [p[k], p[k+1]) y.lo and y.hi fall
    lo = searchsortedlast(B.p, y.lo)
    hi = searchsortedlast(B.p, y.hi)

    # they may be n+1 if y.hi==1
    lo = clamp(lo, 1, length(B))
    hi = clamp(hi, 1, length(B))

    return (lo, hi)
end

"""
	relative_measure((a,b)::Tuple{<:Interval,<:Interval}, (c,d)::Tuple{<:Interval,<:Interval})

Relative measure of the intersection of (a,b) wrt the whole interval (c,d)
Assumes that a,b and c,d are sorted correctly
"""
function relative_measure(
    (a, b)::Tuple{<:Interval,<:Interval},
    (c, d)::Tuple{<:Interval,<:Interval},
)
    # this is a bit lazy because we could compute orientations and tell which is which
    # but it won't matter unless a,b are computed with very low precision
    a, b = min(a, b), max(a, b)

    lower = max(a, c)
    upper = min(b, d)
    intersection = max(upper - lower, 0) / (d - c)
    return intersection
end

"""
Given a preimage of an interval ```I_i```, this iterator returns
its relative intersection with all the elements of the Ulam basis that
have nonzero intersection with it
"""
function Base.iterate(S::ProjectDualElement{BT,DT}, state = S.j_min) where {BT<:Ulam,DT}
    if state == S.j_max + 1
        return nothing
    end
    j = state
    x = relative_measure(S.dual_element, (Interval(S.basis.p[j]), Interval(S.basis.p[j+1])))
    return (j, x), state + 1
end
Base.eltype(f::ProjectDualElement{<:Ulam,DT}) where {DT} = Tuple{Int64,Interval{Float64}}

evaluate(B::Ulam{T}, i, x) where {T} =
    (x > (i - 1) / length(B)) && (x < i / length(B)) ? 1 : 0

evaluate_integral(B::Ulam{S}, i, T::Type) where {S} = T(1) / length(B)

"""
	iterate(S::AverageZero{Ulam{T}}, state = 1) where{T}

Return the elements of the basis of average ``0`` functions in the Ulam 
Basis
"""
function Base.iterate(S::AverageZero{Ulam{T}}, state = 1) where {T}
    n = length(S.basis)
    if state == n
        return nothing
    end
    v = zeros(Float64, n)
    v[1] = 1
    v[state+1] = -1
    return (v, state + 1)
end

"""
	Base.length(S::AverageZero{Ulam})

Return the size of the Average Zero space
"""
Base.length(S::AverageZero{Ulam{T}}) where {T} = length(S.basis) - 1

is_refinement(Bf::Ulam, Bc::Ulam) = Bc.p ⊆ Bf.p

function integral_covector(B::Ulam{T}) where {T}
    n = length(B)
    return 1 / n * ones(Interval{Float64}, n)'
end

is_integral_preserving(B::Ulam{T}) where {T} = true

function one_vector(B::Ulam)
    return ones(length(B))
end

strong_norm(B::Ulam) = TotalVariation
weak_norm(B::Ulam) = L1
aux_norm(B::Ulam) = L1

# See BasisDefinition for docs on these constants
# These must be rounded up correctly!

weak_projection_error(B::Ulam) = 0.5 ⊘₊ Float64(length(B), RoundDown)
aux_normalized_projection_error(B::Ulam) = 0.0
strong_weak_bound(B::Ulam) = Float64(length(B), RoundUp)
aux_weak_bound(B::Ulam) = 1.0
weak_by_strong_and_aux_bound(B::Ulam) = (0.0, 1.0)
bound_weak_norm_from_linalg_norm(B::Ulam) = (1.0, 0.0)
bound_linalg_norm_L1_from_weak(B::Ulam) = 1.0
bound_linalg_norm_L∞_from_weak(B::Ulam) = Float64(length(B), RoundUp)
bound_weak_norm_abstract(B::Ulam, D = nothing; dfly_coefficients = nothing) = 1.0

opnormbound(B::Ulam{T}, N::Type{L1}, A::AbstractVecOrMat{S}) where {T,S} = opnormbound(N, A)
#opnormbound(B::Ulam{T}, N::Type{L1}, Q::IntegralPreservingDiscretizedOperator) where {T} = opnormbound(N, Q.L)
normbound(B::Ulam{T}, N::Type{L1}, v) where {T} = normbound(N, v)

function invariant_measure_strong_norm_bound(
    B::Ulam,
    D::Dynamic;
    dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D),
)
    A, B = dfly_coefficients
    @assert A < 1.0
    return B ⊘₊ (1.0 ⊖₋ A)
end

struct UlamDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    lastpoint::Interval
end
Dual(B::Ulam, D; ϵ, max_iter) =
    UlamDual(preimages(B.p, D, 1:length(B.p)-1; ϵ, max_iter)..., domain(D)[end])

Base.length(dual::UlamDual) = length(dual.x)
Base.eltype(dual::UlamDual) =
    Tuple{eltype(dual.xlabel),Tuple{eltype(dual.x),eltype(dual.x)}}
function Base.iterate(dual::UlamDual, state = 1)
    n = length(dual.x)
    if state < n
        return (dual.xlabel[state], (dual.x[state], dual.x[state+1])), state + 1
    elseif state == n
        return (dual.xlabel[n], (dual.x[n], dual.lastpoint)), state + 1
    else
        return nothing
    end
end
