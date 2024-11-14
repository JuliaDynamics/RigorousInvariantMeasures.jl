"""
Hat basis on the Torus [0,1]
"""

struct HatNP{T<:AbstractVector} <: Basis
    p::T
    # TODO: check in constructor that p is sorted, starts with 0 and ends with 1
end
HatNP(n::Integer) = HatNP(LinRange(0.0, 1.0, n + 1))

"""
	Base.length(B::HatNP)

Return the size of the HatNP basis
"""
Base.length(B::HatNP{T}) where {T} = Base.length(B.p)


"""
Hat function (on the reals)

This is a piecewise linear function such that:
	f(x) = 0 if x <= lo
	f(mi) = 1
	f(hi)
"""
struct HatFunction{T<:Number}
    lo::T
    mi::T
    hi::T
    function HatFunction{T}(lo, mi, hi) where {T<:Number}
        @assert lo <= mi <= hi
        new{T}(lo, mi, hi)
    end
end
HatFunction(lo::T, mi::T, hi::T) where {T<:Number} = HatFunction{T}(lo, mi, hi);
HatFunction(lo::Number, mi::Number, hi::Number) = HatFunction(promote(lo, mi, hi)...)

"""
Evaluate a HatFunction (on the real line).

Must work correctly when S is an interval.
"""
function (f::HatFunction{T})(x::S) where {T} where {S<:Number}
    lo = convert(S, f.lo)
    mi = convert(S, f.mi)
    hi = convert(S, f.hi)

    left_branch = (x - lo) / (mi - lo)
    right_branch = (hi - x) / (hi - mi)
    # 1 is not necessary in practice, but it avoids spurious results from rounding
    return max(min(left_branch, right_branch, 1), 0)
end

"""
	Base.getindex(B::HatNP, i::Int)

makes so that B[j] returns a HatFunction with the j-th basis element
"""
function Base.getindex(B::HatNP, i::Int)
    n = length(B) - 1
    @boundscheck 1 <= i <= n + 1 || throw(BoundsError(B, i))
    if i == 1
        x_minus_one = -1.0 / n
        return HatFunction(x_minus_one, B.p[1], B.p[2])
    elseif i == n + 1
        x_plus_one = 1.0 + 1.0 / n
        return HatFunction(B.p[end-1], B.p[end], x_plus_one)
    end
    return HatFunction(B.p[i-1], B.p[i], B.p[i+1])
end


function is_dual_element_empty(::HatNP, d)
    # TODO: the preim() may indeed be empty, so there could be an additional check here
    return false
end

is_refinement(Bf::HatNP, Bc::HatNP) = Bc.p ⊆ Bf.p

function integral_covector(B::HatNP)
    n = length(B)
    return 1 / n * ones(Interval{Float64}, n)'
end

function one_vector(B::HatNP)
    return ones(length(B))
end


"""
Return the range of indices of the elements of the basis whose support intersects
with the given dual element (i.e., a pair (y, absT')).
The range may end with length(B)+1; this must be interpreted "mod length(B)":
it means that it intersects with the hat function peaked in 0 as well
(think for instance y = 0.9999).
"""
function nonzero_on(B::HatNP, dual_element)
    #@info typeof(dual_element)
    y, absT′ = dual_element
    # Note that this cannot rely on arithmetic unless it is verified

    y = y ∩ Interval(0.0, 1.0) # we assume it's bona-fide interval in [0,1]
    # this should work for preims(), since they are supposed to return
    # a number in [0,1]

    # finds in which semi-open interval [p[k], p[k+1]) y.lo and y.hi fall
    lo = searchsortedlast(B.p, y.lo)
    hi = searchsortedlast(B.p, y.hi)
    lo = min(lo, length(B)) # lo may be n+1 if y.lo==1
    hi = min(hi, length(B)) # hi may be n+1 if y.hi==1
    hi = min(hi + 1, length(B)) # because the hat centered in p[k] is also nonzero in the interval before

    if lo == 1 # 1:N+1 does not make sense and would mean that the first interval is counted twice
        hi = min(hi, length(B))
    end
    return (lo, hi)
end

"""
Given a preimage ```y``` of a point ```x```, this iterator returns
```\\phi_j(y)/T'(y) ```
"""
function Base.iterate(S::ProjectDualElement{T,DT}, state = S.j_min) where {T<:HatNP,DT}
    if state == S.j_max + 1
        return nothing
    end
    y, absT′ = S.dual_element
    j = state
    n = length(S.basis)

    return ((j, S.basis[mod(j, 1:n)](y) / absT′), state + 1)
end

strong_norm(B::HatNP) = Lipschitz
weak_norm(B::HatNP) = Linf
aux_norm(B::HatNP) = L1

evaluate_integral(B::HatNP, i, T) = T(1) / (length(B) - 1)

function Base.iterate(S::AverageZero{HatNP{T}}, state = 1) where {T}
    n = length(S.basis)
    if state == n
        return nothing
    end
    v = zeros(Float64, n)
    v[1] = 1
    v[state+1] = -1
    return (v, state + 1)
end

Base.length(S::AverageZero{HatNP{T}}) where {T} = length(S.basis) - 1

weak_projection_error(B::HatNP) = 0.5 ⊘₊ Float64(length(B) - 1, RoundDown)
aux_normalized_projection_error(B::HatNP) = 0.5 ⊘₊ Float64(length(B) - 1, RoundDown)
strong_weak_bound(B::HatNP) = 2.0 ⊗₊ Float64(length(B) - 1, RoundDown)
aux_weak_bound(B::HatNP) = 1.0
weak_by_strong_and_aux_bound(B::HatNP) = (1.0, 1.0)
bound_weak_norm_from_linalg_norm(B::HatNP) = @error "TODO"
bound_linalg_norm_L1_from_weak(B::HatNP) = @error "TODO"
bound_linalg_norm_L∞_from_weak(B::HatNP) = @error "TODO"
opnormbound(B::HatNP{T}, N::Type{Linf}, A::AbstractVecOrMat{S}) where {S,T} =
    opnormbound(N, A)
normbound(B::HatNP{T}, N::Type{Linf}, v) where {T} = normbound(N, v)

function invariant_measure_strong_norm_bound(B::HatNP, D::Dynamic)
    A, B = dfly(strong_norm(B), aux_norm(B), D)
    @assert A < 1.0
    return B ⊘₊ (1.0 ⊖₋ A)
end

struct HatNPDual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    x′::Vector{Interval}
end

function HatNPDualBranch(y, br::MonotonicBranch, ylabel = 1:length(y); ϵ, max_iter)
    if is_increasing(br)
        endpoint_X = br.X[2]
        der = derivative(br.f)(endpoint_X)
        preim_der = preimages_and_derivatives(y, br, ylabel; ϵ, max_iter)
        return [preim_der[1]; endpoint_X],
        [preim_der[2]; length(preim_der[2]) + 1],
        [preim_der[3]; der]
    else
        endpoint_X = br.X[2]
        der = derivative(br.f)(endpoint_X)
        preim_der = preimages_and_derivatives(B.p, D, 1:length(B.p)-1; ϵ, max_iter)
        return [preim_der[1]; endpoint_X],
        [preim_der[2]; length(preim_with_der[2]) + 1],
        [preim_der[3]; der]
    end
end

function Dual(B::HatNP, D; ϵ = 0.0, max_iter = 100)
    @assert is_full_branch(D)
    results =
        collect(HatNPDualBranch(B.p, b, 1:length(B.p)-1; ϵ, max_iter) for b in branches(D))
    x = vcat((result[1] for result in results)...)
    xlabel = vcat((result[2] for result in results)...)
    x′ = vcat((result[3] for result in results)...)
    return HatNPDual(x, xlabel, x′)
end

Base.length(dual::HatNPDual) = length(dual.x)
Base.eltype(dual::HatNPDual) =
    Tuple{eltype(dual.xlabel),Tuple{eltype(dual.x),eltype(dual.x′)}}
function Base.iterate(dual::HatNPDual, state = 1)
    if state <= length(dual.x)
        return ((dual.xlabel[state], (dual.x[state], abs(dual.x′[state]))), state + 1)
    else
        return nothing
    end
end

function change_of_basis(Bu::Ulam, Bh::HatNP, v::Vector{T}) where {T}
    @assert length(Bh) == length(Bu) + 1
    w = zeros(T, length(Bh))
    w[1] = v[1]
    for i = 2:length(Bu)
        w[i] = (v[i-1] + v[i]) / 2
    end
    w[end] = v[end]
    return w
end
