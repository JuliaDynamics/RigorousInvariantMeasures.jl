module C2BasisDefinition

"""
C2 basis on [0,1]
"""

using ..BasisDefinition, ..DynamicDefinition, IntervalArithmetic
import Base: iterate, length
import ..BasisDefinition: one_vector, integral_covector, is_integral_preserving
import ...RigorousInvariantMeasures: NormKind, Linf

export C2Basis, dual_val, dual_der, C1, C2

#struct C1 <: NormKind end
#struct C2 <: NormKind end

"""
Equispaced partition of size n of [0,1]
"""
struct EquispacedPartitionInterval{T} <: AbstractVector{T}
    n::Integer
end

function Base.getindex(p::EquispacedPartitionInterval{T}, i::Int) where {T}
    @boundscheck 1 <= i <= p.n + 1 || throw(BoundsError(p, i))
    return convert(T, i - 1) / p.n
end

EquispacedPartitionInterval(i::Int) = @error "The real type must be specified explicitly"

Base.size(p::EquispacedPartitionInterval) = (p.n + 1,)
Base.IndexStyle(::Type{<:EquispacedPartitionInterval}) = IndexLinear()
Base.issorted(p::EquispacedPartitionInterval) = true

struct C2Basis{T<:AbstractVector} <: Basis
    p::T
    # TODO: check in constructor that p is sorted and starts with 0
end
C2Basis(n::Integer) = C2Basis(EquispacedPartitionInterval{Float64}(n))

"""
Return the size of the C2 basisBase.length(S::AverageZero) = length(S.basis)-1
"""
Base.length(b::C2Basis) = 2 * length(b.p)

function ϕ(x::Interval{T}) where {T}
    if x ∩ Interval{T}(-1, 1) == ∅
        return zero(x)
    else
        x₋ = x ∩ Interval{T}(-1, 0)
        val₋ = evalpoly(x₋, (1, 0, 0, 10, 15, 6))
        x₊ = x ∩ Interval{T}(0, 1)
        val₊ = evalpoly(x₊, (1, 0, 0, -10, 15, -6))
        return val₋ ∪ val₊
    end
end

function ϕprime(x::Interval{T}) where {T}#Derivative of ϕ
    if x ∩ Interval{T}(-1, 1) == ∅
        return zero(x)
    else
        x₋ = x ∩ Interval{T}(-1, 0)
        val₋ = evalpoly(x₋, (0, 0, 30, 60, 30))
        x₊ = x ∩ Interval{T}(0, 1)
        val₊ = evalpoly(x₊, (0, 0, -30, 60, -30))
        return val₋ ∪ val₊
    end
end

function ν(x::Interval{T}) where {T}
    if x ∩ Interval{T}(-1, 1) == ∅
        return zero(x)
    else
        x₋ = x ∩ Interval{T}(-1, 0)
        val₋ = evalpoly(x₋, (0, 1, 0, -6, -8, -3))
        x₊ = x ∩ Interval{T}(0, 1)
        val₊ = evalpoly(x₊, (0, 1, 0, -6, 8, -3))
        return val₋ ∪ val₊
    end
end

function νprime(x::Interval{T}) where {T} #Derivative of ν

    if x ∩ Interval{T}(-1, 1) == ∅
        return zero(x)
    else
        x₋ = x ∩ Interval{T}(-1, 0)
        val₋ = evalpoly(x₋, (1, 0, -18, -32, -15))
        x₊ = x ∩ Interval{T}(0, 1)
        val₊ = evalpoly(x₊, (1, 0, -18, 32, -15))
        return val₋ ∪ val₊
    end
end

κ(x) = 6 * x * (1 - x)

function Base.getindex(B::C2Basis, i::Int)
    n = length(B.p)
    @boundscheck 1 <= i <= 2 * n || throw(BoundsError(B, i))
    if i <= n
        #println("val")
        return x -> ϕ((n - 1) * x - i + 1)
    elseif n < i <= 2 * n
        return x -> ν((n - 1) * x - i + n + 1) / (n - 1) #we have to subtract n to i
    end
    #return x->κ(x) 
end


function basis_element_with_der(B::C2Basis, i::Int)
    n = length(B.p)
    @boundscheck 1 <= i <= 2 * n || throw(BoundsError(B, i))
    if i <= n
        #println("val")
        return x -> ϕ((n - 1) * x - i + 1), x -> (n - 1) * ϕprime((n - 1) * x - i + 1)
    elseif n < i <= 2 * n
        return x -> ν((n - 1) * x - i + n + 1) / (n - 1),
        x -> νprime((n - 1) * x - i + n + 1) #we have to subtract n to i
    end
    #return x->κ(x) for i in 1:N+1 
end

function BasisDefinition.is_dual_element_empty(::C2Basis, d)
    # TODO: the preim() may indeed be empty, so there could be an additional check here
    return false
end

# Check the absolute values in the formula
dual_val(f::Function, fprime::Function, x, der, derder) = f(x) / der
dual_der(f::Function, fprime::Function, x, der, derder) =
    fprime(x) / der^2 - f(x) * derder / der^3
"""
Return (in an iterator) the pairs (i, (x, |T'(x)|)) where x is a preimage of p[i], which
describe the "dual" L* evaluation(p[i])
"""
function Base.iterate(
    S::DualComposedWithDynamic{T,Dynamic},
    state = (1, 1),
) where {T<:C2Basis}
    i, k = state

    if i == length(S.basis) + 1
        return nothing
    end

    # remark that this version supposes that for each i there exists a preimage
    # another more specific version should be implemented for maps with
    # incomplete branches

    x = preim(S.dynamic, k, S.basis.p[i], S.ϵ)
    absT′ = abs(der(S.dynamic, x))
    derder = der_der(S.dynamic, x)

    if i <= n
        ret = (x, (f, fprime) -> dual_val(f, fprime, x, der, derder))
    else
        ret = (x, (f, fprime) -> dual_der(f, fprime, x, der, derder))
    end

    #if k == nbranches(S.dynamic)
    #	return ((i, ret), (i+1, 1))
    #else
    #	return ((i, ret), (i, k+1))
    #end
    #if k == nbranches(S.dynamic)
    #	return ((i, (x, absT′, derder)), (i+1, 1))
    #else
    #	return ((i, (x, absT′, derder)), (i, k+1))
    #end
end

BasisDefinition.is_refinement(Bf::C2Basis, Bc::C2Basis) = Bc.p ⊆ Bf.p

function integral_covector(B::C2Basis)
    n = length(B.p)
    return 1 / (n - 1) *
           [
        @interval 0.5
        ones(Interval{Float64}, n - 2)
        @interval 0.5
        @interval 1 // 10
        zeros(n - 2)
        @interval -1 // 10
    ]'
end

function one_vector(B::C2Basis)
    return [ones(length(B.p)); zeros(length(B.p))]
end

"""
Return the range of indices of the elements of the basis whose support intersects
with the given dual element (i.e., a pair (y, absT', derder)).
"""
function BasisDefinition.nonzero_on(B::C2Basis, dual_element)
    y = dual_element[1]
    # Note that this cannot rely on arithmetic unless it is verified
    # searchsortedfirst(a, x) return the index of the first value in a greater than or equal to x
    lo = max(searchsortedfirst(B.p, y.lo) - 1, 1)
    hi = searchsortedfirst(B.p, y.hi)
    if lo == 1 # 1:N+1 does not make sense and becomes 1:N
        hi = min(hi, length(B))
    end
    return (lo, hi)
end

"""
Given a preimage ```y``` of a point ```x```, this iterator returns
```\\phi_j(y)/T'(y) ```
"""
function Base.iterate(S::ProjectDualElement{T}, state = (S.j_min, :val)) where {T<:C2Basis}
    dual = S.dual_element[2]
    j = state[1]
    n = length(S.basis.p)

    if state[2] == :val              # we first fill up the coordinates of the ϕ 
        f, fprime = basis_element_with_der(S.basis, j)
        if state[1] == S.j_max
            return ((j, dual(f, fprime)), (S.j_min, :der))
        end
        return ((j, dual(f, fprime)), (j + 1, :val))
    end
    if state[2] == :der             # we then fill up the coordinates of the ν
        if state[1] == S.j_max + 1
            return nothing
        end
        f, fprime = basis_element_with_der(S.basis, j + n)
        return ((j + n, dual(f, fprime)), (j + 1, :der))
    end
end

using IntervalOptimisation
function infnormoffunction(B::C2Basis, v)
    n = length(B.p)
    maximum = -Inf
    for i = 1:length(B.p)-1
        coeff = v[i] * [1, 0, 0, -10, 15, -6] #coeff for ϕ
        coeff += (v[i+n] / (n - 1)) * [0, 1, 0, -6, 8, -3]  #coeff for ν
        # coefficients from the right endpoint
        coeff += v[i+1] * [0, 0, 0, 10, -15, 6]
        coeff += (v[i+n+1] / (n - 1)) * [0, 0, 0, -4, +7, -3]

        dom = Interval(0, 1)
        f(x) = abs(evalpoly(x, coeff))
        maximum = max(maximum, maximise(f, dom)[1])
    end
    return maximum
end

function infnormofderivative(B::C2Basis, v)
    n = length(B.p)
    maximum = -Inf
    for i = 1:length(B.p)-1
        coeff = (n - 1) * v[i] * [0, 0, -30, 60, -30] #coeff for ϕ
        coeff += v[i+n] * [1, 0, -18, 32, -15]  #coeff for ν
        # coefficients from the right endpoint
        coeff += (n - 1) * v[i+1] * [0, 0, 30, -60, 30]
        coeff += v[i+n+1] * [0, 0, -12, +28, -15]

        dom = Interval(0, 1)
        f(x) = abs(evalpoly(x, coeff))
        maximum = max(maximum, maximise(f, dom)[1])
    end
    return maximum
end

C1Norm(B::C2Basis, v) = infnormoffunction(B, v) + infnormofderivative(B, v)
rescaling_factor(B::C2Basis) = 3 * length(B.p)

Base.length(S::AverageZero{T}) where {T<:C2Basis} = length(S.basis) - 1

function Base.iterate(S::AverageZero{T}, state = 1) where {T<:C2Basis}
    n = length(S.basis) ÷ 2
    v = zeros(2 * n)

    i = state
    if i == 2 * n
        return nothing
    elseif 1 <= i < n
        v[i+1] = 1
        v[1] = -2
    elseif i == n
        v[i+1] = 1
        v[1] = -1 / 5
    elseif n < i < 2 * n - 1
        v[i+1] = 1
    elseif i == 2 * n - 1
        v[i+1] = 1
        v[1] = 1 / 5
    end
    return v, state + 1
end


end

using RecipesBase

@userplot PlotC2
@recipe function f(h::PlotC2)
    if length(h.args) != 2 ||
       (typeof(h.args[1]) != C2) ||
       !(typeof(h.args[2]) <: AbstractVector)
        error("Plot C2 needs as an input a C2 Basis and a vector")
    end

    B = h.args[1]
    w = h.args[2]

    seriestype := :path
    collect(B), mid.(w)
end
