using ValidatedNumerics
using .DynamicDefinition, .PwDynamicDefinition
using .Contractors

using .DynamicDefinition: derivative, plottable

using TaylorSeries: Taylor1

struct Iterate <: Dynamic
    D::PwMap
    n::Int
end

Base.show(io::IO, D::Iterate) = print(io, "$(D.n)-times iterate of: ", D.D)

function (D::Iterate)(x::Taylor1)
    y = x
    for i = 1:D.n
        y = (D.D)(y)
    end
    return y
end

DynamicDefinition.nbranches(D::Iterate) = nbranches(D.D)^D.n
DynamicDefinition.is_full_branch(D::Iterate) = is_full_branch(D.D)
DynamicDefinition.domain(D::Iterate) = domain(D.D)

"""
compose_endpoints(D, endpointsf)

Find discontinuity points of f ∘ D, where f has discontinuity points endpointsf
and D is a given PwMap.
"endpoints" always include the extrema of the given domain.
"""
function compose_endpoints(D, endpointsf)
    v = [D.endpoints[1]]
    for k = 1:nbranches(D)
        preims = [preim(D, k, x) for x in endpointsf[2:end-1]]
        if orientation(D, k) > 0
            # TODO: still trouble if the branches are not full
            append!(v, preims)
        else
            append!(v, preims[end:-1:1])
        end
        append!(v, [D.endpoints[k+1]])
    end
    return v
end

function DynamicDefinition.endpoints(D::Iterate)
    @assert D.n >= 1
    endpoints = D.D.endpoints
    for k = 2:D.n
        endpoints = compose_endpoints(D.D, endpoints)
    end
    return endpoints
end

"""
unpack(k, b, D)

Convert an integer k∈[1,b^n] into a tuple ∈[1,b]^n bijectively

This is used to index preimages: the k'th of the b^n preimages of an Iterate
corresponds to choosing the v[i]'th branch when choosing the i'th preimage, for i = 1..k,
where v = unpack(k, b, n)
"""
function unpack(k, b, n)
    @assert 1 ≤ k ≤ b^n
    v = fill(0, n)
    k = k-1
    for i = 1:n
        (k, v[n+1-i]) = divrem(k, b)
    end
    return v .+ 1
end

function evaluate_branch(D::Iterate, k, x)
    @assert 1 ≤ k ≤ nbranches(D)
    n = D.n
    v = unpack(k, nbranches(D.D), n)
    for i = 1:n
        x = branch(D.D, v[i])(x)
    end
    return x
end

DynamicDefinition.branch(D::Iterate, k) = x -> evaluate_branch(D, k, x)

using LinearAlgebra: Bidiagonal

function Jac(fprime, v::Vector{T}) where {T}
    dv = fprime.(v)
    ev = -ones(T, length(v)-1)
    return Bidiagonal{T}(dv, ev, :U)
end

"""
Compute the preimage of an Iterate D in its k'th branch
"""
function DynamicDefinition.preim(D::Iterate, k, y, ϵ=1e-15; max_iter = 100)
    @assert 1 <= k <= nbranches(D)

    n = D.n

    v = unpack(k, nbranches(D.D), n)

    fs = D.D.Ts[v]
    S = [hull(D.D.endpoints[v[i]], D.D.endpoints[v[i]+1]) for i in 1:n]
    return nthpreimage!(y, fs, S)[1]
end

function DynamicDefinition.orientation(D::Iterate, k)
    v = unpack(k, nbranches(D.D), D.n)
    return prod([orientation(D.D, k) for k in v])
end

function DynamicDefinition.plottable(D::Iterate, x)
	@assert 0 <= x <= 1
    for k = 1:D.n
        x = DynamicDefinition.plottable(D.D, x)
        x = clamp(x, 0, 1)
    end
    return x
end

using RecipesBase
@recipe f(::Type{Iterate}, D::Iterate) = x -> plottable(D, x)
