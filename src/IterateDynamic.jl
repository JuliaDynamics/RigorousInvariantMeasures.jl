using ValidatedNumerics
using .DynamicDefinition, .PwDynamicDefinition
using .Contractors

using .DynamicDefinition: derivative

using TaylorSeries: Taylor1

struct Iterate <: Dynamic
    D::PwMap
    n::Int
end

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

using LinearAlgebra: Bidiagonal

function Jac(fprime, v::Vector{T}) where {T}
    dv = fprime.(v)
    ev = -ones(T, length(v)-1)
    return Bidiagonal{T}(dv, ev, :U)
end

"""
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

"""
Compute the preimage of an Iterate D in its k'th branch
"""
function DynamicDefinition.preim(D::Iterate, k, y, ϵ=1e-15; max_iter = 100)
    @assert 1 <= k <= nbranches(D)

    n = D.n

    v = unpack(k, nbranches(D.D), n)

    # Uses the "shooting method", i.e., the multivariate interval Newton method on
    # [f(x_1)-x_2, f(x_2)-x_3, ..., f(x_{n})-y] == 0

    f = x -> [D.D.Ts[v[i]](x[i]) for i in 1:n] - [x[2:end]; y]
    f′ = x -> Bidiagonal([derivative(D.D.Ts[v[i]], x[i]) for i in 1:n],  fill(-Interval(1.), n-1), :U)
    S = IntervalBox(hull(D.D.endpoints[v[i]], D.D.endpoints[v[i]+1]) for i in 1:n)
    return root(f, f′, S, ϵ; max_iter = max_iter)[1]

end
