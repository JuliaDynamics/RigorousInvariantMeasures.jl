export SkewProductMap


# we expect the (T(x), G(x,y)) to be injective
# when restricted to the monotonicity domain of
# T×[0,1]
struct SkewProductMap
    T::PwMap
    G::Array{Function, 1}
end

Base.getindex(D::SkewProductMap, i::Integer) = (D.T.branches[i], D.G[i])

nbranches(D) = nbranches(D.T)

function FiberMap(D::SkewProductMap, x, y)
    intersection_x = intersect_domain_bool(D.T, x)
    @assert sum(intersection_x) == 1 "Intersects many branches, ambiguous" 
    return (D.G[intersection_x][1])(x, y)
end
