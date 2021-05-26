using IntervalArithmetic
using .DynamicDefinition

"""
Return Branches for a given dynamic, in an iterable
"""
function branches(D::PwMap)
    return [Branch(D.Ts[k], (D.endpoints[k], D.endpoints[k+1]), D.is_full[k] ? (Interval(0),Interval(1)) : nothing, D.orientations[k]==1) for k in 1:length(D.Ts)]
end

"""
Return preimages of a certain sequence on all branches of a dynamic
"""
function preimages(seq, D::Dynamic, ϵ = 0.0)
    return [preimages(seq, branch, ϵ) for branch in branches(D)]
end
# TODO: overload for IterateDynamic

"""
Constructs associated "Ulam duals" for a branch b,
i.e., the sequence (k, (T⁻¹(p[k]), T⁻¹(p[k+1]))),
handling endpoints correctly. The duals are then put! to a channel.

This should eventually replace `DualComposedWithDynamic`.
"""
put_duals!(ch, B::Ulam, branch::Branch, preims=nothing) = put_duals!(ch, B, branch, preims)
function put_duals!(ch::Channel, B::Ulam, branch::Branch, preims=nothing)
    if preims === nothing
        preims = preimages(PointSequence(B.p), b)
    end
    n = length(preims.v)

    # duals = Tuple{Int, Tuple{eltype(preims), eltype(preims)}}[]
    # sizehint!(duals, n+2)

    if n==0 # special case: only one preimage
        first_endpoint = preims.increasing ? branch.X[1] : branch.X[2]
        last_endpoint = preims.increasing ? branch.X[2] : branch.X[1]

        put!(ch, (preims.skip, (first_endpoint, last_endpoint)))
        return nothing
    end

    if preims.skip > 0
        first_endpoint = preims.increasing ? branch.X[1] : branch.X[2]
        if first_endpoint != preims.v[1]
            put!(ch, (preims.skip, (first_endpoint, preims.v[1])))
        end
    end
    for k in 1:n-1
        put!(ch, (preims.skip+k, (preims.v[k], preims.v[k+1])))
    end
    if preims.skip + n < length(B.p)  # if there are skipped entries at the end of preims
        last_endpoint = preims.increasing ? b.X[2] : b.X[1]
        if last_endpoint != preims.v[n]
            put!(ch, (preims.skip+n, (preims.v[n], last_endpoint)))
        end
    end
    return nothing
end

function put_duals!(ch, B::Basis, D::Dynamic)
    for (branch, preim) in zip(branches(D), preimages(PointSequence(B.p), D))
        put_duals!(ch, B, branch, preim)
    end
end

"""
Generator that returns duals
"""
function duals(B::Basis, D)
    T = Tuple{Int, Tuple{Interval{Float64}, Interval{Float64}}} # TODO: replace for more genericity?
    Channel{T}() do ch
        put_duals!(ch, B, D)
    end
end
