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

"""
Given an already-computed sequence of preimages, constructs associated "Ulam duals",
i.e., the sequence (k, (T⁻¹(p[k]), T⁻¹(p[k+1]))),
handling endpoints correctly.
This should eventually replace `DualComposedWithDynamic`.
"""
ulam_duals(b, partition, preims=nothing) = ulam_duals(b, PointSequence(partition), preims)
function ulam_duals(b::Branch, partition::PointSequence, preims=nothing)
    if preims === nothing
        preims = preimages(partition, b)
    end
    n = length(preims.v)

    duals = Tuple{Int, Tuple{eltype(preims), eltype(preims)}}[]

    if n==0 # special case: only one preimage
        first_endpoint = preims.increasing ? b.X[1] : b.X[2]
        last_endpoint = preims.increasing ? b.X[2] : b.X[1]

        push!(duals, (preims.skip, (first_endpoint, last_endpoint)))
        return duals
    end

    sizehint!(duals, n+2)

    if preims.skip > 0
        first_endpoint = preims.increasing ? b.X[1] : b.X[2]
        if first_endpoint != preims.v[1]
            push!(duals, (preims.skip, (first_endpoint, preims.v[1])))
        end
    end
    for k in 1:n-1
        push!(duals, (preims.skip+k, (preims.v[k], preims.v[k+1])))
    end
    if preims.skip + n < partition.skip+length(partition.v)  # skipped entries in preims
        last_endpoint = preims.increasing ? b.X[2] : b.X[1]
        if last_endpoint != preims.v[n]
            push!(duals, (preims.skip+n, (preims.v[n], last_endpoint)))
        end
    end
    return duals
end
