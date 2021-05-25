using IntervalArithmetic
using .DynamicDefinition

"""
Return Branches for a given dynamic, in an iterable
"""
function branches(D::PwMap)
    return [Branch(D.Ts[k], hull(D.endpoints[k], D.endpoints[k+1]), D.is_full[k] ? @interval(0,1) : nothing, D.orientations[k]==1) for k in 1:length(D.Ts)]
end

"""
Return preimages of a certain sequence on all branches of a dynamic
"""
function preimages(seq, D::Dynamic, ϵ = 0.0)
    return [preimages(seq, branch, ϵ) for branch in branches(D)]
end

"""
Given a sequence of preimages, constructs associated "duals" (k, (T⁻¹(p[k]), T⁻¹(p[k+1]))),
handling endpoints correctly.
"""
function duals(B::Ulam, seq)
    TODO
end