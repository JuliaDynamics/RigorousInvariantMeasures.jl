using IntervalArithmetic
using .DynamicDefinition

"""
Return MonotonicBranches for a given dynamic, in an iterable
"""
function branches(D::PwMap)
    return [
        MonotonicBranch(
            D.Ts[k],
            (D.endpoints[k], D.endpoints[k+1]),
            D.is_full[k] ? (Interval(0), Interval(1)) :
            (D.Ts[k](Interval(D.endpoints[k]), D.Ts[k](Interval(D.endpoints[k+1])))),
            D.orientations[k] == 1,
        ) for k = 1:length(D.Ts)
    ]
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
handling endpoints correctly. A callback `f(k, (a, b))` is called on each dual.

This should eventually replace `DualComposedWithDynamic`.
"""
function callback_duals(f, B::Ulam, branch::MonotonicBranch, preims = nothing)
    if preims === nothing
        preims = preimages(PointSequence(B.p), b)
    end
    n = length(preims.v)

    # duals = Tuple{Int, Tuple{eltype(preims), eltype(preims)}}[]
    # sizehint!(duals, n+2)

    if n == 0 # special case: only one preimage
        first_endpoint = preims.increasing ? branch.X[1] : branch.X[2]
        last_endpoint = preims.increasing ? branch.X[2] : branch.X[1]

        f(preims.skip, (first_endpoint, last_endpoint))
        return nothing
    end

    if preims.skip > 0
        first_endpoint = preims.increasing ? branch.X[1] : branch.X[2]
        if first_endpoint != preims.v[1]
            f(preims.skip, (first_endpoint, preims.v[1]))
        end
    end
    for k = 1:n-1
        f(preims.skip + k, (preims.v[k], preims.v[k+1]))
    end
    if preims.skip + n < length(B.p)  # if there are skipped entries at the end of preims
        last_endpoint = preims.increasing ? b.X[2] : b.X[1]
        if last_endpoint != preims.v[n]
            f(preims.skip + n, (preims.v[n], last_endpoint))
        end
    end
    return nothing
end

function callback_duals(f, B::Hat, branch::MonotonicBranch, preims = nothing)
    if preims === nothing
        preims = preimages(PointSequence(B.p), b)
    end
    @assert preims.skip == 0
    @assert length(preims.v) == length(B.p) # The Hat basis estimation works only for full-branch maps
    for i = 1:length(B) # this skips p[n+1]==1
        x = preims.v[i]
        absT′ = abs(derivative(branch.f, x))
        f(i, (x, absT′))
    end
end

function callback_duals(f, B::Basis, D::Dynamic)
    for (branch, preim) in zip(branches(D), preimages(PointSequence(B.p), D))
        callback_duals(f, B, branch, preim)
    end
end

function assemble2(B::Basis, D::Dynamic, ϵ = 0.0; T = Float64)
    # putting types here in hope to improve callback inference
    I::Vector{Int64} = Int64[]
    J::Vector{Int64} = Int64[]
    nzvals::Vector{Interval{T}} = Interval{T}[]
    n::Int64 = length(B)

    # TODO: reasonable size hint?

    function assemble_callback(i, dual_element)
        if !is_dual_element_empty(B, dual_element)
            for (j, x) in ProjectDualElement(B, dual_element)
                push!(I, i)
                push!(J, mod(j, 1:n))
                push!(nzvals, x)
            end
        end
    end
    callback_duals(assemble_callback, B, D)

    return sparse(I, J, nzvals, n, n)
end
