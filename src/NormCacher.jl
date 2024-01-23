"""
Types to compute norms iteratively by "adding a column at a time".
"""
abstract type NormCacher{T} end

mutable struct NormCacherL1 <: NormCacher{L1}
    B::Basis
    C::Float64
    function NormCacherL1(B, n)
        new(B, 0)
    end
end
"""
Create a new NormCacher to compute the normbound of the empty matrix with n rows
"""
NormCacher{L1}(B, n) = NormCacherL1(B, n)

mutable struct NormCacherLinf <: NormCacher{Linf}
    B::Basis
    C::Vector{Float64}
    function NormCacherLinf(B, n)
        new(B, zeros(n))
    end
end
NormCacher{Linf}(B, n) = NormCacherLinf(B, n)

"""
Update a NormCacher to add one column to the matrix it is computing a norm of.
This column may be affected by an error ε (in the same norm).
"""
function add_column!(Cacher::NormCacherL1, v::AbstractVector, ε::Float64)
    Cacher.C = max(Cacher.C, opnormbound(Cacher.B, L1, v) ⊕₊ ε)
end

function add_column!(Cacher::NormCacherLinf, v::AbstractVector, ε::Float64)
    Cacher.C = Cacher.C .⊕₊ abs.(v) .⊕₊ ε
end

"""
Return the norm of the matrix the NormCacher is working on.
"""
function get_norm(Cacher::NormCacherL1)
    return Cacher.C
end

function get_norm(Cacher::NormCacherLinf)
    return maximum(Cacher.C)
end
