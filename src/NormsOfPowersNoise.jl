"""
Estimates the norms ||RQ||, ||(RQ)^2||, ... ||(RQ)^m|| on U^0.

U is the matrix [ones(1,n-1); -I_(n-1,n-1)]. It is currently assumed that
f*U==0 (i.e., all elements of f are equal).

f must be an interval vector.

The following constants may be specified as keyword arguments:

normQ, normR, normE, normv0, normEF, normIEF, normN

otherwise they are computed (which may be slower).

e and f must be specified in case is_integral_preserving==false
In case is_integral_preserving is true, they may be specified but they are then ignored.

Implementation note: currently we perform this computation one column at a time,
to be able to scale (slowly) to cases with large size; for moderate sizes, it would
indeed be better to do the computation all columns at the same time, in BLAS level 3.
"""
function norms_of_powers_noise(N::Type{<:NormKind}, m::Integer, Q::DiscretizedOperator, NK::NoiseOperator, f::AbstractArray;
        normv0::Real=-1., #used as "missing" value
        normQ::Real=-1.,
        normNK::Real=-1,
        normE::Real=-1.,
        normEF::Real=-1.,
        normIEF::Real=-1.,
        normN::Real=-1.)

    @assert eltype(f) <: Interval
    T = typeof(zero(eltype(Q.L)).hi) # gets "Float64" from Q.L
    n = size(Q.L, 1)
    M = mid.(Q.L)
        
    R = radius.(Q.L)
    δ = opnormbound(N, R)
    
    # this is the operator radius of the noise operator,
    # we use this function since the operator may be defined implictly
    δₖ = opradius(NK)
    
    γz = gamma(T, max_nonzeros_per_row(Q.L))
    γn = gamma(T, n+3) # not n+2 like in the paper, because we wish to allow for f to be the result of rounding
    
    γₖ = gamma(T, max_nonzero_per_row(NK))
    
    ϵ = zero(T)

    nrmM = opnormbound(N, M)

    #this needs to be implemented
    nrm_mid_MK = opnormbound(N, MK)
    normMK = nrm_mid_MK ⊕₊ δₖ

    # precompute norms
    if !is_integral_preserving(Q)
        if normE == -1.
            normE = opnormbound(N, Q.e)
        end
        if normEF == -1.
            normEF = opnormbound(N, Q.e*f)
        end
        if normIEF == -1.
            normIEF =  opnormbound(N, [Matrix(UniformScaling{Float64}(1),n,n) Q.e*f])
        end
        if normN == -1.
            normN = opnormbound(N, Matrix(UniformScaling{Float64}(1),n,n) - Q.e*f)
        end
    end

    if normQ == -1.
        if is_integral_preserving(Q)
            normQ = nrmM ⊕₊ δ
        else
            defect = opnormbound(N, Q.w)
            normQ = nrmM ⊕₊ δ ⊕₊ normE ⊗₊ defect
        end
    end

    # initialize normcachers
    normcachers = [NormCacher{N}(n) for j in 1:m]

    midf = map(mid, f)

    # main loop

    v = zeros(T, n)
    for j = 1:n-1
        v .= zero(T) # TODO: check for type stability in cases with unusual types
        v[1] = one(T) # TODO: in full generality, this should contain entries of f rather than ±1
        v[j+1] = -one(T)
        if normv0 == -1.
            nrmv = opnormbound(N, v)
        else
            nrmv = normv0
        end
        ϵ = 0.
        nrmw = nrmv # we assume that initial vectors are already integral-preserving
        for k = 1:m
            w = M * v
            if is_integral_preserving(Q)
                v = w
                ϵ = (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
            else
                v = w - Q.e * (midf*w)  # TODO: we are currently assuming that f is not too large, to estimate the error (the result of only one floating point operation)
                new_nrmw = opnormbound(N, w)
                ϵ = γn ⊗₊ normIEF ⊗₊ (new_nrmw ⊕₊ normEF ⊗₊ nrmw) ⊕₊ normN ⊗₊ (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
                nrmw = new_nrmw
            end
            nrmv = opnormbound(N, v)
            v = w
            w = MK * v
            if is_integral_preserving(Q)
                v = w
                ϵ = (γₖ ⊗₊nrm_mid_MK ⊕₊δₖ) ⊗₊ nrmv ⊕₊ normMK ⊗₊ ϵ
            else
                v = w - Q.e * (midf*w)  # TODO: we are currently assuming that f is not too large, to estimate the error (the result of only one floating point operation)
                new_nrmw = opnormbound(N, w)
                ϵ = γn ⊗₊ normIEF ⊗₊ (new_nrmw ⊕₊ normEF ⊗₊ nrmw) ⊕₊ normN ⊗₊ (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmv ⊕₊ normQ ⊗₊ ϵ
                nrmw = new_nrmw
            end
            
            
            
            
            
            add_column!(normcachers[k], v, ϵ) #TODO: Could pass and reuse nrmv in the case of norm-1
        end
    end
    return map(get_norm, normcachers)
end