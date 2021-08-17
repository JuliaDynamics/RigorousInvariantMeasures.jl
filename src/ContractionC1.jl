using .C2BasisDefinition: C1, C1Norm



function opnormbound(N::Type{<:C1}, M, B)
    n = size(M, 2)

    est = 0

    for i in 1:n
        v = zeros(n)
        v[i] = 1
        z = C1Norm(B, v)
        est+= C1Norm(B, M*v)/z
    end
    return est.hi
end

using ProgressMeter
function norms_of_powers_basis(B, m::Integer, Q::DiscretizedOperator, f::AbstractArray;
    normv0::Real=-1., #used as "missing" value
    normQ::Real=-1.,
    normE::Real=-1.,
    normEF::Real=-1.,
    normIEF::Real=-1.,
    normN::Real=-1.) 

    @assert eltype(f) <: Interval
    T = typeof(zero(eltype(Q.L)).hi) # gets "Float64" from Q.L
    n = size(Q.L, 1)
    M = mid.(Q.L)
    R = radius.(Q.L)
    #δ = opnormbound(N, R)
    #γz = gamma(T, max_nonzeros_per_row(Q.L))
    #γn = gamma(T, n+3) # not n+2 like in the paper, because we wish to allow for f to be the result of rounding
    #ϵ = zero(T)
    midf = mid.(f)

    # TODO: correct here
    norm(v) = C1Norm(B, v)

    norms = zeros(m)

    #S = zeros((n, m))
    #k = length(B.p)
    factor = rescaling_factor(B)

    @showprogress for (v, norm_0) in AverageZero(B)
        
        #@info v
        v/= norm_0
        #@info v
        #@info norm(v)
        for i in 1:m
            
            w = M*v 
            v = w - Q.e * (midf*w)[1]
            #@info v

            #@info "norm_$i" norm(v)*factor
           # @info infnormoffunction(B, v)
           # @info infnormofderivative(B, v)
            norms[i] = max(norm(v).hi*factor, norms[i])
        end
        #@info norms
    end

    return norms
    
end








