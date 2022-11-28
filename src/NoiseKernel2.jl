using LinearAlgebra

struct NoiseUlam
    B::Ulam
    k::Integer
    BC
    Ext
    w_ext::Vector
    w_conv::Vector
end

"""
    Extension_Operator
Define an operator that extends a vector v
to [zeros(k); v; zeros(k)], as a sparse matrix,
so it can be used in Blas
"""
function ExtensionOperator(n, k)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for i in 1:n
        push!(I, i+k)
        push!(J, i)
        push!(V, 1)
    end
    push!(I, n+2*k)
    push!(J, n)
    push!(V, 0)
    return sparse(I, J, V)
end

function PeriodicBoundaryConditionOperator(n, k)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for i in 1:n
        
        if i<=k
            push!(I, n-k+i)
            push!(J, i)
            push!(V, 1)

            push!(I, i)
            push!(J, n+k+i)
            push!(V, 1)
        end
        push!(I, i)
        push!(J, i+k)
        push!(V, 1)
    end
    return sparse(I, J, V)
end

function ReflectingBoundaryConditionOperator(n, k)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for i in 1:n
        
        if i<=k
            push!(I, k+1-i)
            push!(J, i)
            push!(V, 1)

            push!(I, n-i+1)
            push!(J, n+k+i)
            push!(V, 1)
        end
        push!(I, i)
        push!(J, i+k)
        push!(V, 1)
    end
    return sparse(I, J, V)
end


"""
    UniformNoiseUlam2(B, k, boundary_condition)

Defines a uniform noise kernel of relative size k
with k odd
"""
function UniformNoiseUlam2(B::Ulam, k; boundary_condition = :periodic)
    @assert k%2 == 1

    n = length(B) 

    w_ext = zeros(n+2*k)
    w_conv = zeros(n+2*k)

    if boundary_condition == :periodic 
        BC = PeriodicBoundaryConditionOperator(n, k)
    elseif boundary_condition == :reflecting
        BC = PeriodicBoundaryConditionOperator(n, k)
    else
        @error "Not implemented"
    end
    Ext = ExtensionOperator(n, k)
    
    return NoiseUlam(B, k, BC, Ext, w_ext, w_conv) 
end

import LinearAlgebra: mul!

function LinearAlgebra.mul!(w::Vector{Float64}, N::NoiseUlam, v::Vector{Float64})
    @info typeof(N)

    mul!(N.w_ext, N.Ext, v)

    n = length(N.B)
    k = N.k
    # # we assume it is extended by zero
    for i in 1:k
        N.w_conv[i] = sum(N.w_ext[i:i+k])
    end 

    for i in 1:length(v)
        N.w_conv[i+k] = sum(N.w_ext[i:i+2*k]) 
    end

    for i in 1:k
        N.w_conv[n+k+i] = sum(N.w_ext[n+i-k:n+i])
    end 

    N.w_conv ./= k

    mul!(w, N.BC, N.w_conv)
end 