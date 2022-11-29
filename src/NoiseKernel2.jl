using LinearAlgebra

struct NoiseUlam
    B::Ulam
    k::Integer
    BC::Matrix
    Ext::Matrix
    w_ext::Vector
    w_conv::Vector
    ξ 
end

NoiseUlam(B, k, BC, Ext) = NoiseUlam(B, k, BC, Ext, zeros(length(B)+k-1), zeros(length(B)+k-1), Float64(k)/length(B))


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
    l = (k-1) ÷ 2

    n = length(B) 

    w_ext = zeros(n+2*k)
    w_conv = zeros(n+2*k)

    if boundary_condition == :periodic 
        BC = PeriodicBoundaryConditionOperator(n, l)
    elseif boundary_condition == :reflecting
        BC = PeriodicBoundaryConditionOperator(n, l)
    else
        @error "Not implemented"
    end
    Ext = ExtensionOperator(n, l)
    
    return NoiseUlam(B, k, BC, Ext) 
end

import LinearAlgebra: mul!

function uniform_convolution!(w_conv, w_ext, l)
    n = length(w_ext)-2*l
    
    for i in 1:l
        w_conv[i] = sum(@view w_ext[i:i+l])
    end 

    for i in 1:n
        w_conv[i+l] = sum(@view w_ext[i:i+2*l]) 
    end

    for i in 1:l
        w_conv[n+l+i] = sum(@view w_ext[n+i:n+l+i])
    end 
    w_conv ./= (2*l+1)
end

function LinearAlgebra.mul!(w::Vector{Float64}, N::NoiseUlam, v::Vector{Float64})
    mul!(N.w_ext, N.Ext, v)
    k = N.k
    l = (k-1) ÷ 2

    uniform_convolution!(N.w_conv, N.w_ext, l)

    mul!(w, N.BC, N.w_conv)
end 

function Base.:*(N::NoiseUlam, v::Vector{Float64})
    w = zeros(length(v))
    mul!(w, N, v)
    return w
end