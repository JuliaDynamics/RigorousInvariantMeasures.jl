using LinearAlgebra

export UniformNoiseUlam2

struct NoiseUlam <: NoiseKernel
    B::Ulam
    k::Integer
    BC!
    Ext::Matrix
    w_ext::Vector
    w_conv::Vector
    ξ 
end

NoiseUlam(B, k, BC, Ext) = NoiseUlam(B, k, BC, Ext, zeros(length(B)+k-1), zeros(length(B)+k-1), Interval(k)/(2*length(B)))


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


function PeriodicBoundaryCondition2!(w_out, w_in; n, l)
    w_out .= @view w_in[l+1:n+l] 
    w_out[1:l] += @view w_in[n+l+1:n+2*l]
    w_out[end-l+1:end] += @view w_in[1:l]
    return nothing
end

function ReflectingBoundaryCondition2!(w_out, w_in; n, l)
    w_out .= @view w_in[l+1:n+l] 
    w_out[1:l] += reverse(@view w_in[1:l])
    w_out[end-l+1:end] += reverse(@view w_in[n+l+1:n+2*l])
    return nothing
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
        BC! = PeriodicBoundaryCondition2!
    elseif boundary_condition == :reflecting
        BC! = ReflectingBoundaryCondition2!
    else
        @error "Not implemented"
    end
    Ext = ExtensionOperator(n, l)
    
    return NoiseUlam(B, k, BC!, Ext) 
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

#using TimerOutputs

#const to = TimerOutput()

function LinearAlgebra.mul!(w::Vector{Float64}, N::NoiseUlam, v::Vector{Float64})
    k = N.k
    l = (k-1) ÷ 2
    n = length(v)
    N.w_ext[l+1:l+n] = v
    uniform_convolution!(N.w_conv, N.w_ext, l)
    N.BC!(w, N.w_conv; n = n, l = l)
end 

function Base.:*(N::NoiseUlam, v::Vector{Float64})
    w = zeros(length(v))
    mul!(w, N, v)
    return w
end

BasisDefinition.opnormbound(B::Ulam, ::Type{L1}, M::NoiseUlam) = 1.0
opradius(::Type{L1}, N::NoiseUlam) = N.k*Interval(radius(Interval(1)/N.k)).hi
nonzero_per_row(N::NoiseUlam) = N.k
dfly(::Type{TotalVariation}, ::Type{L1}, N::NoiseUlam) = (0.0, (1/(2*N.ξ)).hi)