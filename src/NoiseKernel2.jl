using LinearAlgebra

export UniformNoiseUlam2

struct NoiseUlam <: NoiseKernel
    B::Ulam
    k::Integer
    BC!
    Ext::SparseMatrixCSC
    BC::SparseMatrixCSC
    w_ext::Vector
    w_conv::Vector
    ξ 
end

NoiseUlam(B, k, BC!, Ext, BC) = NoiseUlam(B, k, BC!, Ext, BC, zeros(length(B)+k-1), zeros(length(B)+k-1), Interval(k)/(2*length(B)))


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
        BC = PeriodicBoundaryConditionOperator(n, l)
    elseif boundary_condition == :reflecting
        BC! = ReflectingBoundaryCondition2!
        BC = ReflectingBoundaryConditionOperator(n, l)
    else
        @error "Not implemented"
    end
    Ext = ExtensionOperator(n, l)
    
    return NoiseUlam(B, k, BC!, Ext, BC) 
end

import LinearAlgebra: mul!

function uniform_convolution!(w_conv, w_ext, l)
    n = length(w_ext)-2*l
    
    for i in 1:l
        @inbounds w_conv[i] = sum(@view w_ext[i:i+l])
    end 

    for i in 1:n
        @inbounds w_conv[i+l] = sum(@view w_ext[i:i+2*l]) 
    end

    for i in 1:l
        @inbounds w_conv[n+l+i] = sum(@view w_ext[n+i:n+i+l])
    end 
    w_conv ./= (2*l+1)
end

#using TimerOutputs

#const to = TimerOutput()

function mul2!(w::Vector{Float64}, N::NoiseUlam, v::Vector{Float64})
    k = N.k
    l = (k-1) ÷ 2
    n = length(v)
    N.w_ext[l+1:l+n] = v
    uniform_convolution!(N.w_conv, N.w_ext, l)
    N.BC!(w, N.w_conv; n = n, l = l)
end 

function LinearAlgebra.mul!(w::Vector{Float64}, N::NoiseUlam, v::Vector{Float64})
    k = N.k
    l = (k-1) ÷ 2
    #n = length(v)
    mul!(N.w_ext, N.Ext, v)
    uniform_convolution!(N.w_conv, N.w_ext, l)
    mul!(w, N.BC, N.w_conv)
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

using CUDA

if has_cuda() && has_cuda_gpu()

    struct NoiseUlamCuda <: NoiseKernel
        B::Ulam
        k::Integer
        BC!
        Ext::CUDA.CUSPARSE.CuSparseMatrixCSC
        BC::CUDA.CUSPARSE.CuSparseMatrixCSC
        w_ext::CuArray{Float64}
        w_conv::CuArray{Float64}
        ξ
        threads
        blocks 
    end
    
    using Adapt

    NoiseUlamCuda(B, k, BC!, Ext, BC, threads, blocks) = NoiseUlamCuda(B, 
                                                    k, 
                                                    BC!, 
                                                    adapt(CUDA.CUSPARSE.CuSparseMatrixCSC,Ext), 
                                                    adapt(CUDA.CUSPARSE.CuSparseMatrixCSC,BC),
                                                    CUDA.zeros(length(B)+k-1), 
                                                    CUDA.zeros(length(B)+k-1), 
                                                    Interval(k)/(2*length(B)),
                                                    threads,
                                                    blocks)
    
    """
    UniformNoiseUlam2(B, k, boundary_condition)

    Defines a uniform noise kernel of relative size k
    with k odd
    """
    function UniformNoiseUlamCuda(B::Ulam, k; boundary_condition = :periodic)
        @assert k%2 == 1
        l = (k-1) ÷ 2

        n = length(B) 
        
        if boundary_condition == :periodic 
            BC! = PeriodicBoundaryCondition2!
            BC = PeriodicBoundaryConditionOperator(n, l)
        elseif boundary_condition == :reflecting
            BC! = ReflectingBoundaryCondition2!
            BC = ReflectingBoundaryConditionOperator(n, l)
        else
            @error "Not implemented"
        end
        Ext = ExtensionOperator(n, l)
        
        w_conv = CUDA.zeros(Float64, n+2*l)
        w_ext = CUDA.zeros(Float64, n+2*l)
        kernel = @cuda launch=false cuda_central_conv!(w_conv, w_ext, l)
        config = launch_configuration(kernel.fun)

        threads = min(length(w_conv), config.threads)
        blocks = cld(length(w_conv), threads)

        return NoiseUlamCuda(B, k, BC!, Ext, BC, threads, blocks) 
    end


    function cuda_left_conv!(w_conv, w_ext, l)
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = gridDim().x * blockDim().x
        for i = index:stride:l
            #w_conv[i] = sum(@view w_ext[i:i+l])
            for j in 0:l
                @inbounds w_conv[i] += w_ext[i+j]
            end
        end
        return nothing 
    end

    function cuda_central_conv!(w_conv, w_ext, l)
        n = length(w_ext)-2*l
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = gridDim().x * blockDim().x
        for i = index:stride:n
            # w_conv[i+l] = sum(@view w_ext[i:i+2*l])
            for j in 0:2*l
                @inbounds w_conv[i+l] += w_ext[i+j]
            end
        end
        return nothing 
    end

    function cuda_right_conv!(w_conv, w_ext, l)
        n = length(w_ext)-2*l
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = gridDim().x * blockDim().x
        for i = index:stride:l
            for j in 0:l
                #w_conv[n+l+i] = sum(@view w_ext[n+i:n+i+l])
                @inbounds w_conv[n+l+i] += w_ext[n+i+j]
            end
        end
        return nothing 
    end

    function cuda_uniform_convolution!(w_conv, w_ext, l, threads, blocks)
        #@info l
        #n = (length(w_ext)-2*l)
        #@info n
        w_conv .= 0.0

        @cuda threads = l cuda_left_conv!(w_conv, w_ext, l)
        @cuda threads = threads blocks = blocks cuda_central_conv!(w_conv, w_ext, l)
        @cuda threads = l cuda_right_conv!(w_conv, w_ext, l)
        CUDA.synchronize()
        w_conv ./= (2*l+1)
    end

    function LinearAlgebra.mul!(w::CuArray{Float64}, N::NoiseUlamCuda, v::CuArray{Float64})
        k = N.k
        l = (k-1) ÷ 2
        #n = length(v)
        mul!(N.w_ext, N.Ext, v)
        cuda_uniform_convolution!(N.w_conv, N.w_ext, l, N.threads, N.blocks)
        mul!(w, N.BC, N.w_conv)
    end 

    function Base.:*(N::NoiseUlamCuda, v::CuArray{Float64})
        w = CUDA.zeros(Float64, length(v))
        mul!(w, N, v)
        return w
    end
    
    BasisDefinition.opnormbound(B::Ulam, ::Type{L1}, M::NoiseUlamCuda) = 1.0
    opradius(::Type{L1}, N::NoiseUlamCuda) = N.k*Interval(radius(Interval(1)/N.k)).hi
    nonzero_per_row(N::NoiseUlamCuda) = N.k
    dfly(::Type{TotalVariation}, ::Type{L1}, N::NoiseUlamCuda) = (0.0, (1/(2*N.ξ)).hi)


end

###
# Gave a shot at KernelAbstractions, but it is not working for me
###
# using KernelAbstractions

# @kernel function KA_left_conv!(w_conv, @Const(w_ext), l::Int64)
#     I = @index(Global) 
#     for j in 0:l
#         @inbounds w_conv[I] += w_ext[I+j]
#     end
# end

# @kernel function KA_central_conv!(w_conv, @Const(w_ext), l::Int64)
#     I = @index(Global) 
#     for j in 0:2*l
#         @inbounds w_conv[I+l] += w_ext[I+j]
#     end
# end

# @kernel function KA_right_conv!(w_conv, @Const(w_ext), l::Int64)
#     n = length(w_ext)-2*l
#     I = @index(Global) 
#     for j in 0:l
#         #w_conv[n+l+i] = sum(@view w_ext[n+i:n+i+l])
#         @inbounds w_conv[n+l+I] += w_ext[n+I+j]
#     end
# end

# function KA_uniform_convolution!(w_conv, w_ext, l; device = CPU(), wgsize=64)
#     ev1 = KA_left_conv!(device, wgsize)(w_conv, w_ext, 9, ndrange = 9)
#     ev2 = KA_central_conv!(device, wgsize)(w_conv, w_ext, 9, ndrange = 9)
#     ev3 = KA_right_conv!(device, wgsize)(w_conv, w_ext, 9, ndrange = 9)
#     wait(ev1)
#     wait(ev2)
#     wait(ev3)
#     w_conv ./= (2*l+1)
#     return nothing
# end
