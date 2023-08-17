@testset "Noise 2" begin
    
    import RigorousInvariantMeasures: ExtensionOperator, 
                                      PeriodicBoundaryConditionOperator, 
                                      ReflectingBoundaryConditionOperator,
                                      UniformNoiseUlam2,
                                      uniform_convolution!,
                                      PeriodicBoundaryCondition2!,
                                      ReflectingBoundaryCondition2!, 
                                      mul2!,
                                      cuda_uniform_convolution!,
                                      cuda_uniform_convolution_shared!

    using LinearAlgebra
    
    Ext = ExtensionOperator(1024, 1024)
    v = rand(1024)
    w = Ext*v
    @test all(w .== [zeros(1024); v; zeros(1024)])

    E = Matrix(I, 256, 256)
    Per = [E E E]

    BC = PeriodicBoundaryConditionOperator(256, 16)
    @test all(Per[:, 256-16+1:512+16].== BC)

    v = zeros(256+16*2)

    v[256+16+1] = 1 
    w = BC * v

    @test w[1] == 1

    z = zeros(256)
    PeriodicBoundaryCondition2!(z, v, n = 256, l = 16)
    @test all(z .== w)

    v = rand(256+16*2)
    w = BC * v
    z = zeros(256)
    PeriodicBoundaryCondition2!(z, v, n = 256, l = 16)
    @test all(z .== w)

    BC = PeriodicBoundaryConditionOperator(256, 128)
    @test all(Per[:, 256-128+1:512+128].== BC)

    Erev = reverse(E, dims = 1)
    Ref = [Erev E Erev]
    
    RC = ReflectingBoundaryConditionOperator(256, 16)
    @test all(Ref[:, 256-16+1:512+16].== RC)

    v = zeros(256+16*2)
    v[256+16+1] = 1 
    w = RC * v
    @test w[end] == 1
    z = zeros(256)
    ReflectingBoundaryCondition2!(z, v, n = 256, l = 16)
    @test all(z .== w)

    v = rand(256+16*2)
    w = RC * v
    z = zeros(256)
    ReflectingBoundaryCondition2!(z, v, n = 256, l = 16)
    @test all(z .== w)


    RC = ReflectingBoundaryConditionOperator(256, 128)
    @test all(Ref[:, 256-128+1:512+128].== RC)

    l = 9
    w_ext = zeros(1042)
    w_ext[10] = 1
    w_conv = zeros(1042)
    uniform_convolution!(w_conv, w_ext, l)

    @test all(w_conv[1:19] .== Float64(1)/19)  

    w_ext = zeros(1042)
    w_ext[1024+l] = 1
    w_conv = zeros(1042)
    uniform_convolution!(w_conv, w_ext, l)
    @test all(w_conv[1024:end] .== Float64(1)/19)  

    n = 1024
    B = Ulam(n)

    k = 19
    l = (k-1)÷2
    N = UniformNoiseUlam2(B, k)

    v = zeros(length(B))
    v[1] = 1
    w = N*v
    @test all(w[1:l+1].==1/k) && all(w[end-l+1:end].== 1/k)

    N = UniformNoiseUlam2(B, k, boundary_condition = :reflecting)

    v = zeros(length(B))
    v[1] = 1
    w = N*v
    @test all(w[1:l].==2/k) && w[l+1] == 1/k

    z = zeros(length(B))
    mul2!(z, N, v)
    @test all(z .== w)

    using CUDA

    if has_cuda() && has_cuda_gpu()
        l = 9
        w_ext = rand(1042)
        w_conv = zeros(1042)
        uniform_convolution!(w_conv, w_ext, l)
        
        w_ext_d = CuArray{Float64}(w_ext)
        w_conv_d = CUDA.zeros(Float64, 1042)
        cuda_uniform_convolution!(w_conv_d, w_ext_d, l, 1024, 1)
        
        z = w_conv-Array(w_conv_d)

        @test all(z .== 0.0)
    
        # cuda_uniform_convolution_shared!(w_conv_d, w_ext_d, l, 1024, 1)
        
        # z = w_conv-Array(w_conv_d)

        # @test all(z .== 0.0)
    
    end

end