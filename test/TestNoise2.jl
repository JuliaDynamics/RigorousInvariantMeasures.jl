@testset "Noise 2" begin
    
    import RigorousInvariantMeasures: ExtensionOperator, 
                                      PeriodicBoundaryConditionOperator, 
                                      ReflectingBoundaryConditionOperator,
                                      UniformNoiseUlam2

    using LinearAlgebra
    
    Ext = ExtensionOperator(1024, 1024)
    v = rand(1024)
    w = Ext*v
    @test all(w .== [zeros(1024); v; zeros(1024)])

    E = Matrix(I, 256, 256)
    Per = [E E E]

    BC = PeriodicBoundaryConditionOperator(256, 16)
    @test all(Per[:, 256-16+1:512+16].== BC)

    BC = PeriodicBoundaryConditionOperator(256, 128)
    @test all(Per[:, 256-128+1:512+128].== BC)

    Erev = reverse(E, dims = 1)
    Ref = [Erev E Erev]
    
    RC = ReflectingBoundaryConditionOperator(256, 16)
    @test all(Ref[:, 256-16+1:512+16].== RC)

    RC = ReflectingBoundaryConditionOperator(256, 128)
    @test all(Ref[:, 256-128+1:512+128].== RC)

    N = UniformNoiseUlam(Ulam(32), 5)

    



end