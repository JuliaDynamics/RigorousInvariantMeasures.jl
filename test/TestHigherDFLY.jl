using Symbolics, RigorousInvariantMeasures

@testset "Test HigherDFLY" begin
    
    #import Symbolics: @variables
    @variables z
    
    P = RigorousInvariantMeasures.SymbL(3, z)

    @test P.n == 3
    @test P.f === z

    ∂ = Differential(z)

    @variables x
    @variables (f(x))[1:3]  # f[0:k](x)
    @variables (DD(x))[1:3] # if DD(x)[k] = (1/T')^(k)

    # this corresponds to the transfer operator Lf = Σ f(y)/(T'(y))
    v = RigorousInvariantMeasures.compute_dfly_k_fi_DDi(0)
    
    #v = RigorousInvariantMeasures.compute_dlfy_k_fi_DDi(0)
    @test v.n == 1
    @test isequal(v.f, f[1])
    
    # its derivative (Lf)' = Σ f'(y)/(T'(y))^2 + Σ f(y)/(T'(y)) * (1/T')'
    v = RigorousInvariantMeasures.compute_dfly_k_fi_DDi(1)
    @test v[1].n == 2
    @test isequal(v[1].f, f[2])
    @test v[2].n == 1
    @test isequal(v[2].f, DD[1]*f[1])

end