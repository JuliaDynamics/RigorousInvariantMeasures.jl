@testset "Abstract Convergence Rates" begin
    ρ, v = RigorousInvariantMeasures.eig_costants_small_matrix([2 0; 0 1])
    @test 2 ∈ ρ

    D = mod1_dynamic(x->2*x)
    B = Ulam(1024)
    Q = DiscretizedOperator(B, D)
    norms = powernormbounds(B, D, Q = Q)
    normsL1, normsBV = convergencerateabstract(B, D, norms)
    @test normsL1[1:10] == [ones(9); 0.0]
    @test normsBV[1:9] == [0.5^i for i in 1:9]

end
