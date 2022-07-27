@testset "Abstract Convergence Rates" begin
    ρ, v = RigorousInvariantMeasures.eig_costants_small_matrix([2 0; 0 1])
    @test 2 ∈ ρ
        
end