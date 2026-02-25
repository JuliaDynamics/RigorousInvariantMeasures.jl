@testset "Full run FourierAnalytic + Aη" begin
    using RigorousInvariantMeasures

    D = mod1_dynamic(x -> 2 * x)

    B = FourierAnalytic(64, 256; η = 0.1)
    Q = DiscretizedOperator(B, D)

    norms = powernormbounds(B, D; Q = Q)
    @test norms[end] < 1e-3

    w = invariant_vector(B, Q)

    error = distance_from_invariant(B, D, Q, w, norms)
    @test error < 1.0  # Lebesgue measure is invariant for doubling map
end

@testset "Full run FourierAnalytic + W{1,1}" begin
    using RigorousInvariantMeasures

    D = mod1_dynamic(x -> 2 * x)

    B = FourierAnalytic(64, 256, W{1,1})
    Q = DiscretizedOperator(B, D)

    norms = powernormbounds(B, D; Q = Q)
    @test norms[end] < 1e-3

    w = invariant_vector(B, Q)

    error = distance_from_invariant(B, D, Q, w, norms)
    @test error < 1.0
end
