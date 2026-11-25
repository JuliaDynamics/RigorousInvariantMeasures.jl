@testset "Full run FourierAnalytic + Aη" begin
    using RigorousInvariantMeasures

    D = mod1_dynamic(x -> 2 * x)

    B = FourierAnalytic(64, 256; η = 0.1)
    Q = DiscretizedOperator(B, D)

    norms = powernormbounds(B, D; Q = Q)
    @test norms[end] < 1e-3

    w = invariant_vector(B, Q)

    error = distance_from_invariant(B, D, Q, w, norms)
    # Guards two bugs that used to cancel: the placeholder Aη dfly returned
    # B = 0, so invariant_measure_strong_norm_bound was 0, which hid a sign error
    # in weak_projection_error (it returned e^{+2πηN} instead of e^{-2πηN}).
    @test error < 1.0  # Lebesgue measure is invariant for doubling map
    @test error < 1e-8
end

@testset "FourierAnalytic + W{1,1}: no projection error" begin
    using RigorousInvariantMeasures

    # The Fourier projection error is derived per unit of ||f^(s)||_1 for s >= 2 only: for
    # s = 1 the aliased sum Σ_{l≠0} |m + ln|^{-1} diverges. Until 2026-11-25 this basis
    # returned (2πK)^{-1} from an inequality that does not hold, and this test ran the full
    # pipeline on it.
    B = FourierAnalytic(64, 256, W{1,1})
    @test_throws ArgumentError RigorousInvariantMeasures.weak_projection_error(B)
end
