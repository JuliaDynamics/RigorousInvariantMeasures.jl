@testset "Residual-based a posteriori estimate" begin
    using RigorousInvariantMeasures
    using RigorousInvariantMeasures: projection_defect_coefficients
    using IntervalArithmetic: interval, mid

    # (1) soundness against a known invariant density: 3x mod 1, h ≡ 1
    B = Ulam(256)
    D = mod1_dynamic(x -> 3x)
    Q = DiscretizedOperator(B, D)
    norms = powernormbounds(B, D; Q = Q)
    w = invariant_vector(B, Q)
    err_res = distance_from_invariant_residual(B, D, Q, w, norms)
    true_err = sum(abs.(w .- 1.0)) / length(B)
    @test err_res >= true_err
    @test err_res < 1e-4

    # (2) known density with a DECREASING branch: tent map, h ≡ 1
    Dt = PwMap(
        [x -> 2x, x -> 2 - 2x],
        [interval(0), interval(0.5), interval(1)],
        [0 1; 1 0],
    )
    Qt = DiscretizedOperator(B, Dt)
    norms_t = powernormbounds(B, Dt; Q = Qt)
    wt = invariant_vector(B, Qt)
    err_res_t = distance_from_invariant_residual(B, Dt, Qt, wt, norms_t)
    true_err_t = sum(abs.(wt .- 1.0)) / length(B)
    @test err_res_t >= true_err_t

    # (3) non-full-branch map (B_dfly > 0): finite, valid, comparable to a priori bound.
    # Needs a fine enough grid: the closed-form denominator 1 − Kh·ΣC·(c_s B/(1−A) + c_w)
    # must be positive (at n=256 the function correctly refuses with an error).
    B3 = Ulam(4096)
    Dn = mod1_dynamic(x -> 2.5x; full_branch = false)
    Qn = DiscretizedOperator(B3, Dn)
    norms_n = powernormbounds(B3, Dn; Q = Qn)
    wn = invariant_vector(B3, Qn)
    err_res_n = distance_from_invariant_residual(B3, Dn, Qn, wn, norms_n)
    err_std_n = distance_from_invariant(B3, Dn, Qn, wn, norms_n)
    @test isfinite(err_res_n) && err_res_n > 0
    # the a posteriori bound should not be drastically worse than the a priori one
    @test err_res_n < 10 * err_std_n

    # (4) Ulam defect coefficients are the sharp ones (1+A, B)
    dfc = RigorousInvariantMeasures.dfly(
        RigorousInvariantMeasures.strong_norm(B),
        RigorousInvariantMeasures.aux_norm(B),
        Dn,
    )
    cs, cw = projection_defect_coefficients(B, Dn; dfly_coefficients = dfc)
    @test mid(cs) ≈ 1 + dfc[1]
    @test mid(cw) ≈ dfc[2]
end
