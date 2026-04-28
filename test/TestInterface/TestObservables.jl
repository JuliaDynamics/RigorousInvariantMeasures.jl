using RigorousInvariantMeasures
using IntervalArithmetic
using TaylorModels
using Test

# Symbols live in the weakdep extension module; fetch a handle.
const TMExt = Base.get_extension(RigorousInvariantMeasures, :TaylorModelsExt)

@testset "Observables" begin
    B = Ulam(4)

    Obs = TMExt.Observable(B, x -> x)

    @test sup(Obs.inf_bound) >= 1
    @test in_interval(0.125, Obs.v[1])

    D = mod1_dynamic(x -> 2 * x)

    B = Ulam(1024)
    logder = TMExt.discretizationlogder(B, D)

    @test in_interval(log(2), logder.v[5])

    # integrateobservable touches ϕ.inf_bound — guard against the
    # infbound/inf_bound rename drifting again.
    B = Ulam(4)
    Obs = TMExt.Observable(B, x -> x)
    f = [interval(1.0), interval(1.0), interval(1.0), interval(1.0)]
    val = TMExt.integrateobservable(B, Obs, f, 0.0)
    @test in_interval(0.5, val)
end

@testset "VariationBound" begin
    # f(x) = x: |f'| ≡ 1, ∫₀¹ 1 dx = 1
    @test in_interval(1.0, TMExt.VariationBound(x -> x))

    # f(x) = x^2: f' = 2x ≥ 0 on [0,1], ∫₀¹ 2x dx = 1
    @test in_interval(1.0, TMExt.VariationBound(x -> x^2))

    # f(x) = 3x - x^2: f' = 3 - 2x ∈ [1,3] on [0,1], ∫₀¹ (3-2x) dx = 2
    @test in_interval(2.0, TMExt.VariationBound(x -> 3x - x^2))

    # Finer partition should still contain the true value
    v = TMExt.VariationBound(x -> x^2; steps = 4096)
    @test in_interval(1.0, v)
end

@testset "ProjectedFunction" begin
    B = Ulam(4)

    # f(x) = x: cell averages are 1/8, 3/8, 5/8, 7/8.
    # With VarBound = 1 (true total variation of x ↦ x), err_bound = 1/4.
    pf = TMExt.ProjectedFunction(B, x -> x; VarBound = interval(1.0))

    @test in_interval(0.125, pf.v[1])
    @test in_interval(0.375, pf.v[2])
    @test in_interval(0.625, pf.v[3])
    @test in_interval(0.875, pf.v[4])
    @test in_interval(0.25, pf.err_bound)

    # Default VarBound delegates to VariationBound(f; tol = tol).
    # For f(x) = x the total variation is exactly 1, so err_bound ≈ 1/length(B).
    pf_default = TMExt.ProjectedFunction(B, x -> x)
    @test in_interval(0.25, pf_default.err_bound)
    @test in_interval(0.125, pf_default.v[1])

    # tol propagates: tighter tol shouldn't break the enclosure.
    pf_tight = TMExt.ProjectedFunction(B, x -> x; tol = 2^-14)
    @test in_interval(0.125, pf_tight.v[1])

    # `projection` is the public entry point; extension routes it to
    # ProjectedFunction. Make sure dispatch actually reaches the extension.
    pf_via_api = projection(B, x -> x; VarBound = interval(1.0))
    @test pf_via_api isa RigorousInvariantMeasures.ProjectedFunction
    @test in_interval(0.125, pf_via_api.v[1])
end

@testset "Fourier Observable / ProjectedFunction (W^{k,1})" begin
    # cos(2πx) has Fourier series ½(e^{2πix} + e^{-2πix}); only modes ±1 are
    # nonzero, with f̂_{±1} = 1/2. ‖f''‖_{L¹} = 4π² · ∫₀¹|cos(2πx)|dx = 4π²·(2/π) = 8π.
    B = FourierAnalytic(4, 9, RigorousInvariantMeasures.W{2,1})
    f(x) = cos(2 * π * x)
    pf = projection(B, f; Wk1_seminorm = 8π)
    @test pf isa RigorousInvariantMeasures.ProjectedFunction
    @test in_interval(0.0, real(pf.v[1]))            # DC = 0
    @test in_interval(0.5, real(pf.v[2]))            # n=1 → ½
    @test in_interval(0.5, real(pf.v[end]))          # n=-1 → ½
    @test pf.err_bound > 0                            # nonzero L² tail bound

    # Observable interface — same coefficients, plus user-supplied
    # `inf_bound` (here ‖cos(2πx)‖_{L²} = 1/√2 ≤ 1) and the new
    # `proj_error` field (L² weak-norm projection error of ϕ itself).
    obs = RigorousInvariantMeasures.Observable(
        B,
        f;
        inf_bound = 1.0,
        Wk1_seminorm = 8π,
    )
    @test obs isa RigorousInvariantMeasures.Observable
    @test in_interval(0.5, real(obs.v[2]))
    @test obs.inf_bound == 1.0
    @test obs.proj_error > 0
    @test obs.proj_error == pf.err_bound  # same formula as ProjectedFunction

    # FFT grid size: a larger `M` divides the per-coefficient aliasing
    # inflation by (M_new/M_old)^k, so the resulting interval widths
    # should shrink (the truncation tail / err_bound is unaffected since
    # it depends only on k_freq, not on the FFT grid size).
    pf_default_M = projection(B, f; Wk1_seminorm = 8π)  # M = 4*B.k = 16
    pf_small_M = projection(B, f; Wk1_seminorm = 8π, M = length(B))
    @test diam(real(pf_default_M.v[2])) < diam(real(pf_small_M.v[2]))
    @test pf_default_M.err_bound == pf_small_M.err_bound

    # Refusing M < length(B) (would mean fewer samples than basis modes).
    @test_throws ArgumentError projection(B, f; Wk1_seminorm = 8π, M = 3)

    # Legacy 3-arg Observable constructor still works (proj_error = nothing).
    legacy = RigorousInvariantMeasures.Observable(B, pf.v, 1.0)
    @test legacy.proj_error === nothing

    # Both Fourier subtypes dispatch.
    Bj = RigorousInvariantMeasures.FourierAdjoint(
        FourierPoints(9, Float64),
        4,
        RigorousInvariantMeasures.W{2,1}(),
        RigorousInvariantMeasures.L2(),
    )
    pf_adj = projection(Bj, f; Wk1_seminorm = 8π)
    @test pf_adj isa RigorousInvariantMeasures.ProjectedFunction
    @test in_interval(0.5, real(pf_adj.v[2]))

    # Refusing k = 1 (logarithmic-divergence regime, not in this commit).
    B1 = FourierAnalytic(4, 9, RigorousInvariantMeasures.W{1,1})
    @test_throws ArgumentError projection(B1, f; Wk1_seminorm = 8π)
end
