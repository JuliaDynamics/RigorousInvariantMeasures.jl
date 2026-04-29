using RigorousInvariantMeasures
using IntervalArithmetic
using TaylorModels
using Test

# Symbols live in the weakdep extension module; fetch a handle.
const TMExt = Base.get_extension(RigorousInvariantMeasures, :TaylorModelsExt)

@testset "Observables" begin
    B = Ulam(4)

    Obs = TMExt.Observable(B, x -> x)

    @test sup(Obs.weak_dual_bound) >= 1
    @test in_interval(0.125, Obs.v[1])

    D = mod1_dynamic(x -> 2 * x)

    B = Ulam(1024)
    logder = TMExt.discretizationlogder(B, D)

    @test in_interval(log(2), logder.v[5])

    # integrateobservable touches ϕ.weak_dual_bound (was inf_bound) — guard
    # against the field rename drifting again.
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
    # With var_bound = 1 (true total variation of x ↦ x), proj_error = 1/4.
    pf = TMExt.ProjectedFunction(B, x -> x; var_bound = interval(1.0))

    @test in_interval(0.125, pf.v[1])
    @test in_interval(0.375, pf.v[2])
    @test in_interval(0.625, pf.v[3])
    @test in_interval(0.875, pf.v[4])
    @test in_interval(0.25, pf.proj_error)

    # Default var_bound delegates to VariationBound(f; tol = tol).
    # For f(x) = x the total variation is exactly 1, so proj_error ≈ 1/length(B).
    pf_default = TMExt.ProjectedFunction(B, x -> x)
    @test in_interval(0.25, pf_default.proj_error)
    @test in_interval(0.125, pf_default.v[1])

    # tol propagates: tighter tol shouldn't break the enclosure.
    pf_tight = TMExt.ProjectedFunction(B, x -> x; tol = 2^-14)
    @test in_interval(0.125, pf_tight.v[1])

    # `projection` is the public entry point; extension routes it to
    # ProjectedFunction. Make sure dispatch actually reaches the extension.
    pf_via_api = projection(B, x -> x; var_bound = interval(1.0))
    @test pf_via_api isa RigorousInvariantMeasures.ProjectedFunction
    @test in_interval(0.125, pf_via_api.v[1])

    # The unified struct also auto-computes weak_dual_bound (= ‖f‖_{L^∞} for
    # weak-L¹ Ulam). For f(x) = x on [0,1], that's 1.
    @test sup(pf.weak_dual_bound) ≥ 1 - 1e-9

    # ----- Algebra: sum, subtraction, multiplication -----
    pf_x = TMExt.ProjectedFunction(B, x -> x)
    pf_y = TMExt.ProjectedFunction(B, x -> 1 - x)

    sum_pf = pf_x + pf_y
    @test length(sum_pf.v) == length(pf_x.v)
    @test in_interval(1.0, sum_pf.v[1] / 1)        # cell averages of x + (1-x) = 1
    @test in_interval(1.0, sum_pf.v[end] / 1)
    # Triangle inequality: combined bounds are p1.bound + p2.bound (interval
    # arithmetic preserves the upper-bound property).
    @test sup(sum_pf.weak_dual_bound) ≥
          sup(pf_x.weak_dual_bound) + sup(pf_y.weak_dual_bound) - 1e-9
    @test sup(sum_pf.proj_error) ≥
          sup(pf_x.proj_error) + sup(pf_y.proj_error) - 1e-9

    diff_pf = pf_x - pf_y
    @test in_interval(-0.75, diff_pf.v[1] / 1)     # 1/8 - 7/8 = -0.75
    @test sup(diff_pf.weak_dual_bound) ≥
          sup(pf_x.weak_dual_bound) + sup(pf_y.weak_dual_bound) - 1e-9

    # Componentwise multiplication: cell values of x · (1 - x), so v[i] = c_i * (1 - c_i)
    # where c_i is the cell value of x ∈ {1/8, 3/8, 5/8, 7/8}.
    prod_pf = pf_x * pf_y
    @test in_interval((1 / 8) * (7 / 8), prod_pf.v[1])
    @test in_interval((7 / 8) * (1 / 8), prod_pf.v[end])
    # Hölder bound: ‖fg‖_{L^∞} ≤ ‖f‖_∞ · ‖g‖_∞ → weak_dual = product
    @test sup(prod_pf.weak_dual_bound) ≥
          sup(pf_x.weak_dual_bound) * sup(pf_y.weak_dual_bound) * (1 - 1e-9)
end

@testset "Fourier ProjectedFunction (W^{k,1})" begin
    # cos(2πx) has Fourier series ½(e^{2πix} + e^{-2πix}); only modes ±1 are
    # nonzero, with f̂_{±1} = 1/2. ‖f''‖_{L¹} = 4π² · ∫₀¹|cos(2πx)|dx = 4π²·(2/π) = 8π.
    # ‖cos(2πx)‖_{L²} = 1/√2 ≤ 1.
    B = FourierAnalytic(4, 9, RigorousInvariantMeasures.W{2,1})
    f(x) = cos(2 * π * x)
    pf = projection(B, f; weak_dual_bound = 1.0, Wk1_seminorm = 8π)
    @test pf isa RigorousInvariantMeasures.ProjectedFunction
    @test in_interval(0.0, real(pf.v[1]))            # DC = 0
    @test in_interval(0.5, real(pf.v[2]))            # n=1 → ½
    @test in_interval(0.5, real(pf.v[end]))          # n=-1 → ½
    @test pf.proj_error > 0                          # nonzero L² tail bound
    @test pf.weak_dual_bound == 1.0

    # `Observable` is now an alias for `ProjectedFunction`.
    @test RigorousInvariantMeasures.Observable === RigorousInvariantMeasures.ProjectedFunction

    # FFT grid size: a larger `M` shrinks the L²-aliasing tail T₂, so
    # `proj_error = √(T₁² + T₂²)` decreases.
    pf_default_M = projection(B, f; weak_dual_bound = 1.0, Wk1_seminorm = 8π)
    pf_small_M =
        projection(B, f; weak_dual_bound = 1.0, Wk1_seminorm = 8π, M = length(B))
    @test pf_default_M.proj_error < pf_small_M.proj_error
    @test pf_default_M.proj_error > 0
    # Sanity: small-M case has T₂ = T₁, so its proj_error = √2 · T₁; default
    # should be between T₁ and √2·T₁.
    T1_eq_small = pf_small_M.proj_error / sqrt(2)
    @test pf_default_M.proj_error ≥ T1_eq_small * 0.99
    @test pf_default_M.proj_error ≤ pf_small_M.proj_error * 1.01

    # Refusing M < length(B) (would mean fewer samples than basis modes).
    @test_throws ArgumentError projection(
        B,
        f;
        weak_dual_bound = 1.0,
        Wk1_seminorm = 8π,
        M = 3,
    )

    # Legacy 3-arg ProjectedFunction constructor still works
    # (proj_error = nothing).
    legacy = RigorousInvariantMeasures.ProjectedFunction(B, pf.v, 1.0)
    @test legacy.proj_error === nothing

    # Both Fourier subtypes dispatch.
    Bj = RigorousInvariantMeasures.FourierAdjoint(
        FourierPoints(9, Float64),
        4,
        RigorousInvariantMeasures.W{2,1}(),
        RigorousInvariantMeasures.L2(),
    )
    pf_adj = projection(Bj, f; weak_dual_bound = 1.0, Wk1_seminorm = 8π)
    @test pf_adj isa RigorousInvariantMeasures.ProjectedFunction
    @test in_interval(0.5, real(pf_adj.v[2]))

    # k = 1 (BV / W^{1,1}) — total variation of cos(2πx) is 4.
    B1 = FourierAnalytic(4, 9, RigorousInvariantMeasures.W{1,1})
    pf1 = projection(B1, f; weak_dual_bound = 1.0, Wk1_seminorm = 4)
    @test pf1 isa RigorousInvariantMeasures.ProjectedFunction
    @test in_interval(0.5, real(pf1.v[2]))
    @test pf1.proj_error > 0

    # ----- integral_pairing -----
    # Pair cos(2πx) against itself; ∫₀¹ cos²(2πx) dx = 1/2.
    val = integral_pairing(pf, pf)
    @test in_interval(0.5, val)

    # weak_dual_norm_bound of the discrete ρ_N (ℓ² of coefficients).
    @test weak_dual_norm_bound(B, pf.v) ≥ 1 / sqrt(2) - 1e-3

    # Vector form with explicit ρ_w_error and ρ_dual_weak_bound override.
    val2 = integral_pairing(pf, pf.v, pf.proj_error; ρ_dual_weak_bound = 1.0)
    @test in_interval(0.5, val2)
end
