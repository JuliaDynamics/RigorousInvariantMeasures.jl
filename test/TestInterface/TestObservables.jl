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
    @test 0.125 ∈ Obs.v[1]

    D = mod1_dynamic(x -> 2 * x)

    B = Ulam(1024)
    logder = TMExt.discretizationlogder(B, D)

    @test log(2) ∈ logder.v[5]

    # integrateobservable touches ϕ.inf_bound — guard against the
    # infbound/inf_bound rename drifting again.
    B = Ulam(4)
    Obs = TMExt.Observable(B, x -> x)
    f = [interval(1.0), interval(1.0), interval(1.0), interval(1.0)]
    val = TMExt.integrateobservable(B, Obs, f, 0.0)
    @test 0.5 ∈ val
end

@testset "VariationBound" begin
    # f(x) = x: |f'| ≡ 1, ∫₀¹ 1 dx = 1
    @test 1.0 ∈ TMExt.VariationBound(x -> x)

    # f(x) = x^2: f' = 2x ≥ 0 on [0,1], ∫₀¹ 2x dx = 1
    @test 1.0 ∈ TMExt.VariationBound(x -> x^2)

    # f(x) = 3x - x^2: f' = 3 - 2x ∈ [1,3] on [0,1], ∫₀¹ (3-2x) dx = 2
    @test 2.0 ∈ TMExt.VariationBound(x -> 3x - x^2)

    # Finer partition should still contain the true value
    v = TMExt.VariationBound(x -> x^2; steps = 4096)
    @test 1.0 ∈ v
end

@testset "ProjectedFunction" begin
    B = Ulam(4)

    # f(x) = x: cell averages are 1/8, 3/8, 5/8, 7/8.
    # With VarBound = 1 (true total variation of x ↦ x), err_bound = 1/4.
    pf = TMExt.ProjectedFunction(B, x -> x; VarBound = interval(1.0))

    @test 0.125 ∈ pf.v[1]
    @test 0.375 ∈ pf.v[2]
    @test 0.625 ∈ pf.v[3]
    @test 0.875 ∈ pf.v[4]
    @test 0.25 ∈ pf.err_bound

    # Default VarBound delegates to VariationBound(f; tol = tol).
    # For f(x) = x the total variation is exactly 1, so err_bound ≈ 1/length(B).
    pf_default = TMExt.ProjectedFunction(B, x -> x)
    @test 0.25 ∈ pf_default.err_bound
    @test 0.125 ∈ pf_default.v[1]

    # tol propagates: tighter tol shouldn't break the enclosure.
    pf_tight = TMExt.ProjectedFunction(B, x -> x; tol = 2^-14)
    @test 0.125 ∈ pf_tight.v[1]

    # `projection` is the public entry point; extension routes it to
    # ProjectedFunction. Make sure dispatch actually reaches the extension.
    pf_via_api = projection(B, x -> x; VarBound = interval(1.0))
    @test pf_via_api isa TMExt.ProjectedFunction
    @test 0.125 ∈ pf_via_api.v[1]
end
