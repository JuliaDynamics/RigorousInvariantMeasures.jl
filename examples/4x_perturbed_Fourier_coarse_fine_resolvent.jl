# Comparison of two coarse-fine certification strategies for the smooth
# circle map T(x) = 2x + 0.5·x·(1-x)  (mod 1).
#
# 1. **Power-norm strategy** (existing): `finepowernormbounds` propagates a
#    bound on ‖L^n|_{U⁰}‖_w from a coarse discretization to a fine one.
# 2. **Resolvent strategy** (new, from the DFLY appendix of
#    [Nisoli, "Certified spectral approximation of transfer operators and
#    the Gauss map"]):
#      a. Certify the L²-resolvent of the coarse discretization on a circle
#         centered at 0 that numerically separates the eigenvalue 1 from
#         the rest of the spectrum, using `BallArithmetic.CertifScripts`.
#      b. Lift to a strong-resolvent bound K(z) for the infinite-dimensional
#         L via `strong_resolvent_lift` (Proposition A.7).
#      c. Bound the L²-operator-norm of the *fine* discretization with
#         `BallArithmetic.upper_bound_L2_opnorm`.
#      d. Apply `coarse_fine_weak_resolvent` (Proposition A.14) to obtain
#         the fine L²-resolvent bound on the same contour.
#
# We use `FourierAdjoint` — the adjoint representation gives the same
# basis-level constants as `FourierAnalytic` (the basis is the same; only
# the matrix is transposed), and the W^{k,1} methods are now wired up for
# both.

using RigorousInvariantMeasures
using FFTW                         # triggers FFTWExt for Fourier assembly
using BallArithmetic
using IntervalArithmetic
using LinearAlgebra: eigvals

const CS = BallArithmetic.CertifScripts

function runComparison(; n_grid = 16384, k_coarse = 128, k_fine = 512,
                        samples = 256,
                        N_candidates = [5, 8, 10, 12, 15, 20, 25, 30, 40, 50],
                        ρ_candidates = [0.85, 0.9, 0.95],
                        # Small circle around eigenvalue 1 for the Riesz
                        # projector / fixed-point distance bound.
                        projector_radii = [0.05, 0.1, 0.15],
                        projector_N = 15)
    D = mod1_dynamic(x -> 2 * x + 0.5 * x * (1 - x))

    # --- Bases (W^{1,1} strong, L² weak, L¹ aux) -------------------------
    B_coarse = FourierAdjoint(FourierPoints(n_grid, Float64), k_coarse,
                                W{1,1}(), L2())
    B_fine = FourierAdjoint(FourierPoints(n_grid, Float64), k_fine,
                              W{1,1}(), L2())

    a, b = dfly(strong_norm(B_coarse), aux_norm(B_coarse), D)
    @info "DFLY constants  (‖Lf‖_s ≤ a·‖f‖_s + b·‖f‖_aux)" a b

    # --- Discretize ------------------------------------------------------
    time_coarse = @elapsed Q_coarse = DiscretizedOperator(B_coarse, D)
    time_fine = @elapsed Q_fine = DiscretizedOperator(B_fine, D)
    @info "Assembled" time_coarse time_fine

    # --- Strategy 1: existing power-norm coarse-fine ---------------------
    time_norms_coarse = @elapsed norms_coarse = powernormbounds(B_coarse, D; Q = Q_coarse)
    normQ_fine_pow = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
    time_norms_fine_pow = @elapsed norms_fine_pow =
        finepowernormbounds(B_coarse, B_fine, D, norms_coarse;
                            normQ_fine = normQ_fine_pow)

    # Empirical mixing rate from the bound sequence: ρ̂ ≈ ‖L_fine^n‖^{1/n}
    n_pow = length(norms_fine_pow)
    rho_powernorm = norms_fine_pow[end]^(1 / n_pow)
    @info "Power-norm coarse-fine" first_term = norms_fine_pow[1] last_term =
        norms_fine_pow[end] sequence_length = n_pow rho_powernorm time_norms_fine_pow

    # --- Strategy 2: resolvent coarse-fine -------------------------------
    # (a) Numerical eigenvalues of the coarse operator. Sanity check on the
    # spectral structure: largest non-trivial modulus tells us how close to 1
    # we can place a contour before R_w_coarse explodes.
    BM_coarse = BallMatrix(Q_coarse.L)
    ev_coarse = eigvals(BM_coarse.c)
    idx_one = argmin(abs.(ev_coarse .- 1.0))
    second_largest = maximum(abs(ev_coarse[i])
                              for i in eachindex(ev_coarse) if i != idx_one)
    @info "Coarse spectral data" second_largest

    # (b) Coarse-level Schur build (shared across contour candidates).
    time_schur = @elapsed schur_data = CS.compute_schur_and_error(BM_coarse)

    # (c) Bound ‖L_fine‖_{L²} via BallArithmetic (independent of contour).
    BM_fine = BallMatrix(Q_fine.L)
    time_opnorm = @elapsed M_L2_fine = upper_bound_L2_opnorm(BM_fine)
    @info "Fine ‖L_fine‖_{L²} (BallArithmetic)" M_L2_fine time_opnorm

    # (d) Sweep contour radius × N. Prop A.14 closure is sensitive to both:
    # ρ closer to 1 shrinks the strong-resolvent K(z) inflation, and N balances
    # a^N · E_{k,w→s} against |z|^{-N} b_N.
    contour_results = Vector{NamedTuple}()
    for ρ in ρ_candidates
        if ρ <= second_largest
            @info "Skipping ρ ≤ second-largest" ρ second_largest
            continue
        end
        circle = CS.CertificationCircle(0.0 + 0.0im, ρ; samples = samples)
        cert = try
            CS.run_certification(BM_coarse, circle; schur_data = schur_data)
        catch err
            @warn "CertifScripts failed at this ρ" ρ err
            nothing
        end
        cert === nothing && continue
        R_w_coarse = cert.resolvent_original
        K_z = strong_resolvent_lift(B_coarse, D, ρ, R_w_coarse)
        # Sweep N. Honest constants — M is the package default
        # `bound_weak_norm_abstract(B, D)` = b + 1, not a manually-asserted
        # tighter value. Since M > 1 here, b_N = b·(M^N − a^N)/(M − a) grows
        # exponentially as O(M^N), and there is at most a small window of N
        # where (a/|z|)^N · E_{k,w→s} balances (M/|z|)^N · b/(M-a).
        for N in N_candidates
            R_w_fine = coarse_fine_weak_resolvent(
                B_fine, D, ρ, K_z;
                N = N,
                M_k_weak = M_L2_fine,
            )
            push!(contour_results,
                  (ρ = ρ, N = N, R_w_coarse = R_w_coarse, K_z = K_z,
                   R_w_fine = R_w_fine))
        end
        @info "Coarse cert at ρ" ρ R_w_coarse K_z
    end

    # Pick the contour with the smallest finite R_w_fine, if any.
    finite_results = filter(r -> isfinite(r.R_w_fine), contour_results)
    best = isempty(finite_results) ? nothing : argmin(r -> r.R_w_fine, finite_results)
    if best === nothing
        @warn "Prop A.14 did not close at any (ρ, N) — see b_N analysis below."
    end

    # --- b_N diagnostic --------------------------------------------------
    # Print the structure of the perturbation factor at each (ρ, N) so the
    # role of b_N is explicit. The package's DFLY is on (s, aux), so the
    # M inside b_N is the aux-norm contractivity (M_aux = 1 for the L¹-
    # preserving transfer operator) and the b is rescaled by M_2 =
    # aux_weak_bound(B). With M_aux = 1, b_N saturates at b·M_2/(1−a).
    a, b = dfly(strong_norm(B_coarse), aux_norm(B_coarse), D)
    M_aux = 1.0
    M_2 = RigorousInvariantMeasures.aux_weak_bound(B_fine)
    b_eff = b * M_2
    E_kws = RigorousInvariantMeasures.strong_weak_bound(B_fine)
    Δ_k = RigorousInvariantMeasures.weak_projection_error(B_fine)
    @info "Diagnostic constants (package convention)" a b M_aux M_2 b_eff E_kws Δ_k

    # --- Summary ---------------------------------------------------------
    println()
    println("─"^72)
    println("Standing constants (basis interface, k_fine = $(k_fine)):")
    println("  a = $(round(a, sigdigits=4))   b = $(round(b, sigdigits=4))   ",
            "M_aux = $(round(M_aux, sigdigits=4))   M_2 = $(round(M_2, sigdigits=4))")
    println("  effective b in b_N: b·M_2 = $(round(b_eff, sigdigits=4))")
    println("  E_{k,w→s} = strong_weak_bound = $(round(E_kws, sigdigits=4))")
    println("  Δ_k       = weak_projection_error = $(round(Δ_k, sigdigits=4))")
    println("─"^72)
    println("Sweep of (ρ, N) for Prop A.14:")
    println("  ρ      N    b_N           a^N·E_kws    R_w_fine")
    for r in contour_results
        bN_val = bN_constant(a, b_eff, M_aux, r.N)
        aN_Ekws = a^r.N * E_kws
        rs = isfinite(r.R_w_fine) ?
            string(round(r.R_w_fine, sigdigits = 3)) : "Inf"
        println("  $(r.ρ)   $(lpad(r.N, 2))   ",
                "$(rpad(round(bN_val, sigdigits=3), 13)) ",
                "$(rpad(round(aN_Ekws, sigdigits=3), 12)) $rs")
    end
    println("─"^72)
    println("Note: package DFLY is on (strong, aux); iterating in aux uses")
    println("M_aux = $(M_aux) (L¹-preserving), so b_N = b·M_2·(1−a^N)/(1−a)")
    println("saturates at $(round(b_eff/(1-a), sigdigits=4)). a^N·E_kws decays")
    println("as a^N; β̃_N = |z|^{-N}·Δ_k·K_z·(a^N·E_kws + b_N) closes at finite N.")
    println("─"^72)
    println("Comparison of certified spectral information on L_fine (k = $(k_fine)):")
    println("─"^72)
    println("Power-norm strategy:")
    println("  ρ̂ = ‖L_fine^n‖^{1/n} = $(round(rho_powernorm, sigdigits=4))  (n = $(n_pow))")
    println("  $(rho_powernorm < 1 ? "certifies spectral radius < 1 on U⁰." :
                                       "does NOT certify a spectral gap (ρ̂ ≥ 1).")")
    println()
    println("Resolvent strategy (Prop A.7 + A.14):")
    if best === nothing
        println("  Prop A.14 did not close at any (ρ, N) candidate at honest M.")
        println("  Cor. A.17 / Thm A.18: closure as k_fine → ∞ since")
        println("  δ_k · E_{k,w→s}^q → 0 polynomially for W^{1,1}/L² (q < 1).")
        println("  At this k_fine = $(k_fine) we are not yet in the closure regime.")
    else
        println("  Best: ρ = $(best.ρ),  N = $(best.N),  ",
                "‖(zI - L_fine)^{-1}‖_{L²} ≤ $(round(best.R_w_fine, sigdigits=3))")
        println("  Certifies that $(best.ρ) is in ρ(L_fine), hence")
        println("  spectral radius of L_fine on U⁰ ≤ $(best.ρ).")
    end
    println("─"^72)

    # --- Riesz projector / fixed-point distance --------------------------
    # Small circle Γ = {|z-1| = r}. CertifScripts is run on the coarse 128
    # matrix only (small enough to Schur cheaply); Prop A.7 lifts the
    # coarse weak certificate to K(z) = R_s(z, L), and Prop A.14 propagates
    # to R_w(z, L_fine) on Γ. Then Lemma A.11 + contour integral gives
    #   ‖P_L − P_{L_fine}‖_{s→w} ≤ r · sup_Γ R_w(z, L_fine) · δ_fine · K(z).
    println("Riesz projector distance on Γ = {|z-1| = r} (CertifScripts on coarse,")
    println("Prop A.7+A.14 propagation, Lemma A.11 + contour):")
    println("  r       R_w_coarse  K_z        R_w_fine     ‖P_L − P_{L_fine}‖_{s→w}")
    projector_results = NamedTuple[]
    for r in projector_radii
        pres = projector_distance_bound(
            B_coarse, B_fine, D, Q_coarse, Q_fine, r;
            samples = samples, N = projector_N,
        )
        push!(projector_results, pres)
        pd_str = isfinite(pres.projector_distance) ?
                 string(round(pres.projector_distance, sigdigits = 3)) : "Inf"
        rw_str = isfinite(pres.R_w_fine) ?
                 string(round(pres.R_w_fine, sigdigits = 3)) : "Inf"
        println("  $(r)    $(rpad(round(pres.R_w_coarse, sigdigits=3), 11)) ",
                "$(rpad(round(pres.K_z, sigdigits=3), 10)) ",
                "$(rpad(rw_str, 11)) $pd_str")
    end
    println("─"^72)
    println("For an L¹-preserving transfer operator with simple eigenvalue 1,")
    println("‖P_L − P_{L_fine}‖_{s→w} = ‖v − v_fine‖_w · ‖∫‖_{s → ℝ},")
    println("so the column above is (up to ‖∫‖, typically ≤ 1) an upper bound")
    println("on the L² distance between the true invariant measure v and the")
    println("discrete invariant measure v_fine of L_fine.")
    println("─"^72)
    return (
        rho_powernorm = rho_powernorm,
        norms_fine_pow = norms_fine_pow,
        contour_results = contour_results,
        best = best,
        M_L2_fine = M_L2_fine,
        projector_results = projector_results,
    )
end

runComparison()
