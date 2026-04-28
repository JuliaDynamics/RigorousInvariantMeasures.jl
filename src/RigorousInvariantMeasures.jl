module RigorousInvariantMeasures

import IntervalArithmetic: Interval, mid, radius, @interval
using IntervalArithmetic: setdisplay
# Suppress IA 1.0 decoration suffixes (_com, _trv_NG, …) from interval display.
# The NG flag is raised for scalar * Interval (e.g. 2*x) because IA 1.0 does not
# recognise plain Julia scalar multiplication as directed-rounding-safe; this is a
# known conservatism that is being discussed upstream.  The numerical enclosures are
# still correct; we suppress the cosmetic noise until IA resolves the issue.
#
# `setdisplay` mutates global state in IntervalArithmetic. When run at the top level
# of this module the call happens during precompilation but is not replayed when
# the cached module is loaded — so doctests (and any other downstream user)
# would see the decorated form. Calling it from `__init__` runs the side effect
# every time the module is brought into a fresh session.
function __init__()
    setdisplay(:infsup; decorations = false, ng_flag = false)
end
using BallArithmetic: BallMatrix, BallVector, upper_bound_L2_opnorm, upper_bound_norm,
    compute_spectral_projector_schur, SchurSpectralProjectorResult
using BallArithmetic.CertifScripts: CertifScripts

const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

include("Norms.jl")
include("differentiation_interface.jl")

include("Contractors.jl")

include("AbstractDynamicDefinition.jl")
export Dynamic, endpoints, nbranches, branch, max_inverse_derivative, max_distortion

include("Basis/BasisDefinition.jl")
export opnormbound, weak_norm, strong_norm, aux_norm, integral_covector
include("NormBounds.jl")
include("NormCacher.jl")


include("GenericAssembler.jl")
export DiscretizedOperator,
    IntegralPreservingDiscretizedOperator, NonIntegralPreservingDiscretizedOperator

include("GenericEstimate.jl")
export invariant_vector, finepowernormbounds, powernormbounds, distance_from_invariant

include("PwDynamic.jl")
export PwMap, mod1_dynamic
include("DFLY.jl")

include("Basis/UlamBasis.jl")
export Ulam
include("Basis/CircleHatBasis.jl")
export Hat
include("Basis/IntervalHatBasis.jl")
export HatNP

include("pitrig.jl")
include("NormsOfPowers.jl")

include("Preimages.jl")
include("Basis/Fourier/FourierIndex.jl")
export Fourier, FourierAnalytic, FourierAdjoint
include("Basis/NewChebyshev.jl")
export Chebyshev, certify_spectral_gap


include("precompile.jl")

export NormKind, L1, L2, Linf, Lipschitz, TotalVariation, C1, W, Aη, Cω

export PwMap,
    Basis,
    assemble,
    preim,
    Hat,
    Ulam,
    EquispacedPartition,
    norms_of_powers,
    sinpi,
    cospi,
    dfly,
    distance_from_invariant,
    mod1_dynamic,
    is_refinement

import IntervalArithmetic: Interval, @interval, @biginterval, midradius
using IntervalArithmetic
export Interval

include("Basis/C2Basis.jl")
export C2Basis



include("ConvergenceRatesOriginal.jl")
export convergencerateabstract

#include("HigherDFLY.jl")


# a special example, the induced map for the LSV map
include("InducedLSV.jl")
export ApproxInducedLSV
include("PwMapInducedLSV.jl")

include("NoiseKernel.jl")
export UniformNoiseUlam

include("UniformNoiseUlam.jl")
export UniformKernelUlamPeriodic, UniformKernelUlamReflecting, UniformKernelUlam

# Function stub for plotting extension
"""
    plot_noisy_system(D::PwMap, K::UniformKernelUlam, w, L; n_samples=1000, n_points=500)

Plot a 3-panel visualization of a dynamical system with uniform noise.
Requires `using Plots` to load the extension that implements this function.

See the PlotsExt extension for full documentation.
"""
function plot_noisy_system end
export plot_noisy_system

# include("NoiseKernel2.jl")
# export UniformNoiseUlam2

include("NormsOfPowersNoise.jl")
export powernormboundsnoise,
    finepowernormboundsnoise,
    abstractpowernormboundsnoise,
    invariant_vector_noise,
    distance_from_invariant_noise

include("NoiseSpecializedEstimate.jl")
export noise_error_aposteriori,
    noise_error_apriori,
    prepare_Wnorm_estimate,
    prepare_derivative_bounds,
    total_variation


include("Basis/BasisIndex.jl")

#include("Lorenz2DUlam.jl")


#include("SkewProductMapDefinition.jl")
#include("Basis/BasisUlam2DSkewProduct.jl")

include("sample_dynamics.jl")




end
