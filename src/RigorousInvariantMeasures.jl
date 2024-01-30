module RigorousInvariantMeasures

using IntervalArithmetic: range_atan


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
include("FFT.jl")
include("Basis/NewChebyshev.jl")
export Chebyshev


include("precompile.jl")

export NormKind, L1, Linf, Lipschitz, TotalVariation

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

import IntervalArithmetic: Interval, @interval, @biginterval, midpoint_radius
using IntervalArithmetic
export Interval

include("Basis/C2Basis.jl")
export C2Basis



include("ConvergenceRatesOriginal.jl")
export convergencerateabstract

include("HigherDFLY.jl")


# a special example, the induced map for the LSV map
include("InducedLSV.jl")
export ApproxInducedLSV
include("PwMapInducedLSV.jl")

include("NoiseKernel.jl")
export UniformNoiseUlam

include("NoiseKernel2.jl")
export UniformNoiseUlam2

include("NormsOfPowersNoise.jl")
export powernormboundsnoise,
    finepowernormboundsnoise,
    abstractpowernormboundsnoise,
    invariant_vector_noise,
    distance_from_invariant_noise

include("Observables.jl")
export Observable, discretizationlogder, integrateobservable

include("Basis/BasisIndex.jl")

#include("Lorenz2DUlam.jl")


#include("SkewProductMapDefinition.jl")
#include("Basis/BasisUlam2DSkewProduct.jl")

include("sample_dynamics.jl")


#include("FourierCommon.jl")
#include("FourierBasis.jl")
#export L2, Fourier1D, GaussianNoise
#include("FourierAnalytic.jl")
#include("FourierAdjoint.jl")
#include("FourierBasisBack.jl")


end
