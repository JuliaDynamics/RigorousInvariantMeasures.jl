module InvariantMeasures

using IntervalArithmetic: range_atan
abstract type NormKind end
struct L1 <: NormKind end
struct L2 <: NormKind end
struct Linf <: NormKind end
struct Lipschitz <: NormKind end
struct TotalVariation <: NormKind end
struct ℓ1 <: NormKind end
struct ℓinf <: NormKind end

# the module Contractors does not depend on any submodule
include("Contractors.jl")

include("DynamicDefinition.jl")
include("BasisDefinition.jl")
include("GenericAssembler.jl")
include("GenericEstimate.jl")

include("PwDynamicDefinition.jl")
include("Mod1Dynamic.jl")
include("Mod1PwDynamic.jl")
include("IterateDynamic.jl")

include("UlamBasis.jl")
include("HatBasis.jl")

using .DynamicDefinition, .BasisDefinition, .Mod1DynamicDefinition, .Contractors, .PwDynamicDefinition

include("Norms.jl")
include("pitrig.jl")
include("NormsOfPowers.jl")

include("preimages.jl")

include("precompile.jl")

export NormKind, L1, Linf, Lipschitz, TotalVariation

export Ulam, PwMap, Mod1Dynamic, Basis, Dynamic, assemble, preim, Hat,
	EquispacedPartition, norms_of_powers, sinpi, cospi, dfly,
	DiscretizedOperator, IntegralPreservingDiscretizedOperator, NonIntegralPreservingDiscretizedOperator,
	opnormbound, weak_norm, strong_norm, aux_norm, integral_covector, distance_from_invariant,
	mod1_dynamic, Iterate, derivative, distorsion, endpoints, nbranches, branch,
	expansivity, max_distorsion, is_refinement,
	skip_beginning, last_end, preimages

import ValidatedNumerics: Interval
export Interval

include("C2Basis.jl")
using .C2BasisDefinition
export C2Basis
include("ContractionC1.jl")


include("ConvergenceRatesOriginal.jl")
include("HigherDFLY.jl")
include("Chebyshev.jl")

# a special example, the induced map for the LSV map
include("InducedLSV.jl")
using .InducedLSVMapDefinition
export ApproxInducedLSV

include("NoiseKernel.jl")
export UniformNoiseUlam
include("NormsOfPowersNoise.jl")
export powernormboundsnoise, finepowernormboundsnoise, abstractpowernormboundsnoise, invariant_vector_noise, distance_from_invariant_noise

include("Observables.jl")

end
