# Salvaged from ibis/FourierToby, commit 5308cdd, as src/FourierToby.jl.
# Not part of the package; see attic/README.md.

using ..BasisDefinition, ..DynamicDefinition


struct C1 <: NormKind end
struct W{k, l} <:NormKind end
order(::Type{W{k, l}}) where {k, l} = k
regularity(::Type{W{k, l}}) where {k, l} = l

struct Fourier2{T<:AbstractVector} <:Basis
	p::T
	k::Integer
	# TODO: check in constructor that p is sorted, starts with 0 and ends with 1
end

Fourier2Points(n, T) = [i/(n) for i=0:n-1]

function Fourier2(n::Integer, k::Integer; T = Float64) 
    return Fourier2(Fourier2Points(n, T), k)
end
Base.show(io::IO, B::Fourier2) = print(io, "Fourier basis on $(length(B)) points, highest frequency $((length(B)-1)/2)")

"""
Return the size of the Fourier2 basis
"""
Base.length(B::Fourier2) = length(B.p)

##########################################################
# These are the methods to rigorously enclose the image of a Chebyshev polynomial
# they need refactoring, but not now
#########################################################
"""
	_eval_T

Eval the Fourier basis up to degree n on an array of 
points in [0, 1].

"""
function _eval_T(n, x::Array{T}) where {T}
	k = length(x)
	M = zeros(Complex, k, n+1)
	M[:, 1] = ones(k)
	for i in 1:Int64(n/2)
        M[:, i+1] = 0.5*exp.(i*2im*pi*x)
        M[:, end-i+1] = 0.5*exp.(-i*2im*pi*x)
	end
	return M
end
function _eval_T_derivative(n, x::Array{T}) where {T}
	k = length(x)
	M = zeros(Complex, k, n+1)

	for i in 1:Int64(n/2)
        M[:, i+1] = 0.5*i*2im*pi*exp.(i*2im*pi*x)
        M[:,end-i+1] = -0.5*i*2im*pi*exp.(-i*2im*pi*x)
	end
	return M
end

evalFourier2Vector(coeff, B::Fourier2) = real.(_eval_T(length(coeff)-1, B.p)*coeff)
evalFourier2DerivativeVector(coeff,  B::Fourier2) =  real.(_eval_T_derivative(length(coeff)-1, B.p)*coeff)

evalFourier2(coeff, x) = real.(_eval_T(length(coeff)-1, [x])*coeff)[1]
evalFourier2Derivative(coeff, x) = real.(_eval_T_derivative(length(coeff)-1, [x])*coeff)[1]
"""
Maybe it is better to state that you want the mid of your interval in coding rather than here
"""
evalFourier2(coeff::Interval, x) = evalFourier(mid.(coeff), x)
evalFourier2Derivative(coeff::Interval, x) = evalFourier(mid.(coeff), x)

function evalFourierCentered(coeff, x::Interval)
	m = Interval(mid.(x))
	return evalFourier(coeff, m)+evalFourierDerivative(coeff, x)*(x-m)
end

###############################################################################
###############################################################################

function BasisDefinition.weak_projection_error(B::Fourier2)
	n = Float64(length(B), RoundUp)
	ν = B.k
	νf = Float64(B.k, RoundUp)
	den =  n⊗₋(νf ⊖₋ 2.0)⊗₋ π ⊗₋ reduce(⊗₋, [n-i for i in 2:ν-1])
	return (4.0 ⊗₊(n+1))⊘₊den
end
function BasisDefinition.aux_normalized_projection_error(B::Fourier2)
	n = Float64(length(B), RoundUp)
	ν = B.k
	νf = Float64(B.k, RoundUp)
	den =  π ⊗₋ νf ⊗₋ n ⊗₋ reduce(⊗₋, [n-i for i in 1:ν-1])
	return 2.0 ⊘₊den
end

"""
	BasisDefinition.strong_weak_bound(B::Chebyshev)

V.A. Markov estimate from GRADIMIR MILOVANOVIC EXTREMAL PROBLEMS AND 
INEQUALITIES OF MARKOV-BERNSTEIN TYPE FOR POLYNOMIALS
"""
# TODO: Check the indexes

function BasisDefinition.strong_weak_bound(B::Chebyshev)
	n = length(B)-1
	k = B.k-1
	# we want to estimate the norm of f^(k) by the C1 norm of f, so we use Markov estimate
	# for derivative k-1
	den = reduce(⊗₋,[Float64(2k-1-2*i, RoundDown) for i in 0:k-1])
	num = reduce(⊗₊,[Float64(n^2-i^2, RoundDown) for i in 0:k-1])
	return num ⊘₊ den ⊕₊ 1.0 # the 1.0 is to take into account the L1 norm of f
end
BasisDefinition.aux_weak_bound(B::Chebyshev) = 1.0

# Check this!!!
function BasisDefinition.weak_by_strong_and_aux_bound(B::Chebyshev)
	@error "TODO"
	ν = B.k	
	return ( ν, 1.)
end
BasisDefinition.bound_weak_norm_from_linalg_norm(B::Chebyshev) = @error "TODO"
BasisDefinition.bound_linalg_norm_L1_from_weak(B::Chebyshev) = @error "TODO"
BasisDefinition.bound_linalg_norm_L∞_from_weak(B::Chebyshev) = @error "TODO"
BasisDefinition.weak_norm(B::Chebyshev) = C1
BasisDefinition.aux_norm(B::Chebyshev) = L1
BasisDefinition.strong_norm(B::Chebyshev) = W{B.k, 1}

"""
	Base.getindex(B::Chebyshev, i::Int)

Make so that B[j] returns a HatFunctionOnTorus with the j-th basis element
"""
function Base.getindex(B::Fourier2, i::Int)
	n = length(B)
    v = zeros(n)
    v[i] = 1
	@boundscheck 1 <= i <= n || throw(BoundsError(B, i))
	return x-> exp(2*pi*im*(i-1)*x) #evalFourier2(v, x)
end

BasisDefinition.is_refinement(Bc::Fourier2, Bf::Fourier2) = length(Bc) < length(Bf)
BasisDefinition.integral_covector(B::Fourier2; T= Float64) = [Interval{T}(1); zeros(length(B)-1)]'
BasisDefinition.one_vector(B::Fourier2) = [1; zeros(length(B)-1)]

Base.length(S::AverageZero{T}) where T<:Chebyshev = length(S.basis)-1

function Base.iterate(S::AverageZero{T}, state = 1) where T<:Fourier2
	B = S.basis
	i = state
	if i == length(B)
		return nothing
	end
	v = zeros(length(B))
	v[i+1] = 1
	v[1] = -mid.(integral_covector(B)[i+1])
	return v, state+1
end

struct Fourier2Dual <: Dual
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    xp::Vector{Interval}
end

function Fourier2DualBranch(y, br::Branch, ylabel = 1:length(y), ϵ = 0.0)
	if br.increasing
		endpoint_X = br.X[2]
		der = Contractors.derivative(br.f)(endpoint_X)
		preim_der = preimages_and_derivatives(y, br, ylabel, ϵ)
		return [preim_der[1];],#endpoint_X is calculated in preimages
			[preim_der[2];], 
			[preim_der[3];]#same as previous comment
	else
		endpoint_X = br.X[2]
		der = Contractors.derivative(br.f)(endpoint_X)
		preim_der = preimages_and_derivatives(B.p, D, 1:length(B.p)-1, ϵ)
		return [preim_der[1]; endpoint_X], 
			[preim_der[2]; length(preim_with_der[2])+1], 
			[preim_der[3]; der]
	end
end

function Dual(B::Fourier2, D::PwMap, ϵ)
	@assert is_full_branch(D)
	results = collect(Fourier2DualBranch(B.p, b, 1:length(B.p), ϵ) for b in branches(D))
	x = vcat((result[1] for result in results)...)
    xlabel = vcat((result[2] for result in results)...)
	xp = vcat((result[3] for result in results)...)
    return x, xlabel, xp
end

Base.length(dual::Fourier2Dual) = length(dual.x)
Base.eltype(dual::Fourier2Dual) = Tuple{eltype(dual.xlabel), Tuple{eltype(dual.x), eltype(dual.x')}}
function Base.iterate(dual::Fourier2Dual, state=1)
    if state <= length(dual.x)
        return ((dual.xlabel[state], (dual.x[state], abs(dual.x'[state]))), state+1)
    else
        return nothing
    end
end


function Fourier2transform(w)
	#@info sum(w)
	n = length(w)
	z= fft(w)/n
	#t = real.(z[1:length(w)])
	return z#Interval.(t)
end

using ProgressMeter
function assemble_standard(B::Fourier2, D::Dynamic, ϵ=0.0; T = Float64)
	n = length(B.p)
	M = zeros(Complex{Interval{Float64}}, (n, n))
	x, labels, xp = Dual(B, D, ϵ)
	@showprogress for i in 1:n
		ϕ = B[i]
		w = zeros(Complex{Interval{Float64}}, n)
		for j in 1:length(x)
			w[labels[j]]+=ϕ(x[j])/abs(xp[j])
		end
		#@info w
		M[:, i] = Fourier2transform(mid.(w))
	end
	return M
end

function assemble(B::Fourier2, D::Dynamic, ϵ=0.0; T = Float64)
	n = length(B)
	n_less = Int64((n-1)/2)
	B_expanded = Fourier2(2*n,B.k)
	M_big = assemble_standard(B_expanded,D)
	M_upper = hcat(M_big[1:1+n_less,1:1+n_less],M_big[1:1+n_less,end-n_less+1:end])
	M_lower = hcat(M_big[end-n_less+1:end,1:1+n_less],M_big[end-n_less+1:end,end-n_less+1:end])
	return vcat(M_upper,M_lower)
end

using IntervalOptimisation

function infnormoffunction(B::Chebyshev, v)
	val = 0
	try 
		val = maximize(x->abs(evalChebyschevCentered(v, x)), Interval(0, 1))[1]
	catch
		print("Refining grid")
		f(x) = abs(evalChebyshevCentered(v, x))
		ran = range_estimate(f, Interval(0,1), 5)
		Bval = union(val, ran)
	end
	return val
end

function infnormofderivative(B::Chebyshev, v) 
	val = Interval(0)
	try 
		val = maximize(x->abs(evalChebyshevDerivative(v, x)), Interval(0, 1))[1]
	catch
		print("Refining grid")
		f(x) = abs(evalChebyshevDerivative(v, x))
		ran = range_estimate(f, Interval(0,1), 5)
		val = union(val, ran)
	end
	return val
end


BasisDefinition.is_integral_preserving(B::Chebyshev) = false
function BasisDefinition.opnormbound(B::Chebyshev, N::Type{C1}, v::Vector{S}) where {T, S} 
	return normbound(B, N, v)
end

function BasisDefinition.opnormbound(B::Chebyshev, N::Type{C1}, w::LinearAlgebra.Adjoint) where {T, S} 
	return normbound(B, N, w')
end

function BasisDefinition.opnormbound(B::Chebyshev, N::Type{C1}, A::Matrix{S}) where {T, S} 
	n, m = size(A)
	norm = 0.0
	for i in 1:m
		norm = max(norm, Float64(opnormbound(B, N, A[:,i]), RoundUp) ⊘₊ Float64(1+i^2, RoundDown))
		# since these are the images of the basis elements, we can use the fact that
		# |Tₖ(x)|= 1, Tₖ' = k*U_{k-1}, |U_{k-1}|=k
	end
	return norm ⊗₊ log(m+2)
end

BasisDefinition.normbound(B::Chebyshev{T}, N::Type{C1}, v) where {T} = Float64((infnormoffunction(B,v)+infnormofderivative(B,v)).hi, RoundUp)

mutable struct NormCacherC1 <: NormCacher{C1}
	B::Basis
    C::Float64
    function NormCacherC1(B, n)
        new(B, 0.0)
    end
end
NormCacher{C1}(B, n) = NormCacherC1(B, n)

function add_column!(Cacher::NormCacherC1, v::AbstractVector, ε::Float64)
    Cacher.C = max(Cacher.C, opnormbound(Cacher.B, C1, v) ⊕₊ ε)
end

"""
Return the norm of the matrix the NormCacher is working on.
"""
function get_norm(Cacher::NormCacherC1)
    n = length(Cacher.B)
	return Cacher.C ⊗₊ log(n+2)
end
