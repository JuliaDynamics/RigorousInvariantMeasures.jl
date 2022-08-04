
#export DiscretizedNoiseKernel, UniformNoise

using FFTW, LinearAlgebra

abstract type NoiseKernel end
function Base.:*(M::NoiseKernel, v)
	@error "Not implemented"
end

BasisDefinition.opnormbound(N::NormKind, M::NoiseKernel) = @error "Not Implemented"
opradius(N::NormKind, M::NoiseKernel) = @error "Not Implemented"
nonzero_per_row(M::NoiseKernel) = @error "Not Implemented"

struct DiscretizedNoiseKernelFFT{S<:AbstractVector, T<:AbstractVector} <: NoiseKernel
	v::S
	Mfft::T
	rad
	P
end

import IntervalArithmetic: Interval, mid, radius, @interval
import RigorousInvariantMeasures: opnormbound, Linf
DiscretizedNoiseKernelFFT(v::Vector{Real}) = DiscretizedNoiseKernelFFT(v, fft(v), 0, plan_fft(mid.(v)))
DiscretizedNoiseKernelFFT(v::Vector{Interval{T}}) where {T} = DiscretizedNoiseKernelFFT(v, fft(mid.(v)), opnormbound(L2, radius.(v)), plan_fft(mid.(v)))
Mfft(Q::DiscretizedNoiseKernelFFT) = Q.Mfft

#this is the Ulam discretization of the uniform noise of size ξ, on a partition of size k
function UniformNoiseFFT(ξ, k, boundarycondition = :periodic)
	n = Int64(floor(ξ*k))
	if boundarycondition == :periodic
		v = Interval.([ones(n); zeros(k-n)])*(1/(2*ξ))
		v+= Interval.([zeros(k-n); ones(n)])*(1/(2*ξ))
		v[n+1]+= @interval (ξ-n*1/k)*k/(2*ξ)
		v[k-n-1]+= @interval (ξ-n*1/k)*k/(2*ξ)
	end
	return DiscretizedNoiseKernelFFT(v)
end

function Base.:*(M::DiscretizedNoiseKernelFFT, v)
	P = M.P
	w = P*v
	@info w
	w = M.Mfft.*w
	@info w
	return real.(P\w)/length(v)
end 

struct DiscretizedNoiseKernelUlam{S<:AbstractVector} <: NoiseKernel
	B::Ulam
	ξ
	v::S
	rad
	boundarycondition::Symbol
	w::Vector
	z::Vector
end

function UniformNoiseUlam(ξ, B::Ulam, boundarycondition = :periodic)
	k = length(B)
	n = 2*Int64(ceil(ξ*k))
	v = zeros(Interval{Float64}, n)
	a = 1/(2*Interval(ξ))
	v[2:n-1] = a*ones(Interval{Float64}, n-2)
	v[1] = (k-sum(v))/2
	v[n] = v[1]
	if boundarycondition == :periodic
		return DiscretizedNoiseKernelUlam(B,
										  Interval(ξ),
										  mid.(v), 
										  opnormbound(L1, v)-k, 
										  boundarycondition, 
										  zeros(k+n), 
										  zeros(k)) 
	elseif boundarycondition == :reflecting
		return DiscretizedNoiseKernelUlam(B,
										  Interval(ξ),
										  mid.(v), 
										  opnormbound(L1, v)-k, 
										  boundarycondition, 
										  zeros(k+n+2), 
										  zeros(k))
	end
end

#TODO, but at the moment this is fine, it is a Markov operator
BasisDefinition.opnormbound(B::Ulam, ::Type{L1}, M::DiscretizedNoiseKernelUlam) = 1.0
opradius(::Type{L1}, M::DiscretizedNoiseKernelUlam) = M.rad
nonzero_per_row(M::DiscretizedNoiseKernelUlam) = length(M.v)
dfly(::Type{TotalVariation}, ::Type{L1}, N::DiscretizedNoiseKernelUlam) = (0.0, (1/(2*N.ξ)).hi)


function Base.:*(M::DiscretizedNoiseKernelUlam, v)
	mult(M, v, Val(M.boundarycondition))
end

function mult(M::DiscretizedNoiseKernelUlam, v, ::Val{:periodic})
	n = length(M.v)
	k = length(v)
	l =n÷2

	M.w .= 0
	
	M.w[l+1:l+k] = v
	M.w[1:l] = v[end-l+1:end]
	M.w[end-l+1:end] = v[1:l]
		
	for i in 1:k
		h = @view M.w[i:i+n-1]
		v[i] = sum( M.v.* h)/k
	end

	return v
end

function mult(M::DiscretizedNoiseKernelUlam, v::Vector{Interval{T}}, ::Val{:periodic}) where {T}
	n = length(M.v)
	k = length(v)
	l =n÷2

	nrmv = opnormbound(M.B, L1, v)
	midv = mid.(v)
	radv = radius.(v)
	nrmrad = opnormbound(M.B, L1, radv)

	M.w .= 0
	
	M.w[l+1:l+k] = midv
	M.w[1:l] = midv[end-l+1:end]
	M.w[end-l+1:end] = midv[1:l]
		
	for i in 1:k
		h = @view M.w[i:i+n-1]
		v[i] = sum( M.v.* h)/k
	end
	δₖ = opradius(L1, M)
	γₖ = gamma(T, nonzero_per_row(M)) 
    nrm_MK = opnormbound(M.B, L1, M)
    normMK = nrm_MK ⊕₊ δₖ	
	
	ϵ = (γₖ ⊗₊ normMK) ⊗₊ nrmv ⊕₊ normMK ⊗₊ nrmrad 		

	return v + fill(Interval(-ϵ ,ϵ), length(v)) 
end



function mult(M::DiscretizedNoiseKernelUlam, v, ::Val{:reflecting})
	n = length(M.v)
	k = length(v)
	l = n÷2+1

	@info v
	M.w .= 0
	M.w[l+1:l+k] = v
	v .= 0
	
	for i in 1:l
		h = @view M.w[i:i+n-1]
		@info h, i
		val = sum( M.v.* h)/k 
		v[l+1-i] += val 
	end
	
	for i in 1:k
		h = @view M.w[i+l:i+l+n-1]
		@info h, i+l
		val = sum( M.v.* h)/k 
		v[i]+=val 
	end
	@info "stop"
	for i in 1:l
		h = @view M.w[k+l+i:k+l+i+n-1]
		@info h
		val = sum( M.v.* h)/k 
		v[end-i+1] +=val 
	end


	return v
end


