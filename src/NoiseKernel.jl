
#export DiscretizedNoiseKernel, UniformNoise

using FFTW, LinearAlgebra

abstract type NoiseKernel end
function Base.:*(M::NoiseKernel, v)
	@error "Not implemented"
end

opnormbound(N::NormKind, M::NoiseKernel) = @error "Not Implemented"
opradius(N::NormKind, M::NoiseKernel) = @error "Not Implemented"
nonzero_per_row(M::NoiseKernel) = @error "Not Implemented"

struct DiscretizedNoiseKernelFFT{S<:AbstractVector, T<:AbstractVector} <: NoiseKernel
	v::S
	Mfft::T
	rad
	P
end

import IntervalArithmetic: Interval, mid, radius, @interval
import InvariantMeasures: opnormbound, Linf
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
	v::S
	rad
	boundarycondition
end

function UniformNoiseUlam(ξ, k, boundarycondition = :periodic)
	n = Int64(2*ceil(ξ*k))
	v = zeros(Interval{Float64}, n)
	l = Int64(2*floor(ξ*k))
	a = 1/(2*Interval(ξ)*k)
	v[2:n-1] = a*ones(Interval{Float64}, l)
	v[1] = (1-sum(v))/2
	v[n] = v[1]
	return DiscretizedNoiseKernelUlam(mid.(v), opnormbound(L1, v)-1, boundarycondition) 
end

#TODO, but at the moment this is fine, it is a Markov operator
opnormbound(::Type{L1}, M::DiscretizedNoiseKernelUlam) = 1.0
opradius(::Type{L1}, M::DiscretizedNoiseKernelUlam) = M.rad
nonzero_per_row(M::DiscretizedNoiseKernelUlam) = length(M.v)


function Base.:*(M::DiscretizedNoiseKernelUlam, v)
	mult(M, v, Val(M.boundarycondition))
end

function mult(M::DiscretizedNoiseKernelUlam, v, ::Val{:periodic})
	k = length(M.v)
	n = length(v)
	w = zeros(n+2*k)
	w[k+1:k+n] = v
	w[1:k] = v[end-k+1:end]
	w[n+k+1:end] = v[1:k]
	z = zeros(n)
	for i in 1:n
		h = @view w[i:i+k-1]
		z[i] = sum( M.v.* h)
	end
	return z
end

function mult(M::DiscretizedNoiseKernelUlam, v, ::Val{:reflecting})
	k = length(M.v)
	n = length(v)
	w = zeros(n+2*k)
	w[k+1:k+n] = v
	w[1:k] = reverse(v[1:k])
	w[n+k+1:end] = reverse(v[end-k+1:end])
	z = zeros(n)
	for i in 1:n
		h = @view w[i:i+k-1]
		z[i] = sum( M.v.* h)
	end
	return z
end


