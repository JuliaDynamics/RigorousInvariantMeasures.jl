module NoiseKernelDefinition
export DiscretizedNoiseKernel, UniformNoise


using FFTW

struct DiscretizedNoiseKernel{S<:AbstractVector, T<:AbstractVector}
	v::S
	Mfft::T
	rad
end

import IntervalArithmetic: Interval, mid, radius, @interval
import InvariantMeasures: opnormbound, Linf
DiscretizedNoiseKernel(v::Vector{Real}) = DiscretizedNoiseKernel(v, fft(v), 0)
DiscretizedNoiseKernel(v::Vector{Interval{T}}) where {T} = DiscretizedNoiseKernel(v, fft(mid.(v)), opnormbound(Linf, radius.(v)))
Mfft(Q::DiscretizedNoiseKernel) = Q.Mfft

#this is the Ulam discretization of the uniform noise of size ξ, on a partition of size k
function UniformNoise(ξ, k, boundarycondition = :periodic)
	n = Int64(floor(ξ*k))
	if boundarycondition == :periodic
		v = Interval.([ones(n); zeros(k-n)])*(1/(2*ξ))
		v+= Interval.([zeros(k-n); ones(n)])*(1/(2*ξ))
		v[n+1]+= @interval (ξ-n*1/k)*k/(2*ξ)
		v[k-n-1]+= @interval (ξ-n*1/k)*k/(2*ξ)
	end
	return DiscretizedNoiseKernel(v)
end


end