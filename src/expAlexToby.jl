using InvariantMeasures

T(x) = 3.3*x*(1-x)
ϕ(x::T) where {T} = abs(x)<1 ? exp(-1/(1-x^2)) : T(0)

using Plots 
plot(ϕ, -2, 2)

ρ(x) = ϕ(x)/0.443994
τ(x, ξ) = (1/ξ)*ρ(x/ξ)

#dx = InvariantMeasures.make_dx(1024, Float64).-0.5
#dy = mid.(τ.(dx, 0.02))

B = InvariantMeasures.Fourier1D(256, 1024)
Q = DiscretizedOperator(B, T)
P = Q.L

using IntervalArithmetic
Pmid = mid.(real(P))+im*mid.(imag(P))
using Pseudospectra, Plots
spectralportrait(Pmid)
savefig("specportraitnonoise.png")
eigen(Pmid)

function circleshape()
    θ = LinRange(0, 2*π, 500)
    return sin.(θ), cos.(θ)
end

using LinearAlgebra
val = eigen(Pmid).values
plt = scatter(val)
plt = plot!(circleshape(), fillalpha = 0.2)
savefig(plt, "eigenportraitnonoise.png")

using FFTW
trasf_002 = fft(dy)/1024
D = zeros(ComplexF64, 513)
InvariantMeasures.restrictfft!(D, trasf_002, 256)

using LinearAlgebra
NK = Diagonal(D)
P = NK*Q.L

Pmid = mid.(real(P))+im*mid.(imag(P))

using Pseudospectra
spectralportrait(Pmid)

dx = InvariantMeasures.make_dx(1024, Float64).-0.5
dy = mid.(τ.(dx, 0.01))

B = InvariantMeasures.Fourier1D(256, 1024)
Q = DiscretizedOperator(B, T)

using FFTW
trasf_002 = fft(dy)/1024
D = zeros(ComplexF64, 513)
InvariantMeasures.restrictfft!(D, trasf_002, 256)

using LinearAlgebra
NK = Diagonal(D)
P = NK*Q.L

Pmid = mid.(real(P))+im*mid.(imag(P))

using Pseudospectra
spectralportrait(Pmid)

N = 32768
for ξ in [0.001; 0.0001; 0.025; 0.015]#[0.001; 0.0025; 0.005; 0.0075; 0.01; 0.0125; 0.015; 0.0175; 0.002; 0.0025; 0.003]
    dx = InvariantMeasures.make_dx(N, Float64).-0.5
    dy = mid.(τ.(dx, ξ))
    trasf_002 = fft(dy)/N
    D = zeros(ComplexF64, 513)
    InvariantMeasures.restrictfft!(D, trasf_002, 256)
    NK = Diagonal(D)
    P = NK*Q.L
    Pmid = mid.(real(P))+im*mid.(imag(P))
    plt = spectralportrait(Pmid)
    savefig(plt, "specportrait$ξ.png")
end



for ξ in [0.001; 0.0025; 0.005; 0.0075; 0.01; 0.0125; 0.015; 0.0175; 0.002; 0.0025; 0.003]
    dx = InvariantMeasures.make_dx(N, Float64).-0.5
    dy = mid.(τ.(dx, ξ))
    trasf_002 = fft(dy)/N
    D = zeros(ComplexF64, 513)
    InvariantMeasures.restrictfft!(D, trasf_002, 256)
    NK = Diagonal(D)
    P = NK*Q.L
    Pmid = mid.(real(P))+im*mid.(imag(P))
    val = eigen(Pmid).values
    plt = scatter(val)
    plt = plot!(circleshape(), fillalpha = 0.2)
    savefig(plt, "eigenportrait$ξ.png")
end