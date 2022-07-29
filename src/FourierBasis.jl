using ..BasisDefinition, ..DynamicDefinition, ..Contractors, ..PwDynamicDefinition
IntervalArithmetic, LinearAlgebra
import Base: iterate
import ..BasisDefinition: one_vector, integral_covector, is_integral_preserving, strong_norm, weak_norm, aux_norm

struct Fourier1D <: Basis
    N::Int64
    FFTNx::Int64
end

make_dx(N, T) = [Interval{T}(i)/N for i in 0:N-1]

# converts [0,1,..., Nx,-Nx,...,-1] to [1,...,2*Nx+1]
unidimensional_index(i, Nx) = i>=0 ? i+1 : 2*Nx+2+i 
inverse_unidimensional_index(i, Nx) = 1<=i<=Nx+1 ? i-1 : i-2*Nx-2

@inline ϕ(k, x; L=1)  = exp(im*2*pi*x*k/L) # L is the length of the interval
@inline ψ(k, h, x, y; Lx=1, Ly=1) = (-1)^(k+h)*ϕ(k, x; L = Lx)*ϕ(h, y; L= Ly)

"""
Fourier(N) builds a Fourier basis truncating at frequence N
by default it uses a FFT of size 2^(log2(ceil(N))+2)

Fourier(N, FFTNx) builds a Fourier basis with FFT size FFTNx

The base size is 2*N+1 (0, ..., N, -N, ..., -1)
"""
Fourier1D(N::Integer) = Fourier1D(N, 2^Int64(ceil(log2(N))+2))
Fourier1D(N::Integer, FFTNx::Integer) = Fourier1D(N, FFTNx)
Base.length(B::Fourier1D) = 2*B.N+1

evaluate(B::Fourier1D, i, x) = ϕ(i, x)
strong_norm(B::Fourier1D) = L2
weak_norm(B::Fourier1D) = L1
aux_norm(B::Fourier1D) = L1
BasisDefinition.is_refinement(Bfine::Fourier1D, Bcoarse::Fourier1D) = Bfine.N > Bcoarse.N ? true : false

is_integral_preserving(B::Fourier1D) = false
function integral_covector(B::Fourier1D; T= Float64) 
    v = zeros(length(B))
    v[1] = T(1)
    return v'
end

one_vector(B::Fourier1D) = [1.0; zeros(2*B.N)]

#weak_projection_error(B::Fourier1D) = 

function restrictfft!(new::Vector, orig::Vector, Nx)
    FFTNx = length(orig)
    new[1:Nx+1] = @view orig[1:Nx+1]
    new[Nx+1:2*Nx+1] = @view orig[FFTNx-Nx:FFTNx]
    return new
end 

function extendfft(orig::Vector, FFTNx) 
    N = length(orig)
    Nx = (N-1) ÷ 2 
    return [orig[1:Nx]; zeros(Complex{Float64}, FFTNx-N); orig[Nx+1:end]]
end

function  norm_2_upper_bound(v::Vector{Float64})
    sum = 0.0
    for x in v
        sum = sum ⊕₊ square_round(x, RoundUp)
    end
    return sqrt_round(sum, RoundUp)
end  

function abs_upper_square(x::Complex{Float64}) 
    absval = square_round(real(x), RoundUp)
    absval = absval ⊕₊ square_round(imag(x), RoundUp)
    return absval
end

function  norm_2_upper_bound(v::Vector{Complex{Float64}}) 
    sum = 0.0
    for x in v
        sum = sum ⊕₊ abs_upper_square(x)
    end
    return sqrt_round(sum, RoundUp)
end  





using FFTW, FastTransforms

"""
Assembles the truncated Fourier matrix for a dynamic D, 
with frequencies -Nx≤ i≤ Nx

Input:
    B:: Fourier1D basis
    D: Dynamic
    T: Floating point type, supports Float64, BigFloat

    Output:
    M: matrix of intervals that rigorously contain M 
"""
function assemble(B::Fourier1D, D, ϵ = 0.0; T = Float64) 
    Nx = B.N
    FFTNx = B.FFTNx

    I = Interval{T}
    dx = make_dx(B.FFTNx, T)
    
    Lx = T(1)
    #@info dx
    one = [T(1) for x in dx]
    P = plan_fft(one)
    # we take the adjoint since we are computing the adjoint operator
  
    Dx = [D(x) for x in dx]

    N = (2*Nx+1) # we are taking Nx positive and negative frequencies and the 0 frequency
    
    M = zeros(Complex{T},(2*Nx+1, 2*Nx+1))
    
    observablevalue = zeros(Complex{T}, FFTNx) 
    observablerad = zeros(Float64, FFTNx) 
    
    onedtransform = zeros(Complex{T}, FFTNx)
    new = zeros(Complex{T}, 2*Nx+1)
    l2error = 0.0
    # a priori estimate of the FFT error
    # we use the fact that 
    # ||FFT(v)/FFTNx||₂ ≤ 1/√FFTNX
    # and
    # ||Fl(FFT(v))- FFT(v)||₂ ≤ tη/(1-η)||v||₂ 
    # where η = μ + γ₄(√2+μ), t > log₂(FFTNx)
    # and μ is the absolute error in the computation of the wiggle factors
    # (we assume them to be precise to machine precision)
    # and γ₄ = 4u/(1-4u) and u is the unit roundoff
    # this bound is from Higham N. J. - Accuracy and Stability of Numerical Algorithms
    # Second Edition - SIAM

    u = Float64(eps(T), RoundUp)
    γ₄ = (4.0 ⊗₊ u)⊘₊(1.0 ⊖₋ 4.0 ⊗₋ u) 
    #@info γ₄
    μ = u  
    η = μ ⊕₊ γ₄ ⊗₊ (sqrt_round(2.0, RoundUp) ⊕₊ μ) 
    #@info μ
    t = ceil(log2(FFTNx))
    rel_err_fft = (t ⊗₊ η) ⊘₊(1.0 ⊖₋ η)
    f_FFTNx = Float64(FFTNx)
    norm_FFT_normalized_2 = 1.0 ⊘₊(sqrt(f_FFTNx, RoundUp))
    
    for i in 1:2*Nx+1    
        l = inverse_unidimensional_index(i, Nx) # the index in the form [0, ..., Nx, -Nx, ..., -1]
        
        for (ind, val) in pairs(Dx)
            obsval = ϕ(l, val; L=Lx)
            real_m, real_r = midpoint_radius(real(obsval))
            imag_m, imag_r = midpoint_radius(imag(obsval))
            observablevalue[ind] =  real_m+im*imag_m
            observablerad[ind] = sqrt(real_r^2+imag_r^2)
        end
        norm_obs = Float64(norm_2_upper_bound(observablevalue), RoundUp)
        norm_rad = Float64(norm_2_upper_bound(observablerad), RoundUp)
        err_fft = norm_FFT_normalized_2⊗₊(rel_err_fft ⊗₊ norm_obs)⊕₊ norm_FFT_normalized_2⊗₊norm_rad
        l2error= max(l2error, err_fft)        
        
        mul!(onedtransform, P, observablevalue)
        restrictfft!(new, onedtransform, Nx)
            
        for (ind, val) in pairs(new)
                if abs(val)!= 0
                    # this is the adjoint matrix of the Koopman operator
                    M[i, ind] = conj(val)/FFTNx
                end
        end
    end

    entrywise_error = l2error ⊗₊ sqrt_round(Float64(N), RoundUp)

    return Interval.(M) .+ Interval(-entrywise_error, entrywise_error)
end

"""
noise_matrix(B, σ)

returns the matrix of the (periodic) convolution operator with a Gaussian of average 0
and variance σ in the Fourier basis, truncated at frequence B.N 
"""

struct GaussianNoise <: NoiseKernel
    B::Basis
    σ
    NK::Matrix
end

GaussianNoise(B::Fourier1D, σ) = GaussianNoise(B, σ, Diagonal([[exp(-(k*σ)^2/2) for k in 0:B.N]; [exp(-(k*σ)^2/2) for k in -B.N:-1]]))


"""
This function computes a candidate for the top singular value
of M|_{V_0} square, i.e., it returns λ^2 and the associated singular vector 
"""
function norm_2_estimator_square_nonrig(M; n_it= 10)
    n, m = size(M)
    @assert n==m
    v = ones(m)
 
    for i in 1:n_it
        v = M'*(M*v)
        v/=norm(v,2)
    end
    v = M'*M*v

    return sqrt(norm(v, 2)), v/norm(v, 2)
end    


"""
This function uses the classical bound, cited in

RESIDUAL BOUNDS ON APPROXIMATE EIGENSYSTEMS OF NONNORMAL MATRICES*
W. KAHAN, B. N. PARLETT AND E. JIANGt

Theorem 1

"""
function faster_certify(A, λ, v)
    w = Interval.(v)
    ν = Interval(λ)^2
    ρ = (norm(A'*A*w-ν*w, 2)/norm(w, 2)).hi
    ν = Interval(-ρ, ρ) + λ
    return ν
end

using IntervalArithmetic

upper(x::Interval) = x.hi

"""
rigorous_norm(M; k = 10)

Computes the L² and L∞ norm of the powers of M|_V₀ up to the power k.
The algorithm first computes a numeric estimate of the top Singular Value
and then certifies it by an Interval Newton step

Outputs: v2
- v₂ : it is the vector with the rigorous bounds on the L² norm
       these are obtained by observing that
       ||f||₂ = ||w||₂, 
       where w is the vector of the Fourier coefficients of the trigonometric polynomial f
"""
function rigorous_norm(M::Matrix{Interval}; k = 20)
   A = M
   n, m = size(M)
   norms = zeros(Interval{Float64}, k)
   j = 0
   for i in 1:k
       λ, v = norm_2_estimator_square_nonrig(A)
       norms[i] = sqrt(faster_certify(mid.(A), λ, v))
       A = M*A
       j = j+1
       @info(j)
   end     
   rescale = sqrt(Interval(n-1))
   @info "rescale", rescale
   return norms, rescale*norms #check to make it rigorous
end

"""
This function bounds the L2 operator norm of M (complex floating point matrix)
"""
function opnormbound(B::Fourier1D, N::Type{L2}, M; n_it = 10)
    λ, v = norm_2_estimator_square_nonrig(M, n_it = n_it)
    return abs(faster_certify(M, λ, v)).hi
end

function opnormbound(B::Fourier1D, N::Type{L2}, M::Array{Complex{Interval{T}}, 2}; n_it = 10) where {T}
    Mid = mid.(real(M))+im*mid.(imag(M))
    Rad = radius.(real(M))+im*radius.(imag(M))
    mid_bound = opnormbound(B, N, Mid; n_it = n_it)
    rad_bound = opnormbound(B, N, Rad; n_it = n_it) 
    return mid_bound ⊕₊ rad_bound
end

function norms_of_powers_noise_interval( B::Fourier1D,
    N::Type{L2}, 
    m::Integer, 
    Q::DiscretizedOperator,
    NK::GaussianNoise, 
    f::AbstractArray;
    normv0::Real=-1., #used as "missing" value
    normQ::Real=-1.,
    normE::Real=-1.,
    normEF::Real=-1.,
    normIEF::Real=-1.,
    normN::Real=-1.,
    normρ::Real=-1.)


    P = NK.NK*Q.L

    M = P[2:end, 2:end]
    norms = zeros(m)

    A = copy(M)
    @showprogress for i in 1:m 
        nrmA = opnormbound(B, N, A)
        A = M*A 
        norms[i] = nrmA
    end
    return norms
end


function norms_of_powers_noise( B::Fourier1D,
                                N::Type{L2}, 
                                m::Integer, 
                                Q::DiscretizedOperator,
                                NK::GaussianNoise, 
                                f::AbstractArray;
                                normv0::Real=-1., #used as "missing" value
                                normQ::Real=-1.,
                                normE::Real=-1.,
                                normEF::Real=-1.,
                                normIEF::Real=-1.,
                                normN::Real=-1.,
                                normρ::Real=-1.)



    
    P = NK.NK*Q.L
    
    Prestricted = P[2:end, 2:end]

    M = mid.(real(Prestricted))+im*mid.(imag(Prestricted))
    T = eltype(real.(M))
    n = size(M, 1)

    R = radius.(real(Prestricted))+im*radius.(imag(Prestricted))
    δ = opnormbound(B, N, R)
    
    γz = gamma(T, n)
    
    ϵ = zero(T)

    nrmM = opnormbound(B, N, M)
    
    norms = zeros(Float64, m)

    A = copy(M)
    nrmA = nrmM
    norms[1] = nrmM
    @showprogress for i in 1:m 
        nrmA = opnormbound(B, N, A) 
        norms[i] = nrmA+ϵ
        A = M*A
        ϵ = (γz ⊗₊ nrmM ⊕₊ δ) ⊗₊ nrmA ⊕₊ (nrmM ⊕₊ δ) ⊗₊ ϵ       
    end
    return norms
end

"""
Array of "trivial" bounds for the powers of a DiscretizedOperator (on the whole space)
coming from from ||Q^k|| ≤ ||Q||^k
"""
function norms_of_powers_trivial_noise(B::Fourier1D,
                                       ::Type{L2}, 
                                       Q::DiscretizedOperator, 
                                       MK::NoiseKernel,
                                       m::Integer)

    norms = fill(NaN, m)
    P = MK.MK*Q.L
                                   
    norms[1] = opnormbound(B, N, P)
    for i = 2:m
        norms[i] = norms[i-1] ⊗₊ norms[1]
    end
    return norms
end

function norms_of_powers_trivial_noise(B::Fourier1D,
                                       ::Type{L1}, 
                                       Q::DiscretizedOperator, 
                                       MK::NoiseKernel,
                                       m::Integer)
    n = 2*B.N+1

    norms = fill(NaN, m)
    norms[1] = 1.0 ⊕₊ noise_projection_error(B, MK)
    for i = 2:m
        norms[i] = norms[i-1] ⊗₊ norms[1]
    end
    return norms
end


dfly(::Type{L2}, ::Type{L1}, N::GaussianNoise) = (0, (Interval(1)/(Interval(2)*sqrt(Interval(2)))).hi)

function noise_projection_error(B::Fourier1D, N::GaussianNoise)
    n = B.N
    σ = N.σ
    val =  @interval -2*π*σ^2
    return (exp(n^2*val)/(1-exp(val))).hi
end 

strong_weak_bound(B::Fourier1D) = (2*B.n+1)
aux_weak_bound(B::Fourier1D) = 1
weak_by_strong_and_aux_bound(B::Fourier1D) = (0., 1.)

function norms_of_powers_abstract_noise(Bas::Fourier1D, N::NoiseKernel, m)
    A, B = dfly(strong_norm(Bas), aux_norm(Bas), N)
    Eh = noise_projection_error(Bas, N)
    M₁n = BasisDefinition.strong_weak_bound(Bas)
    M₂ = BasisDefinition.aux_weak_bound(Bas)
    S₁, S₂ = BasisDefinition.weak_by_strong_and_aux_bound(Bas)

    norms = fill(NaN, m)
    strongs = fill(NaN, m)

    v = Array{Float64}([M₁n; M₂])
    # We evaluate [S₁ S₂] * ([1 0; Eh 1]*[A B; 0 1])^k * [M₁n; M₂] (with correct rounding)
    for k = 1:m
        # invariant: v[1] bounds ||Q^kf||_s for ||f||_w=1
        # v[2] bounds |||Q^kf||| for ||f||_w=1
        v[1] = A ⊗₊ v[1] ⊕₊ B ⊗₊ v[2]
        v[2] = Eh ⊗₊ v[1] ⊕₊ v[2]
        strongs[k] = v[1]
        norms[k] = S₁ ⊗₊ v[1] ⊕₊ S₂ ⊗₊ v[2]
    end
    return strongs, norms
end

"""
Estimate norms of powers from those on a coarser grid (see paper for details)
TODO: Check the final formula!!!
"""
function norms_of_powers_from_coarser_grid_noise(fine_basis::Fourier1D, 
                                            coarse_basis::Fourier1D, 
                                            Q::DiscretizedOperator,
                                            NK::NoiseKernel, 
                                            coarse_norms::Vector)
    
    ####Check this!!!!
    #if !(BasisDefinition.is_refinement(fine_basis, coarse_basis))
    #    @error "The fine basis is not a refinement of the coarse basis"
    #end
    m = length(coarse_norms)
    fine_norms = fill(NaN, m+1)
    trivial_norms = norms_of_powers_trivial_noise(fine_basis, weak_norm(fine_basis), Q, NK, m+1)

    A, B = dfly(strong_norm(fine_basis), aux_norm(fine_basis), NK)

    # adds a 0th element to strongs
    trivial_norms0(k::Integer) = k==0 ? 1. : trivial_norms[k]
    coarse_norms0(k::Integer) = k==0 ? 1. : coarse_norms[k]

    Kh =  noise_projection_error(coarse_basis, NK)
    
    fine_norms[1] = trivial_norms0(1)

    for k in 1:m
		temp = 0.
		for i in 0:k-1
			temp = temp ⊕₊ coarse_norms0(i) ⊗₊ (trivial_norms0(k-i) ⊕₊ trivial_norms0(k-i-1))
		end
        fine_norms[k+1] = coarse_norms0(k) ⊕₊ B ⊗₊ Kh ⊗₊ (temp ⊕₊ 1.0)  
	end
    return fine_norms
end

function powernormboundsnoise(B::Fourier1D; Q=DiscretizedOperator(B, D), NK = NK::NoiseKernel)
	m = 8
	computed_norms = []
	while true
		computed_norms = norms_of_powers_noise(B, strong_norm(B), m, Q, NK, integral_covector(B))
		if any(computed_norms .< 0.1)
			break
		end
		m = 2*m
	end
	trivial_norms = norms_of_powers_trivial_noise(B, strong_norm(B), Q, NK, m)
	norms = min.(trivial_norms, computed_norms)

	m_extend = 2*m
	better_norms = []
	while true
		better_norms = refine_norms_of_powers(norms, m_extend)
		if better_norms[end] < 1e-8
			break
		end
		m_extend = 2*m_extend
	end

	return better_norms

end