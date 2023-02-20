

function IntervalArithmetic.midpoint_radius(v::Vector{Complex{Interval{T}}}) where {T}
    n = length(v)
    mid_vector = zeros(Complex{T}, n)
    rad_vector = zeros(Complex{T}, n)
    for i in 1:n
        real_m, real_r = midpoint_radius(real(v[i]))
        imag_m, imag_r = midpoint_radius(imag(v[i]))
        real_r = Float64(real_r, RoundUp)
        imag_r = Float64(imag_r, RoundUp)
        mid_vector[i] = real_m+im*imag_m
        rad_vector[i] = sqrt(square_round(real_r, RoundUp)⊕₊square_round(imag_r, RoundUp), RoundUp)
    end
    return mid_vector, rad_vector
end

using FFTW, FastTransforms, IntervalArithmetic
import AbstractFFTs

function interval_fft(P::AbstractFFTs.Plan{Complex{T}}, v::Vector{Complex{Interval{T}}}) where {T}
    n = Float64(length(v), RoundUp)
    u = Float64(eps(T), RoundUp)
    γ₄ = (4.0 ⊗₊ u)⊘₊(1.0 ⊖₋ 4.0 ⊗₋ u) 
    μ = u  
    η = μ ⊕₊ γ₄ ⊗₊ (sqrt_round(2.0, RoundUp) ⊕₊ μ) 
    t = ceil(log2(n))
    rel_err_fft = (t ⊗₊ η) ⊘₊(1.0 ⊖₋ η)
    
    norm_FFT_normalized_2 = 1.0 ⊘₊(sqrt(n, RoundUp))
    vector_mid, vector_radius = midpoint_radius(v)
    norm_obs = BasisDefinition.opnormbound(L2, vector_mid)
    norm_rad = BasisDefinition.opnormbound(L2, vector_radius)
    err_fft = norm_FFT_normalized_2⊗₊(rel_err_fft ⊗₊ norm_obs)⊕₊ norm_FFT_normalized_2⊗₊norm_rad
    mid_fft = P*vector_mid
    w = [Interval(real(z))+im*Interval(imag(z)) for z in mid_fft]
    return w.+(Interval(-err_fft, err_fft)+im*Interval(-err_fft, err_fft))
end

Base.:*(P::AbstractFFTs.Plan{Complex{T}}, v::Vector{Complex{Interval{T}}}) where {T} = interval_fft(P, v)


function interval_fft(v::Vector{Complex{Interval{T}}}) where {T}
    w = ones(T, length(v))
    P = plan_fft(w)
    return interval_fft(P, v)
end

function interval_fft(v::Vector{Interval{T}}) where {T}
    w = ones(T, length(v))
    v+=im*zeros(T, length(v))
    P = plan_fft(w)
    return interval_fft(P, v)
end