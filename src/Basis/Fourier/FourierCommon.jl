using .RigorousInvariantMeasures: Dual, interval_fft

export Fourier, evalFourier, FourierPoints, assemble_common, eval_on_dual


using IntervalArithmetic

abstract type Fourier <: Basis end


FourierPoints(n, T) = [Interval{T}(i) / (n) for i = 0:n-1]

Base.lastindex(B::Fourier) = length(B)



function evalFourier(coeff, x)
    @assert length(coeff) % 2 == 1

    K = length(coeff)
    N = (K - 1) ÷ 2
    pos_coeff = coeff[1:N+1]
    neg_coeff = [typeof(x)(0); reverse(coeff[N+2:K])]

    return evalpoly(exp(2 * pi * im * x), pos_coeff) +
           evalpoly(exp(-2 * pi * im * x), neg_coeff)
end

"""
Make so that B[j] returns the basis function of coordinate j
"""
function Base.getindex(B::Fourier, i::Int)

    K = length(B)
    @boundscheck 1 <= i <= K || throw(BoundsError(B, i))

    N = (K - 1) ÷ 2

    if i == 1
        return x -> typeof(x)(1.0)
    end

    if 1 < i <= N + 1
        return x -> exp(2 * pi * im * (i - 1) * x)
    end

    if N + 2 <= i <= K
        L = i + (-2 * N - 2)
        return x -> exp(2 * pi * im * L * x)
    end
end

is_refinement(Bc::Fourier, Bf::Fourier) = length(Bc) < length(Bf)
integral_covector(B::Fourier; T = Float64) = [Interval{T}(1); zeros(length(B) - 1)]'
one_vector(B::Fourier) = [1; zeros(length(B) - 1)]

Base.length(S::AverageZero{T}) where {T<:Fourier} = length(S.basis) - 1

function Base.iterate(S::AverageZero{T}, state = 1) where {T<:Fourier}
    B = S.basis
    i = state
    if i == length(B)
        return nothing
    end
    v = zeros(length(B))
    v[i+1] = 1
    return v, state + 1
end

abstract type FourierDual <: Dual end

function eval_on_dual(B::Fourier, computed_dual::FourierDual, ϕ) end

using ProgressMeter
function assemble_common(B::Fourier, D; ϵ = 0.0, max_iter = 100, T = Float64)
    n = length(B)

    @info n

    k = (n - 1) ÷ 2

    @info k

    M = zeros(Complex{Interval{Float64}}, (n, n))
    computed_dual = Dual(B, D; ϵ, max_iter)
    #@showprogress enabled=SHOW_PROGRESS_BARS  
    for i = 1:n
        ϕ = B[i]
        w = eval_on_dual(B, computed_dual, ϕ)
        #@info w

        FFTw = interval_fft(w)

        M[:, i] = [FFTw[1:k+1]; FFTw[end-k+1:end]]
    end
    return M
end
