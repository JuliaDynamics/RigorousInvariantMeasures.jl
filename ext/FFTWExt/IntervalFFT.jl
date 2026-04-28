# Per-element midpoint/radius split for `Vector{Complex{Interval}}`. The
# radius is widened up so the disk it describes contains the rectangular
# Complex{Interval} enclosure.
function _midradius_complex_interval(v::Vector{Complex{Interval{T}}}) where {T}
    n = length(v)
    mid_vec = zeros(Complex{T}, n)
    rad_vec = zeros(T, n)
    for i = 1:n
        re_m, re_r = mid(real(v[i])), radius(real(v[i]))
        im_m, im_r = mid(imag(v[i])), radius(imag(v[i]))
        re_r = Float64(re_r, RoundUp)
        im_r = Float64(im_r, RoundUp)
        mid_vec[i] = re_m + im * im_m
        rad_vec[i] =
            sqrt_round(square_round(re_r, RoundUp) ⊕₊ square_round(im_r, RoundUp), RoundUp)
    end
    return mid_vec, rad_vec
end

@doc raw"""
    interval_fft(v::Vector{Complex{Interval{Float64}}})
    interval_fft(v::Vector{Interval{Float64}})

Rigorously enclosing FFT of `v`, normalized by `length(v)`. Higham 1996
a-priori bound: a relative-roundoff term scaled by ``\log_2(N)/\sqrt{N}``
times the L² norm of input midpoints, plus a ``1/\sqrt{N}`` term picking
up the L² norm of the input radii.
"""
function interval_fft(v::Vector{Complex{Interval{Float64}}})
    n = Float64(length(v), RoundUp)
    u = Float64(eps(Float64), RoundUp)
    γ₄ = (4.0 ⊗₊ u) ⊘₊ (1.0 ⊖₋ 4.0 ⊗₋ u)
    η = u ⊕₊ γ₄ ⊗₊ (sqrt_round(2.0, RoundUp) ⊕₊ u)
    rel_err_fft = (ceil(log2(n)) ⊗₊ η) ⊘₊ (1.0 ⊖₋ η)

    inv_sqrt_n = 1.0 ⊘₊ sqrt_round(n, RoundUp)
    vec_mid, vec_rad = _midradius_complex_interval(v)
    norm_obs = opnormbound(L2, vec_mid)
    norm_rad = opnormbound(L2, vec_rad)
    err_fft =
        inv_sqrt_n ⊗₊ (rel_err_fft ⊗₊ norm_obs) ⊕₊ inv_sqrt_n ⊗₊ norm_rad

    mid_fft = FFTW.fft(vec_mid) ./ length(v)
    err_box = interval(-err_fft, err_fft) + im * interval(-err_fft, err_fft)
    return [interval(real(z)) + im * interval(imag(z)) + err_box for z in mid_fft]
end

interval_fft(v::Vector{Interval{Float64}}) = interval_fft(v .+ 0im)
