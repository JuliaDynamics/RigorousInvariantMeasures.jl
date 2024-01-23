
"""
(A, B) = dfly(strongnorm, auxnorm, dynamic)

Constants (A, B) such that ||Lf||_s ≦ A||f||_s + B||f||_aux
"""
dfly(::Type{<:NormKind}, ::Type{<:NormKind}, ::Dynamic) = @error "Not implemented"

# I don't think this is used in production anymore
function dfly(N1::Type{TotalVariation}, N2::Type{L1}, D::Dynamic)
    dist = max_distortion(D)
    lam = max_inverse_derivative(D)

    if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
        @error "The function is not expanding"
    end

    if is_full_branch(D)
        return lam.hi, dist.hi
    else
        if !(abs(lam) < 0.5)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        endpts = endpoints(D)
        min_width = minimum([endpts[i+1] - endpts[i] for i = 1:length(endpts)-1])
        return lam.hi, dist.hi ⊕₊ (2 / min_width).hi
    end
end

function dfly(N1::Type{TotalVariation}, N2::Type{L1}, D::PwMap)
    if has_infinite_derivative_at_endpoints(D)
        return dfly_inf_der(N1, N2, D, 10^-3)
    end

    dist = max_distortion(D)
    lam = max_inverse_derivative(D)
    vec = endpoints(D)
    disc = maximum(2 / abs(vec[i] - vec[i+1]) for i = 1:nbranches(D))

    if is_full_branch(D)
        if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
            @error "The function is not expanding"
        end
        return lam.hi, dist.hi
    else
        if !(abs(2 * lam) < 1)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        return (2 * lam).hi, (dist + disc).hi
    end
end

function dfly(::Type{Lipschitz}, ::Type{L1}, D::Dynamic)
    # TODO: should assert that D is globally C2 instead, but we don't have that kind of infrastructure yet.
    @assert is_full_branch(D)

    dist = max_distortion(D)
    #@info dist
    lam = max_inverse_derivative(D)
    #@info lam

    return ((lam * (2 * dist + 1)).hi, (dist * (dist + 1)).hi)
end
