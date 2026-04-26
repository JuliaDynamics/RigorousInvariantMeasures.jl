
"""
(A, B) = dfly(strongnorm, auxnorm, dynamic)

Constants (A, B) such that ||Lf||_s ≦ A||f||_s + B||f||_aux
"""
dfly(::Type{<:NormKind}, ::Type{<:NormKind}, ::Dynamic) = @error "Not implemented"
# Instance dispatch fallback: forward to type dispatch for simple norms
dfly(n1::NormKind, n2::Type{<:NormKind}, D::Dynamic) = dfly(typeof(n1), n2, D)

# I don't think this is used in production anymore
# function dfly(N1::Type{TotalVariation}, N2::Type{L1}, D::Dynamic)
#     dist = max_distortion(D)
#     lam = max_inverse_derivative(D)

#     if !(abs(lam) < 1) # these are intervals, so this is *not* equal to abs(lam) >= 1.
#         @error "The function is not expanding"
#     end

#     if is_full_branch(D)
#         return sup(lam), sup(dist)
#     else
#         if !(abs(lam) < 0.5)
#             @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
#         end
#         endpts = endpoints(D)
#         min_width = minimum([endpts[i+1] - endpts[i] for i = 1:length(endpts)-1])
#         return sup(lam), sup(dist) ⊕₊ sup(2 / min_width)
#     end
# end

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
        return sup(lam), sup(dist)
    else
        if !(abs(2 * lam) < 1)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        return sup(2 * lam), sup(dist + disc)
    end
end

# W{1,1} → TotalVariation bridge: W^{1,1} seminorm = total variation
# Accept both type and instance dispatch
dfly(::Type{W{1,1}}, N2::Type{L1}, D::PwMap) = dfly(TotalVariation, L1, D)
dfly(::W{1,1}, N2::Type{L1}, D::PwMap) = dfly(TotalVariation, L1, D)

# Higher-order W{k,1} (k≥2) is provided by the SymbolicsExt extension.
# Without Symbolics loaded, the generic dfly fallback will @error "Not implemented".

# Aη DFLY: analytic strip norm
# ||Lf||_{Aη'} ≤ A · ||f||_{Aη} + B · ||f||_{L¹}
# For expanding full-branch maps: A = max|1/T'|, B = distortion bound
# This is a conservative estimate; tighter bounds require complex interval evaluation
function dfly(norm::Aη, ::Type{L1}, D::PwMap)
    if has_infinite_derivative_at_endpoints(D)
        # Fall back to TotalVariation DFLY for infinite derivative maps
        return dfly(TotalVariation, L1, D)
    end
    dist = max_distortion(D)
    lam = max_inverse_derivative(D)

    if is_full_branch(D)
        if !(abs(lam) < 1)
            @error "The function is not expanding"
        end
        return sup(lam), sup(dist)
    else
        if !(abs(2 * lam) < 1)
            @error "Expansivity is insufficient to prove a DFLY. Try with an iterate."
        end
        vec = endpoints(D)
        disc = maximum(2 / abs(vec[i] - vec[i+1]) for i = 1:nbranches(D))
        return sup(2 * lam), sup(dist + disc)
    end
end

function dfly(::Type{Lipschitz}, ::Type{L1}, D::Dynamic)
    # TODO: should assert that D is globally C2 instead, but we don't have that kind of infrastructure yet.
    @assert is_full_branch(D)

    dist = max_distortion(D)
    #@info dist
    lam = max_inverse_derivative(D)
    #@info lam

    return (sup(lam * (2 * dist + 1)), sup(dist * (dist + 1)))
end
