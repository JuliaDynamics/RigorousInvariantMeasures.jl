# Salvaged from the GitHub branch master, commit 41c661d, as the dfly_inf_der region of
# src/PwDynamicDefinition.jl. Not part of the package; see attic/README.md.
#
# main keeps the loop inline in src/PwDynamic.jl; the differences are the function-barrier
# split below, the two-pass tolerance (rough_tol over the thirteen candidate widths, then
# fine_tol on the best one), the error when no candidate gives A < 1, and the return value,
# which here is the pair (A, B) rather than a DFLYCoefficients. main also dispatches on
# VarInterval where this dispatches on TotalVariation.

"""
max_expansivity, int_distortion, bound_distortion = dfly_inf_der_auxiliary_quantities(br::PwDynamicDefinition.MonotonicBranch, radratio, tol=1e-3)

Return three auxiliary quantities used in dfly_inf_der:

* max_expansivity: maximum of the expansivity outside the critical intervals
* int_distortion: integral of the distortion over the critical intervals
* bound_distortion: lower bound to the distortion inside the critical intervals

The "critical intervals" for a branch defined on `(a,b) = br.X` are:
* [a, a+rad], if the branch has infinite derivative in a
* [b-rad, b], if the branch has infinite derivative in b

If the branch does not have infinite derivatives, int_distortion = 0 and bound_distortion = 0.

Three intervals are returned; it is up to the caller to take their .hi.
"""
function dfly_inf_der_auxiliary_quantities(br::PwDynamicDefinition.MonotonicBranch, i, tol=1e-3)
    left, right = has_infinite_derivative_at_endpoints(br) # this could be cached along calls with different values of `rad`, but it is cheap anyway
    int_distortion = Interval(0.0) # integral of the distortion over the critical intervals
    bound_distortion = Interval(0.0) # lower bound to the distortion inside the critical intervals
    rad = radius(hull(br.X[1], br.X[2]))
    f = inverse_derivative(br.f)
    g = distortion(br.f)
    left_endpoint = br.X[1]
    right_endpoint = br.X[2]
    if left
        left_endpoint += rad/2^i
    end
    if right
        right_endpoint-= rad/2^i
    end
    I = hull(left_endpoint, right_endpoint)
    max_result = maximise(x->abs(f(x)), I, tol=tol)
    @debug "max_result: $max_result"
    max_expansivity = max_result[1]
    
    if left
        # we work on the interval [br.X[1], left_endpoint] =: [a, b].
        # we assume here that f(a) = 0, g(a)=∞, and that g is monotonically decreasing on [a,b]
        l_left = abs(g(left_endpoint))
        # we use the fact that the primitive of the distortion g is f=1/T', and compute
        # int_a^b (distorsion) = f(b) - f(a) = f(b)
        int_distortion += abs(f(left_endpoint))
        bound_distortion = max(bound_distortion, abs(l_left))
    end
    if right
        # same reasoning as above but on the right endpoint
        l_right = abs(g(right_endpoint)) 
        int_distortion += abs(f(right_endpoint))
        bound_distortion = max(bound_distortion, abs(l_right))
    end
    return max_expansivity, int_distortion, bound_distortion
end

"""
    candidate_A, candidate_B = dfly_inf_der_single_estimate(D, i, width_term, tol=1e-3)

Computes dfly coefficients for a map with infinite derivatives, wtih a fixed width exponent i in the width of the problematic intervals.
"""
function dfly_inf_der_single_estimate(D, i; width_term=maximum( 2 / (br.X[2] - br.X[1]) for br in D.branches), tol=1e-3)
    results = [dfly_inf_der_auxiliary_quantities(br, i, tol) for br in D.branches]
    max_expansivity = maximum(r[1] for r in results)
    int_distortion = sum(r[2] for r in results)
    @debug "expansivities: $([r[1] for r in results]), int_distortion: $int_distortion"
    bound_distortion = maximum(r[3] for r in results)

    candidate_A = 2*max_expansivity.hi ⊕₊ (int_distortion/2).hi 
    candidate_B = bound_distortion.hi ⊕₊ width_term.hi
    return candidate_A, candidate_B
end
    
using ProgressMeter
"""
dfly inequality for maps with infinite derivatives. 
    
The strategy to compute it follows a variant of Lemma 9.1 in the GMNP paper: 
* we find a "problematic set" I by taking a small interval of size radius(branch domain)/2^3 around each endpoint with infinite derivative; 
* we find l such that T >= l for each point in I
* we compute the dfly coefficients as in the lemma.
* we repeat the computation replacing 2^3 with 2^4, 2^5, ... 2^15 and take the best estimate among these.
"""
function dfly_inf_der(::Type{TotalVariation}, ::Type{L1}, D::PwDynamicDefinition.PwMap; rough_tol=1e-4, fine_tol=1e-6)
    width_term = maximum( 2 / (br.X[2] - br.X[1]) for br in D.branches)
    est = +∞
    best_i = 0;

    @showprogress 1 "Checking candidate values for infinite-derivative DFLY..."  for i in 3:15
        candidate_A, candidate_B = dfly_inf_der_single_estimate(D, i; width_term, tol=rough_tol)
        @debug "i=$i, A=$candidate_A, B=$candidate_B"
        
        if candidate_A<1.0 && candidate_B ⊘₊ (1.0 ⊖₋ candidate_A) < est
            est = candidate_B ⊘₊ (1.0 ⊖₋ candidate_A)
            best_i = i
            @debug "improving on previous best"
        end
    end
    if best_i == 0
        error("Could not compute a dfly() inequality with A<1")
    end
    # recomputes the optimal value with a finer tolerance
    A, B = dfly_inf_der_single_estimate(D, best_i; width_term, tol=fine_tol)
    @info "Computed infinite-derivative dfly coefficients: A=$A, B=$B"
    return A, B
end
