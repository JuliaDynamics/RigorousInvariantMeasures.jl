import Pkg
Pkg.activate("../")

using RigorousInvariantMeasures, IntervalArithmetic

θ = 109 / 64
α = 51 / 64

D0 = Lorenz(θ, α)
D = D0 ∘ D0 ∘ D0

G(x) = (x - 0.5) * log(x - 0.5) - (x - 0.5)
function discretizationlogderLorenz(B::Ulam, θ, α)
    N = length(B)
    v = [G(Interval(B.p[i+1])) - G(Interval(B.p[i])) for i = (N÷2+2):N]
    w = [G(Interval(B.p[(N÷2+2)])); v]
    z = [reverse(w); w] * (α - 1)
    return z * N .+ (log(α) + log(θ))
end

size_coarse = 2^15
size_fine = 2^20

@info "Size coarse: $(size_coarse) Size fine: $(size_fine)"

using JLD

file = jldopen("output_Lorenz_$(size_coarse)_$(size_fine).jld", "w")


@info "Computing Lorenz"

B = Ulam(size_coarse)
Q = DiscretizedOperator(B, D; ϵ = 10^(-14), max_iter = 100)

dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D)
@info dfly_coefficients

norms = powernormbounds(B, D, Q = Q, m = 40)
@info norms

B_fine = Ulam(size_fine)
Q_fine = DiscretizedOperator(B_fine, D; ϵ = 10^(-14), max_iter = 100)

normQ_fine = opnormbound(B_fine, weak_norm(B_fine), Q_fine)
norms2 = refine_norms_of_powers(norms, 400)
norms_fine = finepowernormbounds(
    B,
    B_fine,
    D,
    norms2;
    normQ_fine = normQ_fine,
    dfly_coefficients = dfly_coefficients,
)
w_fine = invariant_vector(B_fine, Q_fine)
error_fine = distance_from_invariant(
    B_fine,
    D,
    Q_fine,
    w_fine,
    norms_fine;
    dfly_coefficients = dfly_coefficients,
)

@info "The fine error is $error_fine"

Bound = RigorousInvariantMeasures.invariant_measure_strong_norm_bound(
    B,
    D;
    dfly_coefficients = dfly_coefficients,
)
@info "The density is bounded above by $Bound"

observable(x) = log(abs(θ * α * abs(x - 0.5)^(α - 1)))
N = length(B_fine)
BoundObs = observable(Interval(B_fine.p[(N÷2+2)]))
v = discretizationlogderLorenz(B_fine, θ, α)

@info "integral small interval $(v[N÷2]/N)"
@info "BoundObs $BoundObs"

notLinf = (2 * (v[N÷2] / N) * Bound + BoundObs * error_fine).hi
@info "notLinf $notLinf"
lyap = (v' * w_fine) / N + Interval(-notLinf, notLinf)
@info "lyap $lyap"



@info "Computed $lyap, with precision $(diam(lyap))"


write(file, "Lorenz", lyap)

# ϕ = discretizationlogder(B_fine, D, degree = 3)
#     lyap = integrateobservable(B_fine, ϕ, w_fine, error_fine)
#     end
#     
# end 

close(file)
