using Test
using IntervalArithmetic
using RigorousInvariantMeasures: is_full_branch
using RigorousInvariantMeasures.DynamicDefinition: is_increasing

@testset "Dynamics" begin
using IntervalArithmetic

D = mod1_dynamic(x->2*x)

@test D.branches[1].f(0.1) == 0.2
@test derivative(D.branches[1].f)(0.1) == 2.0

@test RigorousInvariantMeasures.dfly(RigorousInvariantMeasures.Lipschitz, RigorousInvariantMeasures.L1, D) == (0.5, 0.0)
@test RigorousInvariantMeasures.dfly(RigorousInvariantMeasures.TotalVariation, RigorousInvariantMeasures.L1, D) == (0.5, 0.0)

D = PwMap([ x-> x^2+0.25, x -> 4*x-2, x -> 4*x-3], [0, 0.5, 0.75, 1])

D = PwMap([x->2*x, x->2-2*x], [@interval(0), @interval(0.5), @interval(1)])

@test is_increasing(D.branches[1]) == true
@test is_increasing(D.branches[2]) == false

D = mod1_dynamic(x -> 3.5x)

@test length(D.branches) == 4
@test [D.branches[i].X[1] for i in 1:4] ≈ [0, 2/7, 4/7, 6/7]
@test [D.branches[i].X[2] for i in 1:4] ≈ [2/7, 4/7, 6/7, 1]
@test  D.branches[2].f(0.5) == 0.75
@test [D.branches[i].Y[1] for i in 1:4] == [0, 0, 0, 0]
@test [D.branches[i].Y[2] for i in 1:4] == [1, 1, 1, 0.5]

@test [is_increasing(D.branches[i]) for i in 1:4] == [1, 1, 1, 1]

@test is_full_branch(D) == false

D = mod1_dynamic(x -> 3.5x + 0.5)
@test [D.branches[i].X[1] for i in 1:4] ≈ [0, 1/7, 3/7, 5/7]
@test [D.branches[i].X[2] for i in 1:4] ≈ [1/7, 3/7, 5/7, 1]
@test  D.branches[end].f(1) == 1
@test [D.branches[i].Y[1] for i in 1:4] == [0.5, 0, 0, 0]
@test [D.branches[i].Y[2] for i in 1:4] == [1, 1, 1, 1]
@test [is_increasing(D.branches[i]) for i in 1:4] == [1, 1, 1, 1]

@test is_full_branch(D) == false

D = mod1_dynamic(x -> -3.5x + 0.5)
@test [D.branches[i].X[1] for i in 1:4] ≈ [0, 1/7, 3/7, 5/7]
@test [D.branches[i].X[2] for i in 1:4] ≈ [1/7, 3/7, 5/7, 1]
@test  D.branches[end].f(5/7) == 1
@test [D.branches[i].Y[1] for i in 1:4] == [0.5, 1, 1, 1]
@test [D.branches[i].Y[2] for i in 1:4] == [0, 0, 0, 0]
@test [is_increasing(D.branches[i]) for i in 1:4] == [0, 0, 0, 0]

D0 = mod1_dynamic(x->2*x)
D = D0 ∘ D0

@test is_full_branch(D0) == true
@test is_full_branch(D) == true
@test RigorousInvariantMeasures.DynamicDefinition.domain(D) == (Interval(0), Interval(1))

A, B, C = RigorousInvariantMeasures.preimages_and_derivatives([0.,0.1], D; ϵ =  1e-13, max_iter = 100)
@test A ≈ [0, 0.025, 0.25, 0.275, 0.5, 0.525, 0.75, 0.775]
@test B == [1, 2, 1, 2, 1, 2, 1, 2]
@test C == fill(4, 8)

T = mod1_dynamic(x-> 1-2*x)

D = D0∘T
@test is_full_branch(D) == true

A, B, C = RigorousInvariantMeasures.preimages_and_derivatives([0.,0.1], D; ϵ =  1e-13, max_iter = 100)
@test A ≈ [0, 0.225, 0.25, 0.475, 0.5, 0.725, 0.75, 0.975]
@test B == [2, 1, 2, 1, 2, 1, 2, 1]
@test C == fill(-4, 8)

D = T∘D0
@test is_full_branch(D) == true
@test mid.(A) ≈ [0, 0.225, 0.25, 0.475, 0.5, 0.725, 0.75, 0.975]
@test B == [2, 1, 2, 1, 2, 1, 2, 1]
@test C == fill(-4, 8)

D = T ∘ T

@test is_full_branch(D) == true

A, B, C = RigorousInvariantMeasures.preimages_and_derivatives([0.,0.1], D; ϵ =  1e-13, max_iter = 100)
@test mid.(A) ≈ [0, 0.025, 0.25, 0.275, 0.5, 0.525, 0.75, 0.775]
@test B == [1, 2, 1, 2, 1, 2, 1, 2]
@test C == fill(4, 8)



# Lanford map, so that we test also something that is not linear
f = x->2*x+0.5*x*(1-x)
D0 = mod1_dynamic(f)
D = D0 ∘ D0 ∘ D0
A, B = RigorousInvariantMeasures.preimages([0., 0.5],D; ϵ =  1e-13, max_iter = 100)
g(x) = f(x) - floor(f(x))
@test g.(g.(g.(A[2:2:end]))) ≈ fill(0.5, 8)

@test endpoints(D) ≈ [0.0, 0.07389363935392047, 0.18200396341631753, 0.28117073603385473, 0.4384471871911697, 0.5287080193012084, 0.6633983486837269, 0.7902622123165944, 1.0]

@test [branch(D.E, k)(Interval(0.2)) for k in 1:nbranches(D)] ≈ [∅, ∅, g(g(g(Interval(0.2)))), ∅, ∅, ∅, ∅, ∅]

D = mod1_dynamic(x->2*x)

@test D.full_branch == true
@test D.branches[1].f(0.125) == 0.25
@test D.branches[2].f(0.5+0.125) == 0.25
@test derivative(D.branches[1].f, 0.1) == 2.0
@test derivative(D.branches[2].f, 0.1) == 2.0

@test 0 ∈ D.branches[1].X[1]
@test 0 ∈ D.branches[1].Y[1]

@test 0.5 ∈ D.branches[1].X[2]
@test 1 ∈ D.branches[1].Y[2]

@test 0.5 ∈ max_inverse_derivative(D)
@test 0 <= max_distortion(D)

# Testing composedPwMap

D1 = mod1_dynamic(x -> 2x)
D2 = mod1_dynamic(x -> 3x)

E = RigorousInvariantMeasures.composedPwMap(D1, D2)

@test all(endpoints(E) .≈ 0:1/6:1)
@test all(b.Y == (0,1) for b in E.branches)
@test E.branches[4].f(0.5) == 0
@test E.branches[4].f(2/3) == 1

D3 = mod1_dynamic(x -> (1-x)^2+(1-x))

E = RigorousInvariantMeasures.composedPwMap(D1, D3)
@test all([b.f(b.X[1]) ≈ 1. for b in E.branches])
@test all([isapprox(b.f(b.X[2]), 0., atol=1e-8) for b in E.branches])

x = 0.3
@test E.branches[1].f(x) == 2*((1-x)^2+(1-x)-1)-1
@test E.branches[2].f(x) == 2*((1-x)^2+(1-x)-1)
@test E.branches[3].f(x) == 2*((1-x)^2+(1-x))-1
@test E.branches[4].f(x) == 2*((1-x)^2+(1-x))
end
