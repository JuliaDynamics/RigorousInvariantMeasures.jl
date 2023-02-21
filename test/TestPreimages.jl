using Test
using RigorousInvariantMeasures: preimages
using RigorousInvariantMeasures.DynamicDefinition: is_increasing
using IntervalArithmetic

@testset "preimages" begin

# Testing skip_beginning / skip_end

a1 = LinRange(0,1,12)
a2 = [0, 0..0.1, 0..0.2, 0.1..0.3, 0.2..0.4, 0.3..0.5, 0.3..0.6, 0.6..0.8, 0.8..1, 1..1] #tricky weakly sorted vector

x1 = 0.25..0.25
x2 = 0.1..0.8 
x3 = 2
x4 = -1
x5 = 0
x6 = 0.3

# test for the for wide intervals and 0 derivative 

@test RigorousInvariantMeasures.range_estimate_monotone(x->x, Interval(0, 1)) == Interval(0, 1)
@test RigorousInvariantMeasures.range_estimate_monotone(x->exp(x), Interval(0, 1)) == Interval(1, exp(Interval(1)).hi)

begin 
    preim_branch = RigorousInvariantMeasures.preimage
    br = RigorousInvariantMeasures.MonotonicBranch(x-> x^2, (@interval(0),@interval(1)))

    ϵ = 1.0/1024
    max_iter = 100
    
    
    # if X is the search interval we test when f(X) ⊂ y
    root = preim_branch(Interval(-5, 5), br, @interval(0.0, 1.0); ϵ, max_iter = 100)
    @test root ⊂ (@interval(0.0, 1.0)+@interval(-ϵ , ϵ))

    # if X is the search interval we test when y ⊂ f(X)
    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.0, 1.0); ϵ, max_iter = 100)
    @test root ⊂ (@interval(0, 0.2)+@interval(-ϵ , ϵ))

    # if X is the search interval we test when y ∩ f(X) ≂̸ ∅, the expected result is
    # f^{-1}(y) ∩ X
    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.1, 1.0); ϵ, max_iter = 100)
    @test root ⊂ (@interval(0.1, 0.2)+@interval(-ϵ, ϵ))

    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.2, 1.0); ϵ, max_iter = 100)
    @test root ⊂ (@interval(0.2)+@interval(-ϵ, ϵ))

    # if X is the search interval we test when y ∩ f(X) = ∅ 
    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.3, 1.0); ϵ, max_iter = 10)
    @test isempty(root)
    
    # we test when the domain of f contains the zero of the derivative

    br = RigorousInvariantMeasures.MonotonicBranch(x-> (x-@interval(0.1))^2, (@interval(0.1),@interval(1)))

    # if X is the search interval we test when y ⊂ f(X)
    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.1, 1.0); ϵ, max_iter = 100)
    @test root ⊂ (@interval(0.1, 0.3)+@interval(-ϵ , ϵ))

    # if X is the search interval we test when y ∩ f(X) ≂̸ ∅, the expected result is
    # f^{-1}(y) ∩ X
    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.2, 1.0); ϵ, max_iter = 100)
    @test root ⊂ (@interval(0.2, 0.3)+@interval(-ϵ, ϵ))

    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.0); ϵ, max_iter = 100)
    @test root ⊂ (@interval(0.0)+@interval(-ϵ, ϵ))

    # if X is the search interval we test when y ∩ f(X) = ∅ 
    root = preim_branch(Interval(0.0, 0.04), br, @interval(0.4, 1.0); ϵ, max_iter = 10)
    @test isempty(root)
    
    # we test the exit rule for Krawczyk, i.e., if 0 ∈ f′(x_mid)
    # return x
    # remark that to bypass the unique increasing test we need to call 
    # the full MonotonicBranch constructor. In general this function would not be allowed
    # by the constructor. MonotonicBranch are monotone in Interval(X[1].hi, X[2].lo)
    br = RigorousInvariantMeasures.MonotonicBranch(x-> x^2, 
                                    (@interval(-1), @interval(1)),
                                    (Interval(1), Interval(1)),
                                    true)
    root = preim_branch(Interval(0.0), br, @interval(0.0); ϵ, max_iter = 10)
    @test root == Interval(0.0)


end

for a = (a1, a2)
    for x = (x1, x2, x3, x4, x5,  x6)
        i = RigorousInvariantMeasures.first_overlapping(a, x)
        @test i==0 || Interval(a[i]).hi <= Interval(x).lo
        @test i==length(a) || !(Interval(a[i+1]).hi <= Interval(x).lo)

        j = RigorousInvariantMeasures.last_overlapping(a, x)
        @test j==0 || !(Interval(a[j]).lo >= Interval(x).hi)
        @test j==length(a) || Interval(a[j+1]).lo >= Interval(x).hi
    end
end

# checking for equality with looser tolerance
function approxintervals(a, b)
    return a.lo ≈ b.lo && a.hi ≈ b.hi
end

for y = (a1, a2)
    for f in (x->x/2, x->1+x/2, x-> 1-x/2)
        b = RigorousInvariantMeasures.MonotonicBranch(f, (@interval(0), @interval(1)))
        ylabel = 1:length(y)
        x, xlabel = RigorousInvariantMeasures.preimages(y, b, ylabel; ϵ = 1e-14, max_iter = 100)
        @test x[1] == b.X[1]
        @test x[end] != b.X[2] # to make sure the last entry isn't there
        @test length(x) == length(xlabel)
        z1 = filter(x->!isempty(x),intersect.(map(Interval,y), hull(b.Y[1], b.Y[2]+1e-15))) # the 1e-15 is there to simulate "disjoint intervals"
        z2 = f.(x)
        if !b.increasing
            z2 = reverse(z2)
        end
        
        @test all(approxintervals.(z1, z2))
    end
end

# Test nonmonotone 
f = x->(x-@interval(0.3))^2
b = RigorousInvariantMeasures.MonotonicBranch(f, (@interval(0.3), @interval(1)))
x = RigorousInvariantMeasures.preimage(0.04, b, Interval(0.2, 0.7), ϵ = 1e-14, max_iter = 100)        
@test 0.5 ∈ x

b = RigorousInvariantMeasures.MonotonicBranch(f, (@interval(0), @interval(0.3)))
x = RigorousInvariantMeasures.preimage(0.04, b, Interval(0, 0.4), ϵ = 1e-14, max_iter = 100)        
@test 0.1 ∈ x

#@test_logs (:debug,"Not contracting, fallback to bisection") RigorousInvariantMeasures.preimage(0.04, b, Interval(0, 0.4), ϵ = 1e-14, max_iter = 100)

D = mod1_dynamic(x -> 2x)
DD = ∘(D, D, D, D)
p = [0, 0.2, 0.4, 0.6, 0.8]

x, xlabel = preimages(p, DD; ϵ =  1e-13, max_iter = 100)
@test x ≈ 0:1/80:79/80
@test xlabel == repeat([1,2,3,4,5],16)

# test composing two different functions, to check that composition is handled in the correct order

f = x-> 2x
g = x -> (x+x^2)/2

D1 = mod1_dynamic(f)
D2 = mod1_dynamic(g)

@test RigorousInvariantMeasures.DynamicDefinition.is_increasing(D1)

y = 0:0.2:1
x, xlabel = preimages(y, D1 ∘ D2; ϵ =  1e-13, max_iter = 100)

@test f.(g.(x)) ≈ 0:0.2:1.8
@test xlabel ≈ repeat(1:5, 2)

D1 = PwMap([x->2x, x->6x-3, x->3x-2], [0, 0.5, @interval(2/3), 1], [0 1; 0 1; 0 1])
D2 = PwMap([x->2x, x->4x-2, x->4x-3], [0, 0.5, 0.75, 1], [0 1; 0 1; 0 1])

@test RigorousInvariantMeasures.DynamicDefinition.is_increasing(D1)

z = 0:0.3:1
y, ylabel, y′ = RigorousInvariantMeasures.preimages_and_derivatives(z, D1; ϵ = 0.0, max_iter = 100)

@test y ≈ [0, 0.15, 0.3, 0.45, 0.5, 0.5+0.05, 0.5+0.1, 0.5+0.15, 2/3, 2/3+0.1, 2/3+0.2, 2/3+0.3]
@test ylabel == repeat(1:4, 3)
@test y′ ≈ [2,2,2,2,6,6,6,6,3,3,3,3]

x, xlabel, x′ = RigorousInvariantMeasures.preimages_and_derivatives(z, D1∘D2; ϵ = 0.0, max_iter = 100)

@test all(x .≈ vcat(y/2, 0.5 .+ y/4, 0.75 .+ y/4))
@test x′ == kron([4,12,6,8,24,12,8,24,12], [1,1,1,1])

f = x -> x^2;
g = x -> 1-x^2
Dinc = PwMap([f], [0..0, 1..1])
Ddec = PwMap([g], [0..0, 1..1])

@test is_increasing(Dinc)
@test !is_increasing(Ddec)
@test is_increasing(Ddec ∘ Ddec)
@test !is_increasing(Dinc ∘ Ddec)
@test !is_increasing(Ddec ∘ Ddec ∘ Ddec)

z = [0..0, 0.5..0.5, 1..1]
x, xlabel, x′ = RigorousInvariantMeasures.preimages_and_derivatives(z, Dinc∘Ddec; ϵ = 0.0, max_iter = 100)

@test all(isapprox(x, [0, sqrt(1 - 1/sqrt(2))]))
@test all(xlabel .== [2,1])
@test all(isapprox.(x′, [0, -2*sqrt(2-sqrt(2))]))

x, xlabel, x′ = RigorousInvariantMeasures.preimages_and_derivatives(z, Ddec∘Ddec; ϵ = 0.0, max_iter = 100)
@test all(isapprox(x, [0, sqrt(1 - 1/sqrt(2))]))
@test all(xlabel .== [1,2])
@test all(isapprox.(x′, [0, 2*sqrt(2-sqrt(2))]))

D = PwMap([x->17*x/5, 
	x->(34*((17*x-5)/17)/25+3)*((17*x-5)/17), 
	x->(34*((17*x-10)/17)/25+3)*((17*x-10)/17), 
	x->17*((17*x-15)/17)/5], 
	[Interval(0), Interval(5)/17, Interval(10)/17, Interval(15)/17, Interval(1)],
	[Interval(0) Interval(1);
	 Interval(0) Interval(1);
	 Interval(0) Interval(1);
	 Interval(0) @interval(0.4)]
	)

# we just check that this doesn't throw, for now
RigorousInvariantMeasures.preimages(0:0.25:1, D; ϵ = 0.0, max_iter = 100)

# test has_infinite_derivative_at_endpoints
using RigorousInvariantMeasures: has_infinite_derivative_at_endpoints

D = PwMap([x->17*x/5, 
	x->(34*((17*x-5)/17)/25+3)*((17*x-5)/17), 
	x->(34*((17*x-10)/17)/25+3)*((17*x-10)/17), 
	x->17*((17*x-15)/17)/5], 
	[Interval(0), Interval(5)/17, Interval(10)/17, Interval(15)/17, Interval(1)],
	[Interval(0) Interval(1);
	 Interval(0) Interval(1);
	 Interval(0) Interval(1);
	 Interval(0) @interval(0.4)]
	)

for i in 1:length(D.branches)
    @test has_infinite_derivative_at_endpoints(D.branches[i]) == (false, false)
end
@test !has_infinite_derivative_at_endpoints(D)

D0 = Lorenz()
D = D0 ∘ D0 ∘ D0

@test has_infinite_derivative_at_endpoints(D.E.branches[1]) == (false, true)
for i = 2:7
    @test has_infinite_derivative_at_endpoints(D.E.branches[i]) == (true, true)
end
@test has_infinite_derivative_at_endpoints(D.E.branches[end]) == (true, false)

@test has_infinite_derivative_at_endpoints(D.E)

D = BZ()
@test has_infinite_derivative_at_endpoints(D)

end #testset