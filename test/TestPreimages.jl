using Test

@testset "preimages" begin

# Testing skip_beginning / skip_end

a1 = LinRange(0,1,12)
a2 = [0..0.1, 0..0.2, 0.1..0.3, 0.2..0.4, 0.3..0.5, 0.3..0.6, 0.6..0.8, 0.8..1, 1..1] #tricky weakly sorted vector
a3 = fill(0..1, 10)

x1 = 0.25..0.25
x2 = 0.1..0.8 
x3 = 2
x4 = -1
x5 = 0
x6 = 0.3

for a = (a1, a2, a3)
    for x = (x1, x2, x3, x4, x5,  x6)
        s = skip_beginning(a, x, true)
        @test s==0 || @interval(a[s]) ≺ @interval(x)
        @test s==length(a) || !(@interval(a[s+1]) ≺ @interval(x))
        e = last_end(a, x, true)
        @test e==0 || !(@interval(x) ≺ @interval(a[e]))
        @test e==length(a) || @interval(x) ≺ @interval(a[e+1])
    end
end

for a = (a1, a2, a3)
    for x = (x1, x2, x3, x4, x5,  x6)
        reva = a[end:-1:1]
        s = skip_beginning(reva, x, false)
        @test s==0 || @interval(x) ≺ @interval(reva[s])
        @test s==length(reva) || !(@interval(x) ≺ @interval(reva[s+1]))
        e = last_end(reva, x, false)
        @test e==0 || !(@interval(reva[e]) ≺ @interval(x))
        @test e==length(reva) || @interval(reva[e+1]) ≺ @interval(x)
    end
end

# checking for equality with looser tolerance
function approxintervals(a, b)
    return a.lo ≈ b.lo && a.hi ≈ b.hi
end

# checks various combinations of increasing and decreasing functions
for a = (a1, a2, a1[end:-1:1], a2[end:-1:1])
    for f in (x->x/2, x->1+x/2, x-> 1-x/2)
        f = x-> 1-x/2
        X = 0..1
        (v, skip) = preimages(a, f, @interval(0,1))
        image = hull(f(@interval(X.lo)),f(@interval(X.hi)))
        @test all(approxintervals.(f.(v), intersect.(map(Interval, a), image)[skip+1:skip+length(v)]))
    end
end

end #testset