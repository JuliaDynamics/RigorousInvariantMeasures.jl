using Test

@testset "preimages" begin

# Testing skip_beginning / skip_end

a1 = LinRange(0,1,12)
a2 = [0..0.1, 0..0.2, 0.1..0.3, 0.2..0.4, 0.3..0.5, 0.3..0.6, 0.6..0.8, 0.8..1, 1..1] #tricky weakly sorted vector

x1 = 0.25..0.25
x2 = 0.1..0.8 
x3 = 2
x4 = -1
x5 = 0
x6 = 0.3

for a = (a1, a2)
    for x = (x1, x2, x3, x4, x5,  x6)
        seq = InvariantMeasures.PointSequence(a, 37)
        (s, e) = InvariantMeasures.skipandlast(seq, x)
        @test s==0 || @interval(a[s]) ≺ @interval(x)
        @test s==length(a) || !(@interval(a[s+1]) ≺ @interval(x))
        @test e==0 || !(@interval(x) ≺ @interval(a[e]))
        @test e==length(a) || @interval(x) ≺ @interval(a[e+1])
    end
end

for a = (a1, a2)
    for x = (x1, x2, x3, x4, x5,  x6)
        reva = a[end:-1:1]
        seq = InvariantMeasures.PointSequence(reva, 37)
        (s, e) = InvariantMeasures.skipandlast(seq, x)
        @test s==0 || @interval(x) ≺ @interval(reva[s])
        @test s==length(reva) || !(@interval(x) ≺ @interval(reva[s+1]))
        @test e==0 || !(@interval(reva[e]) ≺ @interval(x))
        @test e==length(reva) || @interval(reva[e+1]) ≺ @interval(x)
    end
end

# checking for equality with looser tolerance
function approxintervals(a, b)
    return a.lo ≈ b.lo && a.hi ≈ b.hi
end

# checks various combinations of increasing and decreasing branches and sequences
for a = (a1, a2, a1[end:-1:1], a2[end:-1:1])
    for f in (x->x/2, x->1+x/2, x-> 1-x/2)
        b = InvariantMeasures.Branch(f, (@interval(0), @interval(1)))
        seq = InvariantMeasures.PointSequence(a, 3)
        preseq = preimages(seq, b)
        @test all(approxintervals.(b.f.(preseq.v), intersect.(map(Interval, a), hull(b.Y...))[preseq.skip-seq.skip+1:preseq.skip-seq.skip+length(preseq.v)]))
    end
end

for a = (a1, a2, a1[end:-1:1], a2[end:-1:1])
        D = PwMap([x->1/3 + 1/4*x, x->1/3-1/4*x, x->2x], [@interval(0), @interval(0.1), @interval(0.2), @interval(1)])
        for branch in InvariantMeasures.branches(D)
            B = Ulam(a)
            d = InvariantMeasures.duals(B, branch)
            # TODO: tricky to set up a good test, we'll resort to integration later?
        end
end

end #testset