using Test

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

for a = (a1, a2)
    for x = (x1, x2, x3, x4, x5,  x6)
        i = InvariantMeasures.first_overlapping(a, x)
        @test i==0 || Interval(a[i]).hi <= Interval(x).lo
        @test i==length(a) || !(Interval(a[i+1]).hi <= Interval(x).lo)

        j = InvariantMeasures.last_overlapping(a, x)
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
        b = InvariantMeasures.Branch(f, (@interval(0), @interval(1)))
        ylabel = 1:length(y)
        x, xlabel = preimages(y, b, ylabel)
        @test x[1] == b.X[1]
        @test x[end] != b.X[2] # to make sure the last entry isn't there
        @test length(x) == length(xlabel)
        y1 = filter(x->!isempty(x),intersect.(map(Interval,y), hull(b.Y[1], b.Y[2]+1e-15))) # the 1e-15 is there to simulate "disjoint intervals"
        y2 = f.(x)
        if !b.increasing
            y2 = reverse(y2)
        end
        @test all(approxintervals.(y1, y2))
    end
end

end #testset