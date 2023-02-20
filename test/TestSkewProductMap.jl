@testset "Skew Product Map" begin
    using RigorousInvariantMeasures
    
    T = mod1_dynamic(x->2*x)

    D = SkewProductMap(T, [(x,y)->y*(x-0.5)^2/4+0.25;(x,y)->y*(x-0.5)^2/4+0.75])

    @test RigorousInvariantMeasures.FiberMap(D, Interval(0.1), 0.2) ⊂ Interval(0.257999, 0.258001)
    @test_throws "Intersects many branches, ambiguous" RigorousInvariantMeasures.FiberMap(D, Interval(0.4, 0.6), 0.2) 
end