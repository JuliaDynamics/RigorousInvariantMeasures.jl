@testset "Skew Product Map" begin
    using RigorousInvariantMeasures
    
    T = mod1_dynamic(x->2*x)

    D = SkewProductMap(T, [(x,y)->y*(x-0.5)^2/4+0.25;(x,y)->-y*(x-0.5)^2/4+0.75])

    @test RigorousInvariantMeasures.FiberMap(D, Interval(0.1), 0.2) ⊂ Interval(0.257999, 0.258001)
    @test_throws "Intersects many branches, ambiguous" RigorousInvariantMeasures.FiberMap(D, Interval(0.4, 0.6), 0.2) 

    @test RigorousInvariantMeasures.intersect_domain_bool(T, 0.1) == [true; false]    

    @info RigorousInvariantMeasures.preimage_fixed_x(D, 1, 0.1, 0.25, 0.2525; ϵ = 0.01, max_iter = 100)
    @test all((Interval(0), Interval(1)/16) .⊆ RigorousInvariantMeasures.preimage_fixed_x(D, 1, 0.1, 0.25, 0.2525; ϵ = 0.01, max_iter = 100))  
    @test all((Interval(0), Interval(1)) .⊆ RigorousInvariantMeasures.preimage_fixed_x(D, 1, 0.1, 0.24, 0.3; ϵ = 0.01, max_iter = 100)) 
    @test all((Interval(0), Interval(1)/8) .⊆  RigorousInvariantMeasures.preimage_fixed_x(D, 2, 0.6, 0.7496875, 0.75; ϵ = 0.01, max_iter = 100)) 
    @test all((Interval(0), Interval(1)) .⊆ RigorousInvariantMeasures.preimage_fixed_x(D, 2, 0.6, 0.7, 0.76; ϵ = 0.01, max_iter = 100)) 
end