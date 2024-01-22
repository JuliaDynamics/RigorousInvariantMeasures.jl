using Test
using IntervalArithmetic
using RigorousInvariantMeasures

@testset "Differentiation Interface" begin

    f = x -> x^3
    @test value_derivative_and_second_derivative(f, 2.0) == (8, 12, 12)
    @test value_derivative_and_second_derivative(f, 1 .. 2) ==
          (@interval(1, 8), @interval(3, 12), @interval(6, 12))
    @test value_and_derivative(f, 2.0) == (8, 12)
    @test value_and_derivative(f, 1 .. 2) == (@interval(1, 8), @interval(3, 12))

    g = @define_with_derivatives x -> 31 x -> 42 x -> 67
    @test value_derivative_and_second_derivative(g, 2 .. 4) == (31, 42, 67)
    @test value_and_derivative(g, 1 .. 2) == (31, 42)

end
