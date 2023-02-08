using TaylorSeries

# TODO: explore if ForwardDiff(), DualNumbers() or other packages work better

"""

    value_and_derivative(f, x)

Evaluate f(x) and f'(x) at the same time. 
We define a single method to do this because often there are shared subexpressions between a function and its derivatives,
or they are computed together by automatic differentiation packages, so this is often more efficient than having two
separate functions for f and f' and calling them separately.

The generic implementation uses `TaylorSeries` to compute derivatives. It can be overwritten for individual functions:
```
f = x -> x^2
value_and_derivative(f::typeof(f), x) = (x^2, 2x)
```

See @define_with_derivatives to define derivatives with three separate expressions.
"""
function value_and_derivative(f, x)
    y = f(Taylor1([x, one(x)]))
    return y[0], y[1]
end

"""
    value_derivative_and_second_derivative(f, x)

Generic method to evaluate f(x), f'(x) and f''(x) at the same time.

See `value_and_derivative` for more detail.
"""
function value_derivative_and_second_derivative(f, x)
    y = f(Taylor1([x, one(x), zero(x)]))
    return y[0], y[1], 2y[2]
end

function derivative(f, x)
    y = f(Taylor1([x, one(x)]))
    return y[1]
end
derivative(f) = x -> derivative(f, x)

function second_derivative(f, x)
    y = f(Taylor1([x, one(x), zero(x)]))
    return 2y[2]
end
second_derivative(f) = x -> second_derivative(f, x)

function derivative_and_second_derivative(f, x)
    y = f(Taylor1([x, one(x), zero(x)]))
    return y[1], 2y[2]
end
derivative_and_second_derivative(f) = x -> derivative_and_second_derivative(f, x)

"""
The expansivity of f is 1/f'.
"""
function expansivity(f, x)
    return 1 / derivative(f, x)
end
expansivity(f) = x -> expansivity(f, x)

"""
The distortion of f is f'' / (f')^2.
"""
function distortion(f, x)
    dy, ddy = derivative_and_second_derivative(f, x)
    return ddy / dy^2
end
distortion(f) = x -> distortion(f, x)

"""
    f = @define_with_derivatives(ex1, ex2, ex3)

Declares a new function f, and redefines the various functions related to derivatives
so that its derivatives are computed with the given expressions.

```jldoctest
julia> f = @define_with_derivatives x->x^2 x->2x x->1;

julia> value_and_derivative(f, 45)
(2025, 90)
```
Note that the three macro parameters are separated just by spaces (no commas or parentheses)

```jldoctest
julia> g = @define_with_derivatives x->12 x->34 x->56;

julia> value_derivative_and_second_derivative(g, -23.5)
(12, 34, 56)
```

This is provided for convenience, but note that in many cases one can find 
common subexpressions in a function and its derivatives; hence
it is more efficient to define the two functions separately.
"""
macro define_with_derivatives(f, df, ddf)
    return quote
        local g = $f
        RigorousInvariantMeasures.value_and_derivative(::typeof(g), x) = ($f(x), $df(x))
        RigorousInvariantMeasures.value_derivative_and_second_derivative(::typeof(g), x) = ($f(x), $df(x), $ddf(x))
        RigorousInvariantMeasures.derivative(::typeof(g), x) = $df(x)
        RigorousInvariantMeasures.second_derivative(::typeof(g), x) = $ddf(x)
        RigorousInvariantMeasures.derivative_and_second_derivative(::typeof(g), x) = ($df(x), $ddf(x))
        g
    end
end

export value_and_derivative, value_derivative_and_second_derivative, derivative, 
        derivative_and_second_derivative, second_derivative, distortion, expansivity, @define_with_derivatives
