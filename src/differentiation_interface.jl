using TaylorSeries

# TODO: explore if ForwardDiff(), DualNumbers() or other packages work better than TaylorSeries

export value_and_derivative,
    value_derivative_and_second_derivative,
    derivative,
    derivative_and_second_derivative,
    second_derivative,
    distortion,
    inverse_derivative,
    @define_with_derivatives,
    check_derivatives


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

Remark: when specializing derivatives by hand, always specialize the two-argoment functions 
(e.g., `value_and_derivative(f, x)`) rather than the one-parameter ones  (e.g., `value_and_derivative(f)`),
since the latter are defined generically in terms of the former.

Remark: tricky point that can cause subtle bugs: functions created in the same line-of-code
will often have the same type; so for instance
# ```jldoctest
# julia> fs = [x -> k*x for k in 1:3];  # the three elements here have the same type

# julia> typeof(fs[1]) == typeof(fs[3])
# true

# julia> RigorousInvariantMeasures.derivative(f::typeof(fs[1]), x) = 1;

# julia> RigorousInvariantMeasures.derivative(f::typeof(fs[2]), x) = 2;

# julia> RigorousInvariantMeasures.derivative(f::typeof(fs[3]), x) = 3;

# julia> RigorousInvariantMeasures.derivative(fs[1], 0.5)  # this returns 3, not 1
# 3
# ```
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
The inverse_derivative of f is 1/f'.
"""
function inverse_derivative(f, x)
    return 1 / derivative(f, x)
end
inverse_derivative(f) = x -> inverse_derivative(f, x)

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

# ```jldoctest
# julia> f = @define_with_derivatives x->x^2 x->2x x->2;
# ERROR: LoadError: UndefVarError: `@define_with_derivatives` not defined
# in expression starting at none:1

# julia> value_and_derivative(f, 45)
# ERROR: UndefVarError: `value_and_derivative` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1
# ```
Note that the three macro parameters are separated just by spaces (no commas or parentheses)

# ```jldoctest
# julia> g = @define_with_derivatives x->12 x->34 x->56;
# ERROR: LoadError: UndefVarError: `@define_with_derivatives` not defined
# in expression starting at none:1

# julia> value_derivative_and_second_derivative(g, -23.5)
# ERROR: UndefVarError: `value_derivative_and_second_derivative` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1
# ```

This is provided for convenience, but note that in many cases one can find 
common subexpressions in a function and its derivatives; hence
it is more efficient to define the two functions separately.

If you define derivatives in this way, it is recommended to run `check_derivatives` to ensure
that you did not make any mistakes (e.g., forgetting a factor 2).
We do not run it automatically because that would require knowing a valid `x` in the domain of `f`.
"""
macro define_with_derivatives(f, df, ddf)
    return quote
        local g = $f
        RigorousInvariantMeasures.value_and_derivative(::typeof(g), x) = ($f(x), $df(x))
        RigorousInvariantMeasures.value_derivative_and_second_derivative(::typeof(g), x) =
            ($f(x), $df(x), $ddf(x))
        RigorousInvariantMeasures.derivative(::typeof(g), x) = $df(x)
        RigorousInvariantMeasures.second_derivative(::typeof(g), x) = $ddf(x)
        RigorousInvariantMeasures.derivative_and_second_derivative(::typeof(g), x) =
            ($df(x), $ddf(x))
        g
    end
end


using Random
"""
    check_derivatives(f, x=rand())

Checks (using assertions) that the derivatives of f agree (up to the square root of machine precision) 
with those computed by TaylorSeries.

This is a useful sanity check if you redefine derivatives.
    
The derivative are checked in a point `x` (default: rand()), which should be in the domain of `f`.

# ```jldoctest
# julia> f = @define_with_derivatives x->x^2 x->2x x->2;
# ERROR: LoadError: UndefVarError: `@define_with_derivatives` not defined
# in expression starting at none:1

# julia> check_derivatives(f, 0.2)
# ERROR: UndefVarError: `check_derivatives` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1

# julia> g = @define_with_derivatives x->x^2 x->2x x->3;
# ERROR: LoadError: UndefVarError: `@define_with_derivatives` not defined
# in expression starting at none:1

# julia> check_derivatives(g, 0.2)
# ERROR: UndefVarError: `check_derivatives` not defined
# Stacktrace:
#  [1] top-level scope
#    @ none:1
# ```
"""
function check_derivatives(f, x = rand())
    y = f(Taylor1([x, one(x), zero(x)]))
    fx, dfx, ddfx = y[0], y[1], 2y[2]

    @assert all(value_and_derivative(f, x) .≈ (fx, dfx))
    @assert all(value_derivative_and_second_derivative(f, x) .≈ (fx, dfx, ddfx))
    @assert derivative(f, x) ≈ dfx
    @assert second_derivative(f, x) ≈ ddfx
    @assert all(derivative_and_second_derivative(f, x) .≈ (dfx, ddfx))
    @assert inverse_derivative(f, x) ≈ 1 / dfx
    @assert distortion(f, x) ≈ ddfx / dfx^2
end
