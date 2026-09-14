###############################################################################
# Directed-rounding arithmetic
#
# The package computes its error bounds with the FastRounding operators, which
# are defined for `Float32` and `Float64` only, through error-free
# transformations. That left every bound Float64-only: `gamma`, the noise
# kernels and everything downstream stopped with a `MethodError` as soon as the
# bound type was `BigFloat`.
#
# MPFR does honour `setrounding`, which Julia no longer supports for `Float64`,
# so the two cases have to be written differently: `Float32` and `Float64`
# delegate to FastRounding unchanged, and `BigFloat` uses the
# `setrounding ... do` blocks that BallArithmetic uses throughout.
#
# These are the package's own functions rather than new methods on
# FastRounding's, since adding `BigFloat` methods to a function another package
# owns, on a type Base owns, would apply to every package sharing the session.
# Nothing at the call sites changes: the names resolve here because no file
# under `src/` brings FastRounding's operators into scope any more.
###############################################################################

import FastRounding
using FastRounding: sqrt_round, square_round

const _FastFloat = Union{Float32,Float64}

for (up, down, tozero, base) in (
    (:⊕₊, :⊕₋, :⊕₀, :+),
    (:⊖₊, :⊖₋, :⊖₀, :-),
    (:⊗₊, :⊗₋, :⊗₀, :*),
    (:⊘₊, :⊘₋, :⊘₀, :/),
)
    for (op, mode) in ((up, RoundUp), (down, RoundDown), (tozero, RoundToZero))
        @eval begin
            $op(x::T, y::T) where {T<:_FastFloat} = FastRounding.$op(x, y)

            $op(x::BigFloat, y::BigFloat) = setrounding(BigFloat, $mode) do
                $base(x, y)
            end
        end
    end
end
