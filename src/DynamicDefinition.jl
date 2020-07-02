module DynamicDefinition

using ValidatedNumerics
import TaylorSeries
export Dynamic, der, der_der, der_n, preim, nbranches, plottable

abstract type Dynamic end


nbranches(S::Dynamic) = @error "Not implemented"
plottable(S::Dynamic) = @error "Not implemented"
preim(S::Dynamic, k, y, ϵ) = @error "Not implemented"


der(S::Dynamic, x) = der_n(S, x, Val(1))
der_der(S::Dynamic, x) = der_n(S, x, Val(2))

der_n(S::Dynamic, x, n) = der_n(S, x, Val(n))
der_n(S::Dynamic, x, ::Val{n}) where {n} = S.T(TaylorSeries.Taylor1([x, 1], n))[n]/factorial(n)

end