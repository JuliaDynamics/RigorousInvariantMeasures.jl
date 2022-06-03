module Contractors
using ValidatedNumerics

export root, range_estimate, ShootingMethod, nthpreimage!, preimage, unique_sign, unique_increasing

"""
unique_sign(x)

Sign of an interval, but throws an error if it is not unique.
Used by various functions to compute orientations
"""
function unique_sign(x)
	s = sign(x)
	@assert isthin(s)
	return s.hi
end

"""
unique_increasing(a, b)

Given intervals a, b, returns `true` if a < b, `false` if b < a, and raises an error if it is not uniquely determined.
"""
function unique_increasing(a::Interval, b::Interval)
	if a.hi < b.lo
		return true
	elseif b.hi < a.lo
		return false
	else
		error("Insufficient precision to check the sign of this function")
	end
end
function unique_increasing(a, b) # Fallback for Float64
	if a < b
		return true
	elseif a > b
		return false
	else
		error("Could not determine sign")
	end
end

# this seems slower
using TaylorSeries
derivative(f) = x-> f(Taylor1([x,1.],1))[1]

#using DualNumbers
#derivative(f) = x->f(Dual(x, 1..1)).epsilon

"""
Compute a single root with (possibly multivariate) interval Newton

x must be an Interval (univariate) or IntervalBox (multivariate)

Stops when the interval reaches a fixed point, when the diameter is smaller than ε,
or when max_iter iterations are reached (with an error)
"""
root(f, x, ϵ; max_iter = 100) = root(f, derivative(f), x, ϵ; max_iter = max_iter)

function root(f, f′, x, ϵ; max_iter = 100)
	for i in 1:max_iter
		x_old = x
		x_mid = Interval(mid(x))
		x = intersect(x, x_mid - f′(x) \ f(x_mid))
		if x_old == x || isempty(x) || diam(x) < ϵ
			return x
		end
	end
	@info "Maximum iterates reached" max_iter, x, f(x)
	return x
end

preimage(y, f, X, ϵ; max_iter=100) = root(x -> f(x)-y, X, ϵ; max_iter)
preimage(y, f, fprime, X, ϵ; max_iter=100) = root(x -> f(x)-y, fprime, X, ϵ; max_iter)

# superseded by IntervalOptimisation.jl
function range_estimate(f, domain, recstep = 5)
	if recstep == 1
		return f(domain)
	else
		a, b = bisect(domain)
		Iₐ = range_estimate(f, a, recstep-1)
		Iᵦ = range_estimate(f, b, recstep-1)
		return Iₐ ∪ Iᵦ
	end
end

using LinearAlgebra

# this function generates the Jacobian for
# f(x_0)=x_1, f(x_1)=x_2, ..., f(x_{n-1})=y

coeff_interval(x::Array{Interval{T}, 1}) where {T} = T

function Jac(fprime, v::Vector{T}) where {T}
    dv = fprime.(v)
    ev = -ones(T, length(v)-1)
    return Bidiagonal{T}(dv, ev, :U)
end

# Used in InducedLSV; won't touch it to avoid breaking stuff --federico
function ShootingMethod(f, fprime, n, x, y, rigstep = 10)
	F = x->(f.(x)-[x[2:end]; y])

	for i in 1:rigstep
		x_mid = Interval{coeff_interval(x)}.(mid.(x))
		x = intersect.(x, x_mid-Jac(fprime, x)\F(x_mid))
	end
	return x
end

"""
Newer version of the 'shooting method' to compute the kth preimage of a point (or interval y)
fs contains k functions, X contains their domains. This computes a solution of
fk(f_{k-1}( ...  f1(x) ... )) = y.
Overwrites X with [x f1(x) f2(f1(x)) ... f_{k-1}(...)], so the
true solution is X[1].

Tries to avoid allocations and stuff.
"""
function nthpreimage!(y, fs, X; max_iter = 100)
	newX = zero(X)
	Xmid = zero(X)
	n = length(X)
	for i in 1:max_iter
		Xmid .= Interval.(mid.(X))
		newX[end] = (fs[end](Xmid[end]) - y) / derivative(fs[end])(X[end])
		for i = length(X)-1:-1:1
			newX[i] = (fs[i](Xmid[i]) - Xmid[i+1] + newX[i+1]) / derivative(fs[i])(X[i])
		end
		newX .= intersect.(Xmid - newX, X)
		if isempty(newX) || newX == X
			return newX
		end
		X .= newX
	end
	@info "Maximum iterates reached" max_iter
	return X
end

end #module
