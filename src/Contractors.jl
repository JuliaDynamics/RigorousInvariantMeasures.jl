module Contractors
using IntervalArithmetic


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
	#@info a, b
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
Compute the root of a monotonic function with interval Newton + bisection.

This function was originally written to work also with IntervalBox'es (multivariate Newton),
but in the end we do not really need it since switching to recursive preimage computation,
which should be more performant than the "shooting method".

The function first attempts to contract x with a step of interval Newton; if this fails, then
a bisection strategy is used.

Stops when the diameter is smaller than ε, or when `max_iter`` iterations are reached (with a warning),
or when the iterations fail to improve the interval (which may happen even with large `diam(x)` if the enclosures computed by `f`
are particularly poor).

The method should be quite robust to the value of `ϵ`, since it is not used as the main stopping criterion.

"""
root(f, x; ϵ, max_iter) = root(f, derivative(f), x; ϵ, max_iter)

function root(f, f′, x; ϵ, max_iter)
	f_endpoints_computed = false
	for i in 1:max_iter
		# attempt 1: Newton

		x_old = x

		x_mid = mid(x)
		Ix_mid = Interval(x_mid)
		f_mid = f(Ix_mid)

		@debug "Step $i of the Newton method:
			   - the interval $x, 
			   - the derivative $(f′(x)),
			   - the value at the midpoint $(f(x_mid))"
		
		Nx = Ix_mid - f′(x) \ f_mid
		@debug "Candidate", Nx 

		x = intersect(x, Nx)

		# if the interval is small enough, we can call it a day
		# this fixes cases in which Newton gets closer and closer to 0 with little practical improvement
		if diam(x) < ϵ || isempty(x)
			return x
		end

		# if Newton worked and improved the interval, continue iterating
		if !(x_old == x)
			continue
		end

		# Attempt 2: use the computed `fx_mid` to exclude half of the interval
		# in a bisection step

		if !f_endpoints_computed
			f_lo = f(Interval(x.lo))
			f_hi = f(Interval(x.hi))
			f_endpoints_computed = true
		end

		# here we rely on monotonicity
		if 0 ∉ hull(f_lo, f_mid)
			x = Interval(x_mid, x.hi)
			continue
		end
		if 0 ∉ hull(f_mid, f_hi)
			x = Interval(x.lo, x_mid)
			continue
		end

		# worst case: both hulls contain 0 (for instance, because 0 ∈ f_mid).
		# This looks like good news at first, but in practice it may just mean
		# that the computed enclosures are very large. x_mid may not be close to the solution at all.

		# Attempt 3: we compute the midpoints of (x.lo, x_mid) and (x_mid, x.hi)
		# to try to exclude the leftmost or rigthmost quarter of the interval
		#
		# Picture the number line like this:
		# x.lo .... a .... x_mid .... b .... x.hi

		a = (x.lo + x_mid) / 2
		if 0 ∉ hull(f_lo, f(Interval(a)))
			x = Interval(a, x.hi)
			continue
		end
		b = (x_mid + x.hi) / 2
		if 0 ∉ hull(f(Interval(b)), f_hi)
			x = Interval(x.lo, b)
			continue
		end

		# if all these attempts failed, then we cannot exclude that the solution belongs to [x.lo,a] nor [b,x.hi]
		# we can only return the current interval, which we cannot shrink anyway by more than a factor 2 at this point.
		# This will likely happen only when f() cannot be computed with sufficient precision.
		return x

	end
	@info "Maximum iterates reached in Newton contractor" max_iter, x, f(x), diam(x)
	return x
end


preimage(y, f, X; ϵ,  max_iter) = preimage(y, f, derivative(f), X; ϵ, max_iter)
preimage(y, f, fprime, X; ϵ, max_iter) = root(x -> f(x)-y, fprime, X; ϵ, max_iter)


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

function range_estimate_der(f, fprime, domain, recstep = 5)
	if recstep == 1
		m = typeof(domain)(mid(domain))
		return f(m)+fprime(domain)*radius(domain)
	else
		a, b = bisect(domain)
		Iₐ = range_estimate_der(f, fprime, a, recstep-1)
		Iᵦ = range_estimate_der(f, fprime, b, recstep-1)
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

# this is now unused
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
