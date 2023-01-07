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
Compute a single root with (possibly multivariate) interval Newton

x must be an Interval (univariate) or IntervalBox (multivariate)

Stops when the interval reaches a fixed point, when the diameter is smaller than ε,
or when max_iter iterations are reached (with an error)
"""
root(f, x; ϵ, max_iter) = root(f, derivative(f), x; ϵ, max_iter)

function root(f, f′, x; ϵ, max_iter)
	for i in 1:max_iter
		
		x_old = x

		x_mid = Interval(mid(x))
		@debug "Step $i of the Newton method:
			   - the interval $x, 
			   - the derivative $(f′(x)),
			   - the value at the midpoint $(f(x_mid))"
		
		Nx = x_mid - f′(x) \ f(x_mid)
		@debug "Candidate", Nx 

		x = intersect(x, Nx)

		if (x_old == x && diam(x) < ϵ) || isempty(x) 
			return x
		end

		if x_old == x
			# we assume our function is monotone 
			# on x and may only be not monotone 
			# on the interval representation of an 
			# endpoint if it is not a representable number
			
			# Isaia: Sketch of proof
			# Suppose the Newton method did not contract,
			# we do a bisection step; we need to estimate the 
			# range on each one of the two bisection intervals
			# since the map is guaranteed monotone by hypothesis
			# with the exception of endpoints coming from representation,
			# its range is contained in the hull of the enlarged endpoints
			# by interval arithmetic inclusion principle.
			
			x_l, x_r = bisect(x)

			y_l = range_estimate_monotone(f, x_l)
			y_r = range_estimate_monotone(f, x_r)
			
			# does a bisection step

			if !(0 ∈ (y_l)) 
				x = x_r
			elseif !(0 ∈ (y_r))
				x = x_l
			else
				# we need to treat the case in which both 
				# range estimates contain $0$; since the range estimate 
				# is obtained by evaluating on the enlarged endpoints
				# and the function is monotone this means that the zero
				# is contained in the enlarged common endpoint
				x = Interval(prevfloat(x_l.hi), nextfloat(x_r.lo))
			end
			x = intersect(x, x_old)
		end
	end
	@debug "Maximum iterates reached" max_iter, x, f(x), diam(x)
	return x
end


preimage(y, f, X; ϵ,  max_iter) = preimage(y, f, derivative(f), X; ϵ, max_iter)
preimage(y, f, fprime, X; ϵ, max_iter) = root(x -> f(x)-y, fprime, X; ϵ, max_iter)

function preimage(y::Interval, f, fprime, X; ϵ, max_iter)
	# I don't really like this, it is checking again if the function 
	# is increasing or decreasing; probably the best would be 
	# to define a preimage method which dispatches on the branch type
	# now, the question is if it should be in Contractors.jl 
	# or PwDynamicDefinition.jl
	if !(0 ∈ fprime(X))
		x_lo = preimage(y.lo, f, fprime, X; ϵ, max_iter)
		x_hi = preimage(y.hi, f, fprime, X; ϵ, max_iter)
		return hull(x_lo, x_hi)
	else
		@error "Preimage of a wide interval through a non monotone function"
	end
end




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

function range_estimate_monotone(f, x)
	low = Interval(prevfloat(x.lo), nextfloat(x.lo))
	high = Interval(prevfloat(x.hi), nextfloat(x.hi))
	return hull(f(low),f(high))
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
