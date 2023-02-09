module Contractors
using IntervalArithmetic


export preimage_monotonic, range_estimate, ShootingMethod, preimage, unique_increasing

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
Compute the guaranteed interval preimage f⁻¹(y) of a point y ∈ R according to the function f: [a,b] → R. 
This preimage is guaranteed to lie inside the real interval `x = hull(a,b)`.

`f` must be strictly monotonic on [a,b], at least in the sense of MonotonicBranch functions (allowing for uncertainty on the endpoints).
`f` must be differentiable; zero and infinite derivative are allowed.

`y1`, `y2` must be guaranteed to be equal to f(a), f(b) or "more outside", so that  y ∉ hull(y1,f(z)) ⟹ f⁻¹(y) ∉ hull(a, z).

Stops when the interval reaches a fixed point, when the diameter is smaller than ε,
or when max_iter iterations are reached (with an error)
"""
preimage_monotonic(y, f::Function, x::Interval, (y1, y2); ϵ, max_iter) = preimage_monotonic(y, f::Function, derivative(f), x::Interval, (y1, y2); ϵ, max_iter)
function preimage_monotonic(y, f::Function, f′, x::Interval, (y1, y2) = (f(Interval(x.lo)), f(Interval(x.hi))); ϵ, max_iter)
	x_old::typeof(x) = ∅
	for i in 1:max_iter

		if diam(x) < ϵ || isempty(x) || x == x_old
			return x
		end

		x_old = x
		x_mid = Interval(mid(x))

		fm::typeof(x) = ∅  # both Newton and Krawczyk will compute f(x_mid) as a byproduct; we save it here.
		enc_der = f′(x)
		if 0 ∉ enc_der
			# Interval Newton step on x -> f(x) - y
			@debug "Step $i, Newton method:
			- the interval $x, 
			- the derivative $(f′(x)),
			- the value at the midpoint $(f(x_mid))
			- the value $(f(x))"

			fm = f(x_mid)
			Nx = x_mid - enc_der \ (fm - y)  # we use the "backward division" operator for good practice, just in case this will ever applied in a multidimensional setting
			@debug "Newton Candidate", Nx
			
		else
			# Krawczyk step on x -> f(x) - y
			@debug "Step $i, Krawczyk method:
			- the interval $x, 
			- the derivative $(f′(x)),
			- the value at the midpoint $(f(x_mid))
			- the value $(f(x))"

			fm = f(x_mid) # TODO: merge this line and the next into a single call
			der_mid = f′(x_mid)

			Nx = x_mid - der_mid \ (fm - y) + (1 - der_mid \ enc_der)*(x-x_mid)

			@debug "Krawczyk Candidate", Nx

		end
		x = intersect(x, Nx)
		if !(x_old == x)  # Newton / Krawczyk was effective and improved the interval
			continue
		end
		# otherwise, let us switch to a zero-th order method based on bisection
		# fm has already been computed, let us use it

		# Notice that this bisection strategy works also when f is not monotonic outside [a,b]; 
		# see https://github.com/JuliaDynamics/RigorousInvariantMeasures.jl/issues/140#issuecomment-1413541452 for a proof
		if y ∉ hull(y1, fm)
			x = Interval(x_mid.hi, x.hi)
			continue
		end
		if y ∉ hull(fm, y2)
			x = Interval(x.lo, x_mid.lo)
			continue
		end

		# this may still fail if y ≈ fm or if the function is computed in a very inaccurate way (relative to the size of `x`)
		# in this case, we try chopping off smaller and smaller parts of `x` at both ends
		# If we cannot chop off even 1/2^6th of `x`, at the next outer loop iteration `x == x_old` and we will return.
		for k = 2:6
			c = x.lo + diam(x) / 2^k
			if y ∉ hull(y1, f(Interval(c)))
				x = Interval(c, x.hi)
				break
			end
			c = x.hi - diam(x) / 2^k
			if y ∉ hull(f(Interval(c)), y2)
				x = Interval(x.lo, c)
				break
			end
		end
	end
	@warn "Maximum iterates reached:" max_iter, x, f(x), diam(x)
	@warn "This should not happen normally, consider debugging Contractors.preimage_monotonic."
	return x
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

end #module
