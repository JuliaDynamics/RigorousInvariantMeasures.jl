module Contractors
using ValidatedNumerics
using DualNumbers


export N_rig, root, range_estimate, ShootingMethod

"""
Step of rigorous Newton (possibly multivariate)
"""
function N_rig(f, f′, x) where {T}
	x_mid = Interval.(mid.(x))
	return intersect(x, x_mid - f′(x) \ f(x_mid))
end

"""
Compute a single root with (possibly multivariate) interval Newton

x must be an Interval (univariate) or IntervalBox (multivariate)
"""
root(f, x::Interval, ϵ; max_iter = 100) = root(f, x->f(Dual(x, 1)).epsilon, x, ϵ; max_iter = max_iter)

function root(f, f′, x, ϵ; max_iter = 100)
	@debug x
	for i in 1:max_iter
		x = N_rig(f, f′, x)
		@debug diam(x)
		if diam(x)<ϵ
			return x
		end
		if isempty(x)
			return x
		end
	end
	@info "Maximum iterates reached" max_iter
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
		return union(Iₐ, Iᵦ)
	end
end

# a shooting method to find the k-th preimage of a point

#function that gives me the type of the endpoints of the interval


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

end
