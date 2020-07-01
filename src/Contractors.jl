module Contractors
using ValidatedNumerics
using DualNumbers


export N_rig, root, range_estimate

function N_rig(f, f′, x::Interval{T}) where {T}
	x_mid = Interval{T}(mid(x))
	return intersect(x, x_mid-f(x_mid)/f′(x))
end

root(f, x::Interval, ϵ; max_iter = 100) = root(f, x->f(Dual(x, 1)).epsilon, x, ϵ; max_iter = max_iter)

function root(f, f′, x::Interval, ϵ; max_iter = 100)
	for i in 1:max_iter
		x = N_rig(f, f′, x)
		if diam(x)<ϵ
			return x
		end
	end
	@info "Maximum iterates reached" max_iter
	return x
end

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

end
