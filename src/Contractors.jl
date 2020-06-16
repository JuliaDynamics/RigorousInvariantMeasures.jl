using ValidatedNumerics
using DualNumbers

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