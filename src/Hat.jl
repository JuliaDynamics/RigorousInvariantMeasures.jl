"""
Equispaced Ulam basis on [0,1] of size n 
"""
struct Hat <:Basis
	n::Integer #TODO, change to partition
end

"""
Returns the size of the Hat basis
"""
Base.length(b::Hat) = b.n


Base.iterate(b::Hat, state = 1) = state < length(b)+1 ? (b[state-1], state+1) : nothing

"""
Returns the point on which the i-th element of the Hat basis has value 1
"""
Base.getindex(b::Hat, i) = Float64(i)/b.n

"""
This iterator returns the preimages of the endpoints 
of the intervals defining the Ulam basis through the dynamic
"""
function Base.iterate(S::DualComposedWithDynamic{Hat, Mod1Dynamic}, state = (1, 1))
	i, k = state
	
	if i == length(S.basis)+1
			return nothing
	end

	# remark that this version supposes that for each i there exists a preimage
	# another more specific version should be implemented for maps with 
	# incomplete branches


	x = preim(S.dynamic, k, getindex(S.basis, i-1), S.ϵ) 
	der = der(S.dynamic, x)

	if k == nbranches(S.dynamic)
		return ((i, (x, der)), (i+1, 1))
	else
		return ((i, (x, der)), (i, k+1))
	end
end

"""
Returns the indices of the elements of the Ulam basis that intersect with the interval y
"""
nonzero_on(B::Hat, y) = max(floor(Int64, y[1].lo*length(B)), 1), min(ceil(Int64, y[2].hi*length(B)), length(B))

"""
Constructor for the ProjectDualElement iterator
"""
function ProjectDualElement(B::Hat, y)
	j_min, j_max = nonzero_on(B, y)
	return ProjectDualElement(B, j_min, j_max, y)
end


ϕ(B::Hat, j, y) = @error Not Implemented


"""
Given a preimage ```y``` of a point ```x```, this iterator returns
```\\phi_j(y)/T'(y) ```
"""
function Base.iterate(S::ProjectDualElement{Hat}, state = S.j_min)
	if state == S.j_max+1
		return nothing
	end
	y, der = S.dual_element

	return ((state, ϕ(S.basis, state, y)/der), 
		    state+1) 
end

evaluate_integral(B::Hat, i, T) = T(i)/length(B)

function Base.iterate(S::AverageZero{Hat}, state = 1) 
	n = length(S.basis)
	if state == n
		return nothing
	end
	v = zeros(Float64, n)
	v[1] = 1
	v[state+1]=-1
	return (v, state+1)
end

using RecipesBase

@userplot PlotHat
@recipe function f(h::PlotHat)
	if length(h.args)!= 2 || (typeof(h.args[1])!= Ulam) || !(typeof(h.args[2])<:AbstractVector)
		error("Plot Ulam needs as an input a Ulam Basis and a vector")
	end

	B = h.args[1]
	w = h.args[2]

	seriestype := :path
	collect(B), mid.(w)
end
