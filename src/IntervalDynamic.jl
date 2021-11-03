struct IntervalDynamic <: Dynamic
    branches::Array{Branch, 1}
	full_branch::Bool
end

Base.getindex(D::IntervalDynamic, k::Int64) = D.branches[k]

IntervalDynamic(Ts, endpoints::Vector{Interval{T}}; full_branch = false) where {T} = IntervalDynamic(Ts, endpoints, hcat([Ts[k](Interval(endpoints[k])) for k in 1:length(Ts)], [Ts[k](Interval(endpoints[k+1])) for k in 1:length(Ts)]); full_branch = full_branch)

function IntervalDynamic(Ts, endpoints, y_endpoints_in; full_branch = false)
    branches = Branch[]
    for k in 1:length(endpoints)-1
        y_endpoints = (y_endpoints_in[k,1], y_endpoints_in[k,2])
        increasing  = unique_increasing(y_endpoints_in[k,1], y_endpoints_in[k,2])
        push!(branches, Branch(Ts[k], (endpoints[k], endpoints[k+1]), y_endpoints, increasing))
    end
    return IntervalDynamic(branches, full_branch)
end

branches(D::IntervalDynamic) = D.branches

DynamicDefinition.nbranches(D::IntervalDynamic) = length(D.branches)
DynamicDefinition.endpoints(D::IntervalDynamic) = [[br.X[1] for br in branches(D)]; branches(D)[end].X[2]]

DynamicDefinition.orientation(D::IntervalDynamic, k) = D.branches[k].increasing ? 1. : -1.
domain(D::IntervalDynamic) = (D.branches[1].X[1], D.branches[end].X[2])

#DynamicDefinition.is_full_branch(D::IntervalDynamic) = all(r == [0.,1.] || r == [1.,0.] for r in eachrow(D.y_endpoints))  # TODO: this assumes domain == [0,1] unnecessarily

"""
Intersect an Interval or TaylorSeries with I
"""
restrict(I, x) = I ∩ x
restrict(I, x::Taylor1) = Taylor1([I ∩ x[0]; x[1:end]], x.order)

"""
function that evaluates the k-th branch of a dynamic on a point x
	(assuming it's in its domain, otherwise ∅)
"""
function DynamicDefinition.branch(D::IntervalDynamic, k)
	return x -> D[k].f(restrict(hull(D[k].X[1], D[k].X[2]), x))
end

### This needs to be better implemented...
BasisDefinition.is_full_branch(D::IntervalDynamic) =  D.full_branch

# Unused as of now
# """
# hull of an iterable of intervals
# """
# common_hull(S) = interval(minimum(x.lo for x in S), maximum(x.hi for x in S))


# Rather than defining derivatives of an IntervalDynamic, we define Taylor1 expansions directly
# and let the generic functions in DynamicDefinition to the work
"""
Function call, and Taylor expansion, of a IntervalDynamic.
Note that this ignores discontinuities; users are free to shoot themselves
in the foot and call this on a non-smooth piecewise map. No better solutions for now.
"""
function (D::IntervalDynamic)(x::Taylor1)
	fx = fill(∅, x.order+1)
	x_restricted = deepcopy(x)
	for i = 1:length(D.branches)
		x_restricted[0] = x[0] ∩ hull(D[i].X[1],D[i].X[2])
		if !isempty(x_restricted[0])
			fx_restricted = D[i].f(x_restricted)
			fx = fx .∪ fx_restricted.coeffs
		end
	end
	@debug "Piecewise f($(x)) = $(Taylor1(fx, x.order))"
	return Taylor1(fx, x.order)
end

function DynamicDefinition.plottable(D::IntervalDynamic, x)
	@assert 0 <= x <= 1
	for k in 1:nbranches(D)
		domain = hull(D[k].X[1], D[k].X[2])
		if x in domain
			return D[k].f(x)
		end
	end
end

using RecipesBase
@recipe f(::Type{PM}, D::PM) where {PM <: IntervalDynamic} = x -> plottable(D, x)
