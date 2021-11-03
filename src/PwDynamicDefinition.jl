
## the function branches is after the module


module PwDynamicDefinition
using ValidatedNumerics
using ..DynamicDefinition
using ..Contractors
using TaylorSeries: Taylor1

using ..DynamicDefinition: derivative, orientation

export PwMap, preim, nbranches, plottable, branches, Branch

"""
Type used to represent a "branch" of a dynamic. The branch is represented by a monotonic map `f` with domain `X=(a,b)` with a≤b (where typically a,b are intervals). 
`Y=(f(a),f(b))` and `increasing` may be provided (for instance if we know that `Y=(0,1)`), otherwise they are computed automatically.
"""
struct Branch{T,S}
    f::T
    X::Tuple{S, S}
    Y::Tuple{S, S}
    increasing::Bool
end
Branch(f, X, Y=(f(Interval(X[1])), f(Interval(X[2]))), increasing=unique_increasing(Y[1], Y[2])) = Branch{typeof(f), typeof(interval(X[1]))}(f, X, Y, increasing)

"""
Dynamic based on a piecewise monotonic map.

The map is defined as T(x) = Ts[k](x) if x ∈ [endpoints[k], endpoints[k+1]).

`y_endpoints` (kx2 matrix) contains the result of applying Ts to the endpoints of each interval. These can be filled in automatically from `endpoints`,
but sometimes they are known to higher accuracy, for instance for `x -> mod(3x, 1)` we know that it is full-branch exactly.
"""
struct PwMap <: Dynamic
	branches::Array{Branch, 1}
	full_branch::Bool
end

function PwMap(Ts, endpoints, y_endpoints_in; full_branch = false)
    branches = Branch[]
    for k in 1:length(endpoints)-1
        y_endpoints = (y_endpoints_in[k,1], y_endpoints_in[k,2])
        increasing  = unique_increasing(y_endpoints_in[k,1], y_endpoints_in[k,2])
        push!(branches, Branch(Ts[k], (endpoints[k], endpoints[k+1]), y_endpoints, increasing))
    end
    return PwMap(branches, full_branch)
end

function PwMap(Ts, endpoints::Vector{T}; full_branch = false) where {T<:Real}  
	return PwMap(Ts, endpoints, hcat([Ts[k](Interval(endpoints[k])) for k in 1:length(Ts)], [Ts[k](Interval(endpoints[k+1])) for k in 1:length(Ts)]); full_branch = full_branch)
end

Base.show(io::IO, D::PwMap) = print(io, "Piecewise-defined dynamic with $(nbranches(D)) branches")

DynamicDefinition.domain(D::PwMap) =  (D.branches[1].X[1], D.branches[end].X[2])

Base.getindex(D::PwMap, k::Int64) = D.branches[k]

DynamicDefinition.nbranches(D::PwMap) = length(D.branches)
DynamicDefinition.endpoints(D::PwMap) = [[br.X[1] for br in branches(D)]; branches(D)[end].X[2]]
DynamicDefinition.branches(D::PwMap) = D.branches

DynamicDefinition.orientation(D::PwMap, k) = D[k].increasing ? 1. : -1.

DynamicDefinition.is_full_branch(D::PwMap) = D.full_branch

function DynamicDefinition.preim(D::PwMap, k, y, ϵ = 1e-15)
	@assert 1 <= k <= nbranches(D)
	domain = hull(D[k].X[1], D[k].X[2])
	return preimage(y, D[k].f, domain, ϵ)
end

"""
Intersect an Interval or TaylorSeries with I
"""
restrict(I, x) = I ∩ x
restrict(I, x::Taylor1) = Taylor1([I ∩ x[0]; x[1:end]], x.order)

"""
function that evaluates the k-th branch of a dynamic on a point x
	(assuming it's in its domain, otherwise ∅)
"""
function DynamicDefinition.branch(D::PwMap, k)
	return x -> D[k].f(restrict(hull(D[k].X[1], D[k].X[2]), x))
end

# Unused as of now
# """
# hull of an iterable of intervals
# """
# common_hull(S) = interval(minimum(x.lo for x in S), maximum(x.hi for x in S))


# Rather than defining derivatives of a PwMap, we define Taylor1 expansions directly
# and let the generic functions in DynamicDefinition to the work
"""
Function call, and Taylor expansion, of a PwMap.
Note that this ignores discontinuities; users are free to shoot themselves
in the foot and call this on a non-smooth piecewise map. No better solutions for now.
"""
function (D::PwMap)(x::Taylor1)
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


function DynamicDefinition.plottable(D::PwMap, x)
	@assert 0 <= x <= 1
	for k in 1:nbranches(D)
		domain = hull(D[k].X[1], D[k].X[2])
		if x in domain
			return D[k].f(x)
		end
	end
end

using RecipesBase
@recipe f(::Type{PM}, D::PM) where {PM <: PwMap} = x -> plottable(D, x)

end

"""
Utility constructor for dynamics Mod 1.
We assume that f is monotonic and differentiable, for now (this is not restrictive, for our purposes)
"""
function mod1_dynamic(f::Function, X = (0.,1.), ε = 0.0; full_branch = false)
    br = Branch(f, X)

    # check monotonicity
    fprime = x -> derivative(f, x)
    @assert minimise(x -> fprime(x) * (br.increasing ? 1 : -1), hull(Interval.(X)...))[1] > 0

    Yhull = hull(Interval.(br.Y)...)
    possible_integer_parts = floor(Int, Yhull.lo):ceil(Int, Yhull.hi)

    x, integer_parts = preimages(possible_integer_parts, br, possible_integer_parts)

    ep = [x; X[end]]
    Ts = [x->f(x)-k for k in integer_parts]

    n = Base.length(x)
    if br.increasing
        y_endpoints::Matrix{Interval{Float64}} = hcat(fill(0., n), fill(1., n))
    else
        y_endpoints = hcat(fill(1., n), fill(0., n))
    end
    y_endpoints[1, 1] = br.Y[begin] - integer_parts[begin]
    if y_endpoints[1, 1] == 0.
        y_endpoints[1, 1] = 0. # hack to get rid of -0..0 intervals
    end
    y_endpoints[end, end] = br.Y[end] - integer_parts[end]
    if y_endpoints[end, end] == 0.
        y_endpoints[end, end] = 0. # hack to get rid of -0..0 intervals
    end


    return PwMap(Ts, ep, y_endpoints; full_branch = full_branch)
end

