
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
    infinite_derivative::Bool
    function PwMap(branches::Array{Branch, 1}; full_branch = false, infinite_derivative = false)
        new(branches, full_branch, infinite_derivative)
    end
end

function PwMap(Ts, endpoints, y_endpoints_in; full_branch = false, infinite_derivative = false)
    branches = Branch[]
    for k in 1:length(endpoints)-1
        y_endpoints = (y_endpoints_in[k,1], y_endpoints_in[k,2])
        increasing  = unique_increasing(y_endpoints_in[k,1], y_endpoints_in[k,2])
        push!(branches, Branch(Ts[k], (endpoints[k], endpoints[k+1]), y_endpoints, increasing))
    end
    return PwMap(branches; full_branch = full_branch, infinite_derivative = infinite_derivative)
end

function PwMap(Ts, endpoints::Vector{T}; full_branch = false, infinite_derivative = false) where {T<:Real}  
	return PwMap(Ts, endpoints, hcat([Ts[k](Interval(endpoints[k])) for k in 1:length(Ts)], [Ts[k](Interval(endpoints[k+1])) for k in 1:length(Ts)]); full_branch = full_branch, infinite_derivative = infinite_derivative)
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


conv_orientation(x::Bool) = x ? 1 : -1
inv_conv_orientarion(x::Int64) = x > 0

function composedPwMap(D1::PwDynamicDefinition.PwMap, D2::PwDynamicDefinition.PwMap)
    new_branches =  Branch[]
    for br2 in branches(D2)
        if br2.increasing
            for br1 in branches(D1)
                y_range = hull(br1.Y[1], br1.Y[2])
                left  = preimage(br1.X[1], br2.f, hull(br2.X[1], br2.X[2]), 10^-13)
                right  = preimage(br1.X[2], br2.f, hull(br2.X[1], br2.X[2]), 10^-13)
                @info left
                @info right
                F = br1.f∘br2.f
                F_increasing = conv_orientation(br1.increasing)*conv_orientation(br2.increasing)
                if left!=∅ && right!=∅
                    y_endpoints = (F(left) ∩ y_range, F(right) ∩ y_range)
                    push!(new_branches, Branch(F, (left, right), y_endpoints, inv_conv_orientarion(F_increasing)))
                elseif left == ∅ && right != ∅
                    left = br2.X[1]
                    y_endpoints = (F(left) ∩ y_range, F(right) ∩ y_range)
                    push!(new_branches, Branch(F, (left, right), y_endpoints, inv_conv_orientarion(F_increasing)))
                elseif left !=∅ && right == ∅
                    right = br2.X[2]
                    y_endpoints = (F(left) ∩ y_range, F(right) ∩ y_range)
                    push!(new_branches, Branch(F, (left, right), y_endpoints, inv_conv_orientarion(F_increasing)))
                end
            end
        else
            for br1 in Iterators.reverse(branches(D1))
                y_range = hull(br1.Y[1], br1.Y[2])
                left  = preimage(br1.X[2], br2.f, hull(br2.X[1], br2.X[2]), 10^-13)
                right  = preimage(br1.X[1], br2.f, hull(br2.X[1], br2.X[2]), 10^-13)
                F = br1.f∘br2.f
                F_increasing = conv_orientation(br1.increasing)*conv_orientation(br2.increasing)
                if left!=∅ && right!=∅
                    y_endpoints = (F(left) ∩ y_range, F(right) ∩ y_range)
                    push!(new_branches, Branch(F, (left, right), y_endpoints, inv_conv_orientarion(F_increasing)))
                elseif left == ∅ && right != ∅
                    left = br2.X[1]
                    y_endpoints = (F(left) ∩ y_range, F(right) ∩ y_range)
                    push!(new_branches, Branch(F, (left, right), y_endpoints, inv_conv_orientarion(F_increasing)))
                elseif left !=∅ && right == ∅
                    right = br2.X[2]
                    y_endpoints = (F(left) ∩ y_range, F(right) ∩ y_range)
                    push!(new_branches, Branch(F, (left, right), y_endpoints, inv_conv_orientarion(F_increasing)))
                end
            end
        end
    end
    # these are conservative in order to be generic, i.e., due to composition
    # we could have a different behaviour, but the statement belows are conservatively true
    full_branch = D1.full_branch && D2.full_branch
    infinite_derivative = D1.infinite_derivative || D2.infinite_derivative
    return PwMap(new_branches; full_branch = full_branch, infinite_derivative = infinite_derivative)
end

function expansivity(D::PwDynamicDefinition.PwMap, tol=1e-3)
	max_exp = Interval(0.0)
    for br in branches(D)
        val = minimise(x -> abs(1/derivative(br.f, x)), hull(br.X[1], br.X[2]), tol=tol)[1]
        @info val
        max_exp = max(val, max_exp)
    end
    return max_exp
end

function max_distortion(D::PwDynamicDefinition.PwMap, tol=1e-3)
	max_dist = Interval(0.0)
    for br in branches(D)
        val = maximise(x -> abs(distorsion(br.f, x)), hull(br.X[1], br.X[2]), tol=tol)[1]
        max_dist = max(val, max_dist)
    end
    return max_dist
end

#function dfly_inf_der(D::PwDynamicDefinition.PwMap, tol=1e-3)
#    max_exp  = 2*expansivity(D, tol)
#    if max_exp > 1
#        @error "More expansivity needed, take an iterate"
#    end

#    aux_der = (1-max_exp)/2
#  @info "aux_der", aux_der
    
#    for br in branches(D)
#        if br.increasing
#            rr = preimage(aux_der, x -> 1/(2*derivative(br.f, x)), hull(br.X[1], br.X[2]), 10^-13)
#            @info rr
#        else
#            rr = preimage(-aux_der, x -> -1/(2*derivative(br.f, x)), hull(br.X[1], br.X[2]), 10^-13)
#            @info rr
#        end
#    end
#end

#dfly_inf_der(D::ComposedDynamic, tol=1e-3) = dfly_inf_der(D.E, tol)