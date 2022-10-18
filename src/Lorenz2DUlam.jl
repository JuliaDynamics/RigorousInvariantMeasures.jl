module Lorenz2D

const α = 1.75::Float64
const s = 3.0::Float64

"""
    ContractingLorenz1D(; α , s)

Construct a dynamic representing the one-dimensional map for the contracting Lorenz Flow,
with parameters `α` and `s`, i. e.

``T(x) = -α*(-x)^s+1`` if ``-1<x<0`` and ``T(x) = α*(x)^s-1`` for ``0<x<1``
"""

import RigorousInvariantMeasures: PwMap, relative_measure, domain

ContractingLorenz1D(; α , s) = PwMap([x-> -α*(-x)^s+1, x-> α*(x)^s-1], [-1, 0, 1])

# problem, this is  a map on [-1, 1]
T = ContractingLorenz1D(α = α , s = s)

# these are the coordinate change maps
ϕ = PwMap([x->2*x-1, ], [0, 1])
ψ = PwMap([x-> x/2+1/2, ], [-1, 1])

# this is the map on [0, 1], incredibly it works!!!
D = ψ∘T∘ϕ


# to simplify implementation, we use the Julia Polyhedra library
import Polyhedra as PH
import GLPK
lib = PH.DefaultLibrary{Float64}(GLPK.Optimizer)

import IntervalArithmetic
import IntervalArithmetic: Interval, mid

import RigorousInvariantMeasures: preimages

Ginverse(; x, r, c) = x > 0 ?
                v -> ((v-c)*2^r)/(x^r) :
                v -> ((v+c)*2^r)/((-x)^r)


"""
    PreimageRectangleLorenz(; preim_x_left, preim_x_right, y_lower, y_upper, k)

Return a Polyhedra that approximates the preimage of the rectangle 
``[T(preimx_{left}), T(preimx_{right})] × [y_{lower}, y_{upper}]``
through the skew product map ``F(x, y) = (T(x), G(x, y)) and the maximum of the 
error in the computation of the vertices.
The argument `k` is the number of equispaced linear interpolation points in the
`x` direction.
"""
function PreimageRectangleLorenz(;  preim_x_left::Interval, 
                                    preim_x_right::Interval,
                                    y_lower, 
                                    y_upper, 
                                    k,
                                    r,
                                    c)
    preim_x = [preim_x_left+(preim_x_right-preim_x_left)/k*i for i in 0:k]

    err_x = maximum(IntervalArithmetic.radius.(preim_x))
    err_y = 0.0

    y_s = NTuple{2, Interval}[]

    for x in preim_x
        G_inv = Ginverse(x = x, r = r, c = c)
        preim_y_low = min(max(G_inv(y_lower), -1), 1)
        preim_y_up = max(min(G_inv(y_upper), 1), -1)
        err_y = maximum([err_y, IntervalArithmetic.radius(preim_y_low), IntervalArithmetic.radius(preim_y_up)])
        push!(y_s, (preim_y_low, preim_y_up))
    end
    n = length(preim_x)
    A = Matrix{Float64}(undef, 2*n, 2)
    for (i, x) in enumerate(preim_x)
        A[i, :] = [mid(x) mid(y_s[i][1])]
    end
    
    for (i, x) in enumerate(reverse(preim_x))
        A[n+i, :] = [mid(x) mid(y_s[end-i+1][2])]
    end
    
    P = PH.polyhedron(PH.vrep(A), lib)
    return P, max(err_x, err_y)
end

using IntervalArithmetic
function _Lorenz_one_dim_map(x::Interval, α, s)
    x_left = x ∩ @interval -1 0
    x_right = x ∩ @interval 0 1
    return -α*(-x_left)^s+1 ∪ α*(x_right)^s-1
end

_Lorenz_left_one_dim_map(x::Interval, α, s) = -α*(-(x ∩ @interval -1 0))^s+Interval(1)
_Lorenz_right_one_dim_map(x::Interval, α, s) = α*(x ∩ @interval 0 1)^s-Interval(1)
_Lorenz_left_fiber_map(x::Interval, y::Interval, r, c) = 2^(-r)*y*(-x)^r-c
_Lorenz_right_fiber_map(x::Interval, y::Interval, r, c) = 2^(-r)*y*x^r+c


function BoundingRandomAttractor(α, s, r, c, ξ, l, k)
    Dict1to2 = Dict{Int64, NTuple{2, Int64}}()
    Dict2to1 = Dict{NTuple{2, Int64}, Int64}()

    left_image_rectangle_x = _Lorenz_left_one_dim_map(Interval(-1,0), α, s)
    left_image_rectangle_y = _Lorenz_left_fiber_map(Interval(-1), Interval(-1, 1), r, c)


    right_image_rectangle_x = _Lorenz_right_one_dim_map(Interval(0.5, 1), α, s)
    right_image_rectangle_y = _Lorenz_right_fiber_map(Interval(0.5), Interval(-1, 1), r, c)
    

    Numberelements = 0

    @info "left"
    for i in 0:2^l-1

        x = Interval(-1)+Interval(i, i+1)/2^l
        left_image_rectangle_x = _Lorenz_left_one_dim_map(x, α, s)
        @info left_image_rectangle_x
        left_image_rectangle_x += Interval(-ξ, ξ)
        @info left_image_rectangle_x
        left_image_rectangle_y = _Lorenz_left_fiber_map(x, Interval(-1, 1), r, c)+Interval(-ξ, ξ)
        
        @debug left_image_rectangle_x
        @debug left_image_rectangle_y

        left_image_integer_x = k*(left_image_rectangle_x/2+1/2)
        left_image_integer_y = k*(left_image_rectangle_y/2+1/2)

        lo_x = floor(Int64, left_image_integer_x.lo)
        hi_x = ceil(Int64, left_image_integer_x.hi)

        lo_y = floor(Int64, left_image_integer_y.lo)
        hi_y = ceil(Int64, left_image_integer_y.hi)
        @debug "x ", lo_x, hi_x,  " y ", lo_y, hi_y

        for i in lo_x:min(hi_x, k)
            for j in lo_y:hi_y
                @debug "i, j" i, j
                if !haskey(Dict2to1, (i, j))
                    Numberelements += 1
                    @debug "N" Numberelements
                    Dict2to1[ (i, j) ] = Numberelements
                    Dict1to2[Numberelements] = (i, j)  
                end
            end
        end

        if hi_x > k
            break
        end        
    end

    @info "right"
        
    for i in 0:2^l-1
        x = Interval(i, i+1)/2^l
        right_image_rectangle_x = _Lorenz_left_one_dim_map(x, α, s)+Interval(-ξ, ξ)
        right_image_rectangle_y = _Lorenz_left_fiber_map(x, Interval(-1, 1), r, c)+Interval(-ξ, ξ)
        
        @debug right_image_rectangle_x
        @debug right_image_rectangle_y

        right_image_integer_x = k*(right_image_rectangle_x/2+1/2)
        right_image_integer_y = k*(right_image_rectangle_y/2+1/2)

        lo_x = floor(Int64, right_image_integer_x.lo)
        hi_x = ceil(Int64, right_image_integer_x.hi)

        lo_y = floor(Int64, right_image_integer_y.lo)
        hi_y = ceil(Int64, right_image_integer_y.hi)
        @debug "x ", lo_x, hi_x,  " y ", lo_y, hi_y

        for i in max(lo_x,0):hi_x
            for j in lo_y:hi_y
                #@info "i, j" i, j
                if !haskey(Dict2to1, (i, j))
                    Numberelements += 1
                    #@info "N" Numberelements
                    Dict2to1[ (i, j) ] = Numberelements
                    Dict1to2[Numberelements] = (i, j)  
                end
            end
        end

        if hi_x > k
            break
        end        
    end

    return Dict1to2, Dict2to1
end

struct LorenzOperator
    L
    d21
    d12
end


using SparseArrays
function compute_transfer_operator(α, s, r, c, k_x, k_y, interp_step)
    T = ContractingLorenz1D(α = α , s = s)
    ϕ = PwMap([x->2*x-1, ], [0, 1])
    ψ = PwMap([x-> x/2+1/2, ], [-1, 1])
    D = ψ∘T∘ϕ
    
    I = Int64[]
	J = Int64[]
	nzvals = Interval{Float64}[]
	n = k_x*k_y

    preim_x = preimages(LinRange(0., 1., k_x+1), D)
    preim_x_coord_change = (2 .*preim_x[1].-1, preim_x[2]) 
    #@info preim_x_coord_change

    part_x = LinRange(-1., 1., k_x+1)
    part_y = LinRange(-1., 1., k_y+1)

    for ind_x_preim in 1:length(preim_x_coord_change[1][1:end]) 
        a = preim_x_coord_change[1][ind_x_preim]
        if ind_x_preim == length(preim_x_coord_change[1][1:end])
            b = domain(D)[end] 
        else
            b = preim_x_coord_change[1][ind_x_preim+1]
        end
        #TODO add what happens near 1
    
        y = hull(a, b)
        #@info "preimage" y
        image_ind_x = preim_x_coord_change[2][ind_x_preim]
        

        preim_a = searchsortedlast(part_x, y.lo)
        preim_b = min(searchsortedlast(part_x, y.hi), k_x)

        for ind_x in preim_a:preim_b
            rel_meas = relative_measure((a, b), (Interval(part_x[ind_x]), Interval(part_x[ind_x+1])))
            #@info ind_x, part_x[ind_x], part_x[ind_x+1]
            #@info a, b, rel_meas
            dom_ind_x = ind_x

            for ind_y in 1:k_y
                # tot_meas = 0.0 
                dom_ind_y = ind_y

                dom_ind_reshaped = (dom_ind_x-1)*k_y+dom_ind_y

                if ind_x <= k_x/2
                    im_y = _Lorenz_left_fiber_map(y, Interval(part_y[ind_y], part_y[ind_y+1]), r, c)
                    ind_y_lower = searchsortedlast(part_y, im_y.lo)
                    ind_y_upper = searchsortedlast(part_y, im_y.hi)
                    indexes_y = ind_y_lower:ind_y_upper
                else
                    im_y = _Lorenz_right_fiber_map(y, Interval(part_y[ind_y], part_y[ind_y+1]), r, c)
                    ind_y_lower = searchsortedlast(part_x, im_y.lo)
                    ind_y_upper = searchsortedlast(part_x, im_y.hi)
                    indexes_y = ind_y_lower:ind_y_upper
                end

                if length(indexes_y)==1
                    @debug "Only 1"
                        
                    image_ind_y = ind_y_lower

                    image_ind_reshaped = (image_ind_x-1)*k_y+image_ind_y

                    push!(I, dom_ind_reshaped)
				    push!(J, image_ind_reshaped)
				    push!(nzvals, rel_meas)
                    @debug ind_y, im_y, ind_y_lower, ind_y_upper, part_y[ind_y_lower], part_y[ind_y_lower+1]
                else
                    @debug "More than 1"
                    for j in indexes_y
                        image_ind_y = j
                        image_ind_reshaped = (image_ind_x-1)*k_y+image_ind_y

                        P, ϵ = PreimageRectangleLorenz(;  preim_x_left = a , 
                                    preim_x_right = b, 
                                    y_lower = part_y[j], 
                                    y_upper = part_y[j+1], 
                                    k = interp_step,
                                    r = r,
                                    c = c)

                        domain_rectangle =  PH.polyhedron(PH.vrep([
                                            part_x[ind_x] part_y[ind_y]
                                            part_x[ind_x+1] part_y[ind_y]
                                            part_x[ind_x+1] part_y[ind_y+1]
                                            part_x[ind_x] part_y[ind_y+1]
                                                ]), lib)
                         

                        intersection = PH.intersect(P, domain_rectangle)
                        @info PH.vrep(intersection)
                        rel_meas_poly = (PH.volume(intersection)*k_x*k_y)/4
                        # TODO: add error for the computed value!
                        
                        push!(I, dom_ind_reshaped)
				        push!(J, image_ind_reshaped)
				        push!(nzvals, rel_meas_poly)

                        @debug ind_y, im_y, ind_y_lower, ind_y_upper, part_y[j], part_y[j+1]
                    end
                    # @debug "tot" tot_meas
                    # @debug "rel_meas" rel_meas 
                    # if abs(rel_meas-tot_meas)> 0.1
                    #    @info "ind_x" ind_x 
                    #    @info "ind_y" ind_y 
                    #    @info "difference" rel_meas-tot_meas
                    # else
                    #    @info "Good"
                    # end
                end
            end
        end
    end
    return sparse(I, J, nzvals, n, n)
end

end