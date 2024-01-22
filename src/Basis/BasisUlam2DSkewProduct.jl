using IntervalArithmetic, LinearAlgebra

export Ulam2DSP,
    length_x, length_y, square_indexes_to_linear, linear_indexes_to_square, getindex_linear

"""
	Ulam2DSP
Ulam basis on [0,1]×[0,1] associated to the partition
``part_x = \\{x_0 = 0, x_1, \\ldots, x_m=1\\}``
and ``part_y = \\{y_0 = 0, y_1, \\ldots, y_n=1\\}``

This version of the Ulam 2D basis uses the Skew Product structure 
(x, y) → (F(x), G(x,y)) to compute the Ulam approximation
"""
struct Ulam2DSP{T<:AbstractVector} <: Basis
    part_x::T
    part_y::T
end

"""
	Ulam2DSP(n::Integer)
	
Equispaced Ulam basis on [0,1]×[0, 1] of size n
"""
Ulam2DSP(n::Integer) = Ulam2DSP(LinRange(0.0, 1.0, n + 1), LinRange(0.0, 1.0, n + 1))

"""
	Ulam2DSP(m::Integer, n::Integer)
	
Equispaced Ulam basis on [0,1]×[0, 1] of size m in the x direction 
and size n in the y direction 
"""
Ulam2DSP(m::Integer, n::Integer) =
    Ulam2DSP(LinRange(0.0, 1.0, m + 1), LinRange(0.0, 1.0, n + 1))

@doc raw"""
	Base.length(B::Ulam2DSP)
Returns the size of the Ulam basis (the size of the underlying vector -1)
"""
Base.length(B::Ulam2DSP) = (length(B.part_x) - 1) * (length(B.part_y) - 1)

length_x(B::Ulam2DSP) = (length(B.part_x) - 1)
length_y(B::Ulam2DSP) = (length(B.part_y) - 1)

@doc raw"""
 	Base.getindex(B::Ulam2DSP, i::Int, j::Int)
    Returns the i-th element of the Ulam basis as a function.
"""
function Base.getindex(B::Ulam2DSP, i::Int, j::Int)
    return (x, y) ->
        (B.part_x[i] <= x < B.part_x[i+1] ? 1 : 0) *
        (B.part_y[j] <= y < B.part_y[j+1] ? 1 : 0)
end

"""
    rectangular_indexes_to_linear(i, j, size_col, size_row)

Maps the indexes i, j in a rectangle of sizes
size_row and size_col to linear indexes

It follows the same convention as reshape, i.e.,
the rectangle
[(1, 1)  (1, 2)  (1, 3)
(2, 1)  (2, 2)  (2, 3)
(3, 1)  (3, 2)  (3, 3)
]
is mapped to [(1, 1), (2, 1), (3, 1), (1, 2), (2, 2), (3, 2), (1, 3), (2, 3), (3, 3)]
"""
function square_indexes_to_linear(ind_row, ind_col, size_row, size_col)
    @assert ind_row <= size_row
    @assert ind_col <= size_col
    return ind_row + size_row * (ind_col - 1)
end

"""
    square_indexes_to_linear(B::Ulam2DSP, i::Integer, j::Integer)

    Remark that the indexing of the basis corresponds to the following 
    partition of [0,1]×[0,1]:
    [(1, 3)  (2, 3)  (3, 3) (4, 3)
    (1, 2)  (2, 2)  (3, 2) (4, 2)
    (1, 1)  (2, 1)  (3, 1) (4, 1)
    ],
    i.e., we start in the left right corner, the first index is the x 
    index and the second is the y index.
    Remark that this is somewhat adjoint to the matrix representation.
"""
function square_indexes_to_linear(B::Ulam2DSP, i::Integer, j::Integer)
    return square_indexes_to_linear(i, j, length_x(B), length_y(B))
end


"""
    linear_indexes_to_square(k, size_x, size_y)

Transforms a linear index k in an index in a rectangle 
of shape size_row, size_col

It follows the reshape convention as in rectangular_indexes_to_linear
"""
function linear_indexes_to_square(k, size_row, size_col)
    i_y, i_x = divrem(k, size_row)

    # we need to treat the case k = l*size_row, with l integer
    if i_x == 0
        i_x = size_row
        i_y -= 1
    end

    return (i_x, i_y + 1)
end

function linear_indexes_to_square(B::Ulam2DSP, k)
    return linear_indexes_to_square(k, length_x(B), length_y(B))
end

"""
    Return the `i`-th element of the basis when we linearize indexes.
    The linearization of indexes is coherent with the behaviour of reshape.
"""
function getindex_linear(B::Ulam2DSP, i::Int)
    return B[linear_indexes_to_square(B, i)...]
end


# We define now the Ulam2DSP Dual; each element in the Dual is a polygon
# please note that we do not save all polygons, but only the preimages on 
# the x direction
struct UlamDual2DSP <: Dual
    B::Ulam2DSP
    D::SkewProductMap
    x::Vector{Interval} #TODO: a more generic type may be needed in future
    xlabel::Vector{Int}
    lastpoint::Interval
    meshsize::Integer
end
Dual(B::Ulam2DSP, D::SkewProductMap; ϵ, max_iter, meshsize = 8) = UlamDual2DSP(
    B,
    D,
    preimages(B.part_x, D.T, 1:length(B.part_x)-1; ϵ, max_iter)...,
    domain(D.T)[end],
    meshsize,
)


# Due to the implementation details, it is not possible to estimate the length 
# of the iterator before running it
# Base.length(dual::UlamDual2DSP) = length(dual.x)*length_y(dual.B)

import Polyhedra as PH
import GLPK
lib = PH.DefaultLibrary{Float64}(GLPK.Optimizer)

# Iterating on UlamDual2DSP returns an index, a polyhedron and an error on the polyhedron
Base.eltype(dual::UlamDual2DSP) = Tuple{eltype(dual.xlabel),Array{Interval,2}}

# we define a specific assemble method for this basis, 
# to avoid the somewhat confusing Dual+iterate code
function assemble(B::Ulam2DSP, D::SkewProductMap; ϵ, max_iter, type = Float64)
    # we will break the assembly of the operator into 
    # sub-operators, one for each injectivity branch
    BranchOperator = SparseMatrixCSC[]
    for k = 1:length(branches(D.T))
        push!(BranchOperator, _assemble_branch(B, D, k; ϵ, max_iter, type))
    end
    # now we sum all the branch operators
    return sum(BranchOperator)
end


function preimage_fixed_x(D::SkewProductMap, branch_idx, x, y_min, y_max; ϵ, max_iter)
    @debug "Preimage_fixed_x"
    @debug "x", x
    @debug "y", y_min, y_max

    g(y) = D.G[branch_idx](x, y)

    try
        # this try block is used to catch the non monotone branch 
        # exception that occurs when the image collapses to a point

        H = MonotonicBranch(g, (Interval(0), Interval(1)))

        if H.increasing
            bound_Y = H.Y[1], H.Y[2]
        else
            bound_Y = H.Y[2], H.Y[1]
        end

        # disjoint
        if y_max < bound_Y[1] || y_min > bound_Y[2]
            @debug "disjoint"
            return (Interval(∅), Interval(∅))
        end

        if y_min <= bound_Y[1] && y_max >= bound_Y[2]
            @debug "contains"
            return (Interval(0), Interval(1))
        end

        preim_y_min = preimage(y_min, H; ϵ, max_iter)
        preim_y_max = preimage(y_max, H; ϵ, max_iter)

        @debug preim_y_min, preim_y_max

        if !H.increasing
            preim_y_min, preim_y_max = preim_y_max, preim_y_min
        end

        if preim_y_min == ∅
            preim_y_min = Interval(0)
        end
        if preim_y_max == ∅
            preim_y_max = Interval(1)
        end

        return (preim_y_min, preim_y_max)
    catch e
        @debug e
        @debug "Collapse to a point"
        val = g(Interval(0, 1))
        if val ⊆ Interval(y_min, y_max)
            return (Interval(0), Interval(1))
        else
            return (Interval(∅), Interval(∅))
        end
    end
end

import Polyhedra as PH
import GLPK
lib = PH.DefaultLibrary{Float64}(GLPK.Optimizer)

function rectangle_preimage(
    D::SkewProductMap,
    branch_idx,
    x_min,
    x_max,
    y_min,
    y_max,
    k;
    ϵ,
    max_iter,
)
    preim_x = [x_min + i * (x_max - x_min) / k for i = 0:k]
    err_y = 0.0
    err_x = maximum(IntervalArithmetic.radius.(preim_x))
    y_s = NTuple{2,Interval}[]

    for x in preim_x

        preim_y_low, preim_y_up =
            preimage_fixed_x(D, branch_idx, x, y_min, y_max; ϵ, max_iter)

        err_y = maximum([
            err_y,
            IntervalArithmetic.radius(preim_y_low),
            IntervalArithmetic.radius(preim_y_up),
        ])
        push!(y_s, (preim_y_low, preim_y_up))
    end
    n = length(preim_x)
    A = Matrix{Float64}(undef, 2 * n, 2)
    for (i, x) in enumerate(preim_x)
        A[i, :] = [mid(x) mid(y_s[i][1])]
    end

    for (i, x) in enumerate(reverse(preim_x))
        A[n+i, :] = [mid(x) mid(y_s[end-i+1][2])]
    end

    P = PH.polyhedron(PH.vrep(A), lib)
    return P, max(err_x, err_y)
end

# function assemble(B::Ulam2DSP, D::SkewProductMap; ϵ, max_iter, type = Float64)
#     P_branches = [_assemble_branch(B, D, i; ϵ, max_iter, type) for i in 1:nbranches(D)]
#     return sum(P_branches)
# end

function _assemble_branch(B::Ulam2DSP, D::SkewProductMap, branch_idx; ϵ, max_iter, type)
    I = Int64[]
    J = Int64[]
    nzvals = Interval{type}[]

    T = D.T.branches[branch_idx]
    G = D.G[branch_idx]

    # we first compute the preimages in the x direction
    preim_x, label_x = preimages(B.part_x, T; ϵ, max_iter)
    # it is important to remember to include the last endpoint
    if T.increasing
        preim_x = [preim_x; T.X[2]]
    else
        preim_x = [preim_x; T.X[1]]
    end

    # now we have two vectors, call j = label[i],
    # then (preim_x[i], preim_x[i+1]) = T^{-1}(I_j)
    # the vector preim is ordered with respect to <
    # while the label vector may be in the reverse order

    for i = 1:length(label_x)
        x_l, x_r = preim_x[i], preim_x[i+1]
        ind_im_x = label_x[i]

        # we compute a bound in the indexes that are intersected 
        # in the vertical direction 
        # is sent into a unique vertical interval
        # if this is true, we do not need to worry about polygon intersections
        ind_im_y_lo, ind_im_y_hi = check_image(B, G, x_l, x_r)

        if (ind_im_y_hi - ind_im_y_lo) == 1
            @info ind_im_y_lo
            # in this case, the problem reduces to a one dimensional estimate, 
            # since the full vertical stripe is sent into 
            # a single vertical element, i.e., 
            # the stripe F([x_l, x_r]×[0, 1]) ⊂ I_{ind_im_x} × I_{ind_im_y_lo}
            ind_x_lo, ind_x_hi = nonzero_on_x(B, x_l, x_r)
            for i_x = ind_x_lo:ind_x_hi
                meas = relative_measure(
                    (x_l, x_r),
                    (Interval(B.part_x[i_x]), Interval(B.part_x[i_x+1])),
                )
                # we now need to fill in this value into the matrix
                for i_y = 1:length(B.part_y)-1
                    @info i_x, i_y
                    @info ind_im_x, ind_im_y_lo

                    i = square_indexes_to_linear(B, i_x, i_y)
                    j = square_indexes_to_linear(B, ind_im_x, ind_im_y_lo)

                    push!(I, i)
                    push!(J, j)
                    push!(nzvals, meas)
                end
            end
        end

        # we need now to treat the case when we have nontrivial intersection: 
        # this method is not well behaved when the endpoint collapse to a point
        # but this should be treated by the former function
        for ind_im_y = ind_im_y_lo:ind_im_y_hi-1
            @info ind_im_y_lo, ind_im_y_hi
            preimP, err = rectangle_preimage(
                D,
                branch_idx,
                x_l,
                x_r,
                B.part_y[ind_im_y],
                B.part_y[ind_im_y+1],
                10;
                ϵ,
                max_iter,
            )
            ind_x_lo, ind_x_hi = nonzero_on_x(B, x_l, x_r)
            for i_x = ind_x_lo:ind_x_hi
                for i_y = 1:length(B.part_y)-1
                    # we now need to compute the intersection and fill it in into the matrix
                    Q_x_l, Q_x_r = B.part_x[i_x], B.part_x[i_x+1]
                    Q_y_l, Q_y_u = B.part_y[i_y], B.part_y[i_y+1]

                    Q_vertices = [Q_x_l Q_y_l; Q_x_r Q_y_l; Q_x_r Q_y_u; Q_x_l Q_y_u]
                    Q = PH.polyhedron(PH.vrep(Q_vertices), lib)
                    preimPintersectQ = PH.intersect(preimP, Q)
                    meas = PH.volume(preimPintersectQ) / PH.volume(Q)

                    i = square_indexes_to_linear(B, i_x, i_y)
                    j = square_indexes_to_linear(B, ind_im_x, ind_im_y_lo)

                    push!(I, i)
                    push!(J, j)
                    push!(nzvals, meas)
                end
            end
        end
    end
    n = length(B)
    return sparse(I, J, nzvals, n, n)
end

function check_image(B::Ulam2DSP, G, x_l, x_r)
    x = hull(x_l, x_r)
    im_bound = hull(G(x, 0), G(x, 1))
    @debug im_bound
    y = im_bound * length(B.part_y)
    y_ceil = Int64(ceil(y.hi))
    y_floor = Int64(floor(y.lo))
    return y_floor, y_ceil
end

# this function is a copy and paste
# of UlamBasis.jl:BasisDefinition.nonzero_on(B::Ulam, (a, b))
function nonzero_on_x(B::Ulam2DSP, a, b)
    y = hull(a, b)

    # finds in which semi-open interval [p[k], p[k+1]) y.lo and y.hi fall
    lo = searchsortedlast(B.part_x, y.lo)
    hi = searchsortedlast(B.part_x, y.hi)

    # they may be n+1 if y.hi==1
    lo = clamp(lo, 1, length(B.part_x))
    hi = clamp(hi, 1, length(B.part_x))

    return (lo, hi)
end

# the state variable is an implementation detail of the iterator, we use it to pass 
# information to the next iteration; in particular
# state = (index on x, index on y, lower bound index y, upper bound index y)


# function Base.iterate(dual::UlamDual2DSP, state = (1, 0, 0, 0))
#     B = dual.B
#     D = dual.D
#     meshsize = dual.meshsize

#     n = length(dual.x)
#     len_y = length_y(B)

#     i_x = state[1]
#     i_y = state[2]
#     lowerbound_y = state[3]
#     upperbound_y = state[4]

#     a, b = (dual.x[i_x], dual.x[i_x+1])


#     # when i_y is equal to 0 we compute the upper and lower bounds
#     # in the y direction, to avoid computing the preimages of 
#     # uninteresting objects
#     if i_y == 0 

#         im_a = hull(D.G(a, 0), D.G(a, 1)) 
#         im_b = hull(D.G(b, 0), D.G(b, 1))
#         im_y = hull(im_a, im_b)

#         # TODO, here it needs a cast to integer
#         lowerbound_y, upperbound_y = len_y*(im_y)

#         # if the difference between lower_bound and upper_bound 
#         # is smaller or equal to 1, this means that the full 
#         # rectangle [a,b] × [0, 1] is sent into one a thin horizontal strip
#         # this is the simplest case
#         if (upper_bound_y-lower_bound_y)<=1
#             A = [a Interval(0.0);
#                 b Interval(0.0);
#                 b Interval(1.0);
#                 a Interval(1.0)]
#             ind_image = square_indexes_to_linear(dual.B, i_x, lower_bound_y)

#             return (ind_image, A), (state[1]+1, 0, 0, 0)
#         else
#             A = [Interval(0.0) Interval(0.0)]
#             return (0, A), (i_x, lower_bound_y, lowerbound_y, upperbound_y)
#         end
#     else
#         # we start by building a mesh between a and b
#         preim_x = [a+(b-a)/meshsize*i for i in 0:meshsize]

#         A = Matrix{Interval}(Interval(0.0), 2*meshsize, 2)
#         for x in preim_x
#             preim_y_low = min(max(G_inv(y_lower), 0), 1)
#             preim_y_up = max(min(G_inv(y_upper), 1), 0)
#         end

#         for (i, x) in enumerate(preim_x)
#             A[i, :] = [mid(x) mid(y_s[i][1])]
#         end

#         for (i, x) in enumerate(reverse(preim_x))
#             A[n+i, :] = [mid(x) mid(y_s[end-i+1][2])]
#         end
#     end
#     @error "behaving badly"
# end





#function BasisDefinition.is_dual_element_empty(::Ulam2DSP, d)
#	return isempty(d[1])
#end
