using ..BasisDefinition, ..DynamicDefinition, ..Contractors, ..PwDynamicDefinition
IntervalArithmetic, LinearAlgebra

export Ulam2DSP, length_x, length_y, square_indexes_to_linear, linear_indexes_to_square, getindex_linear

"""
	Ulam2DSP
Ulam basis on [0,1]×[0,1] associated to the partition
``p_x = \\{x_0 = 0, x_1, \\ldots, x_m=1\\}``
and ``p_y = \\{y_0 = 0, y_1, \\ldots, y_n=1\\}``

This version of the Ulam 2D basis uses the Skew Product structure 
(x, y) → (F(x), G(x,y)) to compute the Ulam approximation
"""
struct Ulam2DSP{T<:AbstractVector} <:Basis
	p_x::T
	p_y::T
end

"""
	Ulam2DSP(n::Integer)
	
Equispaced Ulam basis on [0,1]×[0, 1] of size n
"""
Ulam2DSP(n::Integer) = Ulam2DSP(LinRange(0., 1., n+1), LinRange(0., 1., n+1))

"""
	Ulam2DSP(m::Integer, n::Integer)
	
Equispaced Ulam basis on [0,1]×[0, 1] of size m in the x direction 
and size n in the y direction 
"""
Ulam2DSP(m::Integer, n::Integer) = Ulam2DSP(LinRange(0., 1., m+1), LinRange(0., 1., n+1))

@doc raw"""
	Base.length(B::Ulam2DSP)
Returns the size of the Ulam basis (the size of the underlying vector -1)
"""
Base.length(B::Ulam2DSP) = (length(B.p_x) - 1)*(length(B.p_y) - 1)

length_x(B::Ulam2DSP) = (length(B.p_x) - 1)
length_y(B::Ulam2DSP) = (length(B.p_y) - 1)

@doc raw"""
 	Base.getindex(B::Ulam2DSP, i::Int, j::Int)
    Returns the i-th element of the Ulam basis as a function.
"""
function Base.getindex(B::Ulam2DSP, i::Int, j::Int)
	return (x, y) -> (B.p_x[i]<= x < B.p_x[i+1] ? 1 : 0)*(B.p_y[j]<= y < B.p_y[j+1] ? 1 : 0)
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
    @assert ind_row<=size_row
    @assert ind_col<=size_col
    return ind_row+size_row*(ind_col-1)
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
        i_y-=1
    end

    return (i_x, i_y+1)
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


function BasisDefinition.is_dual_element_empty(::Ulam2DSP, d)
	return isempty(d[1])
end