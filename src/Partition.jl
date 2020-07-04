
# TODO non-homogeneous partitions
"""
This function creates a non-homogeneous partition, provided:
- A list of higher level endpoints, a₁,..., aₙ
- The number of homogeneous elements in betwee each aᵢ and aᵢ+1
"""
struct Partition
	endpoints::Vector{Real}
	elements::Vector{Integer}
	indices_endpoints::Vector{Integer}
end

function Partition( endpoints::Vector{Real}, elements::Vector{Integer})
	indices_endpoints = cumsum(elements)
	return Partition(endpoints,elements, indices_endpoints)
end

#Base.getindex(P::Partition, i) = 