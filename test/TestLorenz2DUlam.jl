@testset "Lorenz2DUlam" begin

A = [(ind_y, ind_x) for ind_y in 1:20, ind_x in 1:10]

Q = reshape(A, 100)

@test Q[square_indexes_to_linear(5, 10)] = A[10, 5]

end