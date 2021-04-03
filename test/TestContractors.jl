using InvariantMeasures: root, nthpreimage!
using StaticArrays
using ValidatedNumerics

@testset "Contractors" begin

@test root(x -> x^2-2, 1..2, 1e-13) ≈ sqrt(2)
# to ensure quadratic convergence
@test root(x -> x^2-2, 1..2, 1e-13; max_iter = 6) ≈ sqrt(2)


fs = (x -> x^2, x -> x^3)
X = [0..1, 0..1]
y = 0.1
nthpreimage!(X, fs, y; max_iter = 100)
@test X[1]^6 ≈ y
@test X[2]^3 ≈ y
X = @MVector[0..1, 0..1]
nthpreimage!(X, fs, y; max_iter = 9)
@test X[1]^6 ≈ y
@test X[2]^3 ≈ y

fs = (x -> x^2, x -> x^3, x -> x^4)
X = @MVector[0..1, 0..1, 0..1]
nthpreimage!(X, fs, y)
@test X[1]^24 ≈ y

end #testset
