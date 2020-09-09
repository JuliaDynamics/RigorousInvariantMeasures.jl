using InvariantMeasures
using ValidatedNumerics
using LinearAlgebra: I, opnorm
using SparseArrays: sparse
using Test

n = 9
m = 10
M = 0.2*randn(n, n)
R = 1e-13*rand(n, n)

LL = interval_from_midpoint_radius.(M, R)

U = [ones(1,n-1); -Matrix(I, n-1,n-1)]

@test norms_of_powers(Linf, m, sparse(LL), true) ≈ [opnorm(M^k*U,Inf) for k = 1:m]

@test norms_of_powers(L1, m, sparse(LL), true) ≈ [opnorm(M^k*U,1) for k = 1:m]

e = ones(n)
f = 1/n*e'

@test norms_of_powers(Linf, m, sparse(LL), false, e=e, f=f) ≈ [opnorm((M+e*(f-f*M))^k*U,Inf) for k = 1:m]
@test norms_of_powers(L1, m, sparse(LL), false, e=e, f=f) ≈ [opnorm((M+e*(f-f*M))^k*U,1) for k = 1:m]
