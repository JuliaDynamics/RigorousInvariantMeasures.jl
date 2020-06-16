using LinearAlgebra
using Arpack

function PerronVector(B::Basis, P::AbstractMatrix{Interval{T}}) where {T}
	PP = mid.(P)
	F = eigs(PP; nev=4, ritzvec=true, v0=ones((0,)))
	return F[2][:, 1]
end

function RigorousResidual(B::Basis, P::AbstractMatrix{Interval{T}}, w) where {T}
	w_int = Interval{T}(w)
	Pw_int = P*w_int
	return rigorous_weak_norm(B, Pw_int-w_int)
end

"""
This function bounds the weak norm of the abstract discretized operator:
the computed discretized operator is an interval of matrices 
that contains the abstract discretized operator.
Therefore, to get a rigorous bound we need to have an a priori bound on the norms
of the powers of the abstract discretized operator
"""
function boundweaknormabsdiscroperator(B::Basis, D::Dynamic, k)
	A, B = dfly(B, D)
	E = boundweak(B)
	h = 1/length(B) ### TODO, add a function
	M₁ = boundstrongbyweak(B)
	M₂ = boundauxiliarybyweak(B)
	S₁, S₂ = boundweakbystrongauxiliary(B)

	SmallMatrix = [1 0; E*h 1]*[A B; 0 1]
	return [S₁ S₂]*(SmallMatrix^k)*[M₁/h; M₂]
end


"""
This function returns a sequence of Cᵢ, \\tilde{C}ᵢ for a matrix P
on a subspace V such that ||P^i|_V||_1\\leq C_i and
||P^i|_V||_{\\infty}\\leq \\tilde{C}_i, with respect to the
"""
function ContractMatrix(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
	# vector of the Cᵢ for P
	C = zeros(m)
	S = zeros((length(B), m))
	tilde_C = zeros(m)

	PP = mid.(P)

	for v in AverageZero(B)
		λ₁, λ₂ = norm(v, 1), norm(v, Inf)
		for i in 1:m
			v = PP*v
			C[i] = max(C[i], norm(v, 1)/λ₁)
			S[:, i]+=abs.(v)/λ₂
		end	
	end


	for k in 1:m
		tilde_C[k] = maximum(S[:,k])
	end

	# we keep track of the error due to the basis we have chosen
	η₁ = SpaceConstant(B, Val(:L1))
	η₂ = SpaceConstant(B, Val(:L∞))

 	return η₁*C, η₂*tilde_C 
end	







