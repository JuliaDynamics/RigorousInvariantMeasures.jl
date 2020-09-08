module GenericEstimate

using LinearAlgebra, Arpack, FastRounding, ValidatedNumerics
using ..DynamicDefinition, ..BasisDefinition

function perronvector(B::Basis, P::AbstractMatrix{Interval{T}}) where {T}
	PP = mid.(P)
	F = eigs(PP; nev=4, ritzvec=true, v0=ones((0,)))
	return F[2][:, 1]
end

function rigorousresidual(B::Basis, P::AbstractMatrix{Interval{T}}, w) where {T}
	w_int = Interval{T}(w)
	Pw_int = P*w_int
	return BasisDefinition.rigorous_weak_norm(B, Pw_int-w_int)
end


"""
This function returns a sequence of Cᵢ, \\tilde{C}ᵢ for a matrix P
on a subspace V such that ||P^i|_V||_1\\leq C_i and
||P^i|_V||_{\\infty}\\leq \\tilde{C}_i, with respect to the
"""
function contractmatrix(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
	# vector of the Cᵢ for P
	C = zeros(m)
	S = zeros((length(B), m))
	tilde_C = zeros(m)

	PP = mid.(P)

	for v in BasisDefinition.AverageZero(B)
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
	η₁ = BasisDefinition.spaceconstant(B, Val(:L1))
	η₂ = BasisDefinition.spaceconstant(B, Val(:L∞))

 	return η₁*C, η₂*tilde_C
end

@deprecate contractmatrix norms_of_powers

"""
This function returns the bound on the weak norm of the discretized operator
"""
function boundnorm(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
	W₁, W₂ = BasisDefinition.bound_weak_norm_from_linalg_norm(B)
	α₁ = BasisDefinition.bound_linalg_norm_L1_from_weak(B)
	α₂ = BasisDefinition.bound_linalg_norm_L∞_from_weak(B)
	C, tilde_C = contractmatrix(B, P, m)
	return (W₁/α₁)*C+(W₂/α₂)*tilde_C
end


"""
This function bounds the norms of the powers of the abstract discretized operator:
the computed discretized operator is an interval of matrices
that contains the abstract discretized operator.
Therefore, to get a rigorous bound we need to have an a priori bound on the norms
of the powers of the abstract discretized operator
"""
function boundstrongauxnormabsdiscroperator(Bas::Basis, D::Dynamic, k)
	@warn "Must be rewritten making sure that operations are correctly rounded"
	A, B = BasisDefinition.dfly(Bas, D)
	E = BasisDefinition.boundweak(Bas)
	h = 1/length(Bas) ### TODO, add a function
	M₁ = BasisDefinition.boundstrongbyweak(Bas)
	M₂ = BasisDefinition.boundauxiliarybyweak(Bas)

	SmallMatrix = [1 0; E*h 1]*[A B; 0 1]
	return (SmallMatrix^k)*[M₁/h; M₂]
end

function boundweaknormabsdiscroperator(B::Basis, D::Dynamic, k)
	@warn "Must be rewritten making sure that operations are correctly rounded"
	S₁, S₂ = BasisDefinition.boundweakbystrongauxiliary(B)
	return [S₁ S₂]*BasisDefinition.boundstrongauxnormabsdiscroperator(B, D, k)
end


"""
This avoids CoarseFine to be compiled for incompatible basis
"""
sanity_check(Bone::Basis, Btwo::Basis) = Val((typeof(Bone)==typeof(Btwo)) && (length(Bone)<length(Btwo)))

"""
This function bounds the norm of a finer operator by using the norms of a
coarse operator
"""
function coarsefine(Bcoarse::Basis, Bfine::Basis, Pfine, D::Dynamic, C)
	return _coarsefine(Bcoarse, Bfine, sanity_check(Bcoarse, Bfine), Pfine, D, C)
end

_coarsefine(Bcoarse, Bfine, ::Val{false}, D, C) = @error "Not the same basis or Coarse>Fine"

function _coarsefine(Bcoarse, Bfine, ::Val{true}, Pfine, D, C)
	@warn "Must be rewritten making sure that operations are correctly rounded"
	n =length(C)
	# please remark that due to the indexes in julia starting with 1,
	# ```R_{k,h,1} = R[k+1]``` and the following vector has length n+2
	R = [boundstrongauxnormabsdiscroperator(Bfine, D, i)[1] for i in 0:n+1]
	Q = BasisDefinition.matrix_norm_estimate(Bfine, Pfine)
	Cfine = Array{Float64}(undef, n)
	K =  BasisDefinition.normapprox(Bcoarse)
	h = 1/length(Bcoarse)
	for m in 1:n
		temp = 0
		for k in 0:m-1
			temp += C[m-k]*(Q*R[k+1]+R[k+2])
		end
		Cfine[m] = (C[m]+2*K*h*temp).hi
	end
	return Cfine
end

end
