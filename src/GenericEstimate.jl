using LinearAlgebra
using Arpack
using FastRounding

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

"""
This function returns the bound on the weak norm of the discretized operator
"""
function BoundNorm(B::Basis, P::AbstractMatrix{Interval{T}}, m) where {T}
	W₁, W₂ = bound_weak_norm_from_linalg_norm(B)
	α₁ = bound_linalg_norm_L1_from_weak(B)
	α₂ = bound_linalg_norm_L∞_from_weak(B)
	C, tilde_C = ContractMatrix(B, P, m)
	return (W₁/α₁)*C+(W₂/α₂)*tilde_C	 
end


"""
This function bounds the norms of the powers of the abstract discretized operator:
the computed discretized operator is an interval of matrices 
that contains the abstract discretized operator.
Therefore, to get a rigorous bound we need to have an a priori bound on the norms
of the powers of the abstract discretized operator
"""
function boundstrongauxnormabsdiscroperator(B::Basis, D::Dynamic, k)
	A, B = dfly(B, D)
	E = boundweak(B)
	h = 1/length(B) ### TODO, add a function
	M₁ = boundstrongbyweak(B)
	M₂ = boundauxiliarybyweak(B)
	
	SmallMatrix = [1 0; E*h 1]*[A B; 0 1]
	return (SmallMatrix^k)*[M₁/h; M₂]
end

function boundweaknormabsdiscroperator(B::Basis, D::Dynamic, k)
	S₁, S₂ = boundweakbystrongauxiliary(B)
	return [S₁ S₂]*boundstrongauxnormabsdiscroperator(B, D, k)
end


"""
This function bounds the norm of a finer operator by using the norms of a 
coarse operator
"""
sanity_check(Bone::Basis, Btwo::Basis) = Val((typeof(Bone)==typeof(Btwo)) && (lenght(Bone)<length(Btwo)))

function CoarseFine(Bcoarse::Basis, Bfine::Basis, Pfine, D::Dynamic, C)
	return _CoarseFine(Bcoarse, Bfine, sanity_check(Bcoarse, Bfine), Pfine, D, C)
end

_CoarseFine(Bcoarse, Bfine, ::Val{false}, D, C) = @error "Not the same basis or Coarse>Fine"

function _CoarseFine(Bcoarse, Bfine, ::Val{true}, Pfine, D, C)
	n =length(C)
	# please remark that due to the indexes in julia starting with 1,
	# ```R_{k,h,1} = R[k+1]``` and the following vector has length n+2
	R = [boundstrongauxnormabsdiscroperator(Bfine, D, i)[1] for i in 0:n+1]
	Q = matrix_norm_estimate(Bfine, Pfine)
	Cfine = undef(Float64, n)
	K =  normapprox(Bcoarse)
	h = 1/length(Bcoarse)
	for m in 1:n
		temp = 0
		for k in 0:m-1 
			temp += C[m-k]*(Q*R[k+1]+R[k+2])
		end
		Cfine[i] = C[i]+2*K*h*temp
	end 
	return Cfine[i]
end