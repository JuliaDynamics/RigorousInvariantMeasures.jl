# This implements the contractions time estimator for the "platonic" 
# operator L, the original small matrix method in Galatolo-Nisoli-Saussol
# An elementary way to rigorously estimate convergence to equilibrium and escape rates

abstract type TrinormWeakCompare end

struct WeakStrongerTrinorm <: TrinormWeakCompare end
struct WeakTrinormUncomparable <: TrinormWeakCompare end


comparetrinormweak(Bas::Basis) = WeakTrinormUncomparable 
comparetrinormweak(Bas::Ulam) = WeakStrongerTrinorm 

#convergencerateabstract(Bas::Basis, D::Dynamic, norms) = _convergenceratesabstract(Bas, D, norms, comparetrinormweak(Bas))

#function _convergenceratesabstract(Bas::Basis, D::Dynamic, norms, ::WeakStrongerTrinorm)
#    boundL = bound_weak_norm_abstract(Bas)
#    A, B = dfly(strong_norm(Bas), aux_norm(Bas), D)
#    m = length(coarse_norms)
#    
#    Kh =  BasisDefinition.weak_projection_error(coarse_basis)
#    
#    strong_norms = fill(NaN, m+1) 
#    weak_norms = fill(NaN, m+1) 
    
#    strong_norms[1] = 1.
#    weak_norms[1] = 1.
    
#    for i in 1:m

using LinearAlgebra
function eig_costants_small_matrix(A)
    Aint = Interval.(A)
    #@info "1, 1 entry" A[1,1]
    #@info "2, 2 entry" A[1,2]
    ρ = Aint[1, 1]+Aint[2, 2]+sqrt((Aint[1,1]-Aint[2,2])^2- 4*Aint[2,1]*Aint[1, 2])
    
    B = Aint'-ρ*I

    v = [B[2, 1] ; ρ-B[1,1]]
    
    v = v/(v[1]+v[2])
    
    return ρ, v 
end



function convergencerateabstract(Bas::Ulam, D::Dynamic, norms)
    boundL = BasisDefinition.bound_weak_norm_abstract(Bas)
    A, B = dfly(strong_norm(Bas), aux_norm(Bas), D)
    @info A, B
    m = length(norms) 
    Kh =  BasisDefinition.weak_projection_error(Bas)

    C = (A⊕₊ 1.0)/(1.0 ⊖₋A)
    D = B ⊗₊ (A⊕₊2.0)

    small_matrices = [[A^n B; Kh*C Kh*n*D+norms[n]] for n in 1:length(norms)] 

    weak_norms = ones(length(norms)) ### This is specific to Ulam
    strong_norms = [A^i+B for i in 1:length(norms)]

    for i in 1:length(small_matrices)
        ρ, v = eig_costants_small_matrix(small_matrices[i])
        if ρ<1 && ρ != ∅
            @info ρ
            strong_norms_here = [strong_norms[1:i-1];[((1/v[1]+B/v[2])*ρ^floor(j÷i)).hi for j in i:length(norms)]]
            strong_norms = min.(strong_norms_here, strong_norms)
            weak_norms_here = [weak_norms[1: i-1]; [((B/v[2])*ρ^floor(j÷i)).hi for j in i:length(norms)]]
            weak_norms = min.(weak_norms_here, weak_norms)
        end
    end
    weak_norms = refine_norms_of_powers(weak_norms, length(weak_norms))
    strong_norms = refine_norms_of_powers(strong_norms, length(strong_norms))
    return weak_norms , strong_norms
end



