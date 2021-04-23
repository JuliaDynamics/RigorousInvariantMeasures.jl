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

function convergencerateabstract(Bas::Ulam, D::Dynamic, norms)
    boundL = BasisDefinition.bound_weak_norm_abstract(Bas)
    A, B = dfly(strong_norm(Bas), aux_norm(Bas), D)
    m = length(norms) 
    Kh =  BasisDefinition.weak_projection_error(Bas)

    C = (A⊕₊ 1.0)/(1.0 ⊖₋A)
    D = B⊗₊ (A⊕₊2.0)

    small_matrices = [[A^n B; Kh*C Kh*n*D+norms[n]] for n in 1:length(norms)] 

    

    return small_matrices
end



