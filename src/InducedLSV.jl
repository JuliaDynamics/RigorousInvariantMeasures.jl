module InducedLSVMapDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors
import ..Hat
import ..BasisDefinition: DualComposedWithDynamic
import ..C2BasisDefinition: C2Basis, dual_val, dual_der

export ApproxInducedLSV, preim, nbranches, plottable

"""
This class introduces an ApproximatedLSV
the induced map for the Liverani-Saussol-Vaienti maps
on the interval I = [0.5, 1].
The interval I is then mapped to [0,1]
"""
struct ApproxInducedLSV <: Dynamic
	TnListPlottable::Array{Function, 1}
	nbranches::Integer
	domains::Array{Interval, 1}
	α::Real
	T::Function
end

# coordinate change that maps [0.5, 1] to [0, 1]
@inline CoordinateChange(x) = 2*x-1
# inverse coordinate change
@inline InvCoordinateChange(x) = x/2+0.5

function derleft(D::ApproxInducedLSV, x) 
	@assert 0<=x<=0.5
	α = D.α
	return 1+2^α*(α+1)*x^α # check that α must be an interv)al? 
end

function derderleft(D::ApproxInducedLSV, x) 
	@assert 0<=x<=0.5
	α = D.α
	return 2^α*(α+1)*α*x^(α-1) # check that α must be an interval? 
end


function ShootingLSV(n, y, α, rigstep = 10; T = Float64)
	x = [Interval{T}(0.5, 1); Interval{T}(0, 0.5)*ones(Interval{T}, n-1)]
	f(x) = 0<=x<=0.5 ? x*(1+(2*x)^α) : 2x-1
	fprime(x) = 0<=x<=0.5 ? 1+(α+1)*(2*x)^α : 2.
	return ShootingMethod(f, fprime, n, x, y, rigstep)
end

function GetDomains(branches, α; T = Float64)
    domains=Interval{T}[]
    left = Interval{T}(0.5)
    for i in branches:-1:2
    	right = ShootingLSV(i, 0.5, α; T = T)[1]
		push!(domains,union(left,right))
    	left = right
    #for i=2:branches
    #    push!(domains, interval(ShootingLSV(i, 0.5, α)[1].lo, ShootingLSV(i-1, 0.5, α)[1].hi))
    end
    push!(domains, union(left, Interval{T}(0.75)))
    push!(domains, union(Interval{T}(0.75), Interval{T}(1)))

    return domains
end

function DefineTnPlottable(branches, α; T = Float64)
    Tleft(x)=x*(1+(2*x)^α)

    Tnlist = Array{Function, 1}(undef, branches+1)
    Tnlist[branches+1] = x->2*x-1

    for i=branches:-1:2
        Tnlist[i] = x-> Tleft(Tnlist[i+1](x))
    end

    right = mid(ShootingLSV(branches, 0.5, α; T = T)[1])

    Tnlist[1] = x-> 0.5*1/(right-0.5)*(x-0.5)+0.5
    return Tnlist
end

"""
This constructor builds the induced LSV map on [0.5, 1],
truncated with k branches
"""
function ApproxInducedLSV(α, k)
	nbranches = k+1
	domains = GetDomains(k, α)
	TnListPlottable = DefineTnPlottable(k, α)
	return ApproxInducedLSV(TnListPlottable, nbranches, domains, α, x-> 0<=x<=0.5 ? x*(1+(2*x)^α) : 2x-1)
end

domains(D::ApproxInducedLSV)=CoordinateChange.(D.domains)



DynamicDefinition.nbranches(S::ApproxInducedLSV)=S.nbranches

DynamicDefinition.is_full_branch(S::ApproxInducedLSV) = true

function _T(x, domains, TnList)
	@assert 0 ≤ x ≤ 1
	x = InvCoordinateChange(x)
	for (i, I) in enumerate(domains)
		if x in I
			return CoordinateChange(TnList[i](x))
		end
	end
	return 0
end

DynamicDefinition.plottable(D::ApproxInducedLSV, x) = _T(x, D.domains, D.TnListPlottable)

function DynamicDefinition.preim(D::ApproxInducedLSV, k, y, ϵ)
	@assert 1 <= k <= D.nbranches
	_y = InvCoordinateChange(y)
	if k == 1
		right = ShootingLSV(D.nbranches-1, 0.5, D.α)[1]
		_x = (2*_y-1)*(right-0.5)+0.5
		return CoordinateChange(_x)
	elseif k == nbranches
		_x = (_y+1)/2
		return CoordinateChange(_x)
	else
		_x = ShootingLSV(D.nbranches-k+1, _y, D.α)[1]
		return CoordinateChange(_x)
	end
end

#returns the preimage of a point in a branch with the der and the second derivative
function preimwithder_derder(D::ApproxInducedLSV, k, y, ϵ) #::NTuple{Interval, 3}
	@assert 1 <= k <= D.nbranches

	y = Interval(y) #hack, please check
	_y = InvCoordinateChange(y)
	
	if k == 1 # the manufactured branch
		right = ShootingLSV(D.nbranches-1, 0.5, D.α)[1]
		_x = (2*_y-1)*(right-0.5)+0.5  
		return CoordinateChange(_x), 0.5/(right-0.5), 0  

	elseif k == D.nbranches # the linear branch
		_x = (_y+1)/2
		return CoordinateChange(_x), 2, 0
	else		
		orbit_x = ShootingLSV(D.nbranches-k+1, _y, D.α)
		
		der = 2
		derder = 0
		for i in 2:length(orbit_x)
			dx = derleft(D, orbit_x[i])
			ddx = derderleft(D, orbit_x[i])
			derder = ddx*der^2+dx*derder
			der*=dx
		end
			
		return CoordinateChange(orbit_x[1]), der, derder/2
		# this /2 follows from the coordinate change
		# i.e., we are looking at the map F(x) = ψ(f(ϕ(x))) where
		# ϕ maps linearly [0,1]->[0.5, 1] and ψ is the
		# inverse map
		# by a direct computation we get  F''(x) = ψ'(f(ϕ(x)))f''(ϕ(x))(ϕ'(x))^2
	end
end

#returns the preimage of a point in a branch with the der and the second derivative
function preimwithder(D::ApproxInducedLSV, k, y, ϵ) #::NTuple{Interval, 3}
	@assert 1 <= k <= D.nbranches

	y = Interval(y) #hack, please check
	_y = InvCoordinateChange(y)
	
	if k == 1 # the manufactured branch
		right = ShootingLSV(D.nbranches-1, 0.5, D.α)[1]
		_x = (2*_y-1)*(right-0.5)+0.5  
		return CoordinateChange(_x), 0.5/(right-0.5), 0  

	elseif k == D.nbranches # the linear branch
		_x = (_y+1)/2
		return CoordinateChange(_x), 2, 0
	else		
		orbit_x = ShootingLSV(D.nbranches-k+1, _y, D.α)
		
		der = 2
		for i in 2:length(orbit_x)
			dx = derleft(D, orbit_x[i])
			der*=dx
		end
			
		return CoordinateChange(orbit_x[1]), der
	end
end

"""
Return (in an iterator) the pairs (i, (x, |T'(x)|)) where x is a preimage of p[i], which
describe the "dual" L* evaluation(p[i])
"""
function Base.iterate(S::DualComposedWithDynamic{T, ApproxInducedLSV}, state = (1, 1)) where T<:Hat
	i, k = state

	if i == length(S.basis)+1
			return nothing
	end

	n = length(S.basis.p)
	
	x, der, derder = preimwithder_derder(S.dynamic, k, S.basis.p[i], S.ϵ)
	ret = x, abs(der)   
	
	if k == nbranches(S.dynamic)
		return ((i, ret), (i+1, 1))
	else
		return ((i, ret), (i, k+1))
	end
end

"""
Return (in an iterator) the pairs (i, (x, |T'(x)|)) where x is a preimage of p[i], which
describe the "dual" L* evaluation(p[i])
"""
function Base.iterate(S::DualComposedWithDynamic{T, ApproxInducedLSV}, state = (1, 1)) where T<:C2Basis
	i, k = state

	if i == length(S.basis)+1
			return nothing
	end

	n = length(S.basis.p)
	
	if i <= n
		x, der, derder = preimwithder_derder(S.dynamic, k, S.basis.p[i], S.ϵ)
		ret = x, (f, fprime) -> dual_val(f, fprime, x, der, derder)   
	else
		x, der, derder = preimwithder_derder(S.dynamic, k, S.basis.p[i-n], S.ϵ)
		ret = x, (f, fprime) -> dual_der(f, fprime, x, der, derder)
	end
	if k == nbranches(S.dynamic)
		return ((i, ret), (i+1, 1))
	else
		return ((i, ret), (i, k+1))
	end
end

function DynamicDefinition.derivative(D::ApproxInducedLSV, x)
	@error "Not implemented"
end


function iterate_LSV(x, i, α)
	@assert i>0
	x = 2*x-1
	for j in 2:i
		x = x*(1+(2*x)^α)
	end
	return x
end

using RecipesBase
@recipe f(::Type{ApproxInducedLSV}, D::ApproxInducedLSV) where {FT} = x -> plottable(D, x)

end

function InvariantMeasures.Dual(B::Chebyshev, D::InducedLSVMapDefinition.ApproxInducedLSV, ϵ = 0.0; T = Float64)
	labels = Int64[]
	x = Interval{T}[]
	x′ = Interval{T}[]
	for k in 1:D.nbranches
		@info "$k-th branch"
		for i in 1:length(B.p)
		x_i, x_i_prime = InducedLSVMapDefinition.preimwithder(D, k, B.p[i], ϵ)
		append!(labels, i)
		append!(x, x_i)
		append!(x′, x_i_prime)
		end
	end
	return x, labels, x′
end



# this function belongs to the InvariantMeasures namespace,
# this is the reason why we define it outside
using TaylorSeries: Taylor1


function dfly(::Type{TotalVariation}, ::Type{L1}, D::InvariantMeasures.InducedLSVMapDefinition.ApproxInducedLSV)
	dist = @interval(0.)
	lam = @interval(0.)
	for i in 1:D.nbranches
		if i==1
			right = InducedLSVMapDefinition.ShootingLSV(D.nbranches-1, 0.5, D.α)[1]
			lam = max(lam, 2*(right-0.5).hi)
			dist =  max(dist, 0)
		elseif i==D.nbranches
			lam = max(lam, 0.5)
			dist =  max(dist, 0)
		else
			f(x) = InducedLSVMapDefinition.iterate_LSV(x, D.nbranches-i+1, D.α)
			fprime(x) = f(Taylor1([x, 1], 1))[1]
			fsecond(x) = f(Taylor1([x, 1], 2))[2]/2
			distorsion(x)=abs(fsecond(x)/(fprime(x)^2))
			lambda(x) = abs(1/fprime(x))
    		dist = max(dist, maximise(distorsion, D.domains[i])[1].hi)
    		lam = max(lam, maximise(lambda, D.domains[i])[1].hi)
		end
	end
	return lam.hi, dist.hi
end

import TaylorSeries
import TaylorSeries: Taylor1 

function derivatives_D(α, k, l; T = Float64)
	w = zeros(Interval, (l, l))
	right = Interval(1.0)
	for i in 1:k
		@info i
		left = InducedLSVMapDefinition.ShootingLSV(i, 0.5, α; T = T)[1]
		dom = hull(left, right)
		f(x) = InducedLSVMapDefinition.iterate_LSV(x, i, α)
		g(x) = 1/(TaylorSeries.derivative(f(Taylor1([x, 1], l))))
				
		dom = Interval(left.lo, right.hi)
		tol = diam(dom)*2^(-10)
		for i = 0:l-1, j=0:l-1
			h(x) = abs((factorial(i)*g(x)[i])*g(x)[0]^j) # \partial^i (1/T') * (1/T')^j
			val = maximise(h, dom, tol = tol)[1]
			w[i+1, j+1] = max(w[i+1, j+1], val) 
		end
		right = left
		if mod(i,20) == 0
			@info w
		end
	end
	return w
end

#import DualNumbers
#function bound_b_ω(α, k; T = Float64)
#	right = Interval(1.0)
#	for i in 1:k
#		@info i
#		left = InducedLSVMapDefinition.ShootingLSV(i, 0.5, α; T = T)[1]
#		dom = hull(left, right)
#		f(α, x) = InducedLSVMapDefinition.iterate_LSV(x, i, α)
#		f_prime_α(x) = f(DualNumbers.Dual(α, 1), x).epsilon
#		f_prime_x(x) = f(α, DualNumbers.Dual(x, 1)).epsilon
#		h(x) = -f_prime_α(x)/f_prime_x(x)
#		dom = Interval(left.lo, right.hi)
#		tol = diam(dom)*2^(-10)
#		val = maximise(h, dom, tol = tol)[1]
#		@info val
#	end
#	return w
#end