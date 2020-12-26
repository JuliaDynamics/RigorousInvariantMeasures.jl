module InducedLSVMapDefinition
using ValidatedNumerics
using ..DynamicDefinition, ..Contractors

export ApproxInducedLSV, preim, nbranches, plottable

"""
This class introduces an ApproximatedLSV
the induced map for the Liverani-Saussol-Vaienti maps
on the interval I = [0.5, 1].
The interval I is then mapped to [0,1]
"""

struct ApproxInducedLSV <: Dynamic
	TnlistPlottable::Array{Function, 1}
	nbranches::Integer
	domains::Array{Interval, 1}
	α::Real
	T::Function
end

# coordinate change that maps [0.5, 1] to [0, 1]
@inline CoordinateChange(x) = 2*x-1
# inverse coordinate change
@inline InvCoordinateChange(x) = x/2+0.5

function ShootingLSV(n, y, α, rigstep = 10; T = Float64)
	x = [Interval{T}(0.5, 1); Interval{T}(0, 0.5)*ones(Interval{T}, n-1)]
	f(x) = 0<=x<=0.5 ? x*(1+(2*x)^α) : 2x-1
	fprime(x) = 0<=x<=0.5 ? 1+(α+1)*(2*x)^α : 2
	return ShootingMethod(f, fprime, n, x, y, rigstep)
end

function GetDomains(branches, α; T = Float64)
    domains=Interval{T}[]
    left = Interval{T}(0.5)
    for i in branches:-1:2
    	right = ShootingLSV(i, 0.5, α)[1]
		push!(domains,union(left,right))
    	left = right
    #for i=2:branches
    #    push!(domains, interval(ShootingLSV(i, 0.5, α)[1].lo, ShootingLSV(i-1, 0.5, α)[1].hi))
    end
    push!(domains, union(left, Interval{T}(0.75)))
    push!(domains, union(Interval{T}(0.75), Interval{T}(1)))

    return domains
end

function DefineTnPlottable(branches, α)
    Tleft(x)=x*(1+(2*x)^α)

    Tnlist = Array{Function, 1}(undef, branches+1)
    Tnlist[branches+1] = x->2*x-1

    for i=branches:-1:2
        Tnlist[i] = x-> Tleft(Tnlist[i+1](x))
    end

    right = mid(ShootingLSV(branches, 0.5, α)[1])

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


function iterate_LSV(x, i, α)
	@assert i>0
	x = 2*x-1
	for j in 2:i
		x = x*(1+(2*x)^α)
	end
	return x
end



end



# this function belongs to the InvariantMeasures namespace,
# this is the reason why we define it outside
using TaylorSeries: Taylor1





function dfly(::Type{TotalVariation}, ::Type{L1}, D::InvariantMeasures.InducedLSVMapDefinition.ApproxInducedLSV)
	dist = 0
	lam = 0
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
			fprime(x) = f(TaylorSeries.Taylor1([x, 1], 1))[1]
			fsecond(x) = f(TaylorSeries.Taylor1([x, 1], 2))[2]/2
			distorsion(x)=abs(fsecond(x)/(fprime(x)^2))
			lambda(x) = abs(1/fprime(x))
    		dist = max(dist, maximise(distorsion, D.domains[i])[1].hi)
    		lam = max(lam, maximise(lambda, D.domains[i])[1].hi)
		end
	end
	return lam, dist
end
