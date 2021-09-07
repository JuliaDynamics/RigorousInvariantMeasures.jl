
import TaylorModels

function integrate(f, I; steps = 1024, degree = 6)
    lo = I.lo
    hi = I.hi
    l = 2*radius(I)
    int_center = Interval(0.0)
    int_error = Interval(0.0)
    for i in 1:steps
        left = lo+(i-1)*(l/steps)
        right = lo+i*l/steps
        r = (1/2)*(l/steps)
        J = interval(left, right)
        TM = TaylorModels.TaylorModel1(degree, J)
        FM = f(TM)
        #@info FM
        for i in 0:Int64(floor(degree/2))
            int_center+=2*(FM.pol[2*i]*r^(2*i+1))/(2*i+1)
            int_error +=2*FM.rem*r
        end
    end
    return int_center+int_error
end

function adaptive_integration(f, I::Interval; tol = 2^-10, steps = 8, degree = 6) # tol 2^-10, steps = 8 are default values
    lo = I.lo
    hi = I.hi
    l = 2*radius(I)
    int_value = Interval(0.0)
    for i in 1:steps
        left = lo+(i-1)*(l/steps)
        right = lo+i*l/steps
        Istep = Interval(left, right)
        val = integrate(f, Istep)
        if radius(val)<tol
            int_value += val
        else
            I₁, I₂ = bisect(I)
            val₁ = adaptive_integration(f, I₁; tol = tol/2, steps = steps, degree = degree+2)
            val₂ = adaptive_integration(f, I₁; tol = tol/2, steps = steps, degree = degree+2)
            int_value +=val₁+val₂
        end
    end
    return int_value
end

struct Observable
    B
    v::Vector
    infbound
end


### TODO: Actually some assumptions are made, as the fact that
# the Ulam base is equispaced 
function Observable(B::Ulam, ϕ::Function; tol = 2^-10)
    v = zeros(Interval{Float64}, length(B))
    for i in 1:length(B)
        I = Interval(B.p[i], B.p[i+1])
        v[i] = adaptive_integration(ϕ, I; tol = tol, steps =1, degree = 2)*length(B)
    end
    infbound = maximise(x->abs(ϕ(x)), Interval(0,1))[1]
    return Observable(B, v, infbound)
end

import TaylorSeries

function discretizationlogder(B, D::PwMap; degree = 7)
    v = zeros(Interval{Float64}, length(B))
    infbound  = emptyinterval()
    endpoints = D.endpoints
    delicate_indexes = Int64.(floor.(length(B)*[x.lo for x in endpoints])).+1
    for i in 1:(length(delicate_indexes)-1)
        for j in delicate_indexes[i]:delicate_indexes[i+1]-1
            I = Interval(B.p[j], B.p[j+1])
            r = Interval(radius(I))
            Tmid = TaylorSeries.Taylor1([Interval(mid(I)), Interval(1)], degree)
            Tint = TaylorSeries.Taylor1([I, Interval(1)], degree)
            Fmid = log(TaylorModels.derivative(D.Ts[i](Tmid)))
            Fint = log(TaylorModels.derivative(D.Ts[i](Tint)))   
            infbound = infbound ∪ abs(Fint[0])
            ϵ = mag(Fint[degree]-Fmid[degree]) 
            
            for k in 0:Int64(floor(Float64(degree)/2))
                v[j]+=(Fmid[2*k]*r^(2*k+1))/(2*k+1)            
            end
            v[j]+=Interval(-ϵ, ϵ)*r^(degree+1)/(degree+1)
        end
    end

    #correction since the endpoints may be wide intervals
    for i in 2:length(endpoints)-1
        x = endpoints[i]
        Tx = TaylorSeries.Taylor1([x, Interval(1)], 1)

        corr = 2*Interval(radius(x))*(abs(log(TaylorModels.derivative(D.Ts[i-1](Tx))))
                +abs(log(TaylorModels.derivative(D.Ts[i](Tx)))))[0]
        v[delicate_indexes[i]]+=corr
    end
    v*=length(B)
    return Observable(B, v, infbound)
end

function integrateobservable(B::Ulam, ϕ::Observable, f::Vector, error)
    val = (ϕ.v)'*f 
    return (val+(ϕ.infbound.hi)*Interval(-error, error))/length(B)
end