#computing the Ulam scheme for the 2-D Lorenz

# I start by implementing the preimage of a square

const α = 1.75::Float64
const s = 3.0::Float64

using Plots

using RigorousInvariantMeasures

T(; α , s) = PwMap([x-> -α*(-x)^s+1, x-> α*(x)^s-1], [-1, 0, 1])

Ttrue = T(α = α, s =s)

preim_x = RigorousInvariantMeasures.preimages(-1:0.1:1, Ttrue)

Tfloat(x) = x > 0 ? α*(x)^s-1 : -α*(-x)^s+1

plot(Tfloat, -1, 1)

@info preim_x[1][1]

G(; x, r, c) = x > 0 ?  
                RigorousInvariantMeasures.Branch(y-> 2^(-r)*y*x^r+c, y-> 2^(-r)*x^r, (Interval(-1), Interval(1))) :
                RigorousInvariantMeasures.Branch(y-> 2^(-r)*y*(-x)^r-c, y-> 2^(-r)*(-x)^r, (Interval(-1), Interval(1)))

Gtrue = G(x = (preim[1])[10], r = 5.0, c = 0.5)

using IntervalArithmetic
plot(x->mid(Gtrue.f(x)), -1, 1)
 
dx = -1:0.1:-0.1
dy = -1:0.001:1

test_values = NTuple{2, Float64}[]

for x in preim_x[1][20:end]
    #@info x
    Gtrue = G(x = x, r = 5.0, c = 0.5)
    for y in dy
        #@info mid(Tfloat(x)), mid(Gtrue.fprime(y))
        push!(test_values, (mid(Tfloat(x)), mid(Gtrue.f(y))))
    end
end
scatter(test_values)

test_values = NTuple{2, Float64}[]
for x in preim_x[1][20:end]
    #@info x
    Gtrue = G(x = x, r = 5.0, c = 0.5)
    preim_y = RigorousInvariantMeasures.preimages(dy, Gtrue)[1]
    #@info preim_y 
    for y in preim_y
        #@info x, y
        push!(test_values, (mid(x), mid(y)))
    end
end
scatter(test_values)

const r = 5.0
test(x) = 2^(-r)*x^r
plot!(Tfloat.(mid.(preim[1][20:end])), mid.(test.(preim[1][20:end])))


struct Square
    vertices::NTuple{4, NTuple{2, Interval}}
end

Square(;xl = a, xr = b, yl = c, yu =d) = Square((xr,yu), (xr, yl), (xl, yl), (xl, yu))


struct PolygonLorenz
    x::Vector{Interval}
    y::Vector{Interval} 
end

function Square_PL(;x_left, x_right, y_lower, y_upper, k)
    x_lo = [Interval(x_left)+i*Interval(x_right-x_left)/k for i in 0:k]
    y_l = [Interval(y_lower) for i in 0:k]
    x_up = [Interval(x_right)+i*Interval(x_left-x_right)/k for i in 0:k]
    y_u = [Interval(y_upper) for i in 0:k]
    return PolygonLorenz([x_lo; x_up], [y_l; y_u])
end 
