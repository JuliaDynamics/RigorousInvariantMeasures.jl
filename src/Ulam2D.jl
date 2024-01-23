#computing the Ulam scheme for the 2-D Lorenz

# I start by implementing the preimage of a square

const α = 1.75::Float64
const s = 3.0::Float64

using Plots

using RigorousInvariantMeasures

T(; α, s) = PwMap([x -> -α * (-x)^s + 1, x -> α * (x)^s - 1], [-1, 0, 1])

Ttrue = T(α = α, s = s)

ϕ = PwMap([x -> 2 * x - 1], [0, 1])
ψ = PwMap([x -> x / 2 + 1 / 2], [-1, 1])

# it works!!! Incredible!!!
D = ψ ∘ Ttrue ∘ ϕ


preim_x = RigorousInvariantMeasures.preimages(-1:0.1:1, Ttrue)

Tfloat(x) = x > 0 ? α * (x)^s - 1 : -α * (-x)^s + 1

plot(Tfloat, -1, 1)

@info preim_x[1][1]

G(; x, r, c) =
    x > 0 ?
    RigorousInvariantMeasures.MonotonicBranch(
        y -> 2^(-r) * y * x^r + c,
        y -> 2^(-r) * x^r,
        (Interval(-1), Interval(1)),
    ) :
    RigorousInvariantMeasures.MonotonicBranch(
        y -> 2^(-r) * y * (-x)^r - c,
        y -> 2^(-r) * (-x)^r,
        (Interval(-1), Interval(1)),
    )

Gtrue = G(x = (preim_x[1])[10], r = 5.0, c = 0.5)

using IntervalArithmetic
plot(x -> mid(Gtrue.f(x)), -1, 1)

dx = -1:0.1:-0.1
dy = -1:0.001:1

test_values = NTuple{2,Float64}[]

# for x in preim_x[1][20:end]
#     #@info x
#     Gtrue = G(x = x, r = 5.0, c = 0.5)
#     for y in dy
#         #@info mid(Tfloat(x)), mid(Gtrue.fprime(y))
#         #push!(test_valueswidth :end]
#     #@info x
#     Gtrue = G(x = x, r = 5.0, c = 0.5)
#     preim_y = RigorousInvariantMeasures.preimages(dy, Gtrue)[1]
#     #@info preim_y 
#     for y in preim_y
#         #@info x, y
#         push!(test_values, (mid(x), mid(y)))
#     end
# end
# scatter(test_values)

const r = 5.0
# test(x) = 2^(-r)*x^r
# plot!(Tfloat.(mid.(preim[1][20:end])), mid.(test.(preim[1][20:end])))


struct Square
    vertices::NTuple{4,NTuple{2,Interval}}
end

Square(; xl = a, xr = b, yl = c, yu = d) = Square((xr, yu), (xr, yl), (xl, yl), (xl, yu))


using IntervalArithmetic

struct PolygonLorenz
    x::Vector{Interval}
    y::Vector{Interval}
end

function Square_PL(; x_left, x_right, y_lower, y_upper, k)
    x_lo = [Interval(x_left) + i * Interval(x_right - x_left) / k for i = 0:k]
    y_l = [Interval(y_lower) for i = 0:k]
    x_up = [Interval(x_right) + i * Interval(x_left - x_right) / k for i = 0:k]
    y_u = [Interval(y_upper) for i = 0:k]
    return PolygonLorenz([x_lo; x_up], [y_l; y_u])
end

import Polyhedra as PH
import GLPK
lib = PH.DefaultLibrary{Float64}(GLPK.Optimizer)

P1 = PH.polyhedron(PH.vrep([
    -1.9 -1.7
    -1.8 0.5
    1.7 0.7
    1.9 -0.3
    0.9 -1.1
]), lib)

P2 = PH.polyhedron(PH.vrep([
    -2.5 -1.1
    -0.8 0.8
    0.1 0.9
    1.8 -1.2
    1.3 0.1
]), lib)

Pint = PH.intersect(P1, P2)
PH.volume(Pint)

import IntervalArithmetic

Ginverse(; x, r, c) = x > 0 ? v -> ((v - c) * 2^r) / (x^r) : v -> ((v + c) * 2^r) / ((-x)^r)

# function PreimageRectangleLorenz(; x_left, x_right, y_lower, y_upper, k)
#     x = [Interval(-1); [Interval(x_left)+i*Interval(x_right-x_left)/k for i in 0:k]; Interval(1)]
#     preim_x = RigorousInvarireantMeasures.preimages(x, Ttrue.branches[2])
#     @info preim_x
#     err_x = maximum(IntervalArithmetic.radius.(preim_x[1]))
#     y_s = NTuple{2, Interval}[]
#     err_y = 0.0
#     for x in preim_x[1][2:end]
#         G_inv = Ginverse(x = x, r = 5.0, c = 0.5)
#         preim_y_low = min(max(G_inv(y_lower), -1), 1)
#         preim_y_up = max(min(G_inv(y_upper), 1), -1)
#         err_y = maximum([err_y, IntervalArithmetic.radius(preim_y_low), IntervalArithmetic.radius(preim_y_up)])
#         push!(y_s, (preim_y_low, preim_y_up))
#     end
#     n = length(preim_x[1][2:end])
#     A = Matrix{Float64}(undef, 2*n, 2)
#     for (i, x) in enumerate(preim_x[1][2:end])
#         A[i, :] = [mid(x) mid(y_s[i][1])]
#     end

#     for (i, x) in enumerate(reverse(preim_x[1][2:end]))
#         A[n+i, :] = [mid(x) mid(y_s[end-i+1][2])]
#     end

#     P = PH.polyhedron(PH.vrep(A), lib)
#     return P, err_x, err_y
# end

function PreimageRectangleLorenz(; preim_x_left, preim_x_right, y_lower, y_upper, k, r, c)
    x_grid = [
        Interval(preim_x_left) + i * Interval(preim_x_right - preim_x_left) / k for i = 0:k
    ]
    err_x = maximum(IntervalArithmetic.radius.(x))
    y_s = NTuple{2,Interval}[]
    err_y = 0.0
    for x in x_grid
        G_inv = Ginverse(x = x, r = r, c = c)
        preim_y_low = min(max(G_inv(y_lower), -1), 1)
        preim_y_up = max(min(G_inv(y_upper), 1), -1)
        err_y = maximum([
            err_y,
            IntervalArithmetic.radius(preim_y_low),
            IntervalArithmetic.radius(preim_y_up),
        ])
        push!(y_s, (preim_y_low, preim_y_up))
    end
    n = length(x_grid)
    A = Matrix{Float64}(undef, 2 * n, 2)
    for (i, x) in enumerate(x_grid)
        A[i, :] = [mid(x) mid(y_s[i][1])]
    end

    for (i, x) in enumerate(reverse(x_grid))
        A[n+i, :] = [mid(x) mid(y_s[end-i+1][2])]
    end

    P = PH.polyhedron(PH.vrep(A), lib)
    return P, err_x, err_y
end


Gfloat(x, y; r = 5.0, c = 0.5) = x > 0 ? 2^(-r) * y * x^r + c : 2^(-r) * y * (-x)^r - c

F(v) = [Tfloat(v[1]); Gfloat(v[1], v[2])]


dx = 0.01:0.01:1
dy = -1:0.01:1

image = Array{Float64}(undef, length(dx) * length(dy), 2)


for (i, (x, y)) in enumerate(Base.Iterators.product(dx, dy))
    v = [x; y]
    image[i, :] = F(F(F(v)))'
end

using IntervalArithmetic
function Tfloat(x::Interval; α = 1.75, s = 3.0)
    xr = x ∩ @interval 0 1
    @debug xr
    yr = α * (xr)^s - 1
    @debug yr
    xl = x ∩ @interval -1 0
    @debug xl
    yl = -α * (-xl)^s + 1
    @debug yl

    return yl ∪ yr
end

F(v::IntervalBox) = Tfloat(v[1]) × Gfloat(v[1], v[2])

image_ib = Array{IntervalBox}(undef, length(dx) * length(dy))

noise_amplitude = 1 / 16
noise_1d = -noise_amplitude .. noise_amplitude
noise = noise_1d × noise_1d

for (i, (x, y)) in enumerate(Base.Iterators.product(dx, dy))
    v = (x + Interval(0, 0.01)) × (y + Interval(0, 0.01))
    image_ib[i] = F(v) + noise_1d
end

plt = plot(image_ib[1])
for i = 2:20000
    plt = plot!(image_ib[i], color = :blue, label = "")
end

plt = plot!(image_ib[1001])

A = fill(false, (4096, 4096))
