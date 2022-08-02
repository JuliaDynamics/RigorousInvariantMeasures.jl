using RigorousInvariantMeasures
using IntervalArithmetic

# Testing whether the "shooting method" beats taking successive preimages

# f = x -> 2*x + 0.1*RigorousInvariantMeasures.sinpi(2*x)
# preim1 = root(x ->f(x)-1., 0..4, 0)

f = x-> 2.5x*(1-x) + 0.01*x*x + 0.3*x*x*x + 0.003*x*x*x*x
preim1 = @interval(0.5)


mypreim(y) = root(x -> f(x)-y, hull(0..0,preim1), 0)

n = 15

y0 = @interval(0.21)

# method 1
ys = fill(y0, n+1)
ys[1] = y0
for i = 1:n
    ys[i+1] = mypreim(ys[i])
end

# method 2

zs = fill(hull(0..0,preim1), n)
fs = fill(f, n)

nthpreimage!(y0, fs, zs)

z = zs[1]

[zs[end:-1:1] ys[2:end]]
[diam.(zs[end:-1:1]) diam.(ys[2:end])]

diam.(zs[end:-1:1]) ./ diam.(ys[2:end])