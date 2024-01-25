# COV_EXCL_START

function val(x, a, b)
    return 0.5 * (a - b) * (a + b - 2 * x)
end

function g(x, v, B)
    w = [v[i] * val(x, B.p[i], B.p[i+1]) for i = 1:length(B)]
    #assumes that the Ulam partition is homogeneous
    return sum(w)
end

h(x, v, δ, B) = x + δ * g(x, v, B)

D_P = mod1_dynamic(x -> h(x, v, 0.5, B))

# COV_EXCL_STOP
