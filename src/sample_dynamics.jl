using JLD2

export BZ, Lorenz

function BZ()
    BZcoeff = load("./examples/bzmapcoeff.jld2")
    BZA::Interval{Float64} = Interval{Float64}(BZcoeff["a"])
    BZB::Interval{Float64} = Interval{Float64}(BZcoeff["b"])
    BZC::Interval{Float64} = Interval{Float64}(BZcoeff["c"]) 

    T_left_leq_1_8(x) = (BZA - (Interval(1)/8-x)^(Interval(1)/3))*exp(-x) + BZB
    T_left_geq_1_8(x) = (BZA + (x-Interval(1)/8)^(Interval(1)/3))*exp(-x) + BZB
    T_right(x) =BZC*(10x*exp(-Interval(10)/3*x))^(19)+BZB

    return PwMap([T_left_leq_1_8,
    T_left_geq_1_8,
    T_right],
	[Interval(0), Interval(1)/8, Interval(3)/10, Interval(1)],
    [T_left_leq_1_8(Interval(0)) T_left_leq_1_8(Interval(1)/8);
    T_left_geq_1_8(Interval(1)/8) T_left_geq_1_8(Interval(3)/10);
    T_right(Interval(3)/10) T_right(Interval(1))
    ]
    )
end

Lorenz(θ=109/64, α=51/64) = PwMap([x->θ*(0.5-x)^α, x->1-θ*(x-0.5)^α],
                    [@interval(0), @interval(0.5), @interval(1)],
                    [θ*(Interval(0.5))^α Interval(0.0);
                    Interval(1.0)  1-θ*(Interval(0.5))^α]; infinite_derivative=true)

