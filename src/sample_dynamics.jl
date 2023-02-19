export BZ, Lorenz

function BZ()
    A_big = Interval{BigFloat}(BigFloat("0.5060735690368223513195993710530479569801417368282037493809901142182256388277772"))
    BZA = Interval{Float64}(A_big)
    B_big = Interval{BigFloat}(BigFloat("0.02328852830307032296813220750095076307514120085284790788241725176646123060202265"), BigFloat("0.02328852830307032296813220750095076307514120085284790788241849636951680283042265"))
    BZB = Interval{Float64}(B_big)
    C_big = Interval{BigFloat}(BigFloat("0.121205692738975111744666848150620569782497212127938371936404761693002104361654"))
    BZC = Interval{Float64}(C_big)

    T_left_leq_1_8(x) = (BZA - (Interval(1)/8-x)^(1//3))*exp(-x) + BZB
    T_left_geq_1_8(x) = (BZA + (x-Interval(1)/8)^(1//3))*exp(-x) + BZB
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
                    Interval(1.0)  1-θ*(Interval(0.5))^α])

