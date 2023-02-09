using JLD2, IntervalArithmetic
BZcoeff = load("./examples/bzmapcoeff.jld2")

using RigorousInvariantMeasures

BZA::Interval{Float64} = Interval{Float64}(BZcoeff["a"])
BZB::Interval{Float64} = Interval{Float64}(BZcoeff["b"])
BZC::Interval{Float64} = Interval{Float64}(BZcoeff["c"]) 


T_left_leq_1_8(x) = (BZA - (Interval(1)/8-x)^(Interval(1)/3))*exp(-x) + BZB
T_left_geq_1_8(x) = (BZA + (x-Interval(1)/8)^(Interval(1)/3))*exp(-x) + BZB
T_right(x) =BZC*(10*x*exp(-Interval(10)/3*x))^(19)+BZB

### Since I'm lazy and I am sure I'm going to make a mess computing derivative, 
### I'm going to use the symbolic package

#using Symbolics 
#@variables x
#D = Differential(x)

#expr_T_l_leq_1_8 = (BZA-(Interval(1)/8-x)^(Interval(1)/3))*exp(-x)+BZB
#der_expr_T_l_leq_1_8 = expand_derivatives(D(expr_T_l_leq_1_8))
#der_T_l_leq_1_8_func = Symbolics.build_function(der_expr_T_l_leq_1_8, x, expression = Val(false))

#expr_T_l_geq_1_8 = (BZA+(x-Interval(1)/8)^(Interval(1)/3))*exp(-x)+BZB
#der_expr_T_l_geq_1_8 = expand_derivatives(D(expr_T_l_geq_1_8))
#der_T_l_geq_1_8_func = Symbolics.build_function(der_expr_T_l_geq_1_8, x, expression = Val(false))

#expr_T_right = C*(Interval(10)*x*exp(-Interval(10)/3*x))^(19)+B
#der_expr_T_right = expand_derivatives(D(expr_T_right))
#der_T_right_func = Symbolics.build_function(der_expr_T_right, x, expression = Val(false))

der_T_l_leq_1_8_func(x) = -exp(-x)*(6*BZA*(1-8*x)^(Interval(2)/3)+24*x-11)/(6*(1-8*x)^(Interval(2)/3))
der_T_l_geq_1_8_func(x) = -exp(-x)*(6*BZA*(8*x-1)^(Interval(2)/3)+24*x-11)/(6*(8*x-1)^(Interval(2)/3))
der_T_right_func(x) = (19*BZC*(10*x*exp(x*Interval(-10)/3))^18)*(10*exp(x*Interval(-10)/3)+10*x*exp(x*Interval(-10)/3)*(Interval(-10)/3))

BZmap() = PwMap([T_left_leq_1_8,
    T_left_geq_1_8,
    T_right],
    [der_T_l_leq_1_8_func,
    der_T_l_geq_1_8_func,
    der_T_right_func
    ],    
	[Interval(0), Interval(1)/8, Interval(3)/10, Interval(1)],
    [T_left_leq_1_8(Interval(0)) T_left_leq_1_8(Interval(1)/8);
    T_left_geq_1_8(Interval(1)/8) T_left_geq_1_8(Interval(3)/10);
    T_right(Interval(3)/10) T_right(Interval(1))
    ]
    )

D = BZmap()
B = Ulam(2^10)
@time Q = DiscretizedOperator(B, D)
