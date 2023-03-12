#
# This file computes (and saves to a file labelled with the current time) timings and accuracy results for various
# matrix assembly and dfly operations. It is useful to assess if changes introduced during the development cause
# major regressions in terms of CPU time or accuracy.
#
# Expect a run time should be about a couple of minutes, since the dfly for Lorenz3 is *slow*.
#

using RigorousInvariantMeasures
using IntervalArithmetic
using Dates

using Statistics

function get_experiment(prefix)
    if prefix=="Lorenz3"
        f = n-> begin
            D0 = Lorenz()
            D = D0 ∘ D0 ∘ D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q    
        end
    elseif prefix=="Lorenz2"
        f = n-> begin
            D0 = Lorenz()
            D = D0 ∘ D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q    
        end
    elseif prefix=="4x_perturbed_Ulam"
        f = n-> begin
            D = mod1_dynamic(x -> 4x + RigorousInvariantMeasures.sinpi(8x)/100)
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="4x_perturbed_Hat"
        f = n-> begin
            D = mod1_dynamic(x -> 4x + RigorousInvariantMeasures.sinpi(8x)/100)
            B = Hat(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="Lanford_Ulam"
        f = n-> begin
            D = mod1_dynamic(x -> 2x+0.5*x*(1-x))
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="Lanford3_Hat"
        f = n-> begin
            D0 = mod1_dynamic(x->2x+0.5*x*(1-x))
            # Taking an iterate is necessary here to get a DFLY inequality with A < 1
            D = D0∘D0∘D0
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end    
    elseif prefix=="PiecewiseLinear2"
        f = n-> begin
            D = PwMap([x->17x/5, 
                x->(34*((17x-5)/17)/25+3)*((17x-5)/17), 
                x->(34*((17x-10)/17)/25+3)*((17x-10)/17), 
                x->17*((17x-15)/17)/5], 
                [Interval(0), Interval(5)/17, Interval(10)/17, Interval(15)/17, Interval(1)],
                [Interval(0) Interval(1);
                Interval(0) Interval(1);
                Interval(0) Interval(1);
                Interval(0) @interval(0.4)]
                )
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="PiecewiseLinear"
        f = n-> begin
            D = PwMap([x->2.5x, x->4x-1, x->4x-2, x-> 4x-3],
                [@interval(0), @interval(0.25), @interval(0.5), @interval(0.75), @interval(1)])
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="175"
        f = n-> begin
            D = mod1_dynamic(x -> 17//5 * x)
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="175_nonlinear"
        f = n-> begin
        D = PwMap(
            [x -> 17x/5,
             x -> 34(x-5//17)^2/25 + 3(x-5//17),
             x -> 34(x-10//17)^2/25 + 3(x-10//17),
             x -> 17(x-15//17)/5
            ],
            [0, @interval(5/17), @interval(10/17), @interval(15/17), 1],
            [0 1; 0 1; 0 1; 0 @interval(0.4)]
            )
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    elseif prefix=="BZ"
        f = n-> begin
            D = BZ()
            B = Ulam(n)
            Q = DiscretizedOperator(B, D)
            return B, D, Q
        end
    else
        @error "Unknown experiment"
    end
end

experiments = [
    "Lanford_Ulam", 
    "Lanford3_Hat",
    "4x_perturbed_Ulam", 
    "4x_perturbed_Hat", 
    "PiecewiseLinear", 
    "PiecewiseLinear2", 
    "175", 
    "175_nonlinear", 
    "Lorenz2", 
    "Lorenz3", 
    "BZ"
    ]

include("warmup.jl")

T = fill(NaN,length(experiments), 6)
for (index, prefix) in enumerate(experiments)
    f = get_experiment(prefix)
    n = 64
    T[index,1] = @elapsed B, D, Q = f(n)
    T[index,2] = maximum(diam.(Q.L))
    T[index,3] = Statistics.mean(diam.(Q.L))
    try
    T[index,4] = @elapsed dfly_coefficients = dfly(strong_norm(B), aux_norm(B), D)
    T[index,5] = dfly_coefficients[1]
    T[index,6] = dfly_coefficients[2]
    catch e
        # DFLY failed, let's leave NaNs in the table
    end
end

headers = ["Experiment" "Ass.time" "max diam" "mean diam" "dfly time" "dfly A" "dfly B"]

display([headers; experiments T])
print("Commit: $(read(`git rev-parse HEAD`,String)[1:7]); Machine: $(Sys.cpu_info()[1].model); Date: $(now()).\n")

filename = "regressions_$(now()).txt"

open(filename, "w") do io
    show(io, "text/plain", [headers; experiments T])
    show(io, "text/plain", "\n\nCommit: $(read(`git rev-parse HEAD`,String)[1:7]); Machine: $(Sys.cpu_info()[1].model); Date: $(now())")
end

print("Saved in '$(filename)'.")
