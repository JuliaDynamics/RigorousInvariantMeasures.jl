using IntervalArithmetic
using RigorousInvariantMeasures


D = BZ()
B = Ulam(2^10)
@time Q = DiscretizedOperator(B, D)
