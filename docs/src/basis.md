We present some of the basis already implemented in the package

# The Ulam basis
```@docs
Ulam
nonzero_on(B::Ulam, (a, b))
relative_measure
iterate(S::AverageZero{Ulam}, state = 1)
```


# The hat basis on ``S^1``
```@docs
Hat
Hat(n::Integer)
HatFunctionOnTorus
Base.getindex(B::Hat, i::Int)
BasisDefinition.nonzero_on(B::Hat, dual_element)
Base.iterate(S::AverageZero{Hat}, state)
```

# The hat basis on ``[0,1]``