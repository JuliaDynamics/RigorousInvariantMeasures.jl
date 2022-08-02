<<<<<<< HEAD
# Basis

These are some of then basis implemented in the package.

## Ulam Basis
```@autodocs
Modules = [RigorousInvariantMeasures]
Pages =["UlamBasis.jl"]
```


## Hat basis
```@docs
Hat
Hat(::Integer)
Base.length(B::Hat{T}) where {T}
```
### Hat functions types
```@docs
HatFunctionOnTorus
=======
We present some of the basis already implemented in the package

# Abstract basis
```@autodocs
Modules = [Base, 
        RigorousInvariantMeasures, 
        RigorousInvariantMeasures.BasisDefinition]
Pages = ["BasisDefinition.jl"]
```

# The Ulam basis
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures, 
            RigorousInvariantMeasures.BasisDefinition]
Pages = ["UlamBasis.jl"]
```

# The hat basis on ``S^1``
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures, 
            RigorousInvariantMeasures.BasisDefinition]
Pages = ["HatBasis.jl"]
```

# The hat basis on ``[0,1]``
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures, 
            RigorousInvariantMeasures.BasisDefinition]
Pages = ["NonPeriodicHatBasis.jl"]
```

# The C2 basis
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures, 
            RigorousInvariantMeasures.BasisDefinition,
            RigorousInvariantMeasures.C2BasisDefinition]
Pages = ["C2Basis.jl"]
```

# The Chebyshev basis
```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures, 
            RigorousInvariantMeasures.BasisDefinition]
Pages = ["NewChebyshev.jl"]
>>>>>>> origin/master
```
