# Implementing a new basis

In the file [BasisDefinition.jl] are declared the necessary methods 
for the implementation of a new basis to work with the coarse-fine framework.

The new basis type has to be defined as a subtype of the
abstract type `Basis`.

```julia
struct NewBasis <: Basis
    # base specific code
end
```

If possible, implement some utility constructors like 
`NewBasis(n)` to initialize a basis of size $n$.

To compute the operator through the generic assembler
we implement a `NewBasisDual <: Dual` object that contains the data 
necessary to compute the raw data necessary 
to compute the entries of the operator, e.g.,
in the Ulam basis case, it contains three objects.
One is a vector of preimages, which contains 
all the preimages ``T_i^{-1}(x_j)`` of the elements of the partition 
through the different branches ``T_i`` of the map,
a second vector that contains the labels associated to these points, i.e., ``j``.
The layout is like this:

``
[T_1^{-1}(x_1) \ldots T_1^{-1}(x_n) T_2^{-1}(x_1) \ldots ]
``

``
[1 \ldots n 1 \ldots]
``

In the case of the `Hat` basis, there is a third vector 
that contains the derivative at the preimages.

Now it is necessary to define an iterator on the `Dual`
object, by defining 
`Base.iterate(dual::NewBasisDual, state = 1)`.
This iterator returns the data necessary for the computation; in the case of the Ulam basis, it returns
the label and the preimages of the endpoints 
of one connected component of the preimage of an interval.

We take now the elements returned by the iteration 
on `NewBasisDual` and for each one of them we define 
construct an object of type `ProjectDualElement`.
This object contains information on which elements
of the basis we need to compute corresponding to the dual element, computed through the `BasisDefinition.nonzero_on` function.
In the Ulam case, it computes which intervals have nonempty intersection with the interval defined in 
the dual element.

Now, we need to define an iterator on the `ProjectDualElement` object that for all indexes
associated to the dual element computes the 
value of the coefficient of the operator.

Resuming:
- define a `NewBasis` structure
- define a `NewDual` structure and the iterator on it
- define a `BasisDefinition.nonzero_on` function which allows the construction of `ProjectDualElement` objects
- define the iterator on `ProjectDualElement`that returns 
the coefficients

To allow the coarse-fine method to work out of the box,
we need to implement different functions that return
bound for the various constants.
A comprehensive list is contained in the file [BasisDefinition.jl].

The first thing we need to define are the norm types
that correspond to our approximation scheme.
```julia 
struct StrongNormNewBasis <: NormKind end
struct WeakNormNewBasis <: NormKind end
```
They must be specialized to our new basis.

```julia
strong_norm(B::NewBasis) = StrongNormNewBasis
weak_projection_error(B::NewBasis)
```

The needed constants and functions are listed there.
If not defined, the compiler will call the most 
general version that applies, the one in BasisDefinition,
that throws an error.

In principle, it is possible to define a new basis by 
starting with only the struct, running the code and looking at the error codes julia is throwing, that correspond to the non implemented functions.


