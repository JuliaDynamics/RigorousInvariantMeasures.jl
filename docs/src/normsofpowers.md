# Computing norms of powers

One of the main tools in our computer aided proofs is 
obtaining an estimate on the distance in the weak norm
from the ``abstract'' invariant density $f$ and the 
approximated invariant density $f_k$ in the weak norm.
We refer to

Galatolo S., Monge M., Nisoli I., Poloni F.
A general framework for the rigorous computation of invariant densities and the coarse-fine strategy
[arXiv:2212.05017]


To obtain this estimate we need to get upper bounds
on 
$$
||L^n_k|_{U_0}||_w.
$$

The main estimate is 
$$
||f-f_k||_w \leq \sum_{i=0}^n ||L^n_k|_{U_0}||_w ||L-L_k f||_{s\to w}$
$$

```@autodocs
Modules = [Base, 
            RigorousInvariantMeasures]
Pages = ["NormsOfPowers.jl"]
```