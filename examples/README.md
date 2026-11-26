# Examples

The examples have their own environment, `examples/Project.toml`, which develops the
package from the parent directory and adds the plotting and I/O packages the scripts
load. From the repository root:

```
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/<file>.jl
```

Scripts that plot should be run headless on a remote machine, with `GKSwstype=100` in
the environment.
