using Documenter, RigorousInvariantMeasures

makedocs(sitename="RigorousInvariantMeasures.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [RigorousInvariantMeasures]
)

