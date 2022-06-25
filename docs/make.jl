using Documenter, RigorousInvariantMeasures

makedocs(sitename="RigorousInvariantMeasures.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [RigorousInvariantMeasures]
)

deploydocs(
    repo = "github.com/USER_NAME/PACKAGE_NAME.jl.git",
)