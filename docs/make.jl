using Documenter, RigorousInvariantMeasures

makedocs(sitename="RigorousInvariantMeasures.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages=[
        "Home" => "index.md",
        "Background" => "background.md",
        "User guide" => "userguide.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ],
)

deploydocs(
    repo = "github.com/JuliaDynamics/RigorousInvariantMeasures.jl.git",
)