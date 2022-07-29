using RigorousInvariantMeasures
using Documenter

DocMeta.setdocmeta!(RigorousInvariantMeasures, :DocTestSetup, :(using RigorousInvariantMeasures); recursive=true)

makedocs(;
    modules=[RigorousInvariantMeasures],
    authors="Isaia Nisoli nisoli@im.ufrj.br and contributors",
    repo="https://github.com/JuliaDynamics/RigorousInvariantMeasures.jl/blob/{commit}{path}#{line}",
    sitename="RigorousInvariantMeasures.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/JuliaDynamics/RigorousInvariantMeasures.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Background" => "background.md",
        "User guide" => "userguide.md",
        "Basis" => "basis.md",
        "Dynamic" => "dynamic.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaDynamics/RigorousInvariantMeasures.jl",
    devbranch="master",
)