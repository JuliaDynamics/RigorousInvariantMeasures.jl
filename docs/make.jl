using RigorousInvariantMeasures
using Documenter

DocMeta.setdocmeta!(
    RigorousInvariantMeasures,
    :DocTestSetup,
    :(using RigorousInvariantMeasures);
    recursive = true,
)

makedocs(;
    modules = [RigorousInvariantMeasures],
    authors = "Isaia Nisoli <nisoli@im.ufrj.br> and contributors",
    sitename = "RigorousInvariantMeasures.jl",
    format = Documenter.HTML(;
        canonical = "https://juliadynamics.github.io/RigorousInvariantMeasures.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        # "Background" => "background.md",
        # "User guide" => [
        #     "General usage" => "userguide.md",
        #     "Implementing a new basis" => "implementingnewbasis.md",
        # ],
        "Basis" => "Basis.md",
        "Dynamic" => "Dynamic.md",
        "Generic assembler interface" => "GenericAssembler.md",
        "Generic estimate interface" => "GenericEstimate.md",      
        "Norms Of Powers" => "NormsOfPowers.md",
        # "Examples" => "examples.md",
        # "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaDynamics/RigorousInvariantMeasures.jl",
    #versions = ["stable" => "v^", "v#.#", devurl => devurl],
    devbranch = "main",
)
