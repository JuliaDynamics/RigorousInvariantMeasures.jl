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
        prettyurls = true,
    ),
    pages = [
        "Home" => "index.md",
        "Background" => "background.md",
        "User guide" => [
            "General usage" => "userguide.md",
            "Implementing a new basis" => "implementingnewbasis.md",
        ],
        "Interface" => [
            "Basis" => "Basis.md",
            "Dynamic" => "Dynamic.md",
            "Generic assembler interface" => "GenericAssembler.md",
            "Generic estimate interface" => "GenericEstimate.md",
            "Norms Of Powers" => "NormsOfPowers.md",
            "Noise Kernel" => "NoiseKernel.md",
            "Observables" => "Observables.md",
            "Accessories" => "api.md",
        ],
        "Examples" => "examples.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaDynamics/RigorousInvariantMeasures.jl",
    #versions = ["stable" => "v^", "v#.#", devurl => devurl],
    devbranch = "main",
)
