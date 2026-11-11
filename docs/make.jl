cd(@__DIR__)

using LagrangianDescriptors

pages = [
    "Overview" => "index.md",
    "Choosing the descriptor and the time window" => "choosing.md",
    "Tutorials" => [
        "ODEs" => "tutorial_ODEs.md",
        "RODEs" => "tutorial_RODEs.md",
        "Alternative via post-processing" => "postprocessing.md",
    ],
    "Examples" => [
        "Overview" => "examples/index.md",
        "Saddle point" => "examples/saddle.md",
        "Pendulum" => "examples/pendulum.md",
        "Van der Pol oscillator" => "examples/vanderpol.md",
        "Double gyre" => "examples/doublegyre.md",
        "Lorenz system" => "examples/lorenz.md",
        "Double pendulum" => "examples/doublependulum.md",
        "Stochastic Duffing equation" => "examples/stochastic_duffing.md",
    ],
    "API" => "api.md",
    "References" => "references.md",
]

# The documentation style of JuliaDynamics: https://github.com/JuliaDynamics/doctheme
import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

build_docs_with_style(
    pages, LagrangianDescriptors;
    bib,
    authors = "Ricardo Rosa",
    htmlkw = (
        canonical = "https://juliadynamics.github.io/LagrangianDescriptors.jl/stable/",
        edit_link = "main",
    ),
)
