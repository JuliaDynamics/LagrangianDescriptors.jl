cd(@__DIR__)

# Draw the Plots figures without opening windows (GR's headless mode), unless set
# otherwise; this must be set before the first plot
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using LagrangianDescriptors

pages = [
    "Overview" => "index.md",
    "The method of Lagrangian descriptors" => "method.md",
    "Choosing the descriptor and the time window" => "choosing.md",
    "Tutorials" => [
        "ODEs" => "tutorial_ODEs.md",
        "RODEs" => "tutorial_RODEs.md",
        "Alternative via post-processing" => "postprocessing.md",
    ],
    "Plotting" => "plotting.md",
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

bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"); style = :authoryear)

build_docs_with_style(
    pages, LagrangianDescriptors;
    bib,
    # the stable documentation only exists after the first release
    linkcheck_ignore = [r"^https://juliadynamics\.github\.io/LagrangianDescriptors\.jl/stable/?$"],
    # the style's default timeout of 2 seconds is too short for some sites on the CI
    linkcheck_timeout = 10,
    authors = "Ricardo Rosa",
    htmlkw = (
        canonical = "https://juliadynamics.github.io/LagrangianDescriptors.jl/stable/",
        edit_link = "main",
    ),
)
