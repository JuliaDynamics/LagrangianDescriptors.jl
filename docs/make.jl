using Documenter
using LagrangianDescriptors

makedocs(
    sitename = "LagrangianDescriptors.jl",
    pages = [
        "Home" => "index.md"
        "Overview" => "overview.md"
        "Problems" => "problems.md"
        "Solutions" => "solutions.md"
        "API" => "api.md"
    ],
    authors = "Ricardo Rosa",
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/JuliaDynamics/LagrangianDescriptors.jl",
        edit_link="main"
    ),
    modules = [LagrangianDescriptors]
)

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/LagrangianDescriptors.jl.git",
        devbranch = "main"
    )
end

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#