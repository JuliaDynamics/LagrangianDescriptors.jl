using Documenter
using LagrangianDescriptors

makedocs(
    sitename = "LagrangianDescriptors",
    format = Documenter.HTML(),
    modules = [LagrangianDescriptors]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
