module LagrangianDescriptors

using SciMLBase, ComponentArrays
import SciMLBase:solve

export augmentprob

include("augmented_problems.jl")
include("lagrangiandescriptors_problems.jl")
include("plot_recipes.jl")

end # module
