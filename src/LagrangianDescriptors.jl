module LagrangianDescriptors

using SciMLBase, ComponentArrays
import SciMLBase:solve

include("augmented_problems.jl")
include("lagrangiandescriptors_problems.jl")
include("plot_recipes.jl")

end # module
