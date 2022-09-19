module LagrangianDescriptors

using SciMLBase, ComponentArrays
import DiffEqBase: solve, init
export LagrangianDescriptorProblem

include("augmented_problems.jl")
include("lagrangiandescriptors_problems.jl")
include("solve.jl")
include("plot_recipes.jl")

end # module
