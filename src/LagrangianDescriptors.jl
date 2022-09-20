module LagrangianDescriptors

using SciMLBase, ComponentArrays, RecipesBase
import DiffEqBase: solve
export LagrangianDescriptorProblem, LagrangianDescriptorSolution

include("augmented_problems.jl")
include("lagrangiandescriptor_problems.jl")
include("lagrangian_descriptor_solution.jl")
include("solve.jl")
include("plot_recipes.jl")

end # module
