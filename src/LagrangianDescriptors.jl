module LagrangianDescriptors

using SciMLBase, DiffEqNoiseProcess
using ComponentArrays, RecipesBase
using QuadGK: quadgk
import DiffEqBase: solve
export LagrangianDescriptorProblem, LagrangianDescriptorSolution, lagrangian_descriptor

include("augmented_problems.jl")
include("lagrangiandescriptor_problems.jl")
include("lagrangiandescriptor_solution.jl")
include("solve.jl")
include("lagrangian_descriptor_postprocessing.jl")
include("plot_recipes.jl")

end # module
