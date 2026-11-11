module LagrangianDescriptors

using SciMLBase
using ComponentArrays, RecipesBase
using QuadGK: quadgk
import SciMLBase: solve
export LagrangianDescriptorProblem, LagrangianDescriptorSolution, lagrangian_descriptor
export IncrementDescriptor

include("augmented_problems.jl")
include("lagrangiandescriptor_problems.jl")
include("increment_descriptor.jl")
include("lagrangiandescriptor_solution.jl")
include("solve.jl")
include("lagrangian_descriptor_postprocessing.jl")
include("plot_data.jl")
include("plot_recipes.jl")

end # module
