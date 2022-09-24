## Current state

Handling of problems of the type `ODEProblem` is essentially done, including some plot recipes for 1D and 2D problems.

What we have so far:

1. A method `augmentprob(prob::ODEProblem, M; direction=:both)` which takes an `ODEProblem` and a local Lagrangian descriptor `M=M(du, u, p, t)` and creates another `ODEProblem` for an augmented system which contains four components: the original forward ODE problem, the associated backward ODE problem, and the equations for the forward and backward Lagrangian descriptors. If `direction=:forward`, then the augmented system only contains the forward equation and the forward Lagrangian descriptor. If `direction=:backard`, then the augmented system only contains the backward equation and the backward Lagrangian descriptor. The given `ODEProblem` can be in-place or out-of-place, but the augmented system is always in-place, and the components are build using `ComponentArrays`, with componentes `fwd` (forward ODE), `bwd` (backward ODE), `lfwd` (forward Lagrangian descriptor), and `lbwd` (backward Lagrangian descriptor).
1. A type with an inner method constructor `LagrangianDescriptorProblem(prob, M, uu0; direction::Symbol=:both)` that takes the original problem `prob`, the local Lagrangian descriptor `M`, a collection of initial conditions `uu0` (e.g. a vector, array, range or anything with a `length` and accessible by index `uu0[i]`), and the optional keyword `direction`. The arguments `prob`, `M` and `direction` are simply passed to `augmentprob(prob, M; direction)` for the creation of an augmented problem `augprob`. Then, an `EnsembleProblem` is created with `augprob` and including a `prob_func` that is supposed to sweep the initial conditions in `uu0` and an `output_func` that is supposed to only save the last element (at time `final(tspan)`) of the foward and backward Lagrangian descriptors. They are then used to create an instance `lagprob` of the `LagrangianDescriptorProblem` that contains two fields, `lagprob.ensprob` being the associated ensemble problem and `lagprob.uu0` being the collection of initial conditions.
1. A dispatch of `solve(lagprob::LagrangianDescriptorProblem, alg, args...; kwargs...)` that solves the ensemble problem `lagprob.ensprob` with `trajectories=length(uu0)` and whatever algoritm and arguments/keyword arguments is necessary or desired. This creates an instance of `LagrangianDescriptorSolution`, with three fields, one being the solution of the associated ensemble problem and the other two being the collection of initial conditions and the direction.
1. Finally, there is a plot recipe for plotting a `LagrangianDescriptorSolution`, provided the array of initial conditions is either an one-dimensional `AbstractVector{<:Number}` or a two-dimensional `AbstractMatrix{<:AbstractVector{<:Number}}`, with the eltype `AbstractVector` being of length two.

## Roadmap

What is currently missing:
1. Support for other types of problems, e.g. `SDEProblem`, `RODEProblem`, mixed differential-algebraic equations, etc.;
1. A more flexible plot recipe;
1. Improve the documentation with more examples;
1. Maybe an adaptive method to refine the set `uu0` of initial conditions!;
1. I don't know whether/how the idea applies to delay type equations, but we should check that out.
1. Register the package.
