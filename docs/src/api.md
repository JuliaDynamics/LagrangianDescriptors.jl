# API

The typical use of the package consists of three steps:

```julia
lagprob = LagrangianDescriptorProblem(prob, M, uu0)  # 1. set up the problem
lagsol = solve(lagprob, alg)                        # 2. solve it
plot(lagsol)                                        # 3. visualize the descriptors
```

where `prob` is a differential equation problem, `M` is an infinitesimal Lagrangian descriptor, `uu0` is a collection of initial conditions, and `alg` is an algorithm to solve `prob`. The Lagrangian descriptors at each initial condition are then given by `lagsol()`, or `lagsol(:forward)`, `lagsol(:backward)`, and so on.

## Problems

```@docs
LagrangianDescriptorProblem
```

## Lagrangian descriptors

The Lagrangian descriptors are given either by the integral of an infinitesimal descriptor `M` along the solutions, or by the sum of the increments of the solutions over a time grid, with an [`IncrementDescriptor`](@ref).

The infinitesimal descriptor `M` is any function with the arguments below, where `du` is the value of the right-hand side, or drift, of the equation at the state `u`, with parameters `p`, at time `t`:

| Problem        | Infinitesimal descriptor | `du`                      |
|:-------------- |:------------------------ |:------------------------- |
| `ODEProblem`   | `M(du, u, p, t)`         | `f(u, p, t)`              |
| `RODEProblem`  | `M(du, u, p, t, W)`      | `f(u, p, t, W)`, with the noise `W` |
| `SDEProblem`   | `M(du, u, p, t)`         | the drift `f(u, p, t)`    |

A common choice is the norm of the velocity, `M(du, u, p, t) = norm(du)`, for which the Lagrangian descriptor is the length of the solution over the time interval.

```@docs
IncrementDescriptor
```

## Solving

```@docs
LagrangianDescriptors.solve
```

## Solutions

```@docs
LagrangianDescriptorSolution
```

## Plotting

A [`LagrangianDescriptorSolution`](@ref) can be plotted with Makie or Plots, e.g. with `plot(lagsol)` or `plot(lagsol, direction)`; see [Plotting](@ref plotting).

## Lower-level functions

The augmented problems solved for each initial condition, and the computation of the Lagrangian descriptors of a single solution of the original problem, by post-processing (see [Alternative via post-processing](@ref)).

```@docs
LagrangianDescriptors.augmentprob
lagrangian_descriptor
```

## Index

```@index
Pages = ["api.md"]
```
