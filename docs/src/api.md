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

A [`LagrangianDescriptorSolution`](@ref) can be plotted with [Plots.jl](https://github.com/JuliaPlots/Plots.jl), through a plot recipe, or with [Makie.jl](https://docs.makie.org), through a package extension that is loaded together with Makie:

```julia
plot(lagsol)                       # Plots or Makie
plot(lagsol, :forward)             # a given direction
plot(lagsol; idxs = (1, 3))        # components 1 and 3 of the initial conditions
heatmap(lagsol, :backward)         # Makie: any compatible plot function
plot(lagsol; seriestype = :contourf) # Plots: any series type
```

The `direction` can be `:forward`, `:backward`, `:both` or `:total` (the sum of both, the default for a problem with both directions), or `:difference` (forward minus backward). When solved with several `realizations`, the mean over the realizations is plotted. The kind of plot depends on the collection `uu0` of initial conditions:

- numbers (scalar problems), or vectors with a single component: a line with the Lagrangian descriptor at each initial condition;
- a matrix of vectors forming a grid in the components `idxs`, in either layout, `[[x, y] for y in ys, x in xs]` or `[[x, y] for x in xs, y in ys]`: a heatmap;
- any other collection of vectors (scattered points, a curve, a matrix not forming a grid, ...): a scatter plot of the components `idxs` of the initial conditions, colored by the Lagrangian descriptor.

The keyword `idxs`, with default `(1, 2)`, chooses the two components of the initial conditions plotted, e.g. for a two-dimensional slice of initial conditions of a higher-dimensional system.

With Makie, a colorbar is added as usual, e.g. `fig, ax, plt = plot(lagsol); Colorbar(fig[1, 2], plt)`. If both Plots and Makie are loaded, their plotting functions must be qualified, e.g. `Makie.plot(lagsol)`.

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
