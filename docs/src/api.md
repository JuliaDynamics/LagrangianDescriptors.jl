# API
## Problems

```@docs
LagrangianDescriptors.augmentprob
```

```@docs
LagrangianDescriptors.LagrangianDescriptorProblem
```

## Solutions

```@docs
LagrangianDescriptors.LagrangianDescriptorSolution
```

```@docs
LagrangianDescriptors.solve
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

## Descriptors

```@docs
LagrangianDescriptors.IncrementDescriptor
```

## Post-processing

```@docs
LagrangianDescriptors.lagrangian_descriptor
```