# [Plotting](@id plotting)

A [`LagrangianDescriptorSolution`](@ref) can be plotted with [Makie.jl](https://docs.makie.org), through a package extension that is loaded together with Makie, or with [Plots.jl](https://github.com/JuliaPlots/Plots.jl), through a plot recipe. The rest of this documentation uses Makie, as the other packages of [JuliaDynamics](https://github.com/JuliaDynamics) do, but both give the same kinds of plots, with the same arguments.

## What is plotted

```julia
plot(lagsol)                         # the default direction
plot(lagsol, :forward)               # a given direction
plot(lagsol; idxs = (1, 3))          # components 1 and 3 of the initial conditions
```

The `direction` can be `:forward`, `:backward`, `:both` or `:total` (the sum of both, the default for a problem with both directions), or `:difference` (forward minus backward). When solved with several `realizations`, the mean over the realizations is plotted. The kind of plot depends on the collection `uu0` of initial conditions:

- numbers (scalar problems), or vectors with a single component: a line with the Lagrangian descriptor at each initial condition;
- a matrix of vectors forming a grid in the components `idxs`, in either layout, `[[x, y] for y in ys, x in xs]` or `[[x, y] for x in xs, y in ys]`: a heatmap;
- any other collection of vectors (scattered points, a curve, a matrix not forming a grid, ...): a scatter plot of the components `idxs` of the initial conditions, colored by the Lagrangian descriptor.

The keyword `idxs`, with default `(1, 2)`, chooses the two components of the initial conditions plotted, e.g. for a two-dimensional slice of initial conditions of a higher-dimensional system.

The colormap is `:inferno` by default, with both Makie and Plots.

As an example, we use the Lagrangian descriptors of the Duffing equation ``\ddot x = x - x^3``, on a grid of initial conditions.

## With Makie

```@example makie
using OrdinaryDiffEq, CairoMakie
using LinearAlgebra: norm
using LagrangianDescriptors

function f!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = x - x^3
end

prob = ODEProblem(f!, [0.0, 0.0], (0.0, 5.0))
M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1.0, 1.0, length = 201), x in range(-1.6, 1.6, length = 321)]
lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())

fig, ax, plt = heatmap(lagsol; axis = (xlabel = L"x", ylabel = L"\dot{x}"))
Colorbar(fig[1, 2], plt; label = "Lagrangian descriptor")
fig
```

Any compatible plot function can be used, with the usual Makie attributes, e.g. a given direction, with a different colormap, and filled contours:

```@example makie
fig = Figure(size = (1000, 400))
ax = Axis(fig[1, 1]; title = "forward", xlabel = L"x", ylabel = L"\dot{x}")
plt = heatmap!(ax, lagsol, :forward; colormap = :viridis)
Colorbar(fig[1, 2], plt)
ax = Axis(fig[1, 3]; title = "backward", xlabel = L"x", ylabel = L"\dot{x}")
plt = contourf!(ax, lagsol, :backward; levels = 20, colormap = :plasma)
Colorbar(fig[1, 4], plt)
fig
```

The default colormap of the package takes precedence over a colormap set with Makie's `set_theme!`, so another colormap is given with `colormap = ...` in each call.

## With Plots

```@example plots
using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using LagrangianDescriptors

function f!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = x - x^3
end

prob = ODEProblem(f!, [0.0, 0.0], (0.0, 5.0))
M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1.0, 1.0, length = 201), x in range(-1.6, 1.6, length = 321)]
lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())

plot(lagsol; xlabel = "x", ylabel = "ẋ")
```

Any series type and the usual Plots attributes can be used, e.g. a given direction, with a different colormap, and filled contours:

```@example plots
plot(
    plot(lagsol, :forward; title = "forward", color = :viridis),
    plot(lagsol, :backward; title = "backward", seriestype = :contourf, levels = 20, color = :plasma);
    layout = (1, 2), size = (1000, 400), margin = 5Plots.mm, xlabel = "x", ylabel = "ẋ",
)
```

With Plots, the colormap can also be changed for all plots with `Plots.default(color = :viridis)`.

## With both

If both Plots and Makie are loaded, their plotting functions must be qualified, e.g. `Makie.heatmap(lagsol)` or `Plots.plot(lagsol)`.
