# Lorenz system

The Lorenz system

```math
\begin{cases}
  \dot x = \sigma (y - x), \\
  \dot y = x (\rho - z) - y, \\
  \dot z = x y - \beta z,
\end{cases}
```

with ``\sigma = 10``, ``\rho = 28`` and ``\beta = 8/3``, has a chaotic attractor. Its phase space is three-dimensional, so we compute the Lagrangian descriptors on a two-dimensional slice of initial conditions, here the plane ``z = 27``, through the nontrivial fixed points ``(\pm\sqrt{72}, \pm\sqrt{72}, 27)``. The plot recipe recognizes the grid in the components that vary, ``x`` and ``y``.

Backward in time, the volumes of phase space expand and the solutions diverge, so we only compute the forward descriptors, which reveal the stable manifolds.

```@example lorenz
using LagrangianDescriptors, OrdinaryDiffEq, CairoMakie
using LinearAlgebra: norm

function lorenz!(du, u, p, t)
    σ, ρ, β = p
    x, y, z = u
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
end

prob = ODEProblem(lorenz!, zeros(3), (0.0, 1.5), (10.0, 28.0, 8 / 3))

M(du, u, p, t) = norm(du)
uu0 = [[x, y, 27.0] for y in range(-25, 25, length = 201), x in range(-20, 20, length = 201)]

lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0; direction = :forward), Tsit5())

fig, ax, plt = heatmap(lagsol; axis = (title = "slice z = 27", xlabel = "x", ylabel = "y"))
Colorbar(fig[1, 2], plt; label = "forward Lagrangian descriptor")
fig
```

The dark regions around the fixed points are the solutions that stay close to them, and the sharp curves are the intersections of the stable manifold of the origin, which separates the two wings of the attractor, with the plane. With a longer time window, more of its folds appear, but they also become thinner.

Other slices are obtained with other collections of initial conditions, and the keyword `idxs` chooses the components to plot. For instance, on the plane ``y = 0``:

```@example lorenz
uu0_xz = [[x, 0.0, z] for z in range(0, 50, length = 201), x in range(-20, 20, length = 201)]

lagsol_xz = solve(LagrangianDescriptorProblem(prob, M, uu0_xz; direction = :forward), Tsit5())

fig, ax, plt = heatmap(lagsol_xz; idxs = (1, 3), axis = (title = "slice y = 0", xlabel = "x", ylabel = "z"))
Colorbar(fig[1, 2], plt; label = "forward Lagrangian descriptor")
fig
```

The two slices can also be seen together in three dimensions, with three components in `idxs`, and a common `colorrange`, so that their colors are comparable. CairoMakie, however, draws each plot over the previous ones, without hiding the parts that should be behind, so we split each plane at the line where they intersect, here through the middle of both, and draw the four half-planes from the farthest to the nearest, for the default view of `Axis3`. (GLMakie and WGLMakie hide the parts behind by themselves.)

```@example lorenz
solve_forward(uu0) = solve(LagrangianDescriptorProblem(prob, M, uu0; direction = :forward), Tsit5())

xs = range(-20, 20, length = 201)
# the plane z = 27, split at y = 0, and the plane y = 0, split at z = 27
back_of_z27 = solve_forward([[x, y, 27.0] for y in range(0, 27, length = 101), x in xs])
front_of_z27 = solve_forward([[x, y, 27.0] for y in range(-27, 0, length = 101), x in xs])
top_of_y0 = solve_forward([[x, 0.0, z] for z in range(27, 54, length = 101), x in xs])
bottom_of_y0 = solve_forward([[x, 0.0, z] for z in range(0, 27, length = 101), x in xs])

# from the farthest to the nearest
halves = (back_of_z27, bottom_of_y0, top_of_y0, front_of_z27)
colorrange = extrema(reduce(vcat, [lagsol() for lagsol in halves]))

fig = Figure(size = (700, 550))
ax = Axis3(fig[1, 1]; title = "slices z = 27 and y = 0", xlabel = "x", ylabel = "y", zlabel = "z", viewmode = :fit)
plts = [plot!(ax, lagsol; idxs = (1, 2, 3), colorrange) for lagsol in halves]
Colorbar(fig[1, 2], plts[1]; label = "forward Lagrangian descriptor")
fig
```

Any of the half-planes can be hidden, e.g. with `plts[4].visible = false` for the front one, to see the ones behind it. The axis then only spans the visible ones, unless its limits are given, e.g. with `limits = ((-20, 20), (-27, 27), (0, 54))` in `Axis3`.
