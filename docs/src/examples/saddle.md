# Saddle point

The linear system

```math
\begin{cases}
  \dot x = x, \\
  \dot y = -y,
\end{cases}
```

has a saddle point at the origin, with the ``y``-axis as its stable manifold and the ``x``-axis as its unstable manifold. Since the solutions are known explicitly, ``x(t) = x_0 e^{t}`` and ``y(t) = y_0 e^{-t}``, this is a good first check of the method.

```@example saddle
using LagrangianDescriptors, OrdinaryDiffEq, CairoMakie
using LinearAlgebra: norm

saddle(u, p, t) = [u[1], -u[2]]
prob = ODEProblem(saddle, [0.0, 0.0], (0.0, 2.0))

M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1, 1, length = 201), x in range(-1, 1, length = 201)]

lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())

fig = Figure(size = (1200, 400))
for (i, direction) in enumerate((:forward, :backward, :total))
    ax = Axis(fig[1, i]; title = string(direction), xlabel = "x", ylabel = "y")
    heatmap!(ax, lagsol, direction)
end
fig
```

The forward descriptor is smallest along the stable manifold ``x = 0``, where the solutions approach the saddle, and grows with ``|x_0|``, since the solutions leave the saddle exponentially fast. The backward descriptor, likewise, reveals the unstable manifold ``y = 0``, and their sum shows both.
