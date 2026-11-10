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
using LagrangianDescriptors, OrdinaryDiffEq, Plots
using LinearAlgebra: norm

saddle(u, p, t) = [u[1], -u[2]]
prob = ODEProblem(saddle, [0.0, 0.0], (0.0, 2.0))

M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1, 1, length = 201), x in range(-1, 1, length = 201)]

lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())

plot(
    plot(lagsol, :forward, title = "forward"),
    plot(lagsol, :backward, title = "backward"),
    plot(lagsol, title = "total"),
    layout = (1, 3), size = (1200, 360), margin = 4Plots.mm, titlefont = 10,
    colorbar = false, xlabel = "x", ylabel = "y",
)
savefig("saddle.png"); nothing # hide
```

![Lagrangian descriptors of a saddle point](saddle.png)

The forward descriptor is smallest along the stable manifold ``x = 0``, where the solutions approach the saddle, and grows with ``|x_0|``, since the solutions leave the saddle exponentially fast. The backward descriptor, likewise, reveals the unstable manifold ``y = 0``, and their sum shows both.
