# Van der Pol oscillator

The van der Pol oscillator

```math
\begin{cases}
  \dot x = y, \\
  \dot y = \mu (1 - x^2) y - x,
\end{cases}
```

with ``\mu > 0``, has an unstable focus at the origin and an attracting limit cycle. Forward in time, all solutions approach the limit cycle. Backward in time, the solutions inside the limit cycle approach the origin, while those outside grow without bound.

The solutions that diverge backward in time need to be stopped, with the keyword `max_trajectory_value` of [`solve`](@ref LagrangianDescriptors.solve). When a problem has both directions, a solution is stopped as soon as either its forward or its backward part exceeds the given value, which also stops the other part. So here we solve the two directions separately, and only stop the backward solutions:

```@example vanderpol
using LagrangianDescriptors, OrdinaryDiffEq, Plots
using LinearAlgebra: norm

function vanderpol!(du, u, μ, t)
    x, y = u
    du[1] = y
    du[2] = μ * (1 - x^2) * y - x
end

prob = ODEProblem(vanderpol!, [0.0, 0.0], (0.0, 10.0), 1.0)

M(du, u, p, t) = norm(du)
xs = range(-3, 3, length = 201)
ys = range(-4, 4, length = 201)
uu0 = [[x, y] for y in ys, x in xs]

forward = solve(LagrangianDescriptorProblem(prob, M, uu0; direction = :forward), Tsit5())
backward = solve(LagrangianDescriptorProblem(prob, M, uu0; direction = :backward), Tsit5();
    max_trajectory_value = 10.0)

count(backward.terminated)
```

The descriptors of the stopped solutions only measure how fast they diverge, so we leave them out of the backward picture, setting them to `NaN`:

```@example vanderpol
Lbwd = backward(:backward)
Lbwd[backward.terminated] .= NaN

plot(
    plot(forward, title = "forward"),
    heatmap(xs, ys, reshape(Lbwd, size(uu0)), title = "backward, without the diverging solutions"),
    layout = (1, 2), size = (900, 360), margin = 4Plots.mm, titlefont = 10,
    colorbar = false, xlabel = "x", ylabel = "y",
)
savefig("vanderpol.png"); nothing # hide
```

![Lagrangian descriptors of the van der Pol oscillator](vanderpol.png)

Forward in time, all solutions approach the limit cycle, and the descriptor only shows the unstable focus at the origin. Backward in time, the limit cycle is repelling, and it appears as the boundary of the region of solutions that remain bounded.
