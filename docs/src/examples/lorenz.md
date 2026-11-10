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
using LagrangianDescriptors, OrdinaryDiffEq, Plots
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

plot(lagsol, xlabel = "x", ylabel = "y", title = "slice z = 27", titlefont = 10)
savefig("lorenz_xy.png"); nothing # hide
```

![Forward Lagrangian descriptors of the Lorenz system on the plane z = 27](lorenz_xy.png)

The dark regions around the fixed points are the solutions that stay close to them, and the sharp curves are the intersections of the stable manifold of the origin, which separates the two wings of the attractor, with the plane. With a longer time window, more of its folds appear, but they also become thinner.

Other slices are obtained with other collections of initial conditions, and the keyword `idxs` chooses the components to plot. For instance, on the plane ``y = 0``:

```@example lorenz
uu0_xz = [[x, 0.0, z] for z in range(0, 50, length = 201), x in range(-20, 20, length = 201)]

lagsol_xz = solve(LagrangianDescriptorProblem(prob, M, uu0_xz; direction = :forward), Tsit5())

plot(lagsol_xz; idxs = (1, 3), xlabel = "x", ylabel = "z", title = "slice y = 0", titlefont = 10)
savefig("lorenz_xz.png"); nothing # hide
```

![Forward Lagrangian descriptors of the Lorenz system on the plane y = 0](lorenz_xz.png)
