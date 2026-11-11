# Pendulum

The pendulum

```math
\begin{cases}
  \dot \theta = \omega, \\
  \dot \omega = -\sin \theta,
\end{cases}
```

has centers at ``\theta = 2k\pi``, surrounded by oscillations, and saddles at ``\theta = (2k + 1)\pi``, connected by heteroclinic orbits, the separatrices, which separate the oscillations from the rotations of the pendulum over the top.

```@example pendulum
using LagrangianDescriptors, OrdinaryDiffEq, CairoMakie
using LinearAlgebra: norm

function pendulum!(du, u, p, t)
    θ, ω = u
    du[1] = ω
    du[2] = -sin(θ)
end

prob = ODEProblem(pendulum!, [0.0, 0.0], (0.0, 10.0))

M(du, u, p, t) = norm(du)
uu0 = [[θ, ω] for ω in range(-3, 3, length = 201), θ in range(-2π, 2π, length = 401)]

lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())

fig, ax, plt = heatmap(lagsol; axis = (xlabel = L"\theta", ylabel = L"\omega"))
Colorbar(fig[1, 2], plt; label = "Lagrangian descriptor")
fig
```

The separatrices appear as sharp curves through the saddles at ``\theta = \pm\pi``. Inside them, the oscillations have small descriptors near the centers, where the solutions move slowly. Outside them, the rotations, which move faster, have larger descriptors.
