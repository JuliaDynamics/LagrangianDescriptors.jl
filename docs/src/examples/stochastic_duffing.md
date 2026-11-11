# Stochastic Duffing equation

The Duffing equation with additive noise in the velocity,

```math
\begin{cases}
  \mathrm{d}x = y \,\mathrm{d}t, \\
  \mathrm{d}y = (x - x^3) \,\mathrm{d}t + \sigma \,\mathrm{d}W_t,
\end{cases}
```

where ``W_t`` is a Wiener process, is a stochastic perturbation of the Duffing equation of the page on [choosing the descriptor and the time window](../choosing.md). We set it up as an `SDEProblem`, with diagonal noise, which is zero in the first component:

```@example stochastic_duffing
using LagrangianDescriptors, StochasticDiffEq, Plots
using LinearAlgebra: norm

drift(u, p, t) = [u[2], u[1] - u[1]^3]
diffusion(u, σ, t) = [0.0, σ]

M(du, u, p, t) = norm(du) # with `du` the drift
uu0 = [[x, y] for y in range(-1, 1, length = 101), x in range(-1.8, 1.8, length = 101)]
nothing # hide
```

Following [BalibreaIniestaEtAl2016](@citet), the Lagrangian descriptors of a single realization of the noise are computed with the same noise path for all initial conditions, with the keyword `shared_noise = true` of [`solve`](@ref LagrangianDescriptors.solve). By default, each initial condition has independent noise. The `seed` of the `SDEProblem` makes the computation reproducible:

```@example stochastic_duffing
prob = SDEProblem(drift, diffusion, [0.0, 0.0], (0.0, 5.0), 0.1; seed = 123)
lagprob = LagrangianDescriptorProblem(prob, M, uu0)
lagsol_path = solve(lagprob, EM(), dt = 0.01; shared_noise = true)
nothing # hide
```

The pictures of different realizations differ, and the mean over realizations shows the structures they share. With the keyword `realizations`, the problem is solved for several realizations, each with its own noise path, shared by all initial conditions, and the solution gives the mean of the descriptors:

```@example stochastic_duffing
lagsol_mean = solve(lagprob, EM(), dt = 0.01; shared_noise = true, realizations = 10)

plot(
    plot(lagsol_path, title = "a single realization"),
    plot(lagsol_mean, title = "mean over 10 realizations"),
    layout = (1, 2), size = (900, 360), margin = 4Plots.mm, titlefont = 10,
    colorbar = false, xlabel = "x", ylabel = "y",
)
savefig("stochastic_duffing.png"); nothing # hide
```

![Lagrangian descriptors of the stochastic Duffing equation](stochastic_duffing.png)

For a single realization, the structures of the deterministic equation, its saddle point and the homoclinic loops, are deformed by the noise, and they change from one realization to another. Their mean shows the structures shared by the realizations.

Without `shared_noise = true`, each initial condition has independent noise, also in each realization. The mean over realizations approximates the same mean at each initial condition, but, with few realizations, the picture is grainier, since the noise changes from one initial condition to the next.

The descriptor of [BalibreaIniestaEtAl2016](@citet), given by the increments of the solutions, is available as [`IncrementDescriptor`](@ref).
