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

Following [BalibreaIniestaEtAl2016](@citet), the Lagrangian descriptors of a single realization are computed with the same noise path for all initial conditions. This is obtained by fixing the `seed` of the `SDEProblem`, so that every solve uses the same noise:

```@example stochastic_duffing
prob = SDEProblem(drift, diffusion, [0.0, 0.0], (0.0, 5.0), 0.1; seed = 123)
lagsol_path = solve(LagrangianDescriptorProblem(prob, M, uu0), EM(), dt = 0.01)
nothing # hide
```

The pictures of different realizations differ, and their mean over realizations shows the structures they share [BalibreaIniestaEtAl2016](@cite). For pictures, it is best to average realizations that each have a single noise path for all initial conditions, here given by the seeds `1` to `10`:

```@example stochastic_duffing
xs = range(-1.8, 1.8, length = 101)
ys = range(-1, 1, length = 101)

realizations = 10
Lmean = sum(1:realizations) do seed
    lagsol = solve(LagrangianDescriptorProblem(remake(prob; seed), M, uu0), EM(), dt = 0.01)
    lagsol()
end / realizations

plot(
    plot(lagsol_path, title = "a single realization"),
    heatmap(xs, ys, reshape(Lmean, size(uu0)), title = "mean over $realizations realizations"),
    layout = (1, 2), size = (900, 360), margin = 4Plots.mm, titlefont = 10,
    colorbar = false, xlabel = "x", ylabel = "y",
)
savefig("stochastic_duffing.png"); nothing # hide
```

![Lagrangian descriptors of the stochastic Duffing equation](stochastic_duffing.png)

For a single realization, the structures of the deterministic equation, its saddle point and the homoclinic loops, are deformed by the noise, and they change from one realization to another. Their mean shows the structures shared by the realizations.

Alternatively, without a fixed seed, the keyword `realizations` of [`solve`](@ref LagrangianDescriptors.solve) solves each initial condition several times, each with independent noise, and gives the mean of the descriptors. It approximates the same mean at each initial condition, but since the noise is independent from one initial condition to the next, the picture is grainier with few realizations.

The descriptor of [BalibreaIniestaEtAl2016](@citet), given by the increments of the solutions, is available as [`IncrementDescriptor`](@ref).
