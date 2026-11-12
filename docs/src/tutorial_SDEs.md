# Stochastic Differential Equations

We consider, now, an example of applying the Lagrangian descriptor method to equations of the type `SDEProblem`.

## Duffing equation with multiplicative noise

We consider the Duffing equation with a multiplicative noise in the velocity, proportional to the cube root of the position,

```math
\begin{cases}
  \mathrm{d}x = y \,\mathrm{d}t, \\
  \mathrm{d}y = (x - x^3) \,\mathrm{d}t + \sigma x^{1/3} \,\mathrm{d}W_t,
\end{cases}
```

where ``\{W_t\}_t`` is a standard Wiener process. The noise vanishes on the line ``x = 0``, so the origin is still an equilibrium, a saddle point, as in the deterministic equation, while the other structures are perturbed by the noise.

With `LagrangianDescriptors.jl`, we i) start by setting up the equation as an `SDEProblem`; next, we ii) wrap that as a `LagrangianDescriptorProblem`; then we iii) `solve` it; and finally, we iv) plot the result with the built-in plot recipe.

So we first load the relevant packages:

```@example sde-duffing
using StochasticDiffEq
using CairoMakie
using LinearAlgebra: norm
using LagrangianDescriptors
```

Next we set up the `SDEProblem`, with the drift `f!` and the diffusion `g!`, here a diagonal noise, which is zero in the first component, and a `seed` for reproducible results:

```@example sde-duffing
function f!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = x - x^3
end

function g!(du, u, σ, t)
    du[1] = 0.0
    du[2] = σ * cbrt(u[1]) # the real cube root, also for negative values
end

σ = 0.05
u0 = [0.0, 0.0]
tspan = (0.0, 5.0)

prob = SDEProblem(f!, g!, u0, tspan, σ; seed = 123)
```

With the SDE problem setup, we choose an infinitesimal Lagrangian descriptor, a collection of initial conditions on the phase space, which is the region to be "painted", and finally we build the `LagrangianDescriptorProblem`. For an `SDEProblem`, the infinitesimal descriptor has the signature `M(du, u, p, t)`, where `du` is the value of the drift:

```@example sde-duffing
M(du, u, p, t) = norm(du)

uu0 = [[x, y] for y in range(-1.0, 1.0, length = 101), x in range(-1.8, 1.8, length = 101)]

lagprob = LagrangianDescriptorProblem(prob, M, uu0)
```

The Lagrangian descriptors are integrated along with the forward and backward solutions by "solving" the `LagrangianDescriptorProblem`, with any SDE algorithm of [StochasticDiffEq.jl](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/), here the adaptive `SOSRI()`, for diagonal noise:

```@example sde-duffing
lagsol = solve(lagprob, SOSRI())
```

By default, each initial condition is solved with its own realization of the noise, and the forward and backward directions are driven by independent noise, as for a two-sided Wiener process. We plot the forward and the backward Lagrangian descriptors, and their sum:

```@example sde-duffing
fig = Figure(size = (1200, 400))
for (i, direction) in enumerate((:forward, :backward, :total))
    ax = Axis(fig[1, i]; title = string(direction), xlabel = L"x", ylabel = L"y")
    heatmap!(ax, lagsol, direction)
end
fig
```

The forward descriptors show the stable manifold of the saddle point, and the backward ones its unstable manifold, as for the deterministic equation, while the noise, independent from one initial condition to the next, makes the pictures slightly grainy.

## A single realization of the noise

Following [BalibreaIniestaEtAl2016](@citet), the Lagrangian descriptors of a single realization of the noise are computed with the same noise path for all initial conditions, with the keyword `shared_noise = true`. With the keyword `realizations`, the problem is solved for several realizations of the noise, and the solution gives the mean of the descriptors over them:

```@example sde-duffing
lagsol_path = solve(lagprob, SOSRI(); shared_noise = true)
lagsol_mean = solve(lagprob, SOSRI(); shared_noise = true, realizations = 10)

fig = Figure(size = (1000, 400))
ax = Axis(fig[1, 1]; title = "a single realization", xlabel = L"x", ylabel = L"y")
heatmap!(ax, lagsol_path)
ax = Axis(fig[1, 2]; title = "mean over 10 realizations", xlabel = L"x", ylabel = L"y")
heatmap!(ax, lagsol_mean)
fig
```

For a single realization, the saddle point stays at the origin, where the noise vanishes, while the homoclinic loops are slightly deformed by the noise, differently in each realization. Their mean shows the structures shared by the realizations.

See also the [Stochastic Duffing equation](@ref) example, with additive noise, and the descriptor of [BalibreaIniestaEtAl2016](@citet), given by the increments of the solutions, available as [`IncrementDescriptor`](@ref).
