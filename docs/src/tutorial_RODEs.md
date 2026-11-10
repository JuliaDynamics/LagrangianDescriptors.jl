# Random Ordinary Differential Equations

We consider, now, some examples of applying the Lagrangian descriptor method to equations of the type `RODEProblem`.

## Random Duffing equation

We consider, first, the random Duffing equation

```math
\ddot x = x - x^3 + a\sin(W_t),
```

where ``\{W_t\}_t`` is a standard Wiener process.

With `LagrangianDescriptors.jl`, we i) start by setting up the equation as a `RODEProblem`; next, we ii) wrap that as a `LagrangianDescriptorProblem`; then we iii) `solve` it; and finally, we iv) plot the result with the built-in plot recipe. Since the noise path is given in the problem, the same path is used for all initial conditions.

So we first load the relevant packages:

```@example rode-duffing
using StochasticDiffEq, DiffEqNoiseProcess
using Plots
using Random: seed!
using LinearAlgebra: norm
using LagrangianDescriptors
```

Next we set up the `RODEProblem`:

```@example rode-duffing
function f!(du, u, p, t, W)
    x, y = u
    du[1] = y
    du[2] = x - x^3 + p * sin(W)
end

p = 0.025
u0 = [0.5, 2.2]
tspan = (0.0, 8.0)
seed!(123) # for a reproducible noise path
W = WienerProcess(0.0, 0.0)

prob = RODEProblem(f!, u0, tspan, p, noise=W)
```

With the RODE problem setup, we choose an infinitesimal Lagrangian descriptor, a collection of initial conditions on the phase space, which is the region to be "painted", and finally we build the `LagrangianDescriptorProblem`:

```@example rode-duffing
M(du, u, p, t, W) = norm(du)

uu0 = [[x, y] for y in range(-1.0, 1.0, length=301), x in range(-1.8, 1.8, length=301)]

lagprob = LagrangianDescriptorProblem(prob, M, uu0)
```

The Lagrangian descriptors are the time-integration of the infinitesimal descriptor along forward and backward solutions of the equation. They are integrated along with the solutions by "solving" the `LagrangianDescriptorProblem`, with an overload of the `solve` method from the [SciML](https://sciml.ai) ecosystem:

```@example rode-duffing
lagsol = solve(lagprob, RandomHeun(), dt=1/100)
```

With the solution at hand, we plot the Lagrangian descriptors to visualize the dynamics of the system:

```@example rode-duffing
plot(lagsol, title="Lagrangian descriptors for the random Duffing equation \$\\ddot x = x - x^3 + a\\sin(W_t)\$\nwith \$a=\$$p", titlefont=8, xlabel="\$x\$", ylabel="\$\\dot x\$")
savefig("duffing_random.png"); nothing # hide
```

![Lagrangian descriptors of the random Duffing equation](duffing_random.png)
