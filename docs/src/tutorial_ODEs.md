# Ordinary Differential Equations

We considere, here, some examples of applying the Lagragian descriptor method to equations of the type `ODEProblem`.

## Periodically-forced Duffing equation

We start with an application of the method of Lagrangian descriptors to the periodically-forced Duffing equation

```math
\ddot x = x - x^3 + A\sin(\omega t).
```

With `LagrangianDescriptors.jl`, we start by setting up the equation as an `ODEProblem` from [SciML/DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Next, we wrap that as a `LagrangianDescriptorProblem` from `LagrangianDescriptors.jl`, which we can solve as an ensemble problem.

Finally, we plot the result with the built-in plot recipe.

So we first load the relevant packages:

```julia duffing
using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using LagrangianDescriptors
```

Next we set up the `ODEProblem`:

```julia duffing
function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * cos(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 13.0)
A = 0.3; ω = π; p = (A, ω)

prob = ODEProblem(f!, u0, tspan, p)
```

With the ODE problem setup, we choose an infinitesimal Lagrangian descriptor, a collection of initial conditions on the phase space, which is the region to be "painted", and finally we build the `LagrangianDescriptorProblem`:

```julia duffing
M(du, u, p, t) = norm(du)

uu0 = [[x, y] for y in range(-1.0, 1.0, length=301), x in range(-1.8, 1.8, length=301)]

lagprob = LagrangianDescriptorProblem(prob, M, uu0)
```

The Lagrangian descriptors are the time-integration of the infinitesimal descriptor along forward and backward solutions of the equation. They are integrated along with the solutions by "solving" the `LagrangianDescriptorProblem`, with an overload of the `solve` method from the [SciML](https://sciml.ai) ecosystem:

```julia duffing
lagsol = solve(lagprob, Tsit5())
```

With the solution at hand, we plot the Lagrangian descriptors to visualize the dynamics of the system:

```julia duffing
plot(lagsol, title="Lagrangian descriptors for the forced Duffing equation \$\\ddot x = x - x^3 + A\\sin(\\omega t)\$\nwith A=$A and ω=$ω", titlefont=8, xlabel="\$x\$", ylabel="\$\\dot x\$")

savefig("img/duffing.png")
```

![Duffing](img/duffing.png)

We may zoom closer to the origin to find the following "painting":

```julia duffing
uu0 = [[x, y] for y in range(-0.6, 0.3, length=601), x in range(-0.2, 0.4, length=401)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)

lagsol = solve(lagprob, Tsit5());

plot(lagsol, title="Lagrangian descriptors for the forced Duffing equation \$\\ddot x = x - x^3 + A\\sin(\\omega t)\$\nwith A=$A and ω=$ω", titlefont=8, xlabel="\$x\$", ylabel="\$\\dot x\$")

savefig("img/duffing2.png")
```

![Duffing](img/duffing2.png)

If we want to change parameters, we just `remake` the original `ODEProblem` (in the future I should add the option to remake the `LagrangianDescriptorProblem` itself.)

```julia
A = 5.0; ω = 2π; p = (A, ω);
prob = remake(prob, p=p)
uu0 = [[x, y] for y in range(-0.5, 0.2, length=501), x in range(-0.3, 0.2, length=401)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)

lagsol = solve(lagprob, Tsit5());

plot(lagsol, title="Lagrangian descriptors for the forced Duffing equation \$\\ddot x = x - x^3 + A\\sin(\\omega t)\$\nwith A=$A and ω=$ω", titlefont=8, xlabel="\$x\$", ylabel="\$\\dot x\$")

savefig("img/duffing3.png")
```

![Duffing](img/duffing3.png)