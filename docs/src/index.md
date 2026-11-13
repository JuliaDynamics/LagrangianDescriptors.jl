# LagrangianDescriptors.jl

```@docs
LagrangianDescriptors
```

## Getting started

- [The method of Lagrangian descriptors](@ref method) explains the method and how the package computes the descriptors.
- [Choosing the descriptor and the time window](choosing.md) shows the effect of the choices that shape the pictures, with some guidance on them.
- The tutorials on [ordinary](tutorial_ODEs.md), [random](tutorial_RODEs.md) and [stochastic](tutorial_SDEs.md) differential equations go through the use of the package step by step, and [Alternative via post-processing](postprocessing.md) compares it with computing the descriptors from solutions of the original equations.
- The [Examples](examples/index.md) show the method on classical systems, from simple ones to chaotic and stochastic ones.
- The [API](api.md) documents all functions and types.
- The [References](references.md) list the literature on the method.

## A first example

Here is the code that computes the Lagrangian descriptors in the image above, of a periodically forced Duffing equation, revealing some of the structures of its dynamics.

```@example index
using OrdinaryDiffEq, CairoMakie
using LinearAlgebra: norm
using LagrangianDescriptors

function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * sin(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 13.0)
A = 5.0; ω = 2π; p = (A, ω);
prob = ODEProblem(f!, u0, tspan, p)

M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1.5, -0.5, length=301), x in range(-0.4, 1.6, length=601)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)

lagsol = solve(lagprob, Tsit5())

fig, ax, plt = heatmap(lagsol; axis = (xlabel = L"x", ylabel = L"\dot{x}"))
Colorbar(fig[1, 2], plt; label = "Lagrangian descriptor")
fig
```

## Developers

[JuliaDynamics/LagrangianDescriptors.jl](https://github.com/JuliaDynamics/LagrangianDescriptors.jl) is currently being developed by [Ricardo M. S. Rosa](https://rmsrosa.github.io).

## Cite

Just cite the github repo [JuliaDynamics/LagrangianDescriptors.jl](https://github.com/JuliaDynamics/LagrangianDescriptors.jl) for now, while the package is not yet registered.
