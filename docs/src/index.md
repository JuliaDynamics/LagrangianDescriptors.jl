# LagrangianDescriptors.jl documentation

## Overview

The dynamics of evolutionary systems can be quite intricate. The method of **Lagrangian Descriptors** comes to help visualize the complicate behavior of such systems.

In a recent article, Wiggins and and García-Garrido call it *painting the phase portrait* of a dynamical system (see [S. Wiggins and V. J. García-Garrido, Painting the Phase Portrait of a Dynamical System with the Computational Tool of Lagrangian Descriptors (AMS Notices, June/July 2022)](https://www.ams.org/journals/notices/202206/noti2489/noti2489.html?adat=June/July%202022&trk=2489&galt=none&cat=feature&pdfissue=202206&pdffile=rnoti-p936.pdf).

The method is akin to droping colored ink in a fluid flow and seeing the dye being transported to reveal the flow patterns, except that the color don't get diffused and an even clearer image appears, better revealing hidden structures.

The next image, for instance, shows the dynamics of a periodically-forced Duffing equation

```math
\ddot x = x - x^3 + A\sin(\omega t).
```

![Duffing](img/duffing.png)

The image above can be obtained by first setting up the equation as an `ODEProblem` from `SciML/DifferentialEquations`, then wrapping that as a `LagrangianDescriptorProblem` from `LagrangianDescriptors.jl`, and finally solving it and plotting the result:

```@example
using Plots
x = range(0.0, 8π, length=201)
plot(x, sin)
```

```julia
using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using LagrangianDescriptors
```

```julia
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

```julia
M(du, u, p, t) = norm(du)

uu0 = [[x, y] for y in range(-1.0, 1.0, length=301), x in range(-1.8, 1.8, length=301)]

lagprob = LagrangianDescriptorProblem(prob, M, uu0)
```

```julia
@time lagsol = solve(lagprob, Tsit5())
```

```julia
plot(lagsol, title="Lagrangian descriptors for the forced Duffing equation \$\\ddot x = x - x^3 + A\\sin(\\omega t)\$\nwith A=$A and ω=$ω", titlefont=8, xlabel="\$x\$", ylabel="\$\\dot x\$")

#savefig("img/blah.png"); nothing # hide
```

`![](img/blah.png)`

```@example
1+4
```



## Developers

LagrangianDescriptors is currently being developed by [Ricardo M. S. Rosa](https://rmsrosa.github.io), but contributors are welcome.

## Cite
