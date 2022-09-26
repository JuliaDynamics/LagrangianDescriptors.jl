# Lagrangian Descriptors

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliadynamics.github.io/LagrangianDescriptors.jl/dev/) ![Main Tests Workflow Status](https://github.com/JuliaDynamics/LagrangianDescriptors.jl/workflows/CI/badge.svg) [![codecov](https://codecov.io/gh/JuliaDynamics/LagrangianDescriptors.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/LagrangianDescriptors.jl) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![GitHub repo size](https://img.shields.io/github/repo-size/JuliaDynamics/LagrangianDescriptors.jl)

Implementation of the method of **Lagrangian Descriptors** to highlight singular features (e.g. stable or unstable invariant manifolds) of the dynamics of an evolutionary system (such as ordinary or partial differential equations, random equations, or stochastic differential equations).

## Example

Here in an example on a periodically forced Duffing system.

```julia
using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using LagrangianDescriptors

function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * cos(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 13.0)
A = 0.3; ω = π

p = (A, ω)
prob = ODEProblem(f!, u0, tspan, p)

M(du, u, p, t) = norm(du)

uu0 = [[x, y] for y in range(-1.0, 1.0, length=301), x in range(-1.8, 1.8, length=301)]

lagprob = LagrangianDescriptorProblem(prob, M, uu0)
```

With all setup, we may solve it as follows:

```
julia> @time lagsol = solve(lagprob, Tsit5())
  4.241800 seconds (204.88 M allocations: 9.173 GiB, 23.40% gc time, 24.99% compilation time)
```

Then we use the built-in plot recipe to get the heatmap of the Lagrangian descriptors:

```julia
plot(lagsol)
```

![Duffing Lagrangian descriptor](docs/src/img/duffing.png)

## References

* [Painting the Phase Portrait of a Dynamical System with the Computational Tool of Lagrangian Descriptors](https://www.ams.org/journals/notices/202206/noti2489/noti2489.html?adat=June/July%202022&trk=2489&galt=none&cat=feature&pdfissue=202206&pdffile=rnoti-p936.pdf)
* [Lagrangian descriptors: A method for revealing phase space structures of general time dependent dynamical systems](https://www.sciencedirect.com/science/article/abs/pii/S1007570413002037)
* [Lagrangian Descriptors - *Discovery and Quantification of Phase Space Structure and Transport*](https://champsproject.github.io/lagrangian_descriptors/docs/authors.html)
* [Frequently Asked Questions about Lagrangian Descriptors](https://acp.copernicus.org/preprints/acp-2016-633/acp-2016-633-SC2-supplement.pdf)

