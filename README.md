# LagrangianDescriptors.jl

[![docsdev](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://juliadynamics.github.io/LagrangianDescriptors.jl/dev/)
[![docsstable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/LagrangianDescriptors.jl/stable/)
[![CI](https://github.com/JuliaDynamics/LagrangianDescriptors.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/LagrangianDescriptors.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/LagrangianDescriptors.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JuliaDynamics/LagrangianDescriptors.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/JuliaDynamics/LagrangianDescriptors.jl/blob/main/LICENSE)

![Lagrangian descriptors of a periodically forced Duffing equation](https://raw.githubusercontent.com/JuliaDynamics/LagrangianDescriptors.jl/main/docs/src/img/LagrangianDescriptors_socialpreview.png)

LagrangianDescriptors.jl is a Julia package for *painting the phase portrait* of dynamical systems with the method of **Lagrangian descriptors**, which reveals singular features of the dynamics, such as stable and unstable manifolds, by coloring each initial condition according to the solution through it. It provides

- Lagrangian descriptors of ordinary, random and stochastic differential equations, given as `ODEProblem`s, `RODEProblem`s and `SDEProblem`s of the [SciML](https://sciml.ai) ecosystem, forward and backward in time.
- Descriptors given by the integral of any infinitesimal descriptor along the solutions, or by the increments of the solutions, as for stochastic equations.
- Means over realizations of the noise, with independent noise for each initial condition, or with the same noise path for all of them.
- Stopping solutions that diverge, separately in each direction of time.
- Plots with [Plots.jl](https://github.com/JuliaPlots/Plots.jl) or [Makie.jl](https://docs.makie.org), as heatmaps of grids of initial conditions, scatter plots of other collections, two-dimensional slices of higher-dimensional systems, or in three dimensions, e.g. as surfaces of initial conditions.

![Lagrangian descriptors of the Lorenz system on slices sweeping through its attractor](https://juliadynamics.github.io/LagrangianDescriptors.jl/dev/examples/lorenz/lorenz_sweeps.gif)

The package is not registered yet. To install it, run `import Pkg; Pkg.add(url = "https://github.com/JuliaDynamics/LagrangianDescriptors.jl")`.

All further information is provided in the documentation, which you can find [online](https://juliadynamics.github.io/LagrangianDescriptors.jl/dev/) or build locally by running the `docs/make.jl` file.
