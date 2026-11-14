# News

## v0.3.0

First registered version.

- Lagrangian descriptors of `ODEProblem`s, `RODEProblem`s and `SDEProblem`s of the SciML ecosystem, forward and backward in time, given by the integral of any infinitesimal descriptor along the solutions, or by the increments of the solutions (`IncrementDescriptor`).
- Means over `realizations` of the noise, with independent noise for each initial condition, or the same noise path for all of them (`shared_noise`).
- Stopping diverging solutions with `max_trajectory_value`, separately in each direction of time.
- Plots with Plots.jl and Makie.jl, as heatmaps of grids of initial conditions, scatter plots of other collections, and lines for scalar problems, for any two components of higher-dimensional systems, or in three dimensions with three components in `idxs`.
