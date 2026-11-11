"""
    solve(prob::LagrangianDescriptorProblem, alg, args...; max_trajectory_value = Inf, realizations = 1, shared_noise = false, kwargs...)

Solve a [`LagrangianDescriptorProblem`](@ref), returning a
[`LagrangianDescriptorSolution`](@ref).

Each initial condition in `prob.uu0` is solved with the algorithm `alg`, as one trajectory of
the ensemble problem `prob.ensprob`. The remaining positional arguments `args` and keyword
arguments `kwargs` are passed on to the ensemble solve, and so to the solve of each
trajectory. For instance, an ensemble algorithm, such as `EnsembleThreads()` (the default),
`EnsembleSerial()` or `EnsembleDistributed()`, can be given after `alg`, and tolerances
(`abstol`, `reltol`), a time step (`dt`), or a `callback` are given as keyword arguments. The
number of trajectories is set from `prob.uu0` and `realizations`.

## Keyword arguments

- `max_trajectory_value::Real = Inf`: if finite, the forward and the backward solutions are
  each stopped as soon as the absolute value of any of their components exceeds it. This is
  useful for systems whose solutions diverge, e.g. dissipative systems backward in time. The
  two directions are treated separately: a direction that is stopped keeps the Lagrangian
  descriptor accumulated up to then, while the other direction goes on. The stopped
  directions are marked in the fields `terminated_forward` and `terminated_backward` of the
  solution, and `terminated` marks the trajectories with any stopped direction.
- `realizations::Integer = 1`: the number of times each initial condition is solved. For
  random or stochastic equations, each realization has different noise, and the solution
  gives the mean of the Lagrangian descriptors over the realizations, as in Balibrea-Iniesta,
  Lopesino, Wiggins and Mancho (2016). The return codes and the `terminated` flags of the
  solution are kept for each of the `realizations * length(prob.uu0)` solves.
- `shared_noise::Bool = false`: for an `SDEProblem`, whether all initial conditions of a
  realization share the same noise path, as in the definition of the Lagrangian descriptors of
  stochastic equations of Balibrea-Iniesta, Lopesino, Wiggins and Mancho (2016), for a fixed
  realization of the noise. By default, each initial condition has independent noise. If the
  `SDEProblem` has a `seed`, the noise is reproducible in both cases. Shared noise is only
  supported for `SDEProblem`s.

## Examples

```julia
lagsol = solve(lagprob, Tsit5())                          # an ODE problem
lagsol = solve(lagprob, Tsit5(), EnsembleSerial(); abstol = 1e-8, reltol = 1e-8)
lagsol = solve(lagprob, Tsit5(); max_trajectory_value = 1e3)
lagsol = solve(lagprob, RandomHeun(), dt = 1 / 100)       # a RODE problem
lagsol = solve(lagprob, SOSRI(); realizations = 50)       # an SDE problem
lagsol = solve(lagprob, EM(), dt = 0.01; shared_noise = true, realizations = 10)
```
"""
function solve(
        prob::LagrangianDescriptorProblem, alg, args...;
        max_trajectory_value::Real = Inf, realizations::Integer = 1, shared_noise::Bool = false,
        kwargs...,
    )
    realizations ≥ 1 || throw(
        ArgumentError("The number of `realizations` must be positive; got `$realizations`."),
    )
    ensprob = _ensemble_problem(prob.prob, prob.M, prob.uu0, prob.direction;
        vmax = max_trajectory_value, shared_noise)
    trajectories = realizations * length(prob.uu0)
    sol = solve(ensprob, alg, args...; trajectories, kwargs...)
    return LagrangianDescriptorSolution(sol, prob.uu0, prob.direction)
end
