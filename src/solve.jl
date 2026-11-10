"""
    solve(prob::LagrangianDescriptorProblem, alg, args...; max_trajectory_value = Inf, realizations = 1, kwargs...)

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

- `max_trajectory_value::Real = Inf`: if finite, a trajectory is stopped as soon as the
  absolute value of any component of its forward or backward solution exceeds it. This is
  useful for systems whose solutions diverge, e.g. dissipative systems backward in time.
  The Lagrangian descriptors accumulated up to that time are kept, and the trajectory is
  marked in the field `terminated` of the solution.
- `realizations::Integer = 1`: the number of times each initial condition is solved. For
  random or stochastic equations, each realization has independent noise, and the solution
  gives the mean of the Lagrangian descriptors over the realizations, as in Balibrea-Iniesta,
  Lopesino, Wiggins and Mancho (2016). The return codes and the `terminated` flags of the
  solution are kept for each of the `realizations * length(prob.uu0)` solves.

## Examples

```julia
lagsol = solve(lagprob, Tsit5())                          # an ODE problem
lagsol = solve(lagprob, Tsit5(), EnsembleSerial(); abstol = 1e-8, reltol = 1e-8)
lagsol = solve(lagprob, Tsit5(); max_trajectory_value = 1e3)
lagsol = solve(lagprob, RandomHeun(), dt = 1 / 100)       # a RODE problem
lagsol = solve(lagprob, SOSRI(); realizations = 50)       # an SDE problem
```
"""
function solve(prob::LagrangianDescriptorProblem, alg, args...; max_trajectory_value::Real=Inf, realizations::Integer=1, kwargs...)
    realizations ≥ 1 || throw(
        ArgumentError("The number of `realizations` must be positive; got `$realizations`."),
    )
    # This first solve was a hack on 1.7.2 rosetta; otherwise the subsequente solve would hang
    # But it is not needed on mac native and probably not on other systems as well
    # solve(prob.ensprob.prob, alg; kwargs...)

    if isfinite(max_trajectory_value)
        # Build termination callback based on direction
        condition = _make_termination_condition(prob.direction, max_trajectory_value)
        affect!(integrator) = terminate!(integrator)
        cb = DiscreteCallback(condition, affect!)

        # Merge with any existing callbacks
        existing_cb = get(kwargs, :callback, nothing)
        merged_cb = existing_cb === nothing ? cb : CallbackSet(existing_cb, cb)
        kwargs = merge(NamedTuple(kwargs), (callback=merged_cb,))
    end

    trajectories = realizations * length(prob.uu0)
    sol = solve(prob.ensprob, alg, args...; trajectories, kwargs...)
    return LagrangianDescriptorSolution(sol, prob.uu0, prob.direction)
end

# Internal helper: create a termination condition function that checks if
# trajectory components exceed `max_val`, for the given `direction`.
function _make_termination_condition(direction::Symbol, max_val::Real)
    if direction == :both
        return (u, t, integrator) -> maximum(abs, u.fwd) > max_val || maximum(abs, u.bwd) > max_val
    elseif direction == :forward
        return (u, t, integrator) -> maximum(abs, u.fwd) > max_val
    elseif direction == :backward
        return (u, t, integrator) -> maximum(abs, u.bwd) > max_val
    else
        error("Unknown direction: $direction")
    end
end
