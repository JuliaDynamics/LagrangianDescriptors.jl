"""
    solve(prob::LagrangianDescriptorProblem, alg, args...; max_trajectory_value::Real=Inf, kwargs...)

Solve a [`LagrangianDescriptorProblem`](@ref), which amounts to solving the associated `EnsembleProblem` in `prob.ensprob` and returning a [`LagrangianDescriptorSolution`](@ref).

You should provide the necessary `args` and the desired `kwargs` for solving the associated ensemble problem for the underlying Differential Equation problem.

## Keyword Arguments

- `max_trajectory_value::Real=Inf`: Maximum allowed value for trajectory components. When any
  component of the trajectory exceeds this threshold, the integration is terminated early.
  This is useful for dissipative systems where backward integration may diverge.
  The Lagrangian descriptor value accumulated up to termination is preserved.
"""
function solve(prob::LagrangianDescriptorProblem, alg, args...; max_trajectory_value::Real=Inf, kwargs...)
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

    sol = solve(prob.ensprob, alg, args...; trajectories = length(prob.uu0), kwargs...)

    # Extract terminated mask from the output ComponentArrays
    terminated = BitVector([Bool(s.terminated) for s in sol.u])

    return LagrangianDescriptorSolution(sol, prob.uu0, prob.direction, terminated)
end

"""
    _make_termination_condition(direction::Symbol, max_val::Real)

Create a termination condition function that checks if trajectory components exceed `max_val`.
"""
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
