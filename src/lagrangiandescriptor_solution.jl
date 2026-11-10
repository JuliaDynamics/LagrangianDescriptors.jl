"""
struct LagrangianDescriptorSolution{T1, T2, T3}
    enssol::T1
    uu0::T2
    direction::T3
    terminated::BitVector
    retcodes::Vector{SciMLBase.ReturnCode.T}
end

Representation of the solution to a [`LagrangianDescriptorProblem`](@ref).

## Fields

- `enssol`: the `EnsembleSolution` of the associated ensemble problem in [`LagrangianDescriptorProblem`](@ref). Each element of `enssol.u` is a `NamedTuple` with the Lagrangian descriptors `lfwd` and/or `lbwd` at the end of the trajectory, depending on the direction, and the solver's return code `retcode`.
- `uu0`: the collection of initial conditions given in 
the [`LagrangianDescriptorProblem`](@ref).
- `direction:` the direction given in the [`LagrangianDescriptorProblem`](@ref).
- `terminated`: a `BitVector` indicating which trajectories were terminated early
  (e.g., due to divergence when using `max_trajectory_value`).
- `retcodes`: the solver's return code (a `SciMLBase.ReturnCode.T`) of each trajectory.
  Trajectories that did not finish normally, e.g. because they reached `maxiters` or
  became unstable, can be found with `.!SciMLBase.successful_retcode.(lagsol.retcodes)`;
  their Lagrangian descriptors are not reliable. Trajectories stopped by
  `max_trajectory_value` have the return code `ReturnCode.Terminated`.
"""
struct LagrangianDescriptorSolution{T1,T2,T3}
    enssol::T1
    uu0::T2
    direction::T3
    terminated::BitVector
    retcodes::Vector{SciMLBase.ReturnCode.T}
end

function LagrangianDescriptorSolution(enssol, uu0, direction)
    retcodes = [s.retcode for s in enssol.u]
    terminated = BitVector(retcodes .== SciMLBase.ReturnCode.Terminated)
    return LagrangianDescriptorSolution(enssol, uu0, direction, terminated, retcodes)
end

const DESCRIPTOR_DIRECTIONS = (:forward, :backward, :both, :total, :difference)

"""
    (lagsol::LagrangianDescriptorSolution)(direction::Symbol=:total)

Depending on the given `direction`, return an array with the sum of the given solution of a Lagrangian descriptor problem, or the forward or the backward or the difference of them, where each element of the array correponds to a given initial condition in the collection of initial conditions.
"""
function (lagsol::LagrangianDescriptorSolution)(direction::Symbol = lagsol.direction)
    _check_direction(direction, DESCRIPTOR_DIRECTIONS)
    u = lagsol.enssol.u
    direction == :forward && return getindex.(u, :lfwd)
    direction == :backward && return getindex.(u, :lbwd)
    direction == :difference && return getindex.(u, :lfwd) - getindex.(u, :lbwd)
    return getindex.(u, :lfwd) + getindex.(u, :lbwd) # :both or :total
end
