"""
struct LagrangianDescriptorSolution{T1, T2, T3}
    enssol::T1
    uu0::T2
    direction::T3
    terminated::BitVector
end

Representation of the solution to a [`LagrangianDescriptorProblem`](@ref).

## Fields

- `enssol`: the `EnsembleSolution` of the associated ensemble problem in [`LagrangianDescriptorProblem`](@ref).
- `uu0`: the collection of initial conditions given in 
the [`LagrangianDescriptorProblem`](@ref).
- `direction:` the direction given in the [`LagrangianDescriptorProblem`](@ref).
- `terminated`: a `BitVector` indicating which trajectories were terminated early
  (e.g., due to divergence when using `max_trajectory_value`).
"""
struct LagrangianDescriptorSolution{T1,T2,T3}
    enssol::T1
    uu0::T2
    direction::T3
    terminated::BitVector
end

# Convenience constructor for backward compatibility (no termination tracking)
function LagrangianDescriptorSolution(enssol::T1, uu0::T2, direction::T3) where {T1,T2,T3}
    return LagrangianDescriptorSolution(enssol, uu0, direction, falses(length(uu0)))
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
