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
  When solving with several `realizations` (see [`solve`](@ref)), `enssol.u`, `terminated`
  and `retcodes` have one entry per solve, for all initial conditions of the first
  realization, then of the second one, and so on.
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

When solved with several `realizations` (see [`solve`](@ref)), it returns the mean over the realizations.
"""
function (lagsol::LagrangianDescriptorSolution)(direction::Symbol = lagsol.direction)
    _check_direction(direction, DESCRIPTOR_DIRECTIONS)
    u = lagsol.enssol.u
    values = if direction == :forward
        getindex.(u, :lfwd)
    elseif direction == :backward
        getindex.(u, :lbwd)
    elseif direction == :difference
        getindex.(u, :lfwd) - getindex.(u, :lbwd)
    else # :both or :total
        getindex.(u, :lfwd) + getindex.(u, :lbwd)
    end
    return _mean_over_realizations(values, length(lagsol.uu0))
end

# Mean over realizations of values ordered as all `n` initial conditions of each realization
function _mean_over_realizations(values, n)
    length(values) == n && return values
    return vec(sum(reshape(values, n, :), dims = 2)) / (length(values) ÷ n)
end

# Compact display, instead of the default one, which prints all initial conditions and
# all trajectory outputs

function Base.show(io::IO, ::MIME"text/plain", prob::LagrangianDescriptorProblem)
    print(io, "LagrangianDescriptorProblem for a", _with_article(nameof(typeof(prob.ensprob.prob))))
    print(io, "\n  initial conditions: ", _describe_collection(prob.uu0))
    print(io, "\n  direction: :", prob.direction)
end

function Base.show(io::IO, ::MIME"text/plain", lagsol::LagrangianDescriptorSolution)
    n = length(lagsol.uu0)
    realizations = length(lagsol.retcodes) ÷ n
    print(io, "LagrangianDescriptorSolution")
    print(io, "\n  initial conditions: ", _describe_collection(lagsol.uu0))
    print(io, "\n  direction: :", lagsol.direction)
    realizations > 1 && print(io, "\n  realizations: ", realizations)
    failed = count(!SciMLBase.successful_retcode, lagsol.retcodes)
    failed > 0 && print(io, "\n  unsuccessful solves: ", failed, " (see `retcodes`)")
    any(lagsol.terminated) &&
        print(io, "\n  terminated by `max_trajectory_value`: ", count(lagsol.terminated))
end

Base.show(io::IO, prob::LagrangianDescriptorProblem) =
    print(io, "LagrangianDescriptorProblem(", length(prob.uu0), " initial conditions, :", prob.direction, ")")
Base.show(io::IO, lagsol::LagrangianDescriptorSolution) =
    print(io, "LagrangianDescriptorSolution(", length(lagsol.uu0), " initial conditions, :", lagsol.direction, ")")

_with_article(name) = (first(string(name)) in "AEIOU" ? "n " : " ") * string(name)

function _describe_collection(uu0)
    n = length(uu0)
    shape = uu0 isa AbstractArray && ndims(uu0) > 1 ? " (" * join(size(uu0), " × ") * ")" : ""
    return string(n, shape, " of type ", typeof(first(uu0)))
end
