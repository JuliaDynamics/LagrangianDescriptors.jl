"""
    LagrangianDescriptorSolution

The solution of a [`LagrangianDescriptorProblem`](@ref), returned by [`solve`](@ref).

## Lagrangian descriptors

The Lagrangian descriptors are obtained by calling the solution, `lagsol(direction)`; see
below. They are plotted with `plot(lagsol)`; see the section on plotting.

## Fields

- `enssol`: the `EnsembleSolution` of the ensemble problem `prob.ensprob`. Each element of
  `enssol.u` is a `NamedTuple` with the Lagrangian descriptors `lfwd` and/or `lbwd` at the
  end of a trajectory, depending on the direction, and the solver's return code `retcode`.
- `uu0`: the collection of initial conditions of the problem.
- `direction`: the direction of the problem, `:both`, `:forward` or `:backward`.
- `terminated`: a `BitVector` marking the trajectories stopped early by the keyword
  `max_trajectory_value` of [`solve`](@ref).
- `retcodes`: the solver's return code (a `SciMLBase.ReturnCode.T`) of each trajectory.
  Trajectories that did not finish normally, e.g. because they reached `maxiters` or became
  unstable, are found with `.!SciMLBase.successful_retcode.(lagsol.retcodes)`; their
  Lagrangian descriptors are not reliable. Trajectories stopped by `max_trajectory_value`
  have the return code `ReturnCode.Terminated`.

When solved with several `realizations` (see [`solve`](@ref)), `enssol.u`, `terminated` and
`retcodes` have one entry per solve: those of all initial conditions of the first
realization, then of the second one, and so on.
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
    (lagsol::LagrangianDescriptorSolution)(direction::Symbol = lagsol.direction)

The Lagrangian descriptors of the solution, as a vector with one value for each initial
condition, in the order of `vec(lagsol.uu0)`; e.g. `reshape(lagsol(), size(lagsol.uu0))`
arranges them like a matrix of initial conditions. The `direction` is one of

- `:forward`: the forward Lagrangian descriptors;
- `:backward`: the backward Lagrangian descriptors;
- `:both` or `:total`: their sum (the default for a problem with both directions);
- `:difference`: the forward minus the backward Lagrangian descriptors.

The directions other than `:forward` and `:backward` need a problem with both directions.
When solved with several `realizations` (see [`solve`](@ref)), the mean over the
realizations is returned.
"""
function (lagsol::LagrangianDescriptorSolution)(direction::Symbol = lagsol.direction)
    _check_direction(direction, DESCRIPTOR_DIRECTIONS)
    _check_available(lagsol, direction)
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

# Whether the descriptors needed for `direction` were computed, i.e. the solution of a problem
# in a single direction only has that direction
function _check_available(lagsol, direction::Symbol)
    lagsol.direction == :both || direction == lagsol.direction || throw(
        ArgumentError(
            "This solution only has the $(lagsol.direction) Lagrangian descriptors, so " *
            "`:$direction` is not available; solve the problem with `direction = :both` for it.",
        ),
    )
    return nothing
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
