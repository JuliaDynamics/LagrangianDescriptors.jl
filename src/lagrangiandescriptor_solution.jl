"""
struct LagrangianDescriptorSolution{T1, T2, T3}
    enssol::T1
    uu0::T2
    direction::T3
end

Representation of the solution to a [`LagrangianDescriptorProblem`](@ref).

## Fields

- `enssol`: the `EnsembleSolution` of the associated ensemble problem in [`LagrangianDescriptorProblem`](@ref).
- `uu0`: the collection of initial conditions given in 
the [`LagrangianDescriptorProblem`](@ref).
- `direction:` the direction given in the [`LagrangianDescriptorProblem`](@ref).
"""
struct LagrangianDescriptorSolution{T1,T2,T3}
    enssol::T1
    uu0::T2
    direction::T3
end

"""
    (lagsol::LagrangianDescriptorSolution)(direction::Symbol=:total)

Depending on the given `direction`, return an array with the sum of the given solution of a Lagrangian descriptor problem, or the forward or the backward or the difference of them, where each element of the array correponds to a given initial condition in the collection of initial conditions.
"""
function (lagsol::LagrangianDescriptorSolution)(direction::Symbol = lagsol.direction)
    if direction == :total || direction == :both
        return getindex.(lagsol.enssol.u, :lfwd) + getindex.(lagsol.enssol.u, :lbwd)
    elseif direction == :forward
        return getindex.(lagsol.enssol.u, :lfwd)
    elseif direction == :backward
        return getindex.(lagsol.enssol.u, :lbwd)
    elseif direction == :difference
        return getindex.(lagsol.enssol.u, :lfwd) - getindex.(lagsol.enssol.u, :lbwd)
    else
        throw(
            ArgumentError(
                "The given argument should be a valid direction:\n" *
                "`:forward`, `:backward`, `:both` or `:total`, or `:difference`",
            ),
        )
    end
end
