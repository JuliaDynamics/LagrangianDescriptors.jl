"""
    LagrangianDescriptorProblem

Defines a Lagrangian descriptor problem associated with a SciML differential equations problem.

## Constructor

    LagrangianDescriptorProblem(prob, M, uu0; direction=:both)

`LagrangianDescriptorProblem` can be constructed by passsing a differential equation problem (currently only `ODEProblem` works, but more problem types will be added), an infinitesimal Lagrangian descriptor with arguments compatible with the differential equation problem, an array of initial conditions, and, optionally, the direction of the flow.

### Arguments

- `prob`: the differential equation problem (e.g. `ODEProblem`, `SDEProblem`, `RODEProblem`, etc.).
- `M`: infinitesimal Lagrangian descriptor (e.g. `M=M(du, u, t, p)` for an ODEProblem).
- `uu0`: collection of initial conditions.
- `direction`: the direction of the flow, with default `:both`, but also accepting `:forward` and `:backward`.

### Fields

With the given arguments, the constructor for `LagrangianDescriptorProblem` returns a type with the following arguments:

- `ensprob::T1`: a suitable ensemble problem to be solved with the given collection of initial conditions `uu0` for each solve, with a suitable `prob_func` to iterate through the collection and a suitable `output_func` to only collect the Lagrangian descriptors at the end of the time interval.
- `uu0::T2`: the given collection of initial conditions.
- `direction::T3`: the given or the default direction of the flow.

## Example Problems

"""
struct LagrangianDescriptorProblem{T1, T2, T3}
    ensprob::T1
    uu0::T2
    direction::T3

    function LagrangianDescriptorProblem(prob, M, uu0; direction::Symbol=:both)
        augprob = augmentprob(prob, M; direction)
        if direction == :both
            prob_func = function (augprob,i,repeat; uu0=uu0)
                remake(augprob, u0 = ComponentVector(fwd=uu0[i], bwd=uu0[i], lfwd=0.0, lbwd=0.0))
            end
            output_func = function (sol, i)
                (ComponentArray(lfwd = last(sol).lfwd, lbwd=last(sol).lbwd), false)
            end
        elseif direction == :forward
            prob_func = function (augprob,i,repeat; uu0=uu0)
                remake(augprob, u0 = ComponentVector(fwd=uu0[i], lfwd=0.0))
            end
            output_func = function (sol, i)
                (ComponentArray(lfwd = last(sol).lfwd), false)
            end
        elseif direction == :backward
            prob_func = function (augprob,i,repeat; uu0=uu0)
                remake(augprob, u0 = ComponentVector(bwd=uu0[i], lbwd=0.0))
            end
            output_func = function (sol, i)
                (ComponentArray(lbwd=last(sol).lbwd), false)
            end
        else
            throw(
                ArgumentError(
                    "Keyword argument `direction = $direction` not implemented; use either `:forward`, `:backward` or `:both`"
                )
            )
        end

        ensprob = EnsembleProblem(augprob,prob_func=prob_func, output_func=output_func)
        new{typeof(ensprob), typeof(uu0), typeof(direction)}(ensprob, uu0, direction)
    end
end
