"""
LagrangianDescriptorProblem{T1, T2, T3}
    ensprob::T1
    uu0::T2
    direction::T3
end
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
