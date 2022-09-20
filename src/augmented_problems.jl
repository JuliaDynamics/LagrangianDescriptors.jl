# This file is for `aumengtprob` of all types.
# Currently it is only implemented for `ODEProblem`.

"""
    augmentprob(prob::ODEProblem, M; direction::Symbol=:both

Create an augmented `ODEProblem` by aggregating the forward and backward Lagrangian descriptors to the given `ODEProblem`.

More precisely, with given `f = prob.f` and `u0 = prob.u0` and `tspan = prob.tspan`, an augmented `ODEProblem` is that solves the forward and backward ODE and the forward and backward Lagrangian descriptors simultaneously.

The (global) Lagrangian descriptors are based on the provided local descriptor `M=M(du, u, p, t)`. The (forward/backward) Lagrangian descriptor is the time-integration of the local descriptor along the (forward/backward) trajectory.

If `direction=:forward` or `direction=:backward` are given, then only the forward or backard problem is constructed in the augmented system. The default is `direction=:both`, where both directions are considered.
"""
function augmentprob(prob::ODEProblem, M; direction::Symbol=:both)
    if direction == :both
        if isinplace(prob.f)
            faug = function (du, u, p, t)
                fwd, bwd = u.fwd, u.bwd
                prob.f(du.fwd, fwd, p, t)
                prob.f(du.bwd, bwd, p, 2last(prob.tspan)-t)
                du.bwd *= -1
                du.lfwd = M(du.fwd, fwd, p, t)
                du.lbwd = M(du.bwd, bwd, p, 2last(prob.tspan)-t)
                nothing
            end
        else
            faug = function (du, u, p, t)
                fwd, bwd = u.fwd, u.bwd
                du.fwd = prob.f(fwd, p, t)
                du.bwd = -prob.f(bwd, p, 2last(prob.tspan)-t)
                du.lfwd = M(du.fwd, fwd, p, t)
                du.lbwd = M(du.bwd, bwd, p, 2last(prob.tspan)-t)
                nothing
            end
        end

        uaug0 = ComponentVector(fwd = prob.u0, bwd = prob.u0, lfwd = zero(eltype(prob.u0)), lbwd = zero(eltype(prob.u0)))
    elseif direction == :forward
        if isinplace(prob.f)
            faug = function (du, u, p, t)
                fwd = u.fwd
                prob.f(du.fwd, fwd, p, t)
                du.lfwd = M(du.fwd, fwd, p, t)
                nothing
            end
        else
            faug = function (du, u, p, t)
                fwd = u.fwd
                du.fwd = prob.f(fwd, p, t)
                du.lfwd = M(du.fwd, fwd, p, t)
                nothing
            end
        end

        uaug0 = ComponentVector(fwd = prob.u0, lfwd = zero(eltype(prob.u0)))
    elseif direction == :backward
        if isinplace(prob.f)
            faug = function (du, u, p, t)
                bwd = u.bwd
                prob.f(du.bwd, bwd, p, 2last(prob.tspan)-t)
                du.bwd *= -1
                du.lbwd = M(du.bwd, bwd, p, 2last(prob.tspan)-t)
                nothing
            end
        else
            faug = function (du, u, p, t)
                bwd = u.bwd
                du.bwd = -prob.f(bwd, p, 2last(prob.tspan)-t)
                du.lbwd = M(du.bwd, bwd, p, 2last(prob.tspan)-t)
                nothing
            end
        end

        uaug0 = ComponentVector(bwd = prob.u0, lbwd = zero(eltype(prob.u0)))
    else
        throw(
            ArgumentError(
                "Keyword argument `direction = $direction` not implemented; use either `:forward`, `:backward` or `:both`"
            )
        )
    end
    augprob = ODEProblem(faug, uaug0, prob.tspan, prob.p; prob.kwargs...)
    return augprob
end
