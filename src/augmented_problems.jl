# This file is for `augmentprob` of all types.
# Currently it is only implemented for `ODEProblem`.

"""
    augmentprob(prob::ODEProblem, M; direction::Symbol=:both)

Create an augmented `ODEProblem` by aggregating the forward and backward Lagrangian descriptors to the given `ODEProblem`.

More precisely, with given `f = prob.f`, `u0 = prob.u0`, and `tspan = prob.tspan`, an augmented `ODEProblem` is created that solves the forward and backward ODE and the forward and backward Lagrangian descriptors simultaneously.

With `tspan = (t0, tf)`, the forward problem is solved over `[t0, tf]` and the backward problem over `[2t0 - tf, t0]`, both starting from `u0` at time `t0`. The backward problem is integrated via the change of variables `t ↦ 2t0 - t`, so that both run along the same time interval `tspan`.

The (global) Lagrangian descriptors are based on the provided local descriptor `M=M(du, u, p, t)`. The (forward/backward) Lagrangian descriptor is the time-integration of the local descriptor along the (forward/backward) trajectory.

If `direction=:forward` or `direction=:backward` are given, then only the forward or backard problem is constructed in the augmented system. The default is `direction=:both`, where both directions are considered.
"""
function augmentprob(prob::ODEProblem, M; direction::Symbol = :both)
    if direction == :both
        if isinplace(prob.f)
            faug = function (du, u, p, t)
                fwd, bwd = u.fwd, u.bwd
                prob.f(du.fwd, fwd, p, t)
                prob.f(du.bwd, bwd, p, 2first(prob.tspan) - t)
                du.bwd *= -1
                du.lfwd = M(du.fwd, fwd, p, t)
                du.lbwd = M(du.bwd, bwd, p, 2first(prob.tspan) - t)
                nothing
            end
        else
            faug = function (du, u, p, t)
                fwd, bwd = u.fwd, u.bwd
                du.fwd = prob.f(fwd, p, t)
                du.bwd = -prob.f(bwd, p, 2first(prob.tspan) - t)
                du.lfwd = M(du.fwd, fwd, p, t)
                du.lbwd = M(du.bwd, bwd, p, 2first(prob.tspan) - t)
                nothing
            end
        end

        uaug0 = ComponentVector(
            fwd = prob.u0,
            bwd = prob.u0,
            lfwd = zero(eltype(prob.u0)),
            lbwd = zero(eltype(prob.u0)),
        )
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
                prob.f(du.bwd, bwd, p, 2first(prob.tspan) - t)
                du.bwd *= -1
                du.lbwd = M(du.bwd, bwd, p, 2first(prob.tspan) - t)
                nothing
            end
        else
            faug = function (du, u, p, t)
                bwd = u.bwd
                du.bwd = -prob.f(bwd, p, 2first(prob.tspan) - t)
                du.lbwd = M(du.bwd, bwd, p, 2first(prob.tspan) - t)
                nothing
            end
        end

        uaug0 = ComponentVector(bwd = prob.u0, lbwd = zero(eltype(prob.u0)))
    else
        throw(
            ArgumentError(
                "Keyword argument `direction = $direction` not implemented; use either `:forward`, `:backward` or `:both`",
            ),
        )
    end
    augprob = ODEProblem(faug, uaug0, prob.tspan, prob.p; prob.kwargs...)
    return augprob
end

"""
    augmentprob(prob::RODEProblem, M; direction::Symbol=:both)

Create an augmented `RODEProblem` by aggregating the forward and backward Lagrangian descriptors to the given `RODEProblem`.

More precisely, with given `f = prob.f`, `u0 = prob.u0`, `W = prob.noise`, and `tspan = prob.tspan`, an augmented `RODEProblem` is created that solves the forward and backward RODE and the forward and backward Lagrangian descriptors simultaneously.

With `tspan = (t0, tf)`, the forward problem is solved over `[t0, tf]` and the backward problem over `[2t0 - tf, t0]`, both starting from `u0` at time `t0`. The backward problem is integrated via the change of variables `t ↦ 2t0 - t`, so that both run along the same time interval `tspan`. Note that the noise `W` is evaluated at the forward time `t` in both directions.

If `prob` has neither `noise` nor `rand_prototype`, the solver's default Wiener process starting at `t0` is used, with the shape of `u0`, as when solving `prob` directly. For a scalar `u0`, `f` and `M` still receive a scalar `W`, but the noise of the augmented problem is stored as a one-element vector, since the augmented problem is in place.

The (global) Lagrangian descriptors are based on the provided local descriptor `M=M(du, u, p, t, W)`. The (forward/backward) Lagrangian descriptor is the time-integration of the local descriptor along the (forward/backward) trajectory.

If `direction=:forward` or `direction=:backward` are given, then only the forward or backard problem is constructed in the augmented system. The default is `direction=:both`, where both directions are considered.
"""
function augmentprob(prob::RODEProblem, M; direction::Symbol = :both)
    # With neither `noise` nor `rand_prototype`, the solver creates a default noise shaped
    # like the augmented state, so we give it the shape of `u0` instead. The augmented
    # problem is in place, so its noise is in place too, which a scalar cannot be: for a
    # scalar `u0`, use a one-element noise and pass its only entry to `f` and `M`.
    scalar_default_noise =
        prob.noise === nothing && prob.rand_prototype === nothing && prob.u0 isa Number
    f, m = if scalar_default_noise
        (u, p, t, W) -> prob.f(u, p, t, first(W)), (du, u, p, t, W) -> M(du, u, p, t, first(W))
    else
        prob.f, M
    end

    if direction == :both
        if isinplace(prob.f)
            faug = function (du, u, p, t, W)
                fwd, bwd = u.fwd, u.bwd
                f(du.fwd, fwd, p, t, W)
                f(du.bwd, bwd, p, 2first(prob.tspan) - t, W)
                du.bwd *= -1
                du.lfwd = m(du.fwd, fwd, p, t, W)
                du.lbwd = m(du.bwd, bwd, p, 2first(prob.tspan) - t, W)
                nothing
            end
        else
            faug = function (du, u, p, t, W)
                fwd, bwd = u.fwd, u.bwd
                du.fwd = f(fwd, p, t, W)
                du.bwd = -f(bwd, p, 2first(prob.tspan) - t, W)
                du.lfwd = m(du.fwd, fwd, p, t, W)
                du.lbwd = m(du.bwd, bwd, p, 2first(prob.tspan) - t, W)
                nothing
            end
        end

        uaug0 = ComponentVector(
            fwd = prob.u0,
            bwd = prob.u0,
            lfwd = zero(eltype(prob.u0)),
            lbwd = zero(eltype(prob.u0)),
        )
    elseif direction == :forward
        if isinplace(prob.f)
            faug = function (du, u, p, t, W)
                fwd = u.fwd
                f(du.fwd, fwd, p, t, W)
                du.lfwd = m(du.fwd, fwd, p, t, W)
                nothing
            end
        else
            faug = function (du, u, p, t, W)
                fwd = u.fwd
                du.fwd = f(fwd, p, t, W)
                du.lfwd = m(du.fwd, fwd, p, t, W)
                nothing
            end
        end

        uaug0 = ComponentVector(fwd = prob.u0, lfwd = zero(eltype(prob.u0)))
    elseif direction == :backward
        if isinplace(prob.f)
            faug = function (du, u, p, t, W)
                bwd = u.bwd
                f(du.bwd, bwd, p, 2first(prob.tspan) - t, W)
                du.bwd *= -1
                du.lbwd = m(du.bwd, bwd, p, 2first(prob.tspan) - t, W)
                nothing
            end
        else
            faug = function (du, u, p, t, W)
                bwd = u.bwd
                du.bwd = -f(bwd, p, 2first(prob.tspan) - t, W)
                du.lbwd = m(du.bwd, bwd, p, 2first(prob.tspan) - t, W)
                nothing
            end
        end

        uaug0 = ComponentVector(bwd = prob.u0, lbwd = zero(eltype(prob.u0)))
    else
        throw(
            ArgumentError(
                "Keyword argument `direction = $direction` not implemented; use either `:forward`, `:backward` or `:both`",
            ),
        )
    end
    noise_kwargs = if prob.noise !== nothing
        (; noise = prob.noise)
    elseif prob.rand_prototype !== nothing
        (; rand_prototype = prob.rand_prototype)
    elseif scalar_default_noise
        (; rand_prototype = [zero(prob.u0)])
    else
        (; rand_prototype = zero(prob.u0))
    end
    augprob = RODEProblem(faug, uaug0, prob.tspan, prob.p; noise_kwargs..., prob.kwargs...)

    return augprob
end
