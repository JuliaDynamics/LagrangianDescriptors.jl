# This file is for `augmentprob` of all types.
# Currently it is implemented for `ODEProblem` and `RODEProblem`.

const FLOW_DIRECTIONS = (:both, :forward, :backward)

function _check_direction(direction::Symbol, allowed = FLOW_DIRECTIONS)
    direction in allowed || throw(
        ArgumentError(
            "Direction `$direction` not implemented; use one of " *
            join(("`:$d`" for d in allowed), ", ", " or "),
        ),
    )
    return direction
end

# Right-hand side of an augmented problem, built from the right-hand side `f` and the
# infinitesimal Lagrangian descriptor `m` of the original problem. It computes the forward
# part (`FWD`) and/or the backward part (`BWD`) of the augmented system, calling `f` in place
# or out of place according to `IIP`. These are type parameters, so the unused branches are
# compiled away. The backward part is integrated via the change of variables `t ↦ 2t0 - t`.
# Any extra arguments after `t` (e.g. the noise `W` of a `RODEProblem`) are passed on to
# `f` and `m`.
struct AugmentedRHS{IIP,FWD,BWD,F,M,T}
    f::F
    m::M
    t0::T
end

AugmentedRHS{IIP,FWD,BWD}(f::F, m::M, t0::T) where {IIP,FWD,BWD,F,M,T} =
    AugmentedRHS{IIP,FWD,BWD,F,M,T}(f, m, t0)

function AugmentedRHS(f, m, t0, iip::Bool, direction::Symbol)
    _check_direction(direction)
    fwd = direction in (:both, :forward)
    bwd = direction in (:both, :backward)
    return AugmentedRHS{iip,fwd,bwd}(f, m, t0)
end

function (rhs::AugmentedRHS{IIP,FWD,BWD})(du, u, p, t, args...) where {IIP,FWD,BWD}
    if FWD
        if IIP
            rhs.f(du.fwd, u.fwd, p, t, args...)
        else
            du.fwd = rhs.f(u.fwd, p, t, args...)
        end
        du.lfwd = rhs.m(du.fwd, u.fwd, p, t, args...)
    end
    if BWD
        tb = 2rhs.t0 - t
        if IIP
            dbwd = du.bwd # a view; `du.bwd .*= -1` would not infer with ComponentArrays
            rhs.f(dbwd, u.bwd, p, tb, args...)
            dbwd .*= -1
        else
            du.bwd = -rhs.f(u.bwd, p, tb, args...)
        end
        du.lbwd = rhs.m(du.bwd, u.bwd, p, tb, args...)
    end
    return nothing
end

# Initial condition of the augmented system, for the initial condition `u0` of the original one.
function _augmented_u0(u0, direction::Symbol)
    l0 = zero(eltype(u0))
    direction == :both && return ComponentVector(fwd = u0, bwd = u0, lfwd = l0, lbwd = l0)
    direction == :forward && return ComponentVector(fwd = u0, lfwd = l0)
    direction == :backward && return ComponentVector(bwd = u0, lbwd = l0)
    _check_direction(direction)
end

"""
    augmentprob(prob::ODEProblem, M; direction::Symbol=:both)

Create an augmented `ODEProblem` by aggregating the forward and backward Lagrangian descriptors to the given `ODEProblem`.

More precisely, with given `f = prob.f`, `u0 = prob.u0`, and `tspan = prob.tspan`, an augmented `ODEProblem` is created that solves the forward and backward ODE and the forward and backward Lagrangian descriptors simultaneously.

With `tspan = (t0, tf)`, the forward problem is solved over `[t0, tf]` and the backward problem over `[2t0 - tf, t0]`, both starting from `u0` at time `t0`. The backward problem is integrated via the change of variables `t ↦ 2t0 - t`, so that both run along the same time interval `tspan`.

The (global) Lagrangian descriptors are based on the provided local descriptor `M=M(du, u, p, t)`. The (forward/backward) Lagrangian descriptor is the time-integration of the local descriptor along the (forward/backward) trajectory.

If `direction=:forward` or `direction=:backward` are given, then only the forward or backard problem is constructed in the augmented system. The default is `direction=:both`, where both directions are considered.
"""
function augmentprob(prob::ODEProblem, M; direction::Symbol = :both)
    faug = AugmentedRHS(prob.f, M, first(prob.tspan), isinplace(prob), direction)
    uaug0 = _augmented_u0(prob.u0, direction)
    return ODEProblem{true}(faug, uaug0, prob.tspan, prob.p; prob.kwargs...)
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
    noise_kwargs = if prob.noise !== nothing
        (; noise = prob.noise)
    elseif prob.rand_prototype !== nothing
        (; rand_prototype = prob.rand_prototype)
    elseif scalar_default_noise
        (; rand_prototype = [zero(prob.u0)])
    else
        (; rand_prototype = zero(prob.u0))
    end

    faug = AugmentedRHS(f, m, first(prob.tspan), isinplace(prob), direction)
    uaug0 = _augmented_u0(prob.u0, direction)
    return RODEProblem{true}(faug, uaug0, prob.tspan, prob.p; noise_kwargs..., prob.kwargs...)
end
