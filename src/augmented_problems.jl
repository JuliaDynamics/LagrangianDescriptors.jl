# This file is for `augmentprob` of all types.
# Currently it is implemented for `ODEProblem`, `RODEProblem` and `SDEProblem`.

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
#
# Each direction is frozen once the absolute value of a component of its solution exceeds
# `vmax` (see `max_trajectory_value` in `solve`): its solution and its Lagrangian descriptor
# stop changing, while the other direction goes on. Since a frozen solution stays above
# `vmax`, no extra state is needed to remember it. With `vmax = Inf`, nothing is frozen.
struct AugmentedRHS{IIP,FWD,BWD,F,M,T,V}
    f::F
    m::M
    t0::T
    vmax::V
end

AugmentedRHS{IIP,FWD,BWD}(f::F, m::M, t0::T, vmax::V = Inf) where {IIP,FWD,BWD,F,M,T,V} =
    AugmentedRHS{IIP,FWD,BWD,F,M,T,V}(f, m, t0, vmax)

function AugmentedRHS(f, m, t0, iip::Bool, direction::Symbol, vmax = Inf)
    _check_direction(direction)
    fwd = direction in (:both, :forward)
    bwd = direction in (:both, :backward)
    return AugmentedRHS{iip,fwd,bwd}(f, m, t0, vmax)
end

# Whether a solution `x`, a number or an array, exceeds `vmax` in absolute value
_exceeds(x, vmax) = isfinite(vmax) && maximum(abs, x) > vmax

function (rhs::AugmentedRHS{IIP,FWD,BWD})(du, u, p, t, args...) where {IIP,FWD,BWD}
    if FWD
        frozen = _exceeds(u.fwd, rhs.vmax)
        if IIP
            dfwd = du.fwd # a view
            frozen ? fill!(dfwd, zero(eltype(dfwd))) : rhs.f(dfwd, u.fwd, p, t, args...)
        else
            du.fwd = frozen ? zero(u.fwd) : rhs.f(u.fwd, p, t, args...)
        end
        du.lfwd = frozen ? zero(eltype(du)) : rhs.m(du.fwd, u.fwd, p, t, args...)
    end
    if BWD
        frozen = _exceeds(u.bwd, rhs.vmax)
        tb = 2rhs.t0 - t
        if IIP
            dbwd = du.bwd # a view; `du.bwd .*= -1` would not infer with ComponentArrays
            if frozen
                fill!(dbwd, zero(eltype(dbwd)))
            else
                rhs.f(dbwd, u.bwd, p, tb, args...)
                dbwd .*= -1
            end
        else
            du.bwd = frozen ? zero(u.bwd) : -rhs.f(u.bwd, p, tb, args...)
        end
        du.lbwd = frozen ? zero(eltype(du)) : rhs.m(du.bwd, u.bwd, p, tb, args...)
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
    augmentprob(prob::ODEProblem, M; direction::Symbol = :both, max_trajectory_value::Real = Inf)

Create an augmented `ODEProblem` by aggregating the forward and backward Lagrangian descriptors to the given `ODEProblem`.

More precisely, with given `f = prob.f`, `u0 = prob.u0`, and `tspan = prob.tspan`, an augmented `ODEProblem` is created that solves the forward and backward ODE and the forward and backward Lagrangian descriptors simultaneously.

With `tspan = (t0, tf)`, the forward problem is solved over `[t0, tf]` and the backward problem over `[2t0 - tf, t0]`, both starting from `u0` at time `t0`. The backward problem is integrated via the change of variables `t ↦ 2t0 - t`, so that both run along the same time interval `tspan`.

The (global) Lagrangian descriptors are based on the provided local descriptor `M=M(du, u, p, t)`. The (forward/backward) Lagrangian descriptor is the time-integration of the local descriptor along the (forward/backward) trajectory.

If `direction=:forward` or `direction=:backward` are given, then only the forward or backward problem is constructed in the augmented system. The default is `direction=:both`, where both directions are considered.

If `max_trajectory_value` is finite, each direction is frozen once the absolute value of a component of its solution exceeds it: its solution and its Lagrangian descriptor stop changing, while the other direction goes on.
"""
function augmentprob(prob::ODEProblem, M; direction::Symbol = :both, max_trajectory_value::Real = Inf)
    faug = AugmentedRHS(prob.f, M, first(prob.tspan), isinplace(prob), direction, max_trajectory_value)
    uaug0 = _augmented_u0(prob.u0, direction)
    return ODEProblem{true}(faug, uaug0, prob.tspan, prob.p; prob.kwargs...)
end

"""
    augmentprob(prob::RODEProblem, M; direction::Symbol = :both, max_trajectory_value::Real = Inf)

Create an augmented `RODEProblem` by aggregating the forward and backward Lagrangian descriptors to the given `RODEProblem`.

More precisely, with given `f = prob.f`, `u0 = prob.u0`, `W = prob.noise`, and `tspan = prob.tspan`, an augmented `RODEProblem` is created that solves the forward and backward RODE and the forward and backward Lagrangian descriptors simultaneously.

With `tspan = (t0, tf)`, the forward problem is solved over `[t0, tf]` and the backward problem over `[2t0 - tf, t0]`, both starting from `u0` at time `t0`. The backward problem is integrated via the change of variables `t ↦ 2t0 - t`, so that both run along the same time interval `tspan`. Note that the noise `W` is evaluated at the forward time `t` in both directions.

If `prob` has neither `noise` nor `rand_prototype`, the solver's default Wiener process starting at `t0` is used, with the shape of `u0`, as when solving `prob` directly. For a scalar `u0`, `f` and `M` still receive a scalar `W`, but the noise of the augmented problem is stored as a one-element vector, since the augmented problem is in place.

The (global) Lagrangian descriptors are based on the provided local descriptor `M=M(du, u, p, t, W)`. The (forward/backward) Lagrangian descriptor is the time-integration of the local descriptor along the (forward/backward) trajectory.

If `direction=:forward` or `direction=:backward` are given, then only the forward or backward problem is constructed in the augmented system. The default is `direction=:both`, where both directions are considered.

If `max_trajectory_value` is finite, each direction is frozen once the absolute value of a component of its solution exceeds it: its solution and its Lagrangian descriptor stop changing, while the other direction goes on.
"""
function augmentprob(prob::RODEProblem, M; direction::Symbol = :both, max_trajectory_value::Real = Inf)
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

    faug = AugmentedRHS(f, m, first(prob.tspan), isinplace(prob), direction, max_trajectory_value)
    uaug0 = _augmented_u0(prob.u0, direction)
    return RODEProblem{true}(faug, uaug0, prob.tspan, prob.p; noise_kwargs..., prob.kwargs...)
end

# Diffusion of an augmented SDE problem, built from the diffusion `g` of the original problem.
# The `fwd` and `bwd` blocks get the diffusion of the original problem, the latter at the
# reversed time `2t0 - t`, and the Lagrangian descriptors `lfwd` and `lbwd` get no noise.
# With diagonal noise (`DIAG`), each component of the augmented state is driven by its own
# Wiener process, so the backward path gets noise independent from the forward one. With
# non-diagonal noise, the original `n × m` noise matrix goes into the block of rows of
# `fwd` and the first `m` columns, and into the block of rows of `bwd` and the next `m`
# columns, so the backward path again gets its own `m` independent noise channels.
#
# A direction frozen by `AugmentedRHS` (see `vmax`) also gets no noise, so it stays frozen.
struct AugmentedDiffusion{IIP,FWD,BWD,DIAG,G,T,R,V}
    g::G
    t0::T
    rows::R # for non-diagonal noise: (fwd rows, bwd rows, fwd columns, bwd columns)
    vmax::V
end

AugmentedDiffusion{IIP,FWD,BWD,DIAG}(g::G, t0::T, rows::R, vmax::V = Inf) where {IIP,FWD,BWD,DIAG,G,T,R,V} =
    AugmentedDiffusion{IIP,FWD,BWD,DIAG,G,T,R,V}(g, t0, rows, vmax)

function (diff::AugmentedDiffusion{IIP,FWD,BWD,true})(du, u, p, t) where {IIP,FWD,BWD}
    if FWD
        frozen = _exceeds(u.fwd, diff.vmax)
        if IIP
            dfwd = du.fwd # a view
            frozen ? fill!(dfwd, zero(eltype(dfwd))) : diff.g(dfwd, u.fwd, p, t)
        else
            du.fwd = frozen ? zero(u.fwd) : diff.g(u.fwd, p, t)
        end
        du.lfwd = zero(eltype(du))
    end
    if BWD
        frozen = _exceeds(u.bwd, diff.vmax)
        tb = 2diff.t0 - t
        if IIP
            dbwd = du.bwd # a view
            frozen ? fill!(dbwd, zero(eltype(dbwd))) : diff.g(dbwd, u.bwd, p, tb)
        else
            du.bwd = frozen ? zero(u.bwd) : diff.g(u.bwd, p, tb)
        end
        du.lbwd = zero(eltype(du))
    end
    return nothing
end

function (diff::AugmentedDiffusion{IIP,FWD,BWD,false})(dG, u, p, t) where {IIP,FWD,BWD}
    # only the two diagonal blocks are written; the other entries of `dG` stay zero
    rfwd, rbwd, cfwd, cbwd = diff.rows
    if FWD
        block = view(dG, rfwd, cfwd)
        if _exceeds(u.fwd, diff.vmax)
            fill!(block, zero(eltype(block)))
        elseif IIP
            diff.g(block, u.fwd, p, t)
        else
            block .= diff.g(u.fwd, p, t)
        end
    end
    if BWD
        block = view(dG, rbwd, cbwd)
        tb = 2diff.t0 - t
        if _exceeds(u.bwd, diff.vmax)
            fill!(block, zero(eltype(block)))
        elseif IIP
            diff.g(block, u.bwd, p, tb)
        else
            block .= diff.g(u.bwd, p, tb)
        end
    end
    return nothing
end

"""
    augmentprob(prob::SDEProblem, M; direction::Symbol = :both, max_trajectory_value::Real = Inf)

Create an augmented `SDEProblem` by aggregating the forward and backward Lagrangian descriptors to the given `SDEProblem`.

More precisely, with given drift `f = prob.f`, diffusion `g = prob.g`, `u0 = prob.u0`, and `tspan = prob.tspan`, an augmented `SDEProblem` is created that solves the forward and backward SDE and the forward and backward Lagrangian descriptors simultaneously.

With `tspan = (t0, tf)`, the forward problem is solved over `[t0, tf]` and the backward problem over `[2t0 - tf, t0]`, both starting from `u0` at time `t0`. The backward problem is integrated via the change of variables `t ↦ 2t0 - t`, so that both run along the same time interval `tspan`: its drift is `-f` and its diffusion is `g`, both evaluated at the time `2t0 - t`.

The noise driving the backward problem is independent of the noise driving the forward problem, as for a two-sided Wiener process, whose increments over `[2t0 - tf, t0]` and over `[t0, tf]` are independent. With diagonal noise, each component of the augmented system has its own Wiener process. With non-diagonal noise, given by an `n × m` `noise_rate_prototype`, the augmented noise matrix has the diffusion of the forward problem and that of the backward problem in separate columns, so each has its own `m` noise channels.

The (global) Lagrangian descriptors are based on the provided local descriptor `M=M(du, u, p, t)`, where `du` is the drift `f(u, p, t)`. The (forward/backward) Lagrangian descriptor is the time-integration of the local descriptor along the (forward/backward) stochastic trajectory.

If `direction=:forward` or `direction=:backward` are given, then only the forward or backward problem is constructed in the augmented system. The default is `direction=:both`, where both directions are considered.

If `max_trajectory_value` is finite, each direction is frozen once the absolute value of a component of its solution exceeds it: its solution and its Lagrangian descriptor stop changing, while the other direction goes on.

Custom noise processes (the `noise` keyword of `SDEProblem`) are not supported, since the augmented system needs noise with a different dimension.
"""
function augmentprob(prob::SDEProblem, M; direction::Symbol = :both, max_trajectory_value::Real = Inf)
    prob.noise === nothing || throw(
        ArgumentError(
            "Custom noise processes are not supported for `SDEProblem`s: the augmented " *
            "problem needs noise with a different dimension. Use the default white noise.",
        ),
    )
    t0 = first(prob.tspan)
    iip = isinplace(prob)
    faug = AugmentedRHS(prob.f, M, t0, iip, direction, max_trajectory_value)
    uaug0 = _augmented_u0(prob.u0, direction)
    fwd = direction in (:both, :forward)
    bwd = direction in (:both, :backward)

    if prob.noise_rate_prototype === nothing
        gaug = AugmentedDiffusion{iip,fwd,bwd,true}(prob.g, t0, nothing, max_trajectory_value)
        noise_kwargs = (;)
    else
        prob.u0 isa AbstractVector || throw(
            ArgumentError("Non-diagonal noise requires the initial condition to be a vector."),
        )
        n, m = size(prob.noise_rate_prototype)
        # `_augmented_u0` lays out the blocks in the order fwd, bwd, lfwd, lbwd
        rfwd = fwd ? (1:n) : (1:0)
        rbwd = bwd ? (fwd ? (n+1:2n) : (1:n)) : (1:0)
        cfwd = fwd ? (1:m) : (1:0)
        cbwd = bwd ? (fwd ? (m+1:2m) : (1:m)) : (1:0)
        nblocks = fwd + bwd
        prototype = zeros(eltype(prob.noise_rate_prototype), length(uaug0), nblocks * m)
        gaug = AugmentedDiffusion{iip,fwd,bwd,false}(prob.g, t0, (rfwd, rbwd, cfwd, cbwd), max_trajectory_value)
        noise_kwargs = (; noise_rate_prototype = prototype)
    end

    return SDEProblem{true}(
        faug, gaug, uaug0, prob.tspan, prob.p;
        noise_kwargs..., seed = prob.seed, prob.kwargs...,
    )
end
