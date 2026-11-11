# `EnsembleProblem` calls `prob_func(prob, i, repeat)` in SciMLBase < 3 and
# `prob_func(prob, ctx::EnsembleContext)` in SciMLBase >= 3.
_trajectory_index(i::Integer, repeat) = i
_trajectory_index(ctx) = ctx.sim_id

# Output of each trajectory: the Lagrangian descriptors `lfwd` and/or `lbwd`, the solver's
# return code, and whether each direction was stopped by `max_trajectory_value`, from the
# final state `u` of the augmented system
function _trajectory_output(u, direction::Symbol, retcode, vmax, lfwd = nothing, lbwd = nothing)
    if direction == :both
        return (lfwd = something(lfwd, u.lfwd), lbwd = something(lbwd, u.lbwd), retcode = retcode,
            stopped_fwd = _exceeds(u.fwd, vmax), stopped_bwd = _exceeds(u.bwd, vmax))
    elseif direction == :forward
        return (lfwd = something(lfwd, u.lfwd), retcode = retcode, stopped_fwd = _exceeds(u.fwd, vmax))
    else
        return (lbwd = something(lbwd, u.lbwd), retcode = retcode, stopped_bwd = _exceeds(u.bwd, vmax))
    end
end

# The augmented problem and the output of each trajectory of the ensemble, for an
# infinitesimal Lagrangian descriptor `M`; see `increment_descriptor.jl` for another method
function _ensemble_setup(prob, M, direction::Symbol, vmax)
    output_func = function (sol, i)
        (_trajectory_output(sol.u[end], direction, sol.retcode, vmax), false)
    end
    return augmentprob(prob, M; direction, max_trajectory_value = vmax), output_func
end

# Seeds of the noise of each trajectory `i`, or `nothing` to leave them to the solver, which
# then gives independent noise to each trajectory. Only `SDEProblem`s use a seed for their noise.
_problem_seed(prob::SDEProblem) = prob.seed
_problem_seed(prob) = zero(UInt64)

function _derived_seed(base, k)
    seed = UInt64(hash((base, k)))
    return iszero(seed) ? one(UInt64) : seed
end

function _seed_function(prob, n::Integer, shared_noise::Bool)
    base = _problem_seed(prob)
    if shared_noise
        prob isa SDEProblem || throw(
            ArgumentError("`shared_noise = true` is only supported for `SDEProblem`s."),
        )
        b = iszero(base) ? rand(UInt64) : base
        # the same seed for all initial conditions of each realization
        return i -> _derived_seed(b, div(i - 1, n) + 1)
    elseif !iszero(base)
        # independent but reproducible noise for each trajectory
        return i -> _derived_seed(base, i)
    else
        return nothing
    end
end

# The ensemble problem solved for the Lagrangian descriptors; with several realizations (see
# `solve`), trajectory `i` starts at `uu0[mod1(i, length(uu0))]`
function _ensemble_problem(prob, M, uu0, direction::Symbol; vmax = Inf, shared_noise::Bool = false)
    augprob, output_func = _ensemble_setup(prob, M, direction, vmax)
    seed_of = _seed_function(prob, length(uu0), shared_noise)
    prob_func = function (augprob, args...)
        i = _trajectory_index(args...)
        u0 = _augmented_u0(uu0[mod1(i, length(uu0))], direction)
        return seed_of === nothing ? remake(augprob; u0) : remake(augprob; u0, seed = seed_of(i))
    end
    return EnsembleProblem(augprob; prob_func, output_func)
end

"""
    LagrangianDescriptorProblem(prob, M, uu0; direction = :both)

A Lagrangian descriptor problem: the computation of the Lagrangian descriptors of the
differential equation problem `prob`, given by the infinitesimal descriptor `M`, at each
initial condition in the collection `uu0`. Solve it with [`solve`](@ref) to obtain a
[`LagrangianDescriptorSolution`](@ref).

## Arguments

- `prob`: an `ODEProblem`, a `RODEProblem` or an `SDEProblem`. Its time span `(t0, tf)`
  sets the time windows: `[t0, tf]` forward and `[2t0 - tf, t0]` backward. Its initial
  condition is replaced by the elements of `uu0`, and everything else, such as its
  parameters and keyword arguments, is kept.
- `M`: the infinitesimal Lagrangian descriptor, integrated along the forward and backward
  solutions. Its signature depends on the problem:
    - `ODEProblem`: `M(du, u, p, t)`, where `du` is the value of the right-hand side;
    - `RODEProblem`: `M(du, u, p, t, W)`, where `W` is the value of the noise;
    - `SDEProblem`: `M(du, u, p, t)`, where `du` is the value of the drift.
  Alternatively, an [`IncrementDescriptor`](@ref), for the descriptor given by the
  increments of the solution.
- `uu0`: a collection of initial conditions, such as a vector, or a matrix for a grid of
  initial conditions of a two-dimensional problem. Each element is used as an initial
  condition of `prob`, so it is a number for a scalar problem, or a vector otherwise.
- `direction`: `:both` (default), `:forward` or `:backward`; the directions of time in which
  the Lagrangian descriptors are computed.

## Fields

- `ensprob`: the `EnsembleProblem` solved by [`solve`](@ref) with its default keyword
  arguments. Its trajectories solve the augmented problem built by [`augmentprob`](@ref),
  and their outputs are the Lagrangian descriptors at the end of the time interval, along
  with the solver's return codes.
- `uu0`: the given collection of initial conditions.
- `direction`: the given direction.
- `prob` and `M`: the given problem and descriptor, from which [`solve`](@ref) builds the
  ensemble problem for other keyword arguments.

## Example

Lagrangian descriptors of a periodically forced Duffing equation, over a grid of initial
conditions:

```julia
using LagrangianDescriptors, OrdinaryDiffEq, Plots
using LinearAlgebra: norm

function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * cos(ω * t)
end

u0 = [0.0, 0.0]
tspan = (0.0, 13.0)
A = 0.3; ω = π; p = (A, ω)

prob = ODEProblem(f!, u0, tspan, p)
M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1.0, 1.0, length = 301), x in range(-1.8, 1.8, length = 301)]

lagprob = LagrangianDescriptorProblem(prob, M, uu0)
lagsol = solve(lagprob, Tsit5())
plot(lagsol)
```
"""
struct LagrangianDescriptorProblem{T1,T2,T3,P,D}
    ensprob::T1
    uu0::T2
    direction::T3
    prob::P
    M::D

    function LagrangianDescriptorProblem(prob, M, uu0; direction::Symbol = :both)
        _check_direction(direction)
        ensprob = _ensemble_problem(prob, M, uu0, direction)
        new{typeof(ensprob),typeof(uu0),typeof(direction),typeof(prob),typeof(M)}(
            ensprob, uu0, direction, prob, M,
        )
    end
end
