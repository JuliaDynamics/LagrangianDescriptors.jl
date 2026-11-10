# `EnsembleProblem` calls `prob_func(prob, i, repeat)` in SciMLBase < 3 and
# `prob_func(prob, ctx::EnsembleContext)` in SciMLBase >= 3.
_trajectory_index(i::Integer, repeat) = i
_trajectory_index(ctx) = ctx.sim_id

# Output of each trajectory: the Lagrangian descriptors at the end of the trajectory and
# the solver's return code
function _descriptors(u, direction::Symbol, retcode)
    direction == :both && return (lfwd = u.lfwd, lbwd = u.lbwd, retcode = retcode)
    direction == :forward && return (lfwd = u.lfwd, retcode = retcode)
    return (lbwd = u.lbwd, retcode = retcode)
end

# The augmented problem and the output of each trajectory of the ensemble, for an
# infinitesimal Lagrangian descriptor `M`; see `increment_descriptor.jl` for another method
function _ensemble_setup(prob, M, direction::Symbol)
    output_func = function (sol, i)
        (_descriptors(sol.u[end], direction, sol.retcode), false)
    end
    return augmentprob(prob, M; direction), output_func
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

- `ensprob`: the `EnsembleProblem` solved by [`solve`](@ref). Its trajectories solve the
  augmented problem built by [`augmentprob`](@ref), and their outputs are the Lagrangian
  descriptors at the end of the time interval, along with the solver's return codes.
- `uu0`: the given collection of initial conditions.
- `direction`: the given direction.

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
struct LagrangianDescriptorProblem{T1,T2,T3}
    ensprob::T1
    uu0::T2
    direction::T3

    function LagrangianDescriptorProblem(prob, M, uu0; direction::Symbol = :both)
        _check_direction(direction)
        # with several realizations (see `solve`), trajectory `i` starts at `uu0[mod1(i, n)]`
        prob_func = function (augprob, args...)
            i = _trajectory_index(args...)
            remake(augprob, u0 = _augmented_u0(uu0[mod1(i, length(uu0))], direction))
        end
        augprob, output_func = _ensemble_setup(prob, M, direction)
        ensprob = EnsembleProblem(augprob, prob_func = prob_func, output_func = output_func)
        new{typeof(ensprob),typeof(uu0),typeof(direction)}(ensprob, uu0, direction)
    end
end
