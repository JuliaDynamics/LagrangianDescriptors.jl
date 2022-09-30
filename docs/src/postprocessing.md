# Alternative via post-processing

Alternatively, instead of augmenting the system and computing the descriptor along with the solutions, we can just solve the original system as usual and then integrate the infinitesimal descriptor over the solution. This can be done for a single solution, backward and/or forward in time, or for an ensemble of solutions of the original system.

## Drawbacks

For individual solutions, this is much simpler to implement, *but* it has some drawbacks:

1. First, we will need to save each trajectory in full, for later post-processing, instead of only ``L_\mathrm{fwd}(t_f)`` and ``L_\mathrm{bwd}(t_f)``. Hence, this is much more memory demanding. Keep in mind we need to solve for a lot of trajectories.
1. If we want the backward Lagrangian descriptors, we also need to set up and solve the system backward.
1. Some solutions may have some spread out time steps saved by the solver. We either have to force it to save on a fine time mesh for an faster and accurate integration, or we use the interpolation present in the solution, which is usually much slower to compute.

## Example of post-processing a single solution

Here we exemplify the idea by integrating a single solution of a scalar cubic equation using [JuliaMath/QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl).

It boils down to taking a solution `sol` of the `ODEProblem` and integrating the infinitesimal descriptor `M` with the apropriate arguments. In the case of a scalar ODE with an out-of-place right hand side `f=f(u, p, t)` and an infinitesimal descriptor `M=M(u, p, t)`, this is simply

```julia
M(t) = M(sol(t), sol.prob.p, t)
first(quadgk(M, minimum(tspan), maximum(tspan)))
```

For a more generic implementation for `ODESolution` types, in either in-place or out-of-place cases, and for either backward of forward descriptors, we define the function

```julia postprocessing
function lagrangian_descriptor(sol::ODESolution, M)
    t0, tf = extrema(sol.prob.tspan)
    integrand = isinplace(sol.prob) ?
        function (t, du = similar(sol.prob.u0))
            sol.prob.f(du, sol(t), sol.prob.p, t)
            M(du, sol(t), sol.prob.p, t)
        end :
        function (t)
            du = sol.prob.f(sol(t), sol.prob.p, t)
            M(du, sol(t), sol.prob.p, t)
        end
    return first(quadgk(integrand, t0, tf))
end
```

With this definition, we compare the two different approaches, making sure they yield the same result, up to approximation errors:

```julia postprocessing
using OrdinaryDiffEq, Test
using LagrangianDescriptors
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm
using QuadGK: quadgk

f(u) = u - u^3
f!(du, u, p, t) = (du .= f.(u))

t0 = 0.0
tf = 5.0
n = 5
u0 = 0.1 .+ 0.7rand(n)

tspanfwd = (t0, tf)
probfwd = ODEProblem(f!, u0, tspanfwd)

tspanbwd = (tf, t0)
probbwd = remake(probfwd, tspan = tspanbwd)

solfwd = solve(probfwd, Tsit5())
solbwd = solve(probbwd, Tsit5())

M(du, u, p, t) = norm(du)

augprob = augmentprob(probfwd, M)
augsol = solve(augprob, Tsit5())

@testset "Augmented vs post-processing" begin
    @test augsol.u[end].lfwd ≈ lagrangian_descriptor(solfwd, M) rtol = 0.01
    @test augsol.u[end].lbwd ≈ lagrangian_descriptor(solbwd, M) atol = 0.01
end
```

```zsh
Test Summary:                | Pass  Total  Time
Augmented vs post-processing |    2      2  0.2s
Test.DefaultTestSet("Augmented vs post-processing", Any[], 2, false, false, true, 1.664213738359913e9, 1.664213738562828e9)
```

## Benchmark

And here is a benchmark with the above setup, with a single solution of the augmented system versus post-processing:

```julia postprocessing
using BenchmarkTools: @btime

@info "Create forward ODE problem:"
@btime ODEProblem($f!, $u0, $tspanfwd)

@info "Remake for backward ODE problem:"
@btime remake($probfwd, tspan = $tspanbwd)

@info "Solve forward ODE problem:"
@btime solve($probfwd, $(Tsit5()))

@info "Solve backward ODE problem:"
@btime solve($probbwd, $(Tsit5()))
    
@info "Create augmented ODE problem:"
@btime augmentprob($probfwd, $M)
    
@info "Solve Augmented ODE problem:"
@btime solve($augprob, $(Tsit5()))

@info "Postprocessing for forward Lagrangian descriptor:"
@btime $lagrangian_descriptor($solfwd, $M)

@info "Postprocessing for backward Lagrangian descriptor:"
@btime $lagrangian_descriptor($solbwd, $M)
```

```zsh
[ Info: Create forward ODE problem:
  4.327 μs (59 allocations: 2.34 KiB)
[ Info: Remake for backward ODE problem:
  2.084 ns (0 allocations: 0 bytes)
[ Info: Solve forward ODE problem:
  6.492 μs (146 allocations: 15.31 KiB)
[ Info: Solve backward ODE problem:
  6.508 μs (146 allocations: 15.31 KiB)
[ Info: Create augmented ODE problem:
  8.722 μs (130 allocations: 6.38 KiB)
[ Info: Solve Augmented ODE problem:
  22.916 μs (388 allocations: 38.05 KiB)
[ Info: Postprocessing for forward Lagrangian descriptor:
  97.750 μs (1488 allocations: 140.91 KiB)
[ Info: Postprocessing for backward Lagrangian descriptor:
  44.375 μs (677 allocations: 63.64 KiB)
1.0964245933425085
```

We see that solving both forward and backward equations separately is a bit faster than solving the augmented system with both forward and backward evolutions together, but the latter also includes the computations of the Lagrangian descriptors. The number of allocations and memory are about the same.

On the other hand, solving the forward and backward equations separately requires a post-processing step for each forward and backward evolutions to obtain the Lagrangian descriptors, and that takes quite a longer time and substantially more allocations and memory. And this was done for a single trajectory. Imagine for the ensemble of solutions, on top of the memory demand of saving the full solutions.

## Ensemble post-processing

Just for the sake of completeness, we also implemented an ensemble version of the post-processing approach. This uses the same API as before, but with the keyword `method = :postprocess` in `LagrangianDescriptorProblem`.

```julia pp
using OrdinaryDiffEq
using Plots
using LinearAlgebra: norm
using LagrangianDescriptors
```

```julia pp
function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * sin(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 13.0)
A = 0.3; ω = π; p = (A, ω)

prob = ODEProblem(f!, u0, tspan, p)
```

```julia pp
M(du, u, p, t) = norm(du)

uu0 = [[x, y] for y in range(-1.0, 1.0, length=101), x in range(-1.8, 1.8, length=101)]

lagprob = LagrangianDescriptorProblem(prob, M, uu0, method=:postprocessed)
```

```julia pp
@time lagsol = solve(lagprob, Tsit5())
```

```julia pp
plot(lagsol, direction=:both)
```
