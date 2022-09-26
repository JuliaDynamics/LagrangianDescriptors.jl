# Alternative via post-processing

Alternatively, instead of augmenting the system and computing the descriptor along with the solutions, we can just solve the original system as usual and then integrate the infinitesimal descriptor over the solution. This can be done for a single solution, backward and/or forward in time, or for an ensemble of solutions of the original system.

## Drawbacks

This is much simpler to implement, *but* it has some drawbacks:

1. First, we will need to save each trajectory in full, for later post-processing, instead of only ``L_\mathrm{fwd}(t_f)`` and ``L_\mathrm{bwd}(t_f)``. Hence, this is much more memory demanding. Keep in mind we need to solve for a lot of trajectories.
1. If we want the backward Lagrangian descriptors, we also need to set up and solve the system backward.
1. Some solutions may have some spread out time steps saved by the solver. We either have to force it to save on a fine time mesh for an faster and accurate integration, or we use the interpolation present in the solution, which is usually much slower to compute.

## Example

Here we exemplify the idea by integrating a single solution of a scalar cubic equation using [JuliaMath/QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl).

```julia postprocessing
using OrdinaryDiffEq, Test
using LagrangianDescriptors
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm
using QuadGK: quadgk
using BenchmarkTools: @btime

@testset "Augmented vs post-processing scalar ODE" begin
    f = function (u)
        u - u^3
    end
    f! = function (du, u, p, t)
        du .= f.(u)
    end
    t0 = 0.0
    tf = 5.0
    tspan = (t0, tf)
    n = 5
    u0 = 0.1 .+ 0.7rand(n)
    prob = ODEProblem(f!, u0, tspan)
    @info "Create forward ODE problem:"
    @btime ODEProblem($f!, $u0, $tspan)
    tspanbwd = (tf, t0)
    probbwd = remake(prob, tspan = (tf, t0))
    @info "Remake for backward ODE problem:"
    @btime remake($prob, tspan = $((tf, t0)))
    solfwd = @test_nowarn solve(prob, Tsit5())
    solbwd = @test_nowarn solve(probbwd, Tsit5())
    @info "Solve forward ODE problem:"
    @btime solve($prob, $(Tsit5()))
    @info "Solve backward ODE problem:"
    @btime solve($probbwd, $(Tsit5()))

    M = function (du, u, p, t)
        norm(du)
    end
    augprob = @test_nowarn augmentprob(prob, M)
    @info "Create augmented ODE problem:"
    @btime augmentprob($prob, $M)
    augsol = @test_nowarn solve(augprob, Tsit5())
    @info "Solve Augmented ODE problem:"
    @btime solve($augprob, $(Tsit5()))

    postproc = function (sol, f, M, tspan)
        first(quadgk(t -> M(f.(sol(t)), nothing, nothing, nothing), first(tspan), last(tspan)))
    end
    @test augsol.u[end].lfwd ≈ postproc(solfwd, f, M, tspan) rtol = 0.01
    @info "Postprocessing for forward Lagrangian descriptor:"
    @btime $postproc($solfwd, $f, $M, $tspan)
    @test augsol.u[end].lbwd ≈ postproc(solbwd, f, M, tspan) atol = 0.01
    @info "Postprocessing for backward Lagrangian descriptor:"
    @btime $postproc($solbwd, $f, $M, $tspan)
end
```

Here are the results of the tests comparing the results of the package with the results obtained computing the Lagrangian descriptor directly via QuadGK, just to make sure everything is being computed properly:

```zsh
Test Summary:                           | Pass  Total   Time
Augmented vs post-processing scalar ODE |    6      6  53.4s
Test.DefaultTestSet("Augmented vs post-processing scalar ODE", Any[], 6, false, false, true, 1.664197063101466e9, 1.664197116484559e9)
```

## Benchmark

And here is the result of the above benchmark, with a single solution of the augmented system versus post-processing:

```zsh
[ Info: Create forward ODE problem:
  6.308 μs (66 allocations: 3.02 KiB)
[ Info: Remake for backward ODE problem:
  2.125 ns (0 allocations: 0 bytes)
[ Info: Solve forward ODE problem:
  6.408 μs (146 allocations: 15.31 KiB)
[ Info: Solve backward ODE problem:
  6.417 μs (146 allocations: 15.31 KiB)
[ Info: Create augmented ODE problem:
  10.166 μs (130 allocations: 6.38 KiB)
[ Info: Solve Augmented ODE problem:
  15.167 μs (249 allocations: 32.25 KiB)
[ Info: Postprocessing for forward Lagrangian descriptor:
  54.833 μs (933 allocations: 88.88 KiB)
[ Info: Postprocessing for backward Lagrangian descriptor:
  30.083 μs (513 allocations: 49.50 KiB)
```

We see that solving both forward and backward equations separately is a bit faster than solving the augmented system with both forward and backward evolutions together, but the latter also includes the computations of the Lagrangian descriptors. The number of allocations and memory are about the same.

On the other hand, solving the forward and backward equations separately requires a post-processing step for each forward and backward evolutions to obtain the Lagrangian descriptors, and that takes quite a longer time and substantially more allocations and memory. And this was done for a single trajectory. Imagine for the ensemble of solutions, on top of the memory demand of saving the full solutions.
