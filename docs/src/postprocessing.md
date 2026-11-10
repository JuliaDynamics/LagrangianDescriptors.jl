# Alternative via post-processing

Alternatively, instead of augmenting the system and computing the descriptor along with the solutions, we can just solve the original system as usual and then integrate the infinitesimal descriptor over the solution. This can be done for a single solution, backward and/or forward in time, or for an ensemble of solutions of the original system.

## Drawbacks

This is much simpler to implement, *but* it has some drawbacks:

1. First, we will need to save each trajectory in full, for later post-processing, instead of only ``L_\mathrm{fwd}(t_f)`` and ``L_\mathrm{bwd}(t_f)``. Hence, this is much more memory demanding. Keep in mind we need to solve for a lot of trajectories.
1. If we want the backward Lagrangian descriptors, we also need to set up and solve the system backward.
1. Some solutions may have some spread out time steps saved by the solver. We either have to force it to save on a fine time mesh for a faster and accurate integration, or we use the interpolation present in the solution, which is usually much slower to compute.

## Example

Here we exemplify the idea by integrating a single solution of a scalar cubic equation using [JuliaMath/QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl).

It boils down to taking a solution `sol` of the `ODEProblem` and integrating the infinitesimal descriptor `M` with the appropriate arguments. In the case of a scalar ODE with an out-of-place right hand side `f=f(u, p, t)` and an infinitesimal descriptor `M=M(u, p, t)`, this is simply

```julia
M(t) = M(sol(t), sol.prob.p, t)
first(quadgk(M, minimum(tspan), maximum(tspan)))
```

For a more generic implementation for `ODESolution` types, in either in-place or out-of-place cases, and for either backward or forward descriptors, one can define a function such as

```julia
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

A version of this, which also avoids some allocations, is available in the package as [`lagrangian_descriptor`](@ref). With it, we compare the two different approaches, making sure they yield the same result, up to approximation errors:

```@example postprocessing
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

tspanbwd = (t0, 2t0 - tf)
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
nothing # hide
```

## Benchmark

And here is a benchmark with the above setup, with a single solution of the augmented system versus post-processing:

```@example postprocessing
using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1 # hide

benchmarks = [
    "Create forward ODE problem" => @benchmark(ODEProblem($f!, $u0, $tspanfwd)),
    "Remake for backward ODE problem" => @benchmark(remake($probfwd, tspan = $tspanbwd)),
    "Solve forward ODE problem" => @benchmark(solve($probfwd, $(Tsit5()))),
    "Solve backward ODE problem" => @benchmark(solve($probbwd, $(Tsit5()))),
    "Create augmented ODE problem" => @benchmark(augmentprob($probfwd, $M)),
    "Solve augmented ODE problem" => @benchmark(solve($augprob, $(Tsit5()))),
    "Post-processing for forward Lagrangian descriptor" => @benchmark(lagrangian_descriptor($solfwd, $M)),
    "Post-processing for backward Lagrangian descriptor" => @benchmark(lagrangian_descriptor($solbwd, $M)),
]

Text(join((rpad(name, 52) * lpad(BenchmarkTools.prettytime(minimum(b).time), 12) *
    lpad("$(minimum(b).allocs) allocations", 18) * lpad(BenchmarkTools.prettymemory(minimum(b).memory), 12)
    for (name, b) in benchmarks), "\n"))
```


We see that solving the augmented system, with the forward and backward evolutions together, takes about as long as solving the forward and backward equations separately, while it also computes the Lagrangian descriptors.

On the other hand, solving the forward and backward equations separately requires a post-processing step for each of the forward and backward evolutions to obtain the Lagrangian descriptors, and that takes considerably longer, with substantially more allocations and memory. And this was done for a single trajectory. Imagine it for an ensemble of solutions, on top of the memory demand of saving the full solutions.
