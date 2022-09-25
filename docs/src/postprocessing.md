# Alternative via post-processing

Alternatively, instead of augmenting the system and computing the descriptor along the solutions, we can just solve an ensemble of solutions of the original system and then integrate the descriptor over each solution.

This is simpler, *but* it has potential drawbacks:

1. First, we will need to save each trajectory in full, instead of only `Lf[end]` and `Lb[end]`, so this is more memory demanding. Keep in mind we need to solve for a lot of trajectories.

2. Secondly, some solutions may have some spread out time steps. We either have to force it to save on a fine time mesh, for better time integration, or we use the interpolation present in the solution, which might be slower to compute.

Here is the result of a simple benchmark with a single solution of the augmented system versus post-processing.

```zsh
[ Info: Create forward ODE problem:
  10.083 μs (108 allocations: 5.00 KiB)
[ Info: Remake for backward ODE problem:
  1.708 ns (0 allocations: 0 bytes)
[ Info: Solve forward ODE problem:
  9.556 μs (138 allocations: 14.95 KiB)
[ Info: Solve backward ODE problem:
  9.514 μs (138 allocations: 14.95 KiB)
[ Info: Create augmented ODE problem:
  16.209 μs (175 allocations: 8.41 KiB)
[ Info: Solve Augmented ODE problem:
  23.625 μs (249 allocations: 32.25 KiB)
[ Info: Postprocessing for forward Lagrangian descriptor:
  68.125 μs (886 allocations: 85.62 KiB)
[ Info: Postprocessing for backward Lagrangian descriptor:
  31.834 μs (405 allocations: 39.30 KiB)
```

We see that solving both forward and backward equations is a bit faster than solving the augmented system with both forward and backward evolutions together, but the latter also includes the computations of the Lagrangian descriptors. On the other hand, solving the forward and backward equations separately requires a post-processing step for each forward and backward evolutions to obtain the Lagrangian descriptors, and that takes quite a long time. And this was done for a single trajectory. Imagine for the ensemble of solutions, on top of the memory demand. It is not worth it. I will include a separate "post-processing" method for the sake of debugging, development, and comparison, but it will not be part of the main API.
