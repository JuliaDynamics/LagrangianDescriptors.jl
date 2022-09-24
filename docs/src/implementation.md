# Implementation

Here are the two initial ideas for implementing such method. We ended up implementing the first one.

## Augmenting the system

We augment the system and compute the descriptors along with the solution.

1. One builds a problem `prob` of a given type from `SciMLBase.jl`, say an `ODEProblem` for `du/dt = f(u, p, t)`.
1. Then we pass it to `ldprob = LagrangianDescriptorProblem(prob, M, uu0)`, where `M = M(du, u, p, t)` is the (infinitesimal) descriptor, which is a scalar function, e.g. `M(du, u, p, t) = norm(du)`, and `uu0` is some iterator with a collection of initial conditions (e.g. an `Array` for a mesh in phase space or a portion of a sub-manifold of the phase space). 
1.  `LagrangianDescriptorProblem` uses `prob.f`, `prob.p`, `prob.u0`, and `prob.tspan` to create, via `ComponentArrays`, a new `ODEProblem` for an augmented system of the form
```math
\begin{cases}
\displaystyle \frac{du}{dt} = f(u, p, t) \\ \\
\displaystyle \frac{dv}{dt}  = -f(v, p, 2t_0-t) \\ \\
\displaystyle \frac{dL_f}{dt}  = M(u, p, t) \\ \\
\displaystyle \frac{dL_b}{dt}  = M(v, p, t) \\
\end{cases}
```
1. Notice $v$ solves the system backwards. If `tspan = (t0, tf)`, then $v$ solves it backwards in the interval `(2t0 - tf, t0) = (t0 - (tf - t_0), t_0)`. So, we solve the system forwards and backwards at the same time.
1. The forward and backward Lagrangian descriptors `L_f` and `L_b` are the forward and backward integrations of the infinitesimal descriptor `M`.
1. Then, solving a `LagrangianDescriptorProblem` works via an `EnsembleProblem`, where at each new solve, a new initial condition is picked.
1. At the end of each of those solves, we only need to save the values of `Lf[end]` and `Lb[end]`.
1. We can visualize the Lagrangian descriptor using a heatmap of `Lf[end]`, `Lb[end]` or `Lf[end] + Lb[end]`.
1. We can also add a flag to build only the forward descriptor `Lf` or the backward descriptor `Lb`.

## Post-processing 

Alternatively, instead of augmenting the system and computing the descriptor along the solutions, we can just solve an ensemble of solutions of the original system and then integrate the descriptor over each solution.

This is simpler, *but* it has potential drawbacks:

1. First, we will need to save each trajectory in full, instead of only `Lf[end]` and `Lb[end]`, so this is more memory demanding. Keep in mind we need to solve for a lot of trajectories.

2. Secondly, some solutions may have some spread out time steps. We either have to force it to save on a fine time mesh, for better time integration, or we use the interpolation present in the solution, which might be slower to compute.

~~Anyway, the plan is to implement both approaches and do some performance comparisons. Maybe we keep both methods and leave a flag to choose between the two.~~

I have decided against implementing this post-processing approach. As I guessed, this would be computationally demanding, as seen in the following benchmark:
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

## References

* [Painting the Phase Portrait of a Dynamical System with the Computational Tool of Lagrangian Descriptors](https://www.ams.org/journals/notices/202206/noti2489/noti2489.html?adat=June/July%202022&trk=2489&galt=none&cat=feature&pdfissue=202206&pdffile=rnoti-p936.pdf)
* [Lagrangian descriptors: A method for revealing phase space structures of general time dependent dynamical systems](https://www.sciencedirect.com/science/article/abs/pii/S1007570413002037)
* [Lagrangian Descriptors - *Discovery and Quantification of Phase Space Structure and Transport*](https://champsproject.github.io/lagrangian_descriptors/docs/authors.html)
* [Frequently Asked Questions about Lagrangian Descriptors](https://acp.copernicus.org/preprints/acp-2016-633/acp-2016-633-SC2-supplement.pdf)
