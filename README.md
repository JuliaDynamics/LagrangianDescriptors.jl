# Lagrangian Descriptors

Implementation of the method of **Lagrangian Descriptors** to highlight singular features (e.g. stable or unstable invariant manifolds) of the dynamics of an evolutionary system (such as ordinary or partial differential equations, random equations, or stochastic differential equations).

Here are some examples on a periodically forced Duffing system.

![Duffing 1](examples/img/duffing1.png)

![Duffing 2](examples/img/duffing2.png)
## Idea

Here are the two initial ideas.

### Augmenting the system

We augment the system and compute the descriptors along with the solution.

1. One builds a problem `prob` of a given type from `SciMLBase.jl`, say an `ODEProblem` for `du/dt = f(u, p, t)`.
2. Then we pass it to `ldprob = LDProblem(prob, M, uu)`, where `M = M(u, p, t)` is the descriptor, which is a scalar function, e.g. `M(u, p, t) = norm(f(u, p, t))`, and `uu` is some iterator with a collection of initial conditions (e.g. an `Array` for a mesh in phase space or a portion of a sub-manifold of the phase space). 
3.  `LDProblem` uses `prob.f` and `prob.tspan` to create, via `ComponentArrays`,  a new `ODEProblem` for an augmented system of the form

$$
\begin{cases}
\displaystyle \frac{du}{dt} = f(u, p, t) \\ \\
\displaystyle \frac{dv}{dt}  = -f(v, p, 2t_0-t) \\ \\
\displaystyle \frac{dL_f}{dt}  = L(u,p,t) \\ \\
\displaystyle \frac{dL_b}{dt}  = L(v,p,t) \\
\end{cases}
$$

5. Notice $v$ solves the system backwards. If `tspan = (t0, tf)`, then $v$ solves it backwards in the interval `(2t0 - tf, t0)`, since `2t0 - tf = t0 - (tf - t0)`. So, we solve the system forwards and backwards at the same time.
6. $L_f$ and $L_b$ are the forward and backward integrations of the Lagrangian descriptor.
7. Then, solving an `LDProblem` works via an `EnsembleProblem`, where at each new solve, a new initial condition is picked.
8. At the end of each of those solves, we only need to save the values of `Lf[end]` and `Lb[end]`.
9. We can visualize the Lagrangian descriptor using a heatmap of `Lf[end]`, `Lb[end]` or `Lf[end] + Lb[end]`.
10. We can also add a flag to build only the forward descriptor `Lf` or the backward descriptor `Lb`.

### Post-processing 

Alternatively, instead of augmenting the system and computing the descriptor along the solutions, we can just solve an ensemble of solutions of the original system and then integrate the descriptor over each solution.

This is simpler, *but* it has some potential drawbacks:

1. First, we will need to save each trajectory in full, instead of only `Lf[end]` and `Lb[end]`, so this is more memory demanding. Keep in mind we need to solve for a lot of trajectories.

2. Secondly, some solutions may have some spread out time steps. We either have to force it to save on a fine time mesh, for better time integration, or we use the interpolation present in the solution, which might be slower to compute.

Anyway, the plan is to implement both approaches and do some performance comparisons. Maybe we keep both methods and leave a flag to choose between the two.

## References

* [Painting the Phase Portrait of a Dynamical System with the Computational Tool of Lagrangian Descriptors](https://www.ams.org/journals/notices/202206/noti2489/noti2489.html?adat=June/July%202022&trk=2489&galt=none&cat=feature&pdfissue=202206&pdffile=rnoti-p936.pdf)
* [Lagrangian descriptors: A method for revealing phase space structures of general time dependent dynamical systems](https://www.sciencedirect.com/science/article/abs/pii/S1007570413002037)
* [Frequently Asked Questions about Lagrangian Descriptors](https://acp.copernicus.org/preprints/acp-2016-633/acp-2016-633-SC2-supplement.pdf)

