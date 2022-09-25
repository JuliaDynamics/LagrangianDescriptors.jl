# LagrangianDescriptors.jl documentation

## About

The dynamics of evolutionary systems can be quite intricate. The method of **Lagrangian Descriptors** helps to visualize the complicate behavior of such systems and to make sense of it. In a recent article, Wiggins and and García-Garrido call it [a method for] *painting the phase portrait* (of a dynamical system) (see [S. Wiggins and V. J. García-Garrido, Painting the Phase Portrait of a Dynamical System with the Computational Tool of Lagrangian Descriptors (AMS Notices, June/July 2022)](https://www.ams.org/journals/notices/202206/noti2489/noti2489.html?adat=June/July%202022&trk=2489&galt=none&cat=feature&pdfissue=202206&pdffile=rnoti-p936.pdf).

The image below, for instance, shows the dynamics of a periodically-forced Duffing equation, with a particular combination of parameters and near time $t=0$ (see [Tutorial: Periodically-forced Duffing equation](tutorial_ODEs.md#Periodically-forced-Duffing-equation)):

![Duffing](img/duffing.png)

## The method

The method is akin to droping colored ink in a fluid flow, tracking the dye as it is transported by the flow, and revealing the pattern created after a certain period of time. The difference being that the color doesn't get diffused as in a real fluid, so the image doesn't get blurred and one gets a clearer picture.

It is similar to drawing a phase portrait, which displays a collection of orbits, but in this method each orbit is painted according to its dynamic behavior, better revealing the overall picture.

This "coloring" is obtained by integrating a **local,** or **infinitesimal descriptor** $M=M(t, u)$ along a solution $u=u(t)$ of the system. One common choice for the infinitesimal descriptor is the velocity of the trajectory, so that, upon integration along the solution over a given period of time, the Lagrangian descriptor becomes the arc-length of that portion of the solution. Depending on whether the solution is at a fixed point, or belongs to a stable or unstable manifold, and so on, the arc-length would be closer for solutions with similar behavior, revealing the common patterns of the dynamics.

Thus, if $u=u(t)$ is a solution over a time span $(t_0, t_f)$, with $t_f > t_0$, then the **forward Lagrangian descriptor** is given by

```math
L_{\mathrm{fwd}} = \int_{t_0}^{t_f} M(t, u(t)) \;\mathrm{d}t.
```

If $t_f < t_0$, so that $u=u(t)$ is a "backward" solution from $t_0$, then we obtain the **backward Lagrangian descriptor**

```math
L_{\mathrm{bwd}} = \int_{t_f}^{t_0} M(t, u(t)) \;\mathrm{d}t = - \int_{t_0}^{t_f} M(t, u(t)) \;\mathrm{d}t.
```

Since the system might be non-autonomous, these Lagragians descriptors are referred as descriptors **near ``t_0``**.

## The computation

Notice the computation of the Lagrangian descriptors only depend on a given solution and on the infinitesimal descriptor, so that they can be computed *a posteriori*. However, this is not the most efficient way of computing them. For a more efficient implementation, one writes the integrals for the Lagrangian descriptors as differential equations, i.e.

```math
\frac{\mathrm{d}L_{\mathrm{fwd}}}{\mathrm{d}t} = M(t, u(t))
```
and similarly for the backward Lagrangian descriptor. Then, one solves the (partially) coupled system to obtain the integrated descriptor.

More explicitly, when both forward and backward Lagrangians are desired, we may write a four-component system for solving both the backward and forward solutions and the backward and forward Lagrangians at the same time:

```math
\begin{cases}
  \displaystyle \frac{\mathrm{d}u}{\mathrm{d}t} = f(u, t), \\ \\
  \displaystyle \frac{\mathrm{d}v}{\mathrm{d}t} = f(u, 2t_0 - t), \\ \\
  \displaystyle \frac{\mathrm{d}L_{\mathrm{fwd}}}{\mathrm{d}t} = M(u, t), \\ \\
  \displaystyle \frac{\mathrm{d}L_{\mathrm{bwd}}}{\mathrm{d}t} = M(v, 2t_0 - t),
\end{cases}
```
with the set of initial conditions
```math
u(t_0) = v(t_0) = u_0, \quad L_{\mathrm{fwd}}(t_0) = L_{\mathrm{bwd}}(t_0) = 0.
```

For a given time interval ``(t_0, t_f)``, with ``t_f > t_0``, notice that ``u=u(t)`` solves the system forward, over the interval ``(t_0, t_f)``, while ``v=v(t)`` solves the system backward, over the interval ``(2t_0 - t_f, t_0) = (t_0 - (t_f - t_0), t_0)``, both with the same initial condition ``u_0``, i.e. the compute the forward and backward parts of the same trajectory, over the same length of time, but in different directions.

## Implementation

The implementation works by 
1. Taking a *differential equation problem* of a type defined by [SciML/DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), with a given time span ``(t_0, t_f)``;
1. Taking an *infinitesimal descriptor* ``M=M(du, u, p, t)`` (or other form suitable to the given problem type) that will be integrated along a solution ``u(t) = u(t; u_0)``, to yield the forward Lagrangian descriptor ``L_{\mathrm{fwd}}(u_0) = \int_{t_0}^{t_f} M(du(t), u(t), p, t)\;\mathrm{d}t`` and, similarly, the backward Lagrangian descriptor ``L_{\mathrm{bwd}}(u_0) = \int_{t_0}^{t_f} M(du(-t), u(-t), p, 2t_0 - t)\;\mathrm{d}t``, for a given initial condition ``u_0``;
1. Generating an *augmented* problem of the same time and with four components, for solving the original equation forward and backward in time, and for solving the Lagrangian descriptors forward and backward in time. The way the augmented systems is constructed is via [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl), so the augmented system is a `ComponentVector`, with components `fwd`, `bwd`, `lfwd`, and `lbwd`, respectively.

1. Creating a `LagrangianDescriptorProblem` by wrapping an [EnsembleProblem](https://diffeq.sciml.ai/dev/features/ensemble/) for the augmented system and with a given collection ``uu_0`` of initial conditions to be used at each trajectory of the ensemble;
1. Solving the wrapped ensemble problem and returning a `LagrangianDescriptorSolution` containing the associated collection of (forward and backward) Lagrangian descriptor values at the final time of each simulation (which is ``t_f`` for the forward components and corresponds to ``2t_0 - t_f`` for the backward ones);
1. Finally, one can visualize either the forward, or the backward, or the sum, or even the difference, of the forward and backward Lagrangian descriptor with a plot recipe for the `LagrangianDescriptorSolution`.

## Current state

The package is still in an embrionary phase and currently accepts differential equations of the type `ODEProblem`. Problems like `SDEProblem` and `RODEProblem` will be implemented soon. Other problems will come eventually.

The plot recipe works for some types of collections of initial conditions (e.g. a `AbstractVector{<:Number}` for scalar problems and `AbstractMatrix{<:AbstractVector{<:Number}}` for two-dimensional problems). More general and flexible plot recipes will also be implemented.

## Developers

[JuliaDynamics/LagrangianDescriptors.jl](https://github.com/JuliaDynamics/LagrangianDescriptors.jl) is currently being developed by [Ricardo M. S. Rosa](https://rmsrosa.github.io), but contributors are welcome.

## Cite

Just cite the github repo [JuliaDynamics/LagrangianDescriptors.jl](https://github.com/JuliaDynamics/LagrangianDescriptors.jl) for now, while the package is not yet registered.
