# Double gyre

The double gyre is a time-dependent model of a fluid flow, with two counter-rotating vortices in the rectangle ``[0, 2] \times [0, 1]``, separated by a boundary that oscillates periodically in time. It is a classical example in the study of transport in fluids [ShaddenEtAl2005](@cite). The velocity field is given by the stream function ``\psi(x, y, t) = A \sin(\pi f(x, t)) \sin(\pi y)``, with ``f(x, t) = a(t) x^2 + b(t) x``, ``a(t) = \varepsilon \sin(\omega t)`` and ``b(t) = 1 - 2a(t)``:

```math
\begin{cases}
  \dot x = -\dfrac{\partial \psi}{\partial y} = -\pi A \sin(\pi f(x, t)) \cos(\pi y), \\[1ex]
  \dot y = \dfrac{\partial \psi}{\partial x} = \pi A \cos(\pi f(x, t)) \sin(\pi y) \dfrac{\partial f}{\partial x}(x, t).
\end{cases}
```

```@example doublegyre
using LagrangianDescriptors, OrdinaryDiffEq, Plots
using LinearAlgebra: norm

function doublegyre!(du, u, p, t)
    A, ε, ω = p
    x, y = u
    a = ε * sin(ω * t)
    b = 1 - 2a
    f = a * x^2 + b * x
    dfdx = 2a * x + b
    du[1] = -π * A * sin(π * f) * cos(π * y)
    du[2] = π * A * cos(π * f) * sin(π * y) * dfdx
end

p = (0.1, 0.25, 2π / 10) # A, ε, ω
prob = ODEProblem(doublegyre!, [0.0, 0.0], (0.0, 15.0), p)

M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(0, 1, length = 201), x in range(0, 2, length = 401)]

lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())

plot(lagsol, xlabel = "x", ylabel = "y", size = (900, 450))
savefig("doublegyre.png"); nothing # hide
```

![Lagrangian descriptors of the double gyre](doublegyre.png)

Since the flow is time dependent, the picture is that at the initial time ``t_0 = 0``. The sharp curves are the stable and unstable manifolds of the hyperbolic trajectories near the oscillating boundary between the gyres. They intersect in a sequence of lobes, which carry fluid from one gyre to the other in each period.
