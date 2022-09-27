# Social preview

Here is the code to draw the social preview image:

```julia social-preview
using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using LagrangianDescriptors

function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * sin(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 13.0)
A = 5.0; ω = 2π; p = (A, ω);
prob = ODEProblem(f!, u0, tspan, p)

M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1.5, -0.5, length=301), x in range(-0.4, 1.6, length=601)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)

lagsol = solve(lagprob, Tsit5())
```

```julia social-preview
plot(lagsol, :forward, size=(1280, 640), colorbar=false, axes=false, ticks=false)

annotate!([(0.6,-1.1,text("LagrangianDescriptors.jl",68, :white, :center, "Times")),(0.6,-1.25, text("Painting the phase portrait of random and deterministic systems", 26, :white, :center, "Times"))])
```

```julia social-preview
savefig("img/LagrangianDescriptors_socialpreview.png")
```
