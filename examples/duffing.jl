using OrdinaryDiffEq, Plots
using LinearAlgebra: norm

include("../src/LagrangianDescriptors.jl")
using .LagrangianDescriptors

# Periodically-forced Duffing system

function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * cos(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 13.0)
A = 0.3;
ω = π;

p = (A, ω)
prob = ODEProblem(f!, u0, tspan, p)

sol = solve(prob, Tsit5())

plot(sol)

plot(sol, idxs = (1, 2))

# augmented

M(du, u, p, t) = norm(du)

augprob = augmentprob(prob, M)

augsol = solve(augprob, Tsit5())

plot(augsol)

plot(augsol, idxs = (0, 1, 2))
plot!(augsol, idxs = (0, 3, 4))
plot(augsol, idxs = 5)
plot!(augsol, idxs = 6)

uu0 = [[x, y] for y in range(-1.0, 1.0, length = 301), x in range(-1.8, 1.8, length = 301)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)

@time lagsol = solve(lagprob, Tsit5())

plot(
    lagsol,
    title = "Lagrangian descriptors for the forced Duffing equation \$\\ddot x = x - x^3 + A\\sin(\\omega t)\$\nwith A=$A and ω=2π",
    titlefont = 8,
    xlabel = "\$x\$",
    ylabel = "\$\\dot x\$",
)

savefig("duffing.png")

plot(
    lagsol,
    :forward,
    title = "Forward Lagrangian descriptors for the forced Duffing with A=$A and ω=$ω",
    titlefont = 8,
)

savefig("duffing_forward.png")

plot(
    lagsol,
    :backward,
    title = "Backward Lagrangian descriptors for the forced Duffing with A=$A and ω=$ω",
    titlefont = 8,
)

savefig("duffing_backward.png")

A = 12.0;
ω = 2π;
p = (A, ω);
uu0 = [[x, y] for y in range(-0.5, 0.0, length = 501), x in range(-0.6, 0.1, length = 501)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)

@time lagsol = solve(lagprob, Tsit5());

plot(lagsol, title = "Lagrangian descriptors - Duffing with A=$A, ω=2π", titlefont = 8)

savefig("img/duffing2.png")
