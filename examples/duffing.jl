using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using DiffEqBase

include("../src/LagrangianDescriptors.jl")
using .LagrangianDescriptors
using .LagrangianDescriptors: augmentprob

# Periodically-forced Duffing system

function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * cos(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 100.0)
A = 0.3; ω = π

p = (A, ω)
prob = ODEProblem(f!, u0, tspan, p)

sol = solve(prob, Tsit5())

plot(sol)

plot(sol, idxs=(1,2))

# augmented

M(du, u, p, t) = sum(abs2, du)

augprob = augmentprob(prob, M)

augsol = solve(augprob, Tsit5())

plot(augsol)

plot(augsol, idxs=(0,1,2))
plot!(augsol, idxs=(0,3,4))
plot(augsol, idxs=5)
plot!(augsol, idxs=6)

uu0 = [[x, y] for x in range(-1.2, 1.2, length=251), y in range(-1.0, 1.0, length=201)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)
# solve(lagprob.ensprob, Tsit5(), trajectories=length(uu0))
lagsol = solve(lagprob, Tsit5())

lsol = fill(0.0, size(uu0)...)
for i in eachindex(lagsol.enssol.u)
    lsol[i] = lagsol.enssol.u[i].lfwd
end

heatmap(getindex.(uu0[:,1], 1), getindex.(uu0[1,:], 2), lsol)


plot(lagsol)
plot(lagsol, :forward)
plot(lagsol, :backward)
