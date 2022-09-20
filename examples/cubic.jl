# # Lagrangian descriptors for a scalar cubic ordinary differential equation

# We apply the method of Lagrangian descriptors to the simple cubic equation
# $$
#   \frac{\mathrm{d}x}{\mathrm{d}t} = x - x^3
# $$
#
# This equation has two stationary solutions, corresponding to the fixed points at $x=0$ and
# $x = 1$, i.e. $x(t) = 0, \forall t$, and $x(t) = 1, \forall t$. The remaining
# non-stationary solutions can be obtained by separation of variables and partial fractions
# explicit solution obtained by separation of variables, integrated by partial fractions:

using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using DiffEqBase

include("../src/LagrangianDescriptors.jl")
using .LagrangianDescriptors
using .LagrangianDescriptors: augmentprob

# Simple cubic equation u' = u - u^3

f(u, p, t) = u - u^3

u0 = 0.5
tspan = (0.0, 5.0)
prob = ODEProblem(f, u0, tspan)

# sol = solve(prob, Tsit5())

#plot(sol)

# augmented

M(du, u, p, t) = sum(abs2, du) # local descriptor

#augprob = LagrangianDescriptors.augmentprob(prob, M)

#augsol = solve(augprob, Tsit5())

#plot(augsol)
#= 
plot(augsol, idxs=(0,1))
plot!(augsol, idxs=(0,2))
plot(augsol, idxs=(0,3))
plot!(augsol, idxs=(0,4)) =#

#

uu0 = range(-1.0, 1.0, length=101)
lagprob = LagrangianDescriptorProblem(prob, M, uu0)
# solve(lagprob.ensprob, Tsit5(), trajectories=length(uu0))
lagsol = solve(lagprob, Tsit5())

plot(lagsol)
plot(lagsol, :forward)
plot(lagsol, :backward)

#plot(uu0, getindex.(lagsol.enssol.u,:lfwd), label="forward Lagrangian descriptor")
#plot!(uu0, getindex.(lagsol.enssol.u,:lbwd), label="backward Lagrangian descriptor")
#plot!(uu0, sum.(lagsol.enssol.u), label="Lagrangian descriptor")
#plot!(uu0, getindex.(lagsol.enssol.u, :lfwd) - getindex.(lagsol.enssol.u, :lbwd), label="Difference lfwd - lbwd")

plot(uu0, lagsol(:forward), label="forward", title="Lagrangian descriptors", titlefont=10)
plot!(uu0, lagsol(:backward), label="backward")
plot!(uu0, lagsol(), label="total")
plot!(uu0, lagsol(:difference), label="difference")