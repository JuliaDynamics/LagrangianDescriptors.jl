using Pkg
Pkg.activate("examples")
include("../src/LagrangianDescriptors.jl")
using .LagrangianDescriptors
using .LagrangianDescriptors: augmentprob
using OrdinaryDiffEq, Plots
using LinearAlgebra: norm
using DiffEqBase

# Simple cubic "reaction-diffusion" equation u' = u - u^3

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
#= 
plot(uu0, getindex.(lagsol.u,:lfwd), label="forward Lagrangian descriptor")
plot!(uu0, getindex.(lagsol.u,:lbwd), label="backward Lagrangian descriptor")
plot!(uu0, sum.(lagsol.u), label="Lagrangian descriptor")
plot!(uu0, getindex.(lagsol.u, :lfwd) - getindex.(lagsol.u, :lbwd), label="Difference lfwd - lbwd")

function prob_func(prob,i,repeat)
    uu0 = range(-1.0, 1.0, length=101)
    remake(prob,u0 = ComponentVector(fwd=uu0[i], bwd=uu0[i], lfwd=0.0, lbwd=0.0))
end

function output_func(sol, i)
    ComponentArray(lfwd = last(sol).lfwd, lbwd=last(sol).lbwd)
end

ensaugprob = EnsembleProblem(augprob,prob_func=prob_func, output_func=output_func)

ensaugsol = solve(ensaugprob, Tsit5(), trajectories = 100)

plot(1:100, i -> getindex(ensaugsol[i], :lfwd))

plot(1:100, i -> getindex(ensaugsol[i], :lbwd))
 =#